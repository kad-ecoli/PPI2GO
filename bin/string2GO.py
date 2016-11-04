#!/usr/bin/env python
docstring='''string2GO.py P75994
string2GO.py -infmt=fasta seq.fasta

    predict GO from STRING protein-protein interaction network.
    
    entry    confidence_of_interaction  GO_term   confidence_of_GO
    partner1             991           GO:0000001        5
    partner2             356           GO:0000001        2
    partner3             474           GO:0000002        2

    confidence of interaction is in [0,1000]
    confidence of GO is in [1,5]

Confidence Score of GO:0000001
[1] GOfreq (number of partners annotated with GO divided by partner number):
    Cscore(GO:0000001)=2/3
[2] swGOfreq (averaged confidence of interaction):
    Cscore(GO:0000001)=(0.991+0.356)/2
[3] gwGOfreq (averaged confidence of GO annotation):
    Cscore(GO:0000001)=(5/5+2/5)/2

options:
    -excludeGO=GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575
        GO terms to be excluded. (default is root and protein binding)
    -infmt={accession,fasta} input format
        accession: (default) input is a uniprot accession. this accession
               will be mapped to string ID using lookup table from STRING
        fasta: input is a fasta sequence. this sequence will be mapped to
               string ID by blast
'''

import sys,os
import subprocess # this script uses zgrep to accelerate searching

from module import obo2csv
from module.fetch import wget
obo_url="http://geneontology.org/ontology/go-basic.obo"

bindir=os.path.abspath(os.path.dirname(__file__))
blastp_exe=os.path.join(bindir,"blastp")

datdir=os.path.abspath(bindir+"/../dat")

protein_sequences=os.path.join(datdir,"protein.sequences.fa")
full_uniprot_2_string=os.path.join(datdir,"full_uniprot_2_string.tsv.gz")
protein_links=os.path.join(datdir,"protein.links.txt.gz")
all_go_knowledge_explicit=os.path.join(datdir,"go_explicit.tsv.gz")

# only when input is fasta sequence
evalue_cutoff="0.001" # evalue cutoff for blastp
seqID_cutoff=0.6      # seqID cutoff for blastp

def mapFASTA2string(query="seq.fasta",stringID_list_file='',
    db=protein_sequences):
    '''using blastp to search fasta "query" against database "db"
    If "stringID_list_file" is present, read the string ID mapping from that
    file instead
    '''
    stringID_dict=dict()

    if stringID_list_file and os.path.isfile(stringID_list_file):
        fp=open(stringID_list_file,'rU')
        txt=fp.read()
        fp.close()
        for line in txt.splitlines():
            line=line.split()
            seqID=float(line[1])
            if seqID>=seqID_cutoff:
                stringID_dict[line[0]]=seqID
        return stringID_dict

    fp=open(query,'rU')
    sequence=''.join([line.strip() for line in fp.read().splitlines() \
        if not line.startswith('>')])
    fp.close()
    seqlen=len(sequence)
    stdin=">query\n%s\n"%sequence

    cmd="%s -query %s -db %s -evalue %s -outfmt 6"%(
        blastp_exe,query,db,evalue_cutoff)
    p=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,
        stdout=subprocess.PIPE)
    stdout,stderr=p.communicate(input=stdin)

    txt=''
    for line in stdout.splitlines():
        line=line.split()
        seqID=float(line[2])*0.01*float(line[3])/seqlen
        txt+="%s\t%.4f\n"%(line[1],seqID)
        if seqID>=seqID_cutoff:
            stringID_dict[line[1]]=seqID

    if stringID_list_file:
        fp=open(stringID_list_file,'w')
        fp.write(txt)
        fp.close()
    return stringID_dict

def uniprot2string(uniprotID="P75994"):
    '''map uniprot ID to a list of matching string ID'''
    cmd='zgrep "\\b%s\\b" %s'%(uniprotID,full_uniprot_2_string)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if not stdout.strip():
        return ''
    # reformat full_uniprot_2_string so that it only contains two column
    for line in stdout.splitlines():
        stringID_dict[line.split()[1]]=1. # perfect match
    return stringID_dict

def get_PPI_partners(stringID_dict={"b1167":1,"4062745":1,"ECDH10B_1216":1},
    PPI_partners_file=''):
    '''
    get a dict whose key is interaction partners and value is combined_score.
    if "PPI_partners_file" is present, read the PPI partner from that file
    instead.
    '''
    partner_dict=dict()
    if PPI_partners_file and os.path.isfile(PPI_partners_file):
        fp=open(PPI_partners_file,'rU')
        txt=fp.read()
        fp.close()
        for line in txt.splitlines():
            line=line.split()
            PPI_partners_dict[line[0]]=float(line[1])
        return partner_dict

    if not stringID_dict:
        return partner_dict

    cmd='zgrep "\\b%s\\b" %s'%("\\|".join(stringID_dict),protein_links)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if not stdout.strip():
        return partner_dict
    for line in stdout.splitlines():
        line=line.split()
        combined_score=float(line[2])/1000.
        stringID_0=line[0].split('.')[-1]
        stringID_1=line[1].split('.')[-1]
        if stringID_0 in stringID_dict:
            partner_dict[stringID_1]=combined_score*stringID_dict[stringID_0]
        elif stringID_1 in stringID_dict:
            partner_dict[stringID_0]=combined_score*stringID_dict[stringID_1]
    if PPI_partners_file:
        fp=open(PPI_partners_file,'w')
        fp.write(''.join(["%s\t%.4f\n"%(stringID,partner_dict[stringID]
            ) for stringID in partner_dict]))
        fp.close()
    return partner_dict

def get_partner_terms(partner_dict=dict(),excludeGO='',obo_dict=dict()):
    '''get mapping from GO terms to string ID'''
    partner2GO_dict={"F":dict(),"P":dict(),"C":dict()}
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    excludeGO=set(excludeGO)
    cmd='zgrep "\\b%s\\b" %s'%("\\|".join(partner_dict.keys()),
        all_go_knowledge_explicit)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if not stdout.strip():
        return partner2GO_dict

    ## map GO term (and parent terms) to partners ##
    for line in stdout.splitlines():
        line=line.split() # stringID GOterm cscore
        partner_ID=line[0]
        Term=line[1]
        Cscore=(float(line[2])+5)/10
        if not partner_ID in partner_dict or Term in excludeGO:
            continue
        for Aspect in partner2GO_dict:
            if not Term in obo_dict[Aspect]["Term"]:
                continue
            Term_list=list(set(obo_dict.is_a(Term,direct=False,
                name=False,number=False).split())-excludeGO)
            if not partner_ID in partner2GO_dict:
                partner2GO_dict[Aspect][partner_ID]=dict()

            for Term in Term_list:
                if not Term in partner2GO_dict[Aspect][partner_ID] or \
                    partner2GO_dict[Aspect][partner_ID][Term]<Cscore:
                    partner2GO_dict[Aspect][partner_ID][Term]=Cscore
    return partner2GO_dict

def partner2GOfreq(Aspect='F',partner_dict=dict(),partner2GO_dict=dict()):
    '''get final confidence score for a specific Aspect by following scoring:
    GOfreq, swGOfreq, gwGOfreq, dwGOfreq
    '''
    GOfreq_pred_list=[]
    swGOfreq_pred_list=[]
    gwGOfreq_pred_list=[]

    GO_list_cat=[] # list of all GO terms annotated to all partners
    for partner_ID in partner2GO_dict[Aspect].keys():
        GO_list_cat+=partner2GO_dict[Aspect][partner_ID].keys()
    GO_list_cat=list(set(GO_list_cat)) # non-redundant list
    
    for Term in GO_list_cat:
        GOfreq_cscore=sum([(1. \
            if Term in partner2GO_dict[Aspect][partner_ID] else 0.
            ) for partner_ID in partner2GO_dict[Aspect]]
            )/len(partner2GO_dict[Aspect])
        swGOfreq_cscore=sum([(partner_dict[partner_ID] \
            if Term in partner2GO_dict[Aspect][partner_ID] else 0.
            ) for partner_ID in partner2GO_dict[Aspect]]
            )/len(partner2GO_dict[Aspect])
        gwGOfreq_cscore=sum([(partner2GO_dict[Aspect][partner_ID][Term] \
            if Term in partner2GO_dict[Aspect][partner_ID] else 0.
            ) for partner_ID in partner2GO_dict[Aspect]]
            )/len(partner2GO_dict[Aspect])

        GOfreq_pred_list.append((GOfreq_cscore,Term,Aspect))
        swGOfreq_pred_list.append((swGOfreq_cscore,Term,Aspect))
        gwGOfreq_pred_list.append((gwGOfreq_cscore,Term,Aspect))
    
    GOfreq_pred=''.join(["%s\t%s\t%.3f\n"%(Term,Aspect,cscore) for \
        cscore,Term,Aspect in sorted(GOfreq_pred_list,reverse=True)])
    swGOfreq_pred=''.join(["%s\t%s\t%.3f\n"%(Term,Aspect,cscore) for \
        cscore,Term,Aspect in sorted(swGOfreq_pred_list,reverse=True)])
    gwGOfreq_pred=''.join(["%s\t%s\t%.3f\n"%(Term,Aspect,cscore) for \
        cscore,Term,Aspect in sorted(gwGOfreq_pred_list,reverse=True)])
    return GOfreq_pred,swGOfreq_pred,gwGOfreq_pred

if __name__=="__main__":
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575"

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-excludeGO="):
            excludeGO=arg[len("-excludeGO="):]
        elif arg.startswith("-infmt="):
            infmt=arg[len("-infmt="):]
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<1:
        sys.stderr.write(docstring)
        exit()

    #### parse GO hierachy ####
    fp=open(wget(obo_url,show_url=True),'rU')
    obo_txt=fp.read()
    fp.close()
    obo_dict=obo2csv.parse_obo_txt(obo_txt)

    #### parse string data ####
    for entry in argv:
        sys.stdout.write("parsing entry %s\n"%entry)
        if infmt=="fasta":
            stringID_list_file=os.path.join(os.path.dirname(argv[0]),"stringID.list")
            stringID_dict=mapFASTA2string(argv[0],stringID_list_file)
        else:
            stringID_dict=uniprot2string(argv[0])
        sys.stdout.write("mapped to string ID: %s\n"%(','.join(stringID_dict)))

        partner_dict=get_PPI_partners(stringID_dict)
        sys.stdout.write("mapped to %d partners\n"%len(partner_dict))

        partner2GO_dict=get_partner_terms(partner_dict,excludeGO,obo_dict)
        for Aspect in ["MF","BP","CC"]:
            sys.stdout.write("Calculating GO-%s from %d partners\n"%(Aspect,
                len(partner2GO_dict[Aspect[1]])))
            GOfreq_pred,swGOfreq_pred,gwGOfreq_pred=partner2GOfreq(
                Aspect[1],partner_dict,partner2GO_dict)

            fp=open("string_GOfreq_"+Aspect,'w')
            fp.write(GOfreq_pred)
            fp.close()

            fp=open("string_swGOfreq_"+Aspect,'w')
            fp.write(swGOfreq_pred)
            fp.close()

            fp=open("string_gwGOfreq_"+Aspect,'w')
            fp.write(gwGOfreq_pred)
            fp.close()
