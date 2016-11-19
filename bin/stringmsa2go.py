#!/usr/bin/env python
docstring='''stringmsa2GO.py seq.fasta string.msa
    predict GO from STRING protein-protein interaction network.
    
    protein2    confidence_of_interaction  GO_term   confidence_of_GO
    partner1             991           GO:0000001        5
    partner2             356           GO:0000001        2
    partner3             474           GO:0000002        2

    (confidence of interaction is in [0,1000])
    (confidence of GO is in [1,5])

    all partner2 interacts with protein1 whose seqID to query is 0.9

Confidence Score of GO:0000001
[1] GOfreq (number of partners annotated with GO divided by partner number):
    Cscore(GO:0000001)=(0.9+0.9)/3
[2] swGOfreq (averaged confidence of interaction):
    Cscore(GO:0000001)=(0.9*0.991+0.9*0.356)/2
[3] gwGOfreq (averaged confidence of GO annotation):
    Cscore(GO:0000001)=(0.9*5/5+0.9*2/5)/2

options:
    -excludeGO=GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575
        GO terms to be excluded. (default is root and protein binding)
    -seqID_cutoff=0.9
        minimum seqID for blastp hit
'''

import sys,os
import subprocess # this script uses zgrep to accelerate searching
import shelve
import cPickle as pickle

from module import obo2csv
from module.fetch import wget
obo_url="http://geneontology.org/ontology/go-basic.obo"

bindir=os.path.abspath(os.path.dirname(__file__))
blastp_exe=os.path.join(bindir,"blastp")

datdir=os.path.abspath(bindir+"/../dat")

protein_links=os.path.join(datdir,"protein.links.shelve")
all_go_knowledge_explicit=os.path.join(datdir,"all_go_knowledge_explicit.pkl")

def mapFASTA2string(query="seq.fasta",msa_file='string.msa',
    seqID_cutoff=0.9):
    '''parse msa_file, return a dict whose key is stringID and value
    is seqID
    '''
    stringID_dict=dict()

    ## read query sequence ##
    fp=open(query,'rU')
    sequence=''.join([line.strip() for line in fp.read().splitlines() \
        if not line.startswith('>')])
    fp.close()
    seqlen=len(sequence)

    ## read msa ##
    fp=open(msa_file)
    txt=fp.read()
    fp.close()
    for line in txt.splitlines():
        if not line.startswith('>'):
            continue
        line=line[1:].strip().split()
        seqID=line[-1]
        identical_residue_num=float(seqID.split('/')[0])
        aligned_residue_num=float(seqID.split('/')[1])
        globalID=identical_residue_num/seqlen
        if globalID>=seqID_cutoff:
            stringID_dict[line[0]]=globalID

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

def get_PPI_partners(stringID_dict={"b1167":1,"4062745":1}):
    '''
    get a dict whose key is interaction partners and value is combined_score.
    if "PPI_partners_file" is present, read the PPI partner from that file
    instead.
    '''
    partner_dict=dict()
    if not stringID_dict:
        return partner_dict

    sparse_db=shelve.open(protein_links,'r')
    for stringID in stringID_dict:
        species=stringID.split('.')[0]
        if not species in sparse_db:
            sys.stderr.write("skip unshelved species %s\n"%species)
            continue
        if not stringID in sparse_db[species]['stringID_dict']:
            sys.stderr.write("skip unshelved protein1 %s\n"%stringID)
        stringID_idx=sparse_db[species]['stringID_dict'][stringID]
        if not stringID_idx in sparse_db[species]['PPI_matrix']:
            sys.stderr.write("skip protein1 %s without PPI\n"%stringID)
            continue
        partner_for_stringID_dict=sparse_db[species]['PPI_matrix'
            ][stringID_idx]
        stringID_list=sparse_db[species]['stringID_list']
        for partner_idx in partner_for_stringID_dict:
            partnerID=stringID_list[partner_idx]
            partner_dict[partnerID]=partner_for_stringID_dict[partner_idx
                ]/1000.*stringID_dict[stringID]
        #sys.stdout.write("%s partners for string ID: %s\n"%(
            #len(partner_for_stringID_dict),stringID))
    sparse_db.close()
    return partner_dict

def get_partner_terms(partner_dict=dict(),excludeGO='',obo_dict=dict()):
    '''get mapping from GO terms to string ID'''
    partner2GO_dict={"F":dict(),"P":dict(),"C":dict()}
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    excludeGO=set(excludeGO)

    fp=open(all_go_knowledge_explicit,'rb')
    all_go_knowledge_dict=pickle.load(fp)
    fp.close()

    ## map GO term (and parent terms) to partners ##
    for partner_ID in partner_dict:
        if not partner_ID in all_go_knowledge_dict:
            sys.stderr.write("skip protein2 %s without GOterm\n"%partner_ID)
            continue
        for Term in all_go_knowledge_dict[partner_ID]:
            if Term in excludeGO:
                continue
            Cscore=all_go_knowledge_dict[partner_ID][Term]
            Cscore=(float(Cscore)+5)/10
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
    seqID_cutoff=0.9 # seqID cutoff for blastp

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-excludeGO="):
            excludeGO=arg[len("-excludeGO="):]
        elif arg.startswith("-seqID_cutoff="):
            seqID_cutoff=float(arg[len("-seqID_cutoff="):])
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<2:
        sys.stderr.write(docstring)
        exit()

    #### parse GO hierachy ####
    fp=open(wget(obo_url,show_url=True),'rU')
    obo_txt=fp.read()
    fp.close()
    obo_dict=obo2csv.parse_obo_txt(obo_txt)

    #### parse string data ####
    stringID_dict=mapFASTA2string(argv[0],argv[1],seqID_cutoff)
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
