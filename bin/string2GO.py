#!/usr/bin/env python
docstring='''string2GO.py P75994
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
'''

import sys,os
import subprocess # this script uses zgrep to accelerate searching

from module import obo2csv
from module.fetch import wget
obo_url="http://geneontology.org/ontology/go-basic.obo"

datdir=os.path.abspath(os.path.dirname(__file__)+"/../dat")

full_uniprot_2_string=os.path.join(datdir,"full_uniprot_2_string.tsv.gz")
protein_links=os.path.join(datdir,"protein.links.txt.gz")
all_go_knowledge_explicit=os.path.join(datdir,"go_explicit.tsv.gz")

def uniprot2string(uniprotID="P75994"):
    '''map uniprot ID to a list of matching string ID'''
    cmd='zgrep "\\b%s\\b" %s'%(uniprotID,full_uniprot_2_string)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if not stdout.strip():
        return ''
    # reformat full_uniprot_2_string so that it only contains two column
    stringID_list=[line.split()[1] for line in stdout.splitlines()]
    return stringID_list

def get_PPI_partners(stringID_list=["b1167","4062745","ECDH10B_1216"]):
    '''get a dict whose key is interaction partners and value is combined_score'''
    partner_dict=dict()
    if not stringID_list:
        return partner_dict
    cmd='zgrep "\\b%s\\b" %s'%("\\|".join(stringID_list),protein_links)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if not stdout.strip():
        return partner_dict
    for line in stdout.splitlines():
        line=line.split()
        combined_score=float(line[2])/100.
        stringID_0=line[0].split('.')[-1]
        stringID_1=line[1].split('.')[-1]
        if stringID_0 in stringID_list:
            partner_dict[stringID_1]=combined_score
        elif stringID_1 in stringID_list:
            partner_dict[stringID_0]=combined_score
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

    stringID_list=uniprot2string(argv[0])
    sys.stdout.write("mapped to string ID: %s\n"%(','.join(stringID_list)))

    partner_dict=get_PPI_partners(stringID_list)
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
