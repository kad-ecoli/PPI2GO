#!/usr/bin/env python
docstring='''string2GO.py P75994
    predict GO from STRING protein-protein interaction network

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
all_go_knowledge_explicit=os.path.join(datdir,"all_go_knowledge_explicit.tsv.gz")

def uniprot2string(uniprotID="P75994"):
    '''map uniprot ID to a list of matching string ID'''
    cmd='zgrep "\\b%s\\b" %s'%(uniprotID,full_uniprot_2_string)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if not stdout.strip():
        return ''
    stringID_list=[line.split()[2] for line in stdout.splitlines()]
    return stringID_list

def get_PPI_partners(stringID_list=["b1167","4062745","ECDH10B_1216"]):
    '''get a list of interaction partners'''
    partner_list=[]
    cmd='zgrep "\\b%s\\b" %s'%("\\|".join(stringID_list),protein_links)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if not stdout.strip():
        return partner_list
    for line in stdout.splitlines():
        line=[e.split('.')[1] for e in line.split()[:2]]
        if line[0] in stringID_list:
            partner_list.append(line[1])
        elif line[1] in stringID_list:
            partner_list.append(line[0])
    return partner_list

def get_partner_terms(partner_list=[],excludeGO='',obo_dict=dict()):
    '''get mapping from GO terms to string ID'''
    partner2GO_dict={"F":dict(),"P":dict(),"C":dict()}
    if isinstance(excludeGO,str):
        excludeGO=excludeGO.split(',')
    excludeGO=set(excludeGO)
    cmd='zgrep "\\b%s\\b" %s'%("\\|".join(partner_list),all_go_knowledge_explicit)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if not stdout.strip():
        return partner2GO_dict

    ## map GO term (and parent terms) to partners ##
    for line in stdout.splitlines():
        line=line.split()
        partner_ID=line[1]
        Term=line[3]
        if not partner_ID in partner_list or Term in excludeGO:
            continue
        for Aspect in partner2GO_dict:
            if not Term in obo_dict[Aspect]["Term"]:
                continue
            Term_list=obo_dict.is_a(
                Term,direct=False,name=False,number=False).split()
            if not partner_ID in partner2GO_dict:
                partner2GO_dict[Aspect][partner_ID]=Term_list
            else:
                partner2GO_dict[Aspect][partner_ID]+=Term_list

    ## remove repeative GO terms ##
    for Aspect in partner2GO_dict:
        for partner_ID in partner2GO_dict[Aspect]:
             partner2GO_dict[Aspect][partner_ID]=list(
                set(partner2GO_dict[Aspect][partner_ID])-excludeGO)
    return partner2GO_dict

def partner2GOfreq(Aspect='F',partner2GO_dict=dict()):
    '''get final confidence score for a specific Aspect'''
    GOfreq_pred=''

    GO_list_cat=[] # list of all GO terms annotated to all partners
    for partner_ID in partner2GO_dict[Aspect]:
        GO_list_cat+=partner2GO_dict[Aspect][partner_ID]
    
    for Term in set(GO_list_cat):
        cscore=1.*sum([(Term in partner2GO_dict[Aspect][partner_ID]) for \
            partner_ID in partner2GO_dict[Aspect]])/len(partner2GO_dict[Aspect])
        GOfreq_pred+="%s\t%s\t%.2f\n"%(Term,Aspect,cscore)
    return GOfreq_pred

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

    partner_list=get_PPI_partners(stringID_list)
    sys.stdout.write("mapped to %d partners\n"%len(partner_list))

    partner2GO_dict=get_partner_terms(partner_list,excludeGO,obo_dict)
    for Aspect in ["MF","BP","CC"]:
        sys.stdout.write("Calculating GO-%s from %d partners\n"%(Aspect,
            len(partner2GO_dict[Aspect[1]])))
        fp=open("string_GOfreq_"+Aspect,'w')
        fp.write(partner2GOfreq(Aspect[1],partner2GO_dict))
        fp.close()
