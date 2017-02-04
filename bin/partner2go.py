#!/usr/bin/env python
docstring='''
partner2GO.py seq.fasta partner.msa string_partner.pkl
    predict GO from STRING protein-protein interaction network.
    
    protein2    confidence_of_interaction  GO_term   confidence_of_GO
    partner1             991           GO:0000001        5
    partner2             356           GO:0000001        2
    partner3             474           GO:0000002        2

    (confidence of interaction is in [0,1000])
    (confidence of GO is in [1,5])

    all protein2 interacts with protein1 whose seqID to query is 0.9.
    seqID(query,protein1)=N_iden/L_query

stringID2GO.py seq.fasta string.msa string.full
    seqID(query,protein1)=N_iden/max{L_query,L_protein1}
    where "string.full" is the full length sequence of matched string ID.

Confidence Score of GO:0000001
[1] GOfreq (number of partners annotated with GO divided by partner number):
    Cscore(GO:0000001)=0.9*(2/3)
[2] swGOfreq (averaged confidence of interaction):
    Cscore(GO:0000001)=0.9*(0.991+0.356)/(0.991+0.356+0.474)
[3] gwGOfreq (averaged confidence of GO annotation):
    Cscore(GO:0000001)=0.9*(0.991*5/5+0.356*2/5)/(0.991+0.356+0.474)

options:
    -excludeGO=GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575
        GO terms to be excluded. (default is root and protein binding)
    -homoflag=0.3
        maximum seqID for funtional templates
'''

import sys,os
import subprocess # this script uses zgrep to accelerate searching
import shelve
import cPickle as pickle

from module import obo2csv
from module.fetch import wget
obo_url="http://geneontology.org/ontology/go-basic.obo"

bindir=os.path.abspath(os.path.dirname(__file__))
datdir=os.path.abspath(bindir+"/../dat")

#protein_links=os.path.join(datdir,"protein.links.shelve")
protein_links="protein.links.shelve" # copy shelve file to current tmpdir
all_go_knowledge_explicit=os.path.join(datdir,"all_go_knowledge_explicit.pkl")

string_partner_file="string_partner.pkl"

def get_homo_partner(query="seq.fasta",msa_file="partner.msa",
    homoflag=1):
    '''parse msa_file, return a list of entries in "msa_file" with greater
    than "homoflag" seqID
    '''
    homo_partner_list=[]

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
        header=line[0]
        seqID=line[-1]
        identical_residue_num=float(seqID.split('/')[0])
        aligned_residue_num=float(seqID.split('/')[1])
        globalID=identical_residue_num/seqlen
        if globalID>=homoflag:
            homo_partner_list.append(header)
    return homo_partner_list

def get_PPI_partners(stringID_dict={"b1167":1,"4062745":1}):
    '''
    get a dict whose key is interaction partners and value is combined_score.
    if "PPI_partners_file" is present, read the PPI partner from that file
    instead.
    '''
    partner_dict=dict()
    if not stringID_dict:
        return partner_dict

    fp=open(string_partner_file,'rb')
    partner_dict=pickle.load(fp)
    fp.close()
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

def partner2GOfreq(Aspect='F',partner_dict=dict(),partner2GO_dict=dict(),
    globalID=1.):
    '''get final confidence score for a specific Aspect by following scoring:
    GOfreq, swGOfreq

    globalID is the global seqID between query and the matched string entry
    '''
    GOfreq_pred_list=[]
    swGOfreq_pred_list=[]
    gwGOfreq_pred_list=[]

    GO_list_cat=[] # list of all GO terms annotated to all partners
    for partner_ID in partner2GO_dict[Aspect].keys():
        GO_list_cat+=partner2GO_dict[Aspect][partner_ID].keys()
    GO_list_cat=list(set(GO_list_cat)) # non-redundant list

    total_string_score=sum([partner_dict[partner_ID] for partner_ID in \
        partner2GO_dict[Aspect]]) # sum of string score for all partners

    for Term in GO_list_cat:
        GOfreq_cscore=globalID*sum([(1. \
            if Term in partner2GO_dict[Aspect][partner_ID] else 0.
            ) for partner_ID in partner2GO_dict[Aspect]]
            )/len(partner2GO_dict[Aspect])
        swGOfreq_cscore=globalID*sum([(partner_dict[partner_ID] \
            if Term in partner2GO_dict[Aspect][partner_ID] else 0.
            ) for partner_ID in partner2GO_dict[Aspect]]
            )/total_string_score
            #)/len(partner2GO_dict[Aspect])
        gwGOfreq_cscore=globalID*sum([(partner_dict[partner_ID] \
            *partner2GO_dict[Aspect][partner_ID][Term] \
            if Term in partner2GO_dict[Aspect][partner_ID] else 0.
            ) for partner_ID in partner2GO_dict[Aspect]]
            )/total_string_score
            #)/len(partner2GO_dict[Aspect])

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

def filter_partner_dict(partner_dict,homo_partner_list):
    '''remove key in partner_dict that is also in homo_partner_list'''
    for partner in homo_partner_list:
        del partner_dict[partner]
    return partner_dict

if __name__=="__main__":
    excludeGO="GO:0005515,GO:0005488,GO:0003674,GO:0008150,GO:0005575"
    seqID_cutoff=0.9 # seqID cutoff for blastp
    homoflag=1.1     # real - keep all functional templates

    argv=[]
    for arg in sys.argv[1:]:
        if arg.startswith("-excludeGO="):
            excludeGO=arg[len("-excludeGO="):]
        elif arg.startswith("-homoflag="):
            homoflag=arg[len("-homoflag="):]
        elif arg.startswith("-"):
            sys.stderr.write("ERROR! Unknown option %s\n"%arg)
            exit()
        else:
            argv.append(arg)

    if len(argv)<3:
        sys.stderr.write(docstring)
        exit()

    if homoflag=="real":
        homoflag=1.1
    elif homoflag=="benchmark":
        homoflag=0.3
    else:
        try:
            homoflag=float(homoflag)
        except Exception,error:
            sys.stderr.write("ERROR! Unknown homoflag value %s\n%s\n"%(
                str(homoflag),error))
            exit()

    #### parse GO hierachy ####
    fp=open(wget(obo_url,show_url=True),'rU')
    obo_txt=fp.read()
    fp.close()
    obo_dict=obo2csv.parse_obo_txt(obo_txt)

    #### get homologous partners ####
    homo_partner_list=[]
    if homoflag<1:
        homo_partner_list=get_homo_partner(argv[0],argv[1],homoflag)

    fp=open(argv[2],'rb')
    obj=pickle.load(fp)
    fp.close()
    partner_dict,stringID_dict=obj
    sys.stdout.write("mapped to %d partners\n"%len(partner_dict))
    globalID=1
    if stringID_dict:
        globalID=sorted(stringID_dict.values())[0]
    sys.stdout.write("mapped to string ID: %s\n"%(','.join(stringID_dict)))

    partner_dict=filter_partner_dict(partner_dict,homo_partner_list)

    partner2GO_dict=get_partner_terms(partner_dict,excludeGO,obo_dict)
    for Aspect in ["MF","BP","CC"]:
        sys.stdout.write("Calculating GO-%s from %d partners\n"%(Aspect,
            len(partner2GO_dict[Aspect[1]])))
        #GOfreq_pred,swGOfreq_pred=partner2GOfreq(
        GOfreq_pred,swGOfreq_pred,gwGOfreq_pred=partner2GOfreq(
            Aspect[1],partner_dict,partner2GO_dict,globalID)

        fp=open("string_GOfreq_"+Aspect,'w')
        fp.write(GOfreq_pred)
        fp.close()

        fp=open("string_swGOfreq_"+Aspect,'w')
        fp.write(swGOfreq_pred)
        fp.close()

        fp=open("string_gwGOfreq_"+Aspect,'w')
        fp.write(gwGOfreq_pred)
        fp.close()
