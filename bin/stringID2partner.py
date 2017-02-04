#!/usr/bin/env python
docstring='''
stringID2partner.py seq.fasta string.msa string_partner.txt
    map string ID from input sequence "seq.fasta" and blast msa "string.msa"
    to PPI partner "string_partner.txt"

options:
    -seqID_cutoff=0.9
        minimum seqID for blastp hit
'''

import sys,os
import subprocess # this script uses zgrep to accelerate searching
import shelve
import cPickle as pickle

bindir=os.path.abspath(os.path.dirname(__file__))
datdir=os.path.abspath(bindir+"/../dat")

#protein_links=os.path.join(datdir,"protein.links.shelve")
protein_links="protein.links.shelve" # copy shelve file to current tmpdir

string_partner_file="string_partner.pkl"

def mapFASTA2string(query="seq.fasta",msa_file="string.msa",
    seqID_cutoff=0.9,full_fasta=''):
    '''parse msa_file, return a dict whose key is stringID and value
    is seqID. if full_fasta is not empty, paste the full length fasta
    sequence file and use it for sequence length caculation
    '''
    stringID_dict=dict()

    ## read full length sequence file ##
    string_len_dict=dict()
    if os.path.isfile(full_fasta):
        fp=open(full_fasta,'rU')
        for block in fp.read().split('>'):
            if not block.strip():
                continue
            header=block.split()[0]
            sequence=''.join(block.splitlines()[1:])
            string_len_dict[header.split()[0]]=len(sequence)
        fp.close()

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
    max_globalID=0
    for line in txt.splitlines():
        if not line.startswith('>'):
            continue
        line=line[1:].strip().split()
        header=line[0]
        seqID=line[-1]
        identical_residue_num=float(seqID.split('/')[0])
        aligned_residue_num=float(seqID.split('/')[1])
        if header in string_len_dict:
            globalID=identical_residue_num/max(
                [seqlen,string_len_dict[header]])
        else:
            globalID=identical_residue_num/seqlen
        if globalID>=seqID_cutoff and globalID>=max_globalID:
            stringID_dict[header]=globalID # only take the top hit
            max_globalID=globalID
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

    if os.path.isfile(string_partner_file):
        fp=open(string_partner_file,'rb')
        obj=pickle.load(fp)
        fp.close()
        partner_dict,stringID_dict=obj
    else:
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
                partner_dict[partnerID]=partner_for_stringID_dict[
                    partner_idx]/1000. #*stringID_dict[stringID]
            #sys.stdout.write("%s partners for string ID: %s\n"%(
                #len(partner_for_stringID_dict),stringID))
        sparse_db.close()
        
        obj=[partner_dict,stringID_dict]
        fp=open(string_partner_file,'wb')
        pickle.dump(obj,fp)
        fp.close()
    return partner_dict

if __name__=="__main__":
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

    #### parse string data ####
    stringID_dict=mapFASTA2string(argv[0],argv[1],seqID_cutoff)
    globalID=1
    if stringID_dict:
        globalID=sorted(stringID_dict.values())[0]
    sys.stdout.write("mapped to string ID: %s\n"%(','.join(stringID_dict)))

    partner_dict=get_PPI_partners(stringID_dict)
    sys.stdout.write("mapped to %d partners\n"%len(partner_dict))

    txt=''.join([partner+'\n' for partner in partner_dict.keys()])

    if len(argv)>=3:
        fp=open(argv[2],'w')
        fp.write(txt)
        fp.close()
    else:
        sys.stdout.write(txt)
