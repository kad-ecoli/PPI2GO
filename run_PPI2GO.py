#!/usr/bin/env python
docstring='''
run_PPI2GO.py config.py
    run PPI2GO protein function prediction

input file:
    config.py (configuration file in the following format:

        seq="seq.txt" 
            # fasta file for all queries
        outdir="."
            # output directory, default value is directory of
            # "seq.txt". Each prediction target should has its
            # own folder. Under each folder, there should be a
            # fasta format sequence file called "seq.fasta"
        evalue=0.001
            # blastp and psiblast evalue cutoff for function transfer
        Q="default"
            # queue destination
        run="real"
            # if "real", preserve all templates
            # if "benchmark", remove templates sharing >=0.3 seqID
            # if set to a number, remove templates sharing
            # at least specified seqID

output file:
    blastp.xml.gz       (XML format blastp output)
    blastp.msa*         (reorgnized MSA)
    blastp_globalID_*   (CSV format GO prediction file using global seqID)
    blastp_localID_*    (CSV format GO prediction file using local seqID)
    blastp_GOfreq_*     (CSV format GO prediction file using GO frequency)

    psiblast.xml.gz     (XML format PSI-BLAST output)
    psiblast.msa*       (reorgnized MSA)
    psiblast_globalID_* (CSV format GO prediction file using global seqID)
    psiblast_localID_*  (CSV format GO prediction file using local seqID)
    psiblast_GOfreq_*   (CSV format GO prediction file using GO frequency)

    combine_GOfreq_*    (final GO prediction by combining PSIBLAST and BLAST)
'''
import os,sys
from module import jobsubmit
from module.split_fasta import split_fasta
from string import Template

if len(sys.argv)!=2:
    sys.stderr.write(docstring)
    exit()
#### docstring END #################

#### parse input START #############
bindir=os.path.join(os.path.dirname(os.path.abspath(__file__)),"bin")
# path to sequence database "lib.fasta" and GO mapping
datdir=os.path.join(os.path.dirname(os.path.abspath(__file__)),"dat")

seq=''        # fasta file for all sequences
outdir='.'  # input directory
evalue=0.001   # blastp evalue cutoff
run='real'

config=sys.argv[1] # config file
if not os.path.isfile(config):
    sys.stderr.write("ERROR! "+config+": No such file.\n")
    exit()
execfile(config)

if not os.path.isfile(seq) and not os.path.isdir(outdir):
    sys.stderr.write("ERROR! Cannot find input FASTA sequence %s\n"%seq)
    exit()

if not outdir:
    outdir=os.path.dirname(seq)
    sys.stderr.write("'outdir' not specified. Using %s\n"%outdir)
outdir=os.path.abspath(outdir)
if not seq:
    ss=[s for s in os.listdir(outdir) if \
        os.path.isdir(os.path.join(outdir,s,"seq.fasta"))]
else:
    ss=split_fasta(infile=seq,outfile="seq.fasta",outdir=outdir)
#### parse input END ###############

#### job template START ############
recorddir=os.path.join(outdir,'record')
if not os.path.isdir(recorddir):
    os.makedirs(recorddir)
jobmod=Template('''#!/bin/bash
#### prepare tmp directory ####
tmpdir="/scratch/$USER/$SLURM_JOBID"; # on COMET computing node
if [ ! -d "$tmpdir" ] || [ -z "$SLURM_JOBID" ] ; then
    tmpdir="/tmp/$USER"
fi
tmpdir="$tmpdir/$tag"

mkdir -p $tmpdir
rm -rf $tmpdir/*
cd     $tmpdir
cp $outdir/$s/seq.fasta .
cp -rp $bindir/* .

#### get blastp hits ####
if [ -s $outdir/$s/string.msa ];then
    cp $outdir/$s/string.msa .
else
    ./blastp -query seq.fasta             \\
         -db $datdir/protein.sequences.fa \\
         -out string.xml                  \\
         -evalue $evalue                  \\
         -outfmt 5
    ./blast2msa.py seq.fasta string.xml string.msa
    gzip string.xml
fi
grep -ohP '^>\S+' string.msa|sed 's/^>//g'|./blastdbcmd -entry_batch - \\
    -db $datdir/protein.sequences.fa -out string.full

#### convert blast hits to GO prediction ####
if [ -s $outdir/$s/go-basic.obo ];then
    cp $outdir/$s/go-basic.obo .
fi
$bindir/stringmsa2go.py seq.fasta string.msa string.full

#### copy result back ####
cp $tmpdir/string.xml* $outdir/$s/
cp $tmpdir/string.msa  $outdir/$s/
cp $tmpdir/string.full $outdir/$s/
cp $tmpdir/string_*_*  $outdir/$s/

#### clean up ####
rm -rf $tmpdir
''')
#### job template END ##############

#### Job submit START ##############
walltime="walltime=10:00:00,mem=10gb"
JOBS = jobsubmit.JOBS(walltime=walltime, priority=Q)
for s in ss:
    datadir=os.path.join(outdir,s)
    tag="PIG_"+s+'_'+str(run) # uniq name for job
    jobname=os.path.join(recorddir,tag)
    
    jobOutput=[os.path.join(datadir,"string_swGOfreq_MF"),
               os.path.join(datadir,"string_swGOfreq_BP"),
               os.path.join(datadir,"string_swGOfreq_CC")]
    
    mod=jobmod.safe_substitute(dict(
        #tmpdir=os.path.join("/tmp",os.getenv("USER"),tag),
        tag=tag,
        bindir=bindir,
        outdir=outdir,
        datdir=datdir,
        evalue=evalue,
        homoflag=run,
        s=s,
    ))

    fp=open(jobname,'w')
    fp.write(mod)
    fp.close()
    os.chmod(jobname, os.stat(jobname).st_mode|0111)
        
    JOBS.SubmitJob(jobInput=jobname, jobOutput=jobOutput)
