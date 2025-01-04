# Durusdinium Transcriptome Annotation, version September 30, 2024
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (@gmail.com) for use on FAU's HPC (KoKo)


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
# setup

# To install Bioperl in your bin directory, please follow these instructions:
cd bin
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/z0on/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes
rm launcher_creator.py

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/mstudiva/Durusdinium-annotated-transcriptome.git
mv Durusdinium-annotated-transcriptome/* .
rm -rf Durusdinium-annotated-transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
# getting transcriptomes

# Durusdinium (Shoguchi; November 2020)
wget https://marinegenomics.oist.jp/symbd/download/102_symbd_transcriptome_nucl.fa.gz
gzip -d 102_symbd_transcriptome_nucl.fa.gz
mv 102_symbd_transcriptome_nucl.fa Durusdinium.fasta

# Durusdinium (Camp; June 2021)
# this transcriptome has already been converted to protein translations, so skip to line 168
# download from https://osf.io/gsn8p/ and scp to your annotate directory
mv D1a.annotated.fa Durusdinium.fasta
sed -i 's/TRINITY_DN/Durusdinium/g' Durusdinium.fasta

# Durusdinium (Ladner; Feb 2013)
# from https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=GAFP01
gunzip GAFP01.1.fsa_nt.gz
mv GAFP01.1.fsa_nt Durusdinium.fasta

# Test each of the transcriptomes on a subset of our sequenced samples, and pick the one with the best alignment rates
# See repository mstudiva/tag-based_RNAseq/tagSeq_processing_README.txt, L164-168 for details

# transcriptome statistics
conda activate bioperl
echo "seq_stats.pl Durusdinium.fasta > seqstats_Durusdinium.txt" >> seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Durusdinium.fasta (Shoguchi)
-------------------------
122923 sequences.
1097 average length.
11620 maximum length.
201 minimum length.
N50 = 1624
134.8 Mb altogether (134792983 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
# uniprot annotations with blast

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# unzipping
gunzip uniprot_sprot.fasta.gz &

# indexing the fasta database
module load blast-plus-2.11.0-gcc-9.2.0-5tzbbls
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch mdb.slurm

# splitting the transcriptome into 200 chunks, or however many is needed to keep the number of seqs per chunk under 1000
splitFasta.pl Durusdinium.fasta 200

# blasting all 200 chunks to uniprot in parallel, 4 cores per chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l
# you should end up with the same number of queries as sequences from the seq_stats script

# combining all blast results
cat subset*br > myblast.br
rm subset*

# creating isogroup lookups
grep ">" Durusdinium.fasta | perl -pe 's/>Durusdinium(\d+)(\S+)\s.+/Durusdinium$1$2\tDurusdinium$1/' > Durusdinium_seq2iso.tab
cat Durusdinium.fasta | perl -pe 's/>Durusdinium(\d+)(\S+).+/>Durusdinium$1$2 gene=Durusdinium$1/' > Durusdinium_iso.fasta


#-------------------------
# extracting coding sequences and corresponding protein translations:
conda activate bioperl # if not active already
echo "perl ~/bin/CDS_extractor_v2.pl Durusdinium_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch cddd


#------------------------------
# GO annotation
# updated based on Misha Matz's GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# selecting the longest contig per isogroup (also renames using isogroups based on Durusdinium annotations):
fasta2SBH.pl Durusdinium_iso_PRO.fas >Durusdinium_out_PRO.fas

# scp your *_out_PRO.fas file to laptop, submit it to
http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*_out_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# Durusdinium (Shoguchi) status: http://eggnog-mapper.embl.de/job_status?jobname=MM_8i8imc8h

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_8i8imc8h/out.emapper.annotations # Durusdinium (Shoguchi)

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Durusdinium_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Durusdinium_iso2geneName.tab


#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Durusdinium_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Durusdinium_iso2kogClass1.tab > Durusdinium_iso2kogClass.tab


#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
srun fasta2SBH.pl Durusdinium_iso.fasta >Durusdinium_4kegg.fasta

# scp *4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\*4kegg.fasta .
# use web browser to submit _4kegg.fasta file to KEGG's KAAS server (http://www.genome.jp/kegg/kaas/)
# select SBH method, upload nucleotide query
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1727803452&key=E12GTnNw # Durusdinium (Shoguchi)

# once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1727803452/query.ko # Durusdinium (Shoguchi)

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Durusdinium_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
# Converting isogoup lookup tables to orthogroups (for mixed symbiont assemblages)

srun perl ~/bin/rename_isogroup_by_orthogroup.pl orthologs_unique.txt Durusdinium_iso2geneName.tab Durusdinium_ortho2geneName.tab
srun perl ~/bin/rename_isogroup_by_orthogroup.pl orthologs_unique.txt Durusdinium_iso2go.tab Durusdinium_ortho2go.tab
srun perl ~/bin/rename_isogroup_by_orthogroup.pl orthologs_unique.txt Durusdinium_iso2kegg.tab Durusdinium_ortho2kegg.tab
srun perl ~/bin/rename_isogroup_by_orthogroup.pl orthologs_unique.txt Durusdinium_iso2kogClass1.tab Durusdinium_ortho2kogClass1.tab
srun perl ~/bin/rename_isogroup_by_orthogroup.pl orthologs_unique.txt Durusdinium_iso2kogClass.tab Durusdinium_ortho2kogClass.tab
srun perl ~/bin/rename_seq2iso_by_orthogroup.pl Durusdinium_seq2iso.tab orthologs_unique.txt Durusdinium_seq2ortho.tab


#------------------------------
# file transfer

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/\* .
