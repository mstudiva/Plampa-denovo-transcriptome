# Transcriptome annotation for the sponge Pione lampa, version December 13, 2021
# Adapted by Michael Studivan (studivanms@gmail.com) based on a repo by Misha Matz (https://github.com/z0on/annotatingTranscriptomes.git) for use on FAU's HPC (KoKo)


#------------------------------
## BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
## Script and workplace setup

# To install Bioperl in your bin directory, please follow these instructions:
cd bin
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/mstudiva/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/mstudiva/Plampa-denovo-transcriptome.git
mv Plampa-denovo-transcriptome/* .
rm -rf Plampa-denovo-transcriptome/

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
## Getting transcriptome

# Pione lampa transcriptome (December 2021)
wget https://www.dropbox.com/s/qkelwwnm5rjbvvr/Plampa.fasta

# use the stream editor to find and replace all instances of "TRINITY" with "Plampa" in the  transcriptome
sed -i 's/TRINITY_DN/Plampa/g' Plampa.fasta

# transcriptome statistics
echo "seq_stats.pl Plampa.fasta > seqstats_Plampa.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

nano seqstats_Plampa.txt

Plampa.fasta
-------------------------
220471 sequences.
1261 average length.
17411 maximum length.
500 minimum length.
N50 = 1536
277.9 Mb altogether (277944908 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## Uniprot annotations with blast

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# unzipping
gunzip uniprot_sprot.fasta.gz &

# indexing the fasta database
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch mdb.slurm

# splitting the transcriptome into 200 chunks
splitFasta.pl Plampa.fasta 200

# blasting all 200 chunks to uniprot in parallel, 4 cores per chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l
# you should end up with the same number of queries as sequences from the seq_stats script (220471 sequences)

# combining all blast results
cat subset*br > myblast.br
mv subset* ~/annotate/backup/

# generating isogroup designations for each contig/component (there can be multiple contigs per isogroup)
grep ">" Plampa.fasta | perl -pe 's/>Plampa(\d+)(\S+)\s.+/Plampa$1$2\tPlampa$1/'>Plampa_seq2iso.tab
cat Plampa.fasta | perl -pe 's/>Plampa(\d+)(\S+).+/>Plampa$1$2 gene=Plampa$1/'>Plampa_iso.fasta


#-------------------------
## Extracting coding sequences and corresponding protein translations:
echo "perl ~/bin/CDS_extractor_v2.pl Plampa_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch cddd


#------------------------------
## GO annotations

# selecting the longest contig per isogroup (also renames using isogroups based on Plampa  annotations):
fasta2SBH.pl Plampa_iso_PRO.fas >Plampa_out_PRO.fas

# scp your *_out_PRO.fas file to laptop, submit it to http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*_out_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# status: go on web to http://eggnog-mapper.embl.de/job_status?jobname=MM_t9rukmp6

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_t9rukmp6/out.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Plampa_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Plampa_iso2geneName.tab


#------------------------------
## KOG annotations

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Plampa_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Plampa_iso2kogClass1.tab > Plampa_iso2kogClass.tab


#------------------------------
## KEGG annotations

# selecting the longest contig per isogroup:
srun fasta2SBH.pl Plampa_iso.fasta >Plampa_4kegg.fasta

# scp *4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*4kegg.fasta .
# use web browser to submit Plampa_4kegg.fasta file to KEGG's KAAS server http://www.genome.jp/kegg/kaas/
# select SBH method, upload nucleotide query
https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1639448600&key=bZ5qj3k5

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1639448600/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Plampa_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
## Done! Transfer files

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/* .
