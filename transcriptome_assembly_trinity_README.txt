## De novo transcriptome assembly  pipeline, version December 13, 2021
# Adapted by Michael Studivan (studivanms@gmail.com) and John Morris (john.morris@noaa.gov) based on repos by Misha Matz (https://github.com/z0on/annotatingTranscriptomes.git) and Brian Strehlow (https://github.com/bstrehlow/Transcriptome_assembly.git) for use on the FAU KoKo HPC


#------------------------------
## BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions – please make sure to read them before copy-pasting.

# log onto cluster
ssh mstudiva@koko-login.hpc.fau.edu


#------------------------------
## Installing RNA-seq scripts and setting up the workspace

# switch to home directory
cd

# unless you have done it in the past, make directory called bin,
# all your scripts should go in there:
mkdir bin

# switch to bin:
cd bin

# clone github repositories
git clone https://github.com/mstudiva/tag-based_RNAseq.git
git clone https://github.com/mstudiva/Plampa-denovo-transcriptome.git
git clone https://github.com/hrivera28/Oculina_arbuscula_transcriptome.git

# move files from subdirectories to bin/:
mv tag-based_RNAseq/* .
mv Plampa-denovo-transcriptome/* .
mv Oculina_arbuscula_transcriptome/Sarahs_scripts/* .

chmod +x ~/bin/bs

## Installing Trinity
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.13.2/trinityrnaseq-v2.13.2.FULL.tar.gz
tar -vxf trinityrnaseq-v2.13.2.FULL.tar.gz
cd trinityrnaseq-v2.13.2/
make install
make plugins
# test it
cd sample_data/test_Trinity_Assembly/
./runMe.sh

## Installing jellyfish
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
tar -vxf jellyfish-2.3.0.tar.gz
cd jellyfish-2.3.0
./configure
make
make install


#------------------------------
## Downloading files via ftp

sftp morris_6888@dnaseq2.igsp.duke.edu
get -r Morris2_6888_210518B6

## Get reference transcriptomes for de novo assembly (in progress)


#------------------------------
## Unzipping reads with a launcher_creator script

# creating and launching a cluster job to unzip all files:
 ls *.gz | perl -pe 's/(\S+)/gunzip $1/' > gunz
 launcher_creator.py -j gunz -n gunz -q shortq7 -t 6:00:00 -e studivanms@gmail.com
 sbatch --mem=200GB gunz.slurm

# look at the reads:
# head -50 SampleName.fastq
# note that every read has four lines, the ID line starts with @A00

# this little one-liner will show sequence-only in file:
# head -100 SampleName.fastq | grep -E '^[NACGT]+$'

# to count the number of reads in all samples
echo "countreads.pl > countreads_raw.txt" > count
launcher_creator.py -j count -n count -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count.slurm

#------------------------------
## Trimming and quality filtering

# Create conda environment
# Uncomment and run below if you don't have a conda env. set up
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n sctld cutadapt

## Removing adaptors and low quality reads
echo '#!/bin/bash' > trim.sh
echo 'conda activate sctld' >> trim.sh
for F in *.fastq; do
echo "cutadapt $F -b GGGGGGGG -b AGATCGG -q 15 -m 50 -o ${F/.fastq/}.trim" >>trim.sh;
done

# Does not work with launcher_creator, consider breaking up script and running multiple jobs
sbatch -o trim.o%j -e trim.e%j --mem=200GB trim.sh

# how the job is doing?
squeue -u mstudiva

# double check you have the same number of files as samples
ll *.trim | wc -l

# but did the trimming really work?
# Use the same one-liner as before on the trimmed file to see if it is different
# from the raw one that you looked at before:

# head -100 SampleName.fq | grep -E '^[NACGT]+$'

# head -100 SampleName.trim | grep -E '^[NACGT]+$'
# the long runs of base A should be gone

# to save time in case of issues, move the raw fastq files to backup directory
mv *.fastq ~/backup/

## to count the number of reads in trimmed samples
echo "countreads_trim.pl > countreads_trim.txt" > count_trim
launcher_creator.py -j count_trim -n count_trim -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_trim.slurm


#------------------------------
## Re-pairing, deduplicating, and assembling forward and reverse reads

# copy the re-pair.sh script from bin to current directory
cp ~/bin/re-pair.sh .
launcher_creator.py -j re-pair.sh -n re-pair -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch re-pair.slurm

cd ..
mkdir cleanReads
cd cleanReads
mv rawReads/Un_* .
mv rawReads/Un_* .

## De-duplicating paired reads (requires left [R1], right [R2], and unpaired read files)
cp ~/bin/dedup.sh .
launcher_creator.py -j dedup.sh -n dedup -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch dedup.slurm

## to count the number of reads in deduplicated samples
echo "countreads_dedup.pl > countreads_dedup.txt" > count_dedup
launcher_creator.py -j count_dedup -n count_dedup -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_dedup.slurm

## Assembling deduplicated reads (and renames)
cp ~/bin/assemble.sh .
launcher_creator.py -j assemble.sh -n assemble -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch assemble.slurm

# when it's done, move all the .dedup files to backup
mv *.dedup ../backup/


#------------------------------
## Transcriptome assembly using Trinity

# first, concatenate all forward and reverse reads into a single set of fastq files per species
echo "cat R1_Pione* > R1_Pione.fastq" > concat
echo "cat R2_Pione* > R2_Pione.fastq" >> concat
launcher_creator.py -j concat -n concat -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch concat.slurm

# create directories for each of the species in your annotate directory, and move the concatenated fastq files there
mkdir pione

# Trinity
# the python module required by launcher_creator is messing up the steps involving numpy, so we'll do it the old school way
echo '#!/bin/bash' > trinity.sh
echo '#SBATCH --partition=longq7' >> trinity.sh
echo '#SBATCH -N 1' >> trinity.sh
echo '#SBATCH --exclusive' >> trinity.sh
echo '#SBATCH --mem-per-cpu=16000' >> trinity.sh
echo 'module load python3/3.7.7' >> trinity.sh
echo "~/bin/trinityrnaseq-v2.13.2/Trinity --seqType fq --left R1_Pione.fastq --right R2_Pione.fastq --CPU 20 --max_memory 100G --output pione_trinity" >> trinity.sh
sbatch -o trinity.o%j -e trinity.e%j trinity.sh

mv pione_trinity.Trinity.fasta Plampa_Trinity.fasta
echo "seq_stats.pl Plampa_Trinity.fasta > seqstats_Plampa_Trinity.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Plampa_Trinity.fasta
-------------------------
592764 sequences.
668 average length.
17411 maximum length.
183 minimum length.
N50 = 917
396.2 Mb altogether (396185464 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## First cleaning step:
# Assemblies include many small contigs that are unlikely to provide significant matches, so for analyses based on sequence homology we consider only contigs ≥500 bp.

echo "perl ~/bin/noshorts.pl Plampa_Trinity.fasta 500" > noshorts
launcher_creator.py -j noshorts -n noshorts -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch noshorts.slurm

echo "seq_stats.pl noshorts_Plampa_Trinity.fasta > seqstats_noshorts_Plampa_Trinity.txt" > seq_stats2
launcher_creator.py -j seq_stats2 -n seq_stats2 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats2.slurm

noshorts_Plampa_Trinity.fasta
-------------------------
220699 sequences.
1260 average length.
17411 maximum length.
500 minimum length.
N50 = 1536
278.1 Mb altogether (278146740 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## Second cleaning step:
Remove ribosomal/mitochondrial contaminants using SILVA and NCBI databases and available references for the target species

# ribosomal RNA, including microbial sequences
# Go to https://www.arb-silva.de/no_cache/download/archive/current/Exports/ and copy the links for the LSURef_NR99 and SSURef_NR99 compressed fasta files
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
gunzip *.gz
# Concatenate both into a single fasta file for blasting
cat SILVA_138.1_LSURef_NR99_tax_silva.fasta SILVA_138.1_SSURef_NR99_tax_silva.fasta > SILVA_SSU_LSU_combined.fasta

# Go to https://www.arb-silva.de/search/ and search for your species; if it's not available, pick a related species
# Add it to your cart using the checkbox on the left, then select Download in the top right
# Choose FASTA without gaps, and tar.gz
# Cliona varians is available in this GitHub repo as 'arb-silva.de_2021-12-01_id1089989'
# Originally from Redmond et al. (2013) doi: 10.1093/icb/ict078
cp ~/bin/arb-silva.de_2021-12-01_id1089989.tgz .
tar -vxf arb-silva.de_2021-12-01_id1089989.tgz

# mitochondrial RNA
# Closest relative available (Cliona varians) is available in this GitHub repo as 'Cvarians_mitoRNA'
# Originally from Plese et al. (2021) doi: 10.1016/j.ympev.2020.107011
cp ~/bin/Cvarians_mitoRNA.fasta .

# Running a perl script that blasts the transcriptome against inputted contamination sequences to generate a final, clean transcriptome
echo "perl ~/bin/RemoveContamSeq_blast+.pl type=blastn score=45 reads=noshorts_Plampa_Trinity.fasta contam=rRNA,arb-silva.de_2021-12-01_id1089989_tax_silva.fasta contam=Mt,Cvarians_mitoRNA.fasta table=Plampa_contamination.txt passed=Plampa.fasta" > contam
launcher_creator.py -j contam -n contam -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch contam.slurm

Plampa.fasta
-------------------------
220699 sequences in input file
228 sequences look like contaminants
	rRNA	197
	Mt	31
220699 sequences passed all tests
-------------------------

echo "seq_stats.pl noshorts_Plampa_Trinity.fasta > seqstats_Plampa.txt" > seq_stats3
launcher_creator.py -j seq_stats3 -n seq_stats3 -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats3.slurm
