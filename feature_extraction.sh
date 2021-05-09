#/bin/sh
# Start with raw nanopore sequencing data (.fast5).

# Step 1: Base call with Albacore
# Input directory containing .fast5 files, the script will read files recursively. The output is fastq files.
# Note: single fast5 file represents one single read. The output fastq file may contain multiple base-called reads
read_fast5_basecaller.py --input /home/yuxin/home/yuxin/Nanopore/DATA/fast5/pass \
--recursive --worker_threads 5 \
--flowcell FLO-MIN107 --kit SQK-RNA001 \
--save_path /home/yuxin/home/yuxin/Nanopore/DATA/fastq

# Step 2: Alignment with Minimap2
# you may need to specify the path: e.g., /home/yuxin/minimap2-2.17_x64-linux/minimap2
minimap2 -d ref.mmi /data/kunqidir/hg19/hg19.fa # indexing
minimap2 -ax map-ont ref.mmi /home/yuxin/home/yuxin/Nanopore/DATA/fastq/workspace/pass.tar.gz > \
/home/yuxin/home/yuxin/Nanopore/DATA/HEK_WT_Alignment.sam
# Sam to bam, sort bam, index bam
samtools view -b -S /home/yuxin/home/yuxin/Nanopore/DATA/HEK_WT_Alignment.sam > \/home/yuxin/home/yuxin/Nanopore/DATA/HEK_WT_Alignment.bam
samtools faidx /home/share/yuxin/DATA/hg19.fa
samtools sort /home/share/yuxin/Nanopore/DATA/HEK_WT_Alignment.bam -o /home/share/yuxin/Nanopore/DATA/HEK_WT_Alignment_sorted.bam
samtools index /home/share/yuxin/Nanopore/DATA/HEK_WT_Alignment_sorted.bam /home/share/yuxin/Nanopore/DATA/HEK_WT_Alignment.bam.bai

# Step 3A: Tombo feature extraction
# Resquiggle
tombo resquiggle /home/yuxin/home/yuxin/Nanopore/DATA/fast5/pass /data/kunqidir/hg19/hg19.fa \
 --processes 4 --num-most-common-errors 5
# Modification base detection
tombo detect_modifications de_novo --fast5-basedirs /home/yuxin/home/yuxin/Nanopore/DATA/fast5/pass \
--statistics-file-basename /home/yuxin/home/yuxin/Nanopore/DATA/fast5/HEK_WT_de_novo_detection --processes 4
# Output the result
tombo text_output browser_files --fast5-basedirs /home/yuxin/home/yuxin/Nanopore/DATA/fast5/pass \
--statistics-filename /home/yuxin/home/yuxin/Nanopore/DATA/fast5/HEK_WT_de_novo_detection.tombo.stats \
--browser-file-basename /home/yuxin/home/yuxin/Nanopore/DATA/fast5/HEK_WT_de_novo_detect \
--file-types coverage dampened_fraction
# Results are stored in wiggle and bedgraph files, which can be processed as text file

# Step 3B: Base call error extraction with EpiNano
# First, you need to create an sequence dictionary for the reference genome with Picard
java -jar /home/share/yuxin/miniconda3/envs/bioinfo/share/picard-2.23.8-0/picard.jar \
CreateSequenceDictionary -R /home/share/yuxin/DATA/hg19.fa -O /home/share/yuxin/DATA/hg19.dict
# There may be something wrong with the EpiNano software: the result file does not include the seqname (chromosome).
# You may separate the bam file by chromosome, and call the function for splitted file in parallel
# Take chromosome 1 as an example
samtools view -b -h /home/share/yuxin/Nanopore/DATA/HEK_WT_Alignment_sorted.bam chr1 > \
/home/share/yuxin/Nanopore/DATA/HEK_WT_Alignment_sorted_chr1.bam
python3 /home/share/yuxin/EpiNano/Epinano_Variants.py -n 2 \
-R /home/share/yuxin/DATA/hg19.fa -b /home/share/yuxin/Nanopore/DATA/HEK_WT_Alignment_sorted_chr1.bam \
-s /home/share/yuxin/jvarkit/dist/sam2tsv.jar --type g

# Warning: Installation and execuation of those software can very time consuming and troublesome.
