##make directory
mkdir ref_genome

##wget your reference genome and output it as fasta.gz
wget -O ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz

##use the gunzip command to change from fasta.gz to fasta file
gunzip ecoli_rel606.fasta.gz

##to know the name of the genome
head ecoli_rel606.fasta ##CP000819.1 Escherichia coli B str. REL606

##to download my reads and output it as tar.gz file
wget -O sub.tar.gz https://ndownloader.figshare.com/files/14418248

##to untar your downloaded reads
tar xvf sub.tar.gz

##to rename my untar file from sub to trimmed.fastq
mv sub trimmed.fastq

##to create directories for our bam, sam, bcf, and vcf results using a single line code
mkdir -p results/sam results/bam results/bcf results/vcf

##index our reference genome
bwa index ecoli_rel606.fasta

##to align our reads both forward and reverse strands for the three datasets to reference genome and send the results as a sam file in the results directory
bwa mem ecoli_rel606.fasta trimmed.fastq/SRR2584863_1.trim.sub.fastq trimmed.fastq/SRR2584863_2.trim.sub.fastq >../results/sam/SRR2584863_aligned.sam
bwa mem ecoli_rel606.fasta trimmed.fastq/SRR2584866_1.trim.sub.fastq trimmed.fastq/SRR2584866_2.trim.sub.fastq >../results/sam/SRR2584866_aligned.sam
bwa mem ecoli_rel606.fasta trimmed.fastq/SRR2589044_1.trim.sub.fastq trimmed.fastq/SRR2589044_2.trim.sub.fastq >../results/sam/SRR2589044_aligned.sam

##convert the sam files in the results directory to bam files
samtools view -S -b sam/SRR2584863_aligned.sam >bam/SRR2584863_aligned.bam
samtools view -S -b sam/SRR2584866_aligned.sam >bam/SRR2584866_aligned.bam
samtools view -S -b sam/SRR2589044_aligned.sam >bam/SRR2589044_aligned.bam

##to sort the files
samtools sort -o bam/SRR2584863_aligned.sorted.bam bam/SRR2584863_aligned.bam
samtools sort -o bam/SRR2584866_aligned.sorted.bam bam/SRR2584866_aligned.bam
samtools sort -o bam/SRR2589044_aligned.sorted.bam bam/SRR2589044_aligned.bam

##install bcf and vcf tools
conda install -c bioconda bcftools
conda install -c bioconda bcftools

##for variant calling. step1: calculate the read coverage of position in the genome using the reference fasta and the sorted file in bam
bcftools mpileup -O b -o bcf/SRR2584863_raw.bcf -f ../ref_genome/ecoli_rel606.fasta bam/SRR2584863_aligned.sorted.bam
bcftools mpileup -O b -o bcf/SRR2584866_raw.bcf -f ../ref_genome/ecoli_rel606.fasta bam/SRR2584866_aligned.sorted.bam
bcftools mpileup -O b -o bcf/SRR2589044_raw.bcf -f ../ref_genome/ecoli_rel606.fasta bam/SRR2589044_aligned.sorted.bam

##step 2: detect the single nucleotide variants (SNVs)
bcftools call --ploidy 1 -m -v -o vcf/SRR2584863_variants.vcf bcf/SRR2584863_raw.bcf
bcftools call --ploidy 1 -m -v -o vcf/SRR2584866_variants.vcf bcf/SRR2584866_raw.bcf
bcftools call --ploidy 1 -m -v -o vcf/SRR2589044_variants.vcf bcf/SRR2589044_raw.bcf

##step 3: filter and report SNV in variant calling formats (VCF)
vcfutils.pl varFilter vcf/SRR2584866_variants.vcf > vcf/SRR2584866_final_variants.vcf
vcfutils.pl varFilter vcf/SRR2584863_variants.vcf > vcf/SRR2584863_final_variants.vcf
vcfutils.pl varFilter vcf/SRR2589044_variants.vcf > vcf/SRR2589044_final_variants.vcf

##explore the vcf variants 
less -S vcf/SRR2584863_final_variants.vcf\

##visualizing the alignment. step 1: index the bam files
samtools index bam/SRR2584863_aligned.sorted.bam
samtools index bam/SRR2584866_aligned.sorted.bam
samtools index bam/SRR2589044_aligned.sorted.bam

##step 2:visualizing the genome
samtools tview bam/SRR2584863_aligned.sorted.bam ../ref_genome/ecoli_rel606.fasta
samtools tview bam/SRR2584866_aligned.sorted.bam ../ref_genome/ecoli_rel606.fasta
samtools tview bam/SRR2589044_aligned.sorted.bam ../ref_genome/ecoli_rel606.fasta


##reproducing the tutorial
mkdir genome
mkdir gen_results

##downloading the reference genome
wget -O genome/organisms.fasta.gz https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

##unzip the reference 
gunzip organisms.fasta.gz

## downloading the reads
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

##create directory genome_fastq and send all the downloaded reads to the directory
mkdir genome_fastq
mv SLGFSK-N_231335_r1_chr5_12_17.fastq.gz SLGFSK-N_231335_r2_chr5_12_17.fastq.gz SLGFSK-T_231336_r1_chr5_12_17.fastq.gz SLGFSK-T_231336_r2_chr5_12_17.fastq.gz genome_fastq/

##unzipped the reads
gunzip *.gz

