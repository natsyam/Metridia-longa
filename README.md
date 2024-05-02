# Metrirdia Longa
Ref_protein - as the ref protein we decided to take the _Isopenicillin N synthase_ (IPNS_STRCL; interpro - P10621; conserv-e site - IPR002057) from _Streptomyces clavuligerus_ gene _pcbC_. Because it's the one, which was detected by our collegues via comparing metagenomics data from metridia londa and clean sea water. 
Also chosing of IPNS_STRCL is based on:
> Были выдвинуты гипотезы (Oba et al., 2009) о биосинтезе целентеразина:
> механизм нерибосомального синтеза циклизация и дальнейшая модификация аминокислотных остатков FYY,  как части более длинного пептида
> 
> Поиск пептидов с мотивом FYY обнаружил наличие гомологов оксигеназ и изопенициллин-N-синтаз у биолюм. гребневиков и не обнаружил у небиолюм. (Francis et al., 2015)

**Steps:**
1. Find needed proteins 
   - find the ORFs in the assembled transcriptome
   - to pick out the ORFs that end with FYY
   - blast the resulting sequences

2. Differential expression analysis
   - mapping reads to assembled transcriptome
   - calculate Counts table

3. Connect the previous steps
   - see if any of the proteins with a significant difference in expression have the FYY motif

### Step 1. Find needed proteins (ORF finding)

##### Installing
```
conda create -n transdecoder
conda activate transdecoder
conda install -c bioconda transdecoder
/path to/TransDecoder.LongOrfs --version
```
TransDecoder.LongOrfs 5.5.0

```
conda create -n hmm
conda activate hmm
conda install -c biocore hmmer
hmmsearch -h
```
hmmsearch :: search profile(s) against a sequence database
HMMER 3.1b2 (February 2015); http://hmmer.org
Usage: hmmsearch [options] <hmmfile> <seqdb>

 **1. Running for step extract the long open reading frames TransDecoder**
```
TransDecoder.LongOrfs -t soft_filtered_transcripts.fasta
_______________________________________________
we got /.../soft_filtered_transcripts.fasta.transdecoder_dir with bunch of files, including **longest_orfs.pep**
```

**2. Search for peptides with an FYY motif**

First we found the sections that contain the FYY motif and then picked the ones that end with FYY.
```
awk 'BEGIN{RS=">";FS="\n"} {seq=$2; gsub("\n","",$2); if(index(seq, "FYY") || index(seq, "YYF")) {print ">"$1"\n"$2"\n"}}' longest_orfs.pep > filtered_hypothetical_peptides.fasta

grep -B 1 "FYY$" filtered_hypothetical_peptides.fasta > FYYend_hypothetical_peptides.fasta
```
Also created a file that contains sections with the FYY motif before the stop codon.
```
grep -B 1 --no-group-separator  "FYY\*" longest_orfs.pep >fyy_peptides.fasta
```
Then merged two files:
```
cat fyy_peptides.fasta FYYend_hypothetical_peptides.fasta > all_fyy.fasta
```

**3. BLAST the results file**

```
update_blastdb.pl --decompress  "refseq_protein" --force --verbose

blastp -db /local/workdir/dmm2017/blast_protein/refseq_protein -query all_fyy.fasta -out blast_results/blastp_results.txt -evalue 0.001 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -num_threads 4
```

### Step 2. Differential expression analysis

**1. Mapping reads to assembled transcriptome**
It was decided to run mapping in two ways: separately for each sample and also run for groups of luminous (2,4,6) and non-luminous parts (3,5,7).

Examples of running mapping using hisat2:
```
hisat2-build -p 16 soft_filtered_transcripts.fasta tr_index

hisat2 -x tr_index -1 S1_R1.fastq.gz -2 S1_R2.fastq.gz | samtools view -Su - | samtools sort --threads 30 -m 12G -o metridia.1.hisat.bam -
```
```
hisat2 -x tr_index -1 S2_R1.fastq.gz,S4_R1.fastq.gz,S6_R1.fastq.gz -2 S2_R2.fastq.gz,S4_R2.fastq.gz,S6_R2.fastq.gz | samtools view -Su - | samtools sort --threads 30 -m 12G -o metridia.246.hisat.bam -

hisat2 -x tr_index -1 S3_R1.fastq.gz,S5_R1.fastq.gz,S7_R1.fastq.gz -2 S3_R2.fastq.gz,S5_R2.fastq.gz,S7_R2.fastq.gz | samtools view -Su - | samtools sort --threads 30 -m 12G -o metridia.357.hisat.bam -
```

**2. Calculate Counts table**

Calculating transcript abundance, using FeatureCounts and then performing differential expression using DESeq2.
