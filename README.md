# Effect of food source availability in the salivary gland transcriptome of the unique burying beetle *Nicrophorus pustulatus* (Coleoptera: Silphidae)

The objective of this repository is to provide supplemental information about the software packages, and their parameters, that were used during the RNAseq analysis of the *Nicrophorus pustulatus* salivary gland transcriptomes for the paper: 

#### Effect of food source availability in the salivary gland transcriptome of the unique burying beetle Nicrophorus pustulatus (Coleoptera: Silphidae) (https://doi.org/10.1371/journal.pone.0255660)

Read quality control, transcriptome assembly, transcriptome annotation and read pseudo-mapping and abundance determination were done using the resources of the Oklahoma State University High Performance Computing Center. Differential expression analysis and gene set enrichment analysis were done in a local machine running R (R Core Team, 2020).
All the scripts can be found in the `scripts/` directory of this repository. A brief explanation of the steps of the analysis is detailed below

## 1. Raw reads quality control

Raw reads (Accession: ***pending***) quality was assessed using FastQC (Andrews, 2010) with the following parameters:
```
fastqc -o fastqc_out -t 12 *.fastq.gz
```

## 2. Quality trimming and transcriptome assembly

Quality trimming of the reads was done using Trimmomatic (Bolger, et al., 2014) as part of the Trinity (Grabherr, et al., 2011) pipeline. Transcriptome assembly was done pooling together the reads from all males and females fed and starved:
```
Trinity --seqType fq --max_memory 750G --left FPF1_1.fq.gz,FPF2_1.fq.gz,FPF3_1.fq.gz,FPS1_1.fq.gz,FPS2_1.fq.gz,FPS3_1.fq.gz,MPF1_1.fq.gz,MPF2_1.fq.gz,MPF3_1.fq.gz,MPS1_1.fq.gz,MPS2_1.fq.gz,MPS3_1.fq.gz  --right FPF1_2.fq.gz,FPF2_2.fq.gz,FPF3_2.fq.gz,FPS1_2.fq.gz,FPS2_2.fq.gz,FPS3_2.fq.gz,MPF1_2.fq.gz,MPF2_2.fq.gz,MPF3_2.fq.gz,MPS1_2.fq.gz,MPS2_2.fq.gz,MPS3_2.fq.gz --CPU 32 --trimmomatic --output Nicropus_transcriptome
```

## 3. Reduction to unigene data set

The assembled transcriptome was deduplicated to obtain an unigete data set by removing sequences with over 95% similarity using CD-HIT (Fu, et al., 2012; Li and Godzick, 2006:

```
cd-hit-est -i Nicropus_transcriptome.fasta -o Nicropus_transcriptome95.fasta -c 0.95 -n 10 -d 0 -T 32
```

## 4. Transcriptome and unigene set completeness assessment

Completeness of the transcriptome and the unigene set was performed using BUSCO (v4) (Seppey, et al., 2019) and the Arthropoda ortholog data set (arthropoda_odb10):

```
busco -m transcriptome -i Nicropus_transcriptome.fasta -o Nicropus_transcriptome_busco -l arthropoda_odb10
busco -m transcriptome -i Nicropus_transcriptome95.fasta -o Nicropus_transcriptome95_busco -l arthropoda_odb10
```

## 5. Functional annotation

Annotation of the unigene data sets was performed following the Trinotate pipeline (Bryant, et al., 2017), which includes searches in multiple biologically relevant databases.

### 5.1 Obtaining peptide sequences 

Peptide sequences of the unigene data set were obtained using Transdecoder (Bryant, et al., 2017):

```
TransDecoder.LongOrfs -t Nicropus_transcriptome95.fasta
TransDecoder.Predict -t Nicropus_transcriptome95.fasta
```

### 5.2 Homology search against the Uniprot-KB/Swissprot database

```
blastx -query Nicropus_transcriptome95.fasta -db uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -out blastx.out

# Nicropus_transcriptome95.pep are the peptide sequences obtained from transdecoder
blastp -query Nicropus_transcriptome95.pep -db uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -out blastx.out
```

### 5.3 Homology search against the NCBI Ref-Seq database

```
blastx -db refseq_protein -query Nicropus_transcriptome95.fasta -num_threads 32 -max_target_seqs 1 -out ref_seq_blast.out -outfmt '6 qseqid sseqid pident evalue bitscore stitle staxids'
```

### 5.4 Protein domain identification

Identification of protein domains was performed using HMMER (Eddy, 2001)

```
hmmscan --cpu 32 --domtblout HmmerPFAM.out Pfam-A.hmm Nicropus_transcriptome95.fasta > pfam.log
```

### 5.5 Identification of rRNA transcripts

Ribosomal RNA transcripts were identified using Rnammer (Lagesen, et al., 2007) with the script provided with Trinotate:

```
$TRINOTATE_HOME/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Nicropus_transcriptome95.fasta --path_to_rnammer /opt/rnammer/1.2/prebuilt/rnammer
```

## 6. Transcriptome pseudo-alignment and abundace determination

Transcript abundance was determined using the pseudo-aligner Kallisto (Bray, et al., 2016).

### 6.1 Build index from unigene data set
```
kallisto index -i Nicropus_transcriptome95.fasta
```

### 6.2 Transcript abundance

The following command was run per each of the samples:

```
kallisto quant -i Nicropus_kallisto_index -t 32 -b 500 -o kallisto_Nicropus-counts/MPF1/ -t 32 MPF1_1.fq.gz MPF1_2.fq.gz
```

## 7. Differential expression analysis

Differential expression analysis was performed in R using the script `7_nicropus_deg.R` in the `scripts/` directory


### 

## References
- Andrews S. FastQC: a quality control tool for high throughput sequence data 2010. Available from: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.
- Bolger, A.M.; Lohse, M.; Usadel, B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 2014, 30, 2114–2120. doi:10.1093/bioinformatics/btu170.
- Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. Nature Biotechnology. 2016;34(5):525-7. doi: 10.1038/nbt.3519.
- Bryant DM, Johnson K, DiTommaso T, Tickle T, Couger MB, Payzin-Dogru D, Lee TJ, Leigh ND, Kuo TH, Davis FG, Bateman J, Bryant S, Guzikowski AR, Tsai SL, Coyne S, Ye WW, Freeman RM Jr, Peshkin L, Tabin CJ, Regev A, Haas BJ, Whited JL. A Tissue-Mapped Axolotl De Novo Transcriptome Enables Identification of Limb Regeneration Factors. Cell Rep. 2017 Jan 17;18(3):762-776. doi: 10.1016/j.celrep.2016.12.063. PubMed PMID: 28099853; PubMed Central PMCID: PMC5419050.
- Eddy SR. Accelerated profile HMM searches. PLoS computational biology. 2011;7(10):e1002195.
- Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics. 2012;28(23):3150-2. Epub 2012/10/13. doi:10.1093/bioinformatics/bts565. PubMed PMID: 23060610; PubMed Central PMCID: PMCPMC3516142.
- Lagesen K, Hallin P, Rødland E, Stærfeldt H, Rognes T, Ussery D. RNammer: consistent annotation of rRNA genes in genomic sequences. Nucleic Acids Res. 2007;35(9):3100-8.
- Li W, Godzik A. Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics. 2006;22(13):1658-9. doi: 10.1093/bioinformatics/btl158.
- R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria, 2020.
- Seppey M., Manni M., Zdobnov E.M. (2019) BUSCO: Assessing Genome Assembly and Annotation Completeness. In: Kollmar M. (eds) Gene Prediction. Methods in Molecular Biology, vol 1962. Humana, New York, NY. 2019 doi.org/10.1007/978-1-4939-9173-0_14
