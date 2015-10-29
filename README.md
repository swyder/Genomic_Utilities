### University of Zurich
### URPP Evolution in Action
![URPP logo](Logo_URPP_kl2.png)

Stefan Wyder & Heidi Lischer

stefan.wyder@uzh.ch  
heidi.lischer@ieu.uzh.ch



## Frequently used bioinformatics utilities

For many tasks we regularly do in bioinformatics handy tools already exist. By combining them we can achieve many things.
Of course we could code their functionality ourself but it's usually better to use exisiting tools as
a) these tools have more users and therefore bugs are more likely to be identified and b) they are usually optimized and fast. 
  
Many things you could also achieve using R or python. But the power of these command-line tools lies in combining them with Linux/Unix tools.
  
**A word of caution**  
Never trust software - Do sanity checks often!


![Tools](Tools.png)


**Alignments (SAM, BAM)**  
1. samtools  
2. Picard  
**Genomic Intervals (bed, gtf, gff)**  
3. bedtools  
**Variants (vcf)**  
4. vcftools  
5. bcftools  
6. Picard  
  
  
Chapters 1-3 Stefan  
Chapters 4-5 Heidi


### Installation

Under Ubuntu installation is simple

```
sudo apt-get install samtools bedtools vcftools picard-tools igv
```

Under Mac OS check `homebrew` or your favourite package manager.


Note that package managers often do not install the latest software versions. If you regularly use a tool it is a good idea to suscribe to its mailing list or to check its website from time to time.

```
sudo ln -s -t /usr/bin/ ~/software/BEDTOOLS/bedtools2.25.0/bin/bedtools
```

## SAM/BAM files

The Sequence Alignment/Mapping (SAM) format for mapping data (and its binary equivalent, BAM). The SAM/BAM formats are the standard formats for 
storing sequencing read mapped to a reference. Unlike more general formats like BED and GTF/GFF, SAM/BAM are designed specifically for storing 
the position reads align to and other information about the alignment. Importantly, SAM/BAM do **not** contain the reference sequence.

A SAM/BAM file consists of a header and an alignment section. The detailed SAM/BAM file specifications are available [http://samtools.github.io/hts-specs/]here.


### samtools

[samtools](http://www.htslib.org/) are part of almost every NGS workflow. Samtools is a set of utilities that manipulate alignments in 
the BAM format. It imports from and exports to the SAM (Sequence Alignment/Map) format, does sorting, merging and indexing, and allows to 
retrieve reads in any regions swiftly.

Here we use version 1.2. Check your version by typing `samtools` in the command line. Package managers usually install 0.1.18. You can do the 
exercises with an older version, but older versions require some flags that are now unnecessary, like `-S` with `samtools view` to specify input 
is a SAM file (now filetype is autodetected).

The full options are available [here](http://www.htslib.org/doc/samtools.html).

```
samtools 

Program: samtools (Tools for alignments in the SAM format)
Version: 1.2 (using htslib 1.2.1)

Usage:   samtools <command> [options]

Commands:
  -- indexing
         faidx       index/extract FASTA
         index       index alignment
  -- editing
         calmd       recalculate MD/NM tags and '=' bases
         fixmate     fix mate information
         reheader    replace BAM header
         rmdup       remove PCR duplicates
         targetcut   cut fosmid regions (for fosmid pool only)
  -- file operations
         bamshuf     shuffle and group alignments by name
         cat         concatenate BAMs
         merge       merge sorted alignments
         mpileup     multi-way pileup
         sort        sort alignment file
         split       splits a file by read group
         bam2fq      converts a BAM to a FASTQ
  -- stats
         bedcov      read depth per BED region
         depth       compute the depth
         flagstat    simple stats
         idxstats    BAM index stats
         phase       phase heterozygotes
         stats       generate stats (former bamcheck)
  -- viewing
         flags       explain BAM flags
         tview       text alignment viewer
         view        SAM<->BAM<->CRAM conversion
```

Help for individual task is available by typing `samtools COMMAND`, e.g. `samtools view`


Indexes BAM file for fast access.
This index is needed when region arguments are used to limit samtools view and similar commands to particular regions of interest

***Limitations***   
Samtools paired-end `rmdup` does not work for unpaired reads (e.g. orphan reads or ends mapped to different chromosomes). If this is a concern, 
please use Picard's MarkDuplicates which correctly handles these cases

A faster and multithreaded alternative to `samtools` is [sambamba](http://lomereiter.github.io/sambamba/)



### Picard

[Picard](http://broadinstitute.github.io/picard/) is another set of tools to work with SAM/BAM files developed in the Broad institute by the same people who develop the GATK (Genome Analysis Toolkit) variant calling pipeline.

You have to use Picard before running GATK for variant calling. Notably, Picard provides **BAM file preprocessing** which often improves SNP calling accuracy (See GATK best practices). 

However, Picard tools are also interesting for getting information on BAM files if you use a different variant caller.

Picard commands have the form

```
java -Xmx2g -jar picard.jar COMMAND OPTION1=value1 OPTION2=value2...
```


You can get a list of Picard commands by doing

```
java -jar software/PICARD/picard-tools-1.140/picard.jar
```

and help on individual commands e.g. for Markduplicates by doing

```
java -jar software/PICARD/picard-tools-1.140/picard.jar MarkDuplicates -h
```

The documentation lists 58 Picard commands. It provides similar functionality as samtools to generate a BAM index (BuildBamIndex), view (ViewSam) or sort a SAM/BAM file (SortSam), mark PCR/optical duplicates (MarkDuplicates/MarkDuplicatesWithMateCigar), fix mate information (FixMateInformation), concatenates (GatherBamFiles), merges (MergeSamFiles), converts formats (SamFormatConverter).

In addition, Picard has functionality for filtering (FilterSamReads), comparing (CompareSAMs) and downsampling (DownsampleSam) SAM/BAM files.

One of Picard's strength are the many command to generate quality reports (some also producing QC figures) like e.g. MeanQualityByCycle and QualityScoreDistribution. It also provides functions to do Quality control for amplicons or captured sequences.

It provides also functions for working with fasta files (NormalizeFasta, ExtractSequences) or with vcf files (FilterVcf,MergeVcfs,SortVcf, LiftoverVcf,VcfFormatConverter,SplitVcfs,GenotypeConcordance,MakeSitesOnlyVcf).



### bedtools

Many types of genomic data are linked to a specific genomic region and regions can be represented as a range
of consecutive positions on a chromosome. Annotation data and genomic data like gene models, sequence variants 
(SNPs), promoter regions, transposable elements, pairwise diversity and many others can all be representated 
as genomic ranges. Also alignments of sequencing reads (e.g. from genome sequencing or RNA-Seq) can be representated as genomic ranges.

Once our genomic data is represented as ranges on chromosomes, there are numerous range operations at our disposal to tackle tasks like finding and counting overlaps, cal‐ culating coverage, finding nearest ranges, and extracting nucleotide sequences from specific ranges.
Specific problems like finding which SNPs overlap coding sequences, or counting the number of read alignments that overlap an exon have simple  
  

bedtools help in format conversion, counting, filtering, annotating, comparing, ...


```
bedtools: flexible tools for genome arithmetic and DNA sequence analysis.
usage:    bedtools <subcommand> [options]

The bedtools sub-commands include:

[ Genome arithmetic ]
    intersect     Find overlapping intervals in various ways.
    window        Find overlapping intervals within a window around an interval.
    closest       Find the closest, potentially non-overlapping interval.
    coverage      Compute the coverage over defined intervals.
    map           Apply a function to a column for each overlapping interval.
    genomecov     Compute the coverage over an entire genome.
    merge         Combine overlapping/nearby intervals into a single interval.
    cluster       Cluster (but don't merge) overlapping/nearby intervals.
    complement    Extract intervals _not_ represented by an interval file.
    subtract      Remove intervals based on overlaps b/w two files.
    slop          Adjust the size of intervals.
    flank         Create new intervals from the flanks of existing intervals.
    sort          Order the intervals in a file.
    random        Generate random intervals in a genome.
    shuffle       Randomly redistrubute intervals in a genome.
    sample        Sample random records from file using reservoir sampling.
    annotate      Annotate coverage of features from multiple files.

[ Multi-way file comparisons ]
    multiinter    Identifies common intervals among multiple interval files.
    unionbedg     Combines coverage intervals from multiple BEDGRAPH files.

[ Paired-end manipulation ]
    pairtobed     Find pairs that overlap intervals in various ways.
    pairtopair    Find pairs that overlap other pairs in various ways.

[ Format conversion ]
    bamtobed      Convert BAM alignments to BED (& other) formats.
    bedtobam      Convert intervals to BAM records.
    bamtofastq    Convert BAM records to FASTQ records.
    bedpetobam    Convert BEDPE intervals to BAM records.
    bed12tobed6   Breaks BED12 intervals into discrete BED6 intervals.

[ Fasta manipulation ]
    getfasta      Use intervals to extract sequences from a FASTA file.
    maskfasta     Use intervals to mask sequences from a FASTA file.
    nuc           Profile the nucleotide content of intervals in a FASTA file.

[ BAM focused tools ]
    multicov      Counts coverage from multiple BAMs at specific intervals.
    tag           Tag BAM alignments based on overlaps with interval files.

[ Statistical relationships ]
    jaccard       Calculate the Jaccard statistic b/w two sets of intervals.
    reldist       Calculate the distribution of relative distances b/w two files.

[ Miscellaneous tools ]
    overlap       Computes the amount of overlap from two intervals.
    igv           Create an IGV snapshot batch script.
    links         Create a HTML page of links to UCSC locations.
    makewindows   Make interval "windows" across a genome.
    groupby       Group by common cols. & summarize oth. cols. (~ SQL "groupBy")
    expand        Replicate lines based on lists of values in columns.

[ General help ]
    --help        Print this help menu.
    --version     What version of bedtools are you using?.
    --contact     Feature requests, bugs, mailing lists, etc.
```

```
- ranges-qry.bed is a simple BED file containing six ranges. These are the query ranges used in the GenomicRanges findOverlaps examples (except, since these are in BED format, they are 0-indexed).
- ranges-sbj.bed is the counterpart of ranges-sbj; these are the subject ranges used in theGenomicRangesfindOverlapsexamples. Bothranges-sbj.bedandranges- qry.bed are depicted in Figure 9-11 (though, this visualization uses 1-based coor‐ dinates).
- Mus_musculus.GRCm38.75_chr1.gtf.gz are features on chromosome 1 of Mus_mus culus.GRCm38.75.gtf.gz. The latter file is Ensembl’s annotation file for mm10 (which is the same as GRCh38) and was downloaded from Ensembl’s FTP site (ftp:// ftp.ensembl.org/pub/release-75/gtf/mus_musculus).
- Mus_musculus.GRCm38_genome.txt is a tab-delimited file of all chromosome names and lengths from the mm10/GRCm38 genome version.
```


### Other tools




## Exercises

You find the genome and BAM files on  


#### 1. samtools

1. Make an index of the example BAM file  
2. Check from header the number of chromosomes and their length  
3. View the BAM file and check where in the genome the read _5:1:8:11124:15157 maps to  
1. Extract all reads mapped to positions 99-110 (Can also be used to generate a BAM per chromosome)  
1. Make an index of a reference fasta file   
2. Extract nucleotide positions 99-110 from the reference genome  


http://www.htslib.org/workflow/#mapping_to_variant



#### 2. Sliding-window analysis

Sliding-window analysis is an often used simple approach to look at heterogeneity across the genome. Combined with different summary statistics/metrices
it is used for many applications such as looking for regions under selection (Fst), recombination breakpoints, copy number variants.  

1. Define windows on genome
2. Count number of reads per window
3. Plot the coverage across the genome
4. Calculate some metric/statistic
5. Plot the metric/statistic across the genome

#### 3. Identify low-coverage regions of an alignment

Due to technical bias and random sampling the read coverage varies across a genome. Often we have low-coverage regions
with too few reads to call SNPs with confidence. Her we want to know for a example BAM how large these low-coverage regions are and what proportion of the genome they cover.

1. Count read coverage for each base in the genome
2. Filter out the positions below a cut-off
3. Merge consecutive positions
4. Filter out very short regions
 
### Sources

- Respective documentations

- [Vince Buffalo. Bioinformatics Data Skills. O'reilly 2015](http://shop.oreilly.com/product/0636920030157.do)  
  This practical book teaches the skills that scientists need for turning large sequencing datasets into reproducible and robust biological findings.
  Also covers methods on Sequence and Alignment Data as well as Genomic Ranges. 



### Solutions


#### samtools

1. Make an index of the example BAM file  
```
samtools index MiSeq_Ecoli_DH10B_110721_PF_subsample.bam
```
2. Check from header the number of chromosomes and their length  
```
samtools view -H MiSeq_Ecoli_DH10B_110721_PF_subsample.bam
@HD	VN:1.3	SO:coordinate
@PG	ID:Illumina.SecondaryAnalysis.SortedToBamConverter
@SQ	SN:EcoliDH10B.fa	LN:4686137	M5:28d8562f2f99c047d792346835b20031
@RG	ID:_5_1	PL:ILLUMINA	SM:DH10B_Sample1
```
3. View the BAM file and check where in the genome the read _5:1:8:11124:15157 maps to  
```
samtools view MiSeq_Ecoli_DH10B_110721_PF_subsample.bam | grep "_5:1:8:11124:15157"
_5:1:8:11124:15157	99	EcoliDH10B.fa	11646	254	125M1D25M	=	11879	383	AGGAAGTTTCAGCGCCAGATCGTTGGTTTCGTTACGCGGCATTGCAATGGCGCCGAGGAGTTTATGGTCGTTTGCCTGCGCCGTGCAGCACAGCATCAGGCTAATCGCCAGGCTGGCGGAAATCGCAAAACGGATTTCATACGGAATCTC	??@BD?D?<FHHHIIGG>CCGG<BCF)CFHI@D;GD@@6;AHHHHCH377;9933=88&+9ACC:::++)058834:41599B9@<?:A<?C9?CACC>:98:4+++(525(++29?B5.59(:2(25<C(922&&44:43(88<08>44	RG:Z:_5_1	BC:Z:1	XD:Z:125^T$A15A3T4	SM:i:773	NM:i:4	AS:i:1163
_5:1:8:11124:15157	147	EcoliDH10B.fa	11879	254	150M	=	11646	-383	TATCATTTTTTTAGGAGTACGACTGTGCTTGGGTTTAATTCTATAAAAAAATAAAGTTGTTGCAAATTTTCCGTGTTCAGCTGCCATATCGCGAAATTTCTGCGCAAAAGCACAAAAAATTTTTGCATCTCCCCCTTGATGACGTGGTTT	+(+((&&(:(:+(+((8(+)+++(+(+(((((+(34((::4((()&((+((:>(+((((((++((((&(*&&(((4(+((8(?C@83,23;5.(6@>CDEDA4C4@7=F@8('69/<?8?B9?9?*?1))))2+++2,+2++=1FDB8:+	RG:Z:_5_1	BC:Z:1	XD:Z:70GAC2G3T6GA2G3T6CA1A2A5G2T1G2GC2AG1CAT1CA2AA4TTG1A3	SM:i:390	NM:i:31	AS:i:1163
```
1. Extract all reads mapped to positions 99-110 (Can also be used to generate a BAM per chromosome)  
```
samtools view URPP_Tutorials/EcoliDB10/MiSeq_Ecoli_DH10B_110721_PF_subsample.bam EcoliDH10B.fa:99-110
```
1. Make an index of a reference fasta file  
```
samtools faidx EcoliDH10B.fa
```  
2. Extract nucleotide positions 99-110 from the reference genome  
```
samtools faidx URPP_Tutorials/EcoliDB10/EcoliDH10B.fa EcoliDH10B.fa:99-110
>EcoliDH10B.fa:99-110
ATTAAAATTTTA
```


#### 2. Sliding-window analysis

Sliding-window analysis is an often used simple approach to look at heterogeneity across the genome. Combined with different summary statistics/metrices
it is used for many applications such as looking for regions under selection (Fst), recombination breakpoints, copy number variants.  

1. Define windows on genome  
```
more EcoliDH10B.fa.fai | cut -f1-2 > Chr_Length
bedtools makewindows -g Chr_Length -w 5000  > SlidingWindows.bed
```
2. Count number of reads per window  
```
bedtools coverage -counts -abam MiSeq_Ecoli_DH10B_110721_PF_subsample.bam -b SlidingWindows.bed > Counts_SlidingWin5kb
```
3. Plot the coverage across the genome  
```
# Import into R
counts.5kb <- read.table("Counts_SlidingWin5kb", sep="\t")
colnames(counts.5kb) <- c("Chr", "Start", "End", "Coverage")
```
4. Calculate some metric/statistic  
5. Plot the metric/statistic across the genome  
```
# We have to sort the file according start position
counts.5kb.sorted <- counts.5kb[order(decreasing=F, counts.5kb$Start), ] 
head(counts.5kb.sorted)
              Chr Start   End Coverage
287 EcoliDH10B.fa     0  5000     1417
288 EcoliDH10B.fa  5000 10000     1579
289 EcoliDH10B.fa 10000 15000     1604
36  EcoliDH10B.fa 15000 20000     1129
290 EcoliDH10B.fa 20000 25000     1392
291 EcoliDH10B.fa 25000 30000     1582
plot(counts.5kb.sorted$Coverage)
```


#### 3. Identify low-coverage regions of an alignment

Often some regions of the genome are low coverage only (or even without any aligned reads) consequently we cannot tell whether polymorphisms
exist in these regions. We want to identify such regions and find out whether they overlap with genes.
The bedtools utilities are convenient for working with genomic coordinates. For example bedtools allows one to intersect, merge, count, 
complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF. You will 
find the documentation under http://bedtools.readthedocs.org/en/latest/index.html


Use `MiSeq_Ecoli_DH10B_110721_PF_subsample.bam` from the first exercise.  

- Try to find out how many nucleotides in the *E. coli* genome do not reach a minimal read coverage of 5, 10 or 20 reads (Hint: Use samtools 
depth, direct the output into a file and then use awk or R to process the file)

```
samtools depth -q 20 MiSeq_Ecoli_DH10B_110721_PF_subsample.bam | awk '$3<5 {print}' | wc -l
8206

samtools depth -q 20 MiSeq_Ecoli_DH10B_110721_PF_subsample.bam  | awk '$3<10 {print}' | wc -l
15328

samtools depth -q 20 MiSeq_Ecoli_DH10B_110721_PF_subsample.bam  | awk '$3<20 {print}' | wc -l
44441
```

or the same can be done using `bedtools genomecov` or by parsing the output of `samtools stats` or using R:

```
samtools depth -q 20 MiSeq_Ecoli_DH10B_110721_PF_subsample.bam  | cut -f 2-3 > genome coverage

Rscript Make_Cov_Histogram.R 
[1] "minimum, lower-hinge, median, upper-hinge, maximum"
[1]  0 36 41 47 80
[1] "Nr of bases with Coverage = 0: 301"
[1] "Nr of bases with Coverage < 5: 8206"
[1] "Nr of bases with Coverage < 10: 15328"
``` 

It saves the histogram as Rplots.pdf. The Make_Cov_Histogram.R looks like this:

```
# script to Make Histogram from Genome Coverage per Base
f <- read.table("genome.coverage")$V2
print("minimum, lower-hinge, median, upper-hinge, maximum")
print(fivenum(f))
print(paste("Nr of bases with Coverage = 0:", length(which(f==0))))
print(paste("Nr of bases with Coverage < 5:", length(which(f<5))))
print(paste("Nr of bases with Coverage < 10:", length(which(f<10))))
hist(f)
```

- We now counted individual nucleotides but actually we would like to know whether there are some regions with low coverage (Hint: use mergeBed from bedtools)

```
samtools depth MiSeq_Ecoli_DH10B_110721_PF_subsample.bam | awk 'BEGIN{OFS="\t"} $3<10 {print $1,$2-1,$2}' | bedtools merge -i - | head
EcoliDH10B.fa   0       17
EcoliDH10B.fa   15618   15749
EcoliDH10B.fa   16379   16496
EcoliDH10B.fa   19963   19964
EcoliDH10B.fa   19981   20114
```
- Are there any regions with unexpected high coverage? As mapping problems are likely (due to repetitive regions) they are often excluded from the analysis.

```
samtools depth -q 20 MiSeq_Ecoli_DH10B_110721_PF_subsample.bam  | awk 'BEGIN {OFS="\t"} $3>57 {print}' > positions_CoverageOverMean2SD
```

- Do some low coverage nucleotides overlap with annotated genes? (Hint: use the command intersect from bedtools)

```
# Extract positions with coverage < 5
samtools depth -q 20 MiSeq_Ecoli_DH10B_110721_PF_subsample.bam  | awk 'BEGIN {OFS="\t"} $3<5 {print $1,$2,$2}' > positions_CoverageBelow5
# Extract exon positions from gff file
awk 'BEGIN {OFS="\t"} $3=="CDS" || $3=="exon" {print}' EcoliDH10B.gff > EcoliDH10B_exons.gff
# Find overlaps
EcoliDH10B.fa	RefSeq	CDS	15445	16557	.	+	0	ID=cds15;Name=YP_001728999.1;Parent=gene15;Dbxref=ASAP:AEC-000
0088,Genbank:YP_001728999.1,GeneID:6061805;gbkey=CDS;product=IS186%2FIS421 transposase;protein_id=YP_001728999.1;transl_table=11
EcoliDH10B.fa	RefSeq	CDS	19811	20314	.	-	0	ID=cds20;Name=YP_001729004.1;Parent=gene21;Dbxref=ASAP:AEC-000
0092,Genbank:YP_001729004.1,GeneID:6060829;gbkey=CDS;product=IS1 transposase InsAB%27;protein_id=YP_001729004.1;transl_table=11
EcoliDH10B.fa	RefSeq	CDS	20233	20508	.	-	0	ID=cds21;Name=YP_001729005.1;Parent=gene22;Dbxref=ASAP:AEC-000
0093,Genbank:YP_001729005.1,GeneID:6061774;gbkey=CDS;product=IS1 repressor protein InsA;protein_id=YP_001729005.1;transl_table=11
EcoliDH10B.fa	RefSeq	exon	197875	199416	.	+	.	ID=id10;Parent=rna2;Dbxref=ASAP:AEC-0004218,GeneID:6062217;gbk
ey=rRNA;product=16S ribosomal RNA
...
# 66 genes have a coverage < 5 at least at 1 position
intersectBed -u -a EcoliDH10B_exons.gff -b positions_CoverageBelow5 | awk 'BEGIN {FS="\t"} {split($9, elements, ";"); print elements[2]}' | sed 's/Parent=//;s/Name=//' | sort -u |
```


