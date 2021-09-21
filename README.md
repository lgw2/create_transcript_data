### Requirements
* Python 3+
* numpy
* networkx

### About

The script `input_and_truth_from_gtf.py` was originally written by Ariel with
updates by me (Lucy).
It takes an input gene transfer format (GTF) file and writes the associated
splice graphs and paths for each gene on the positive strand (note: I don't really
understand the positive strand stuff, but basically we only look at lines from
the input file with a `+`).

We are concerned with two types of inputs:
1. transcripts from the [human gene annotation](http://www.ensembl.org/Homo_sapiens/Info/Annotation)
with simulated abundances.
2. transcripts and abundances built from real reads from experiments stored in the [Sequence
	Read Archive](https://www.ncbi.nlm.nih.gov/sra) using existing
	state-of-the-art transcript assembly pipelines, e.g., [Hisat2](http://daehwankimlab.github.io/hisat2/) for read
	alignment and [StringTie](https://ccb.jhu.edu/software/stringtie/) for transcript assembly.

For type 1 data,
an input GTF file for humans can be downloaded from http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/
by clicking the `Homo_sapiens.GRCh38.104.gtf.gz`. A smaller test version of this file
is included in the repository as `small_human.gtf`.

For type 2 data, there are a number of steps needed to create a GTF from the
reads. In the future maybe we will automate this process. For now, here is a
summary of the pipeline for reads in SRA accession number SRR307903 along with
commands to run. (Note that the programs would be need to installed before
running the commands.)

1. Download SRA archive file from https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR307903
2. Use `fastq-dump` from the SRA Toolkit to extract two `fastq` files files from the SRA archive file (two files because paired-end reads)
```
fastq-dump -I --split-files SRR307903
```
3. Download pre-built GRCh38 index from https://daehwankimlab.github.io/hisat2/download/#h-sapiens. I used the genome_snp_tran, but maybe I should have used a different one.
4. Use `hisat2` with the two `fastq` files and the pre-built index to produce a sam file.
```
hisat2 -p 6 -x ../grch38_snp_tran/genome_snp_tran -1 SRR307903_1.fastq -2 SRR307903_2.fastq -S aligned_reads.sam
```
5. Use `samtools` to convert the `sam` file to a `bam` file and then to a sorted `bam` file.
```
samtools view -bS aligned_reads.sam > aligned_reads.bam
samtools sort eg1.bam -o aligned_reads_sorted.bam
```
6. Run `stringtie` on the aligned, sorted reads in the sorted `bam` file.
```
stringtie -o assembly.gtf aligned_reads_sorted.bam
```

### How to run

Running

```
input_and_truth_from_gtf.py small_human.gtf
```

will create a directory called `small_human` with two files: `1.sg` and
`1.truth`. The first contains the splice graphs for all transcripts for each
gene in `small_human.gtf` in what we are calling sg format, which is a list of edges, each with a
weight. The nodes are are genomic coordinate ranges, e.g.,
`(1479049,1479108)`. We add an edge from the source `(0,0)` to the first node
in each transcript an from the last node to the sink `(-1,-1)`.
The second file gives true sequence of nodes for each
transcript. See `small_human/1.sg` and `small_human/1.truth` for examples.
For a larger input file, there will be GFA and truth files for each
chromosome. This process will work for both annotation files (type 1 above, for
example the one from
ensemble) and assemblies (type 2 above, for example the one produce by
stringtie).
For example, the following works for a GTF file assembled with
stringtie:

```
python input_and_truth_from_gtf.py small_SRR307903_assembly.gtf
```

Type 2 inputs already have coverage values. We may want to simulate
coverages for type 1 as described in the next paragraph.

`input_and_truth_from_gtf.py` can  also simulate
weights for the transcripts using the [Numpy lognormal
distribution](https://numpy.org/doc/stable/reference/random/generated/numpy.random.lognormal.html) with a mean of -4 and variance of 4, as in the
default setting for the
[RNA-Seq reads
simulator](http://alumni.cs.ucr.edu/~liw/rnaseqreadsimulator.html).
We also multiply these by 1000, because they tend to be close to zero, but the
catfish inputs have weights ranging from 1 to in the tens of thousands (from
what I remember).
Then, the output files will
have weighted edges and ground truth paths based on these simulated weights.
This can be done using the `--simulate_cov` flag. For example:
```
python input_and_truth_from_gtf.py small_human.gtf --simulate_cov
```

### A note on coverage values

Coverage values are rounded to the nearest integer. If that integer is zero, the
path is not included in the ground truth. No zero-weight edges are written to
the `sg` file, either.

### Pre-computed outputs

If you just want outputs from the two data sets described above, they are at
https://drive.google.com/drive/folders/11bRstTOTOTsYbcRgqiZ7BLAf85SlFXI2?usp=sharing.

The commands to create these files are
```
python input_and_truth_from_gtf.py human.gtf --simulate_cov
python input_and_truth_from_gtf.py SRR307903_assembly.gtf
```

