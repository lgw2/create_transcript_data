### Requirements:
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
* transcripts from the [human gene annotation](http://www.ensembl.org/Homo_sapiens/Info/Annotation)
with simulated abundances.
* transcripts and abundances built from real reads from experiments stored in the [Sequence
	Read Archive](https://www.ncbi.nlm.nih.gov/sra) using existing
	state-of-the-art transcript assembly pipelines, e.g., [Hisat2](http://daehwankimlab.github.io/hisat2/) for read
	alignment and [StringTie](https://ccb.jhu.edu/software/stringtie/) for transcript assembly.


The input GTF file for humans can be downloaded from http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/
by clicking the `Homo_sapiens.GRCh38.104.gtf.gz`. A small version of this file
is included in the repository as `small_human.gtf`.


For example, running

```
input_and_truth_from_gtf.py small_human.gtf
```

will create a directory called `small_human` with two files: `1.gfa` and
`1.truth`. The first contains the splice graphs for all transcripts for each
gene in `small_human.gtf` in GFA format, which is a list of edges, each with a
weight. The nodes are are genomic coordinate ranges, e.g.,
`(1479049,1479108)`. The second file gives true sequence of nodes for each
transcript.

`input_and_truth_from_gtf.py` can  also simulate
weights for the transcripts using the [Numpy lognormal
distribution](https://numpy.org/doc/stable/reference/random/generated/numpy.random.lognormal.html) with a mean of -4 and variance of 4, as in the
default setting for the
[RNA-Seq reads simulator](http://alumni.cs.ucr.edu/~liw/rnaseqreadsimulator.html)
Then, the output files will
have weighted edges and ground truth paths based on these simulated weights.
This can be done using the `--simulate_cov` flag. For example:
```
input_and_truth_from_gtf.py small_human.gtf --simulate_cov
```
