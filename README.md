The script `input_and_truth_from_gtf.py` was originally written by Ariel with
updates by me (Lucy).

It takes an input gene transfer format (GTF) file and writes the associated
splice graphs for each gene on the positive strand (note: I don't really
understand the positive strand stuff, but basically we only look at lines from
the input file with a `+`) in GFA format.

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

Not yet implemented: `input_and_truth_from_gtf.py` will take in a parameter
saying how to simulate weights for the transcripts. Then, the output files will
have weighted edges and ground truth paths based on these simulated weights.
