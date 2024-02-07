# floxer: FM-index longread PEX-based aligner 

An exact longread aligner applying the following techniques to be not orders of magnitude slower than approximate tools like [minimap](https://github.com/lh3/minimap2).

* 2-3 error approximate FM-Index search using search schemes
* PEX hierarchical verification
* parallel and vectorized pairwise alignment implementation

## Usage

Basic usage:
```
./floxer --reference hg38.fasta --query reads.fastq --errors 7 --output mapped_reads.bam
```

For more examples about command line options, run:
```
./floxer --help
```
