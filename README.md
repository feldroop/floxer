# floxer: FM-index longread PEX-based aligner

An exact longread aligner applying the following techniques to be not orders of magnitude slower than approximate tools like [minimap](https://github.com/lh3/minimap2).

* 2-3 error approximate FM-Index search using search schemes
* PEX hierarchical verification
* parallel and vectorized pairwise alignment implementation

## Installation on Linux

Requires a recent C++20-capable compiler and CMake.

```
git clone --recurse-submodules https://github.com/feldroop/floxer
mkdir floxer/build && cd floxer/build
cmake .. -DCMAKE_BUILD_TYPE:STRING=Release
make
./bin/floxer --help
```

## Usage

Basic usage:
```
./floxer --reference hg38.fasta --query reads.fastq --errors 7 --output mapped_reads.bam
```

For more examples about command line options, run:
```
./floxer --help
```
