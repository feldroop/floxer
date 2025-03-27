# floxer: FM-Index longread PEX-based aligner

An *exact*\* longread aligner applying the following techniques to not be orders of magnitude slower than more approximate tools like [minimap](https://github.com/lh3/minimap2).

* 2-3 error approximate FM-Index search using optimal search schemes
* heuristic anchor selection\*
* PEX hierarchical verification guided by a novel PEX tree generation strategy
* parallel and vectorized pairwise alignment implementation

*Exactness*\* here means adhering to a specific formal definition. It is of course impossible to exactly solve the biological read mapping problem. This tool is guaranteed to find a similar (representative) alignment for every linear alignment that matches the query with at most a given error ratio in edit distance\*. This means that large indel and structural variant resolution are currently out of scope of this project.

\*The exactness property is not held in highly repetetive regions where seeds produce many anchors/hits/matches. Here the heuristic anchor selection is used to identify possibly non-repetitive anchors.

This is an experimental research prototype (for my master's thesis) and currently not competitive to state-of-the-art tools like minimap2, in most regards. In addition to the limitations above, it is much slower and therefore not well suited for most applications with large amounts of data.

## Installation on Linux

Requires a C++20-capable compiler and CMake.
```
git clone https://github.com/feldroop/floxer
mkdir floxer/build && cd floxer/build
cmake .. -DCMAKE_BUILD_TYPE:STRING=Release
make
```

Execute the following command inside the build directory to run the tests:
```
make check
```

## Usage

Basic usage:
```
./floxer --reference hg38.fasta --queries reads.fastq --error-probability 0.08 --output mapped_reads.bam
```

For a list and descriptions of the basic command line options, run:
```
./floxer --help
```

For all available options, including ones intended exclusively for research and evaluation, run:
```
./floxer --advanced-help
```
