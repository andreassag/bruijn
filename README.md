# De Bruijn Graph Assembler (C++20)

A memory-aware De Bruijn graph assembler implemented in modern C++20.

This implementation focuses on practical memory efficiency by combining:

- POD read storage (`uint8_t` contiguous base buffer + POD spans)
- Two-pass multi-threaded k-mer counting with a counting Bloom filter prefilter
- Optional canonical k-mer filtering (reverse-complement collapsing for abundance thresholding)
- Compact open-addressed hash maps for k-mer and node indexing
- Graph cleaning passes for short tips and bounded BFS bubble resolution
- Contig extraction from maximal non-branching paths

## Features

- K-mer Counting
  - Rolling 2-bit encoding (`A,C,G,T -> 0..3`) into `uint64_t`
  - Two-pass counting:
    - Pass 1: approximate frequency with counting Bloom filters (thread-local then merged)
    - Pass 2: exact counting only for candidate k-mers above threshold
  - Optional `--canonical-kmers` mode:
    - Reverse-complement pairs are collapsed while applying the abundance cutoff
    - Directional k-mer nodes are still counted for graph traversal
  - Multi-threaded with `std::thread`

- Graph Cleaning
  - Tip clipping:
    - Removes short dead-end chains from sources and sinks
  - Bubble popping:
    - Bounded BFS explores candidate paths from each outgoing branch
    - Supports reconvergence across more than two branches
    - Removes lower-coverage alternatives while keeping the strongest path

- Contig Generation
  - Traverses maximal non-branching paths
  - Handles isolated cycles
  - Outputs contigs in FASTA format sorted by length

## Project Layout

- `include/bruijn/dna.hpp`
  - DNA encoding/decoding helpers and rolling k-mer iteration
- `include/bruijn/read_buffer.hpp`, `src/read_buffer.cpp`
  - POD-centric read storage and FASTA/FASTQ parsing
- `include/bruijn/bloom_filter.hpp`
  - Counting Bloom filter
- `include/bruijn/compact_kmer_map.hpp`
  - Compact open-addressed hash map for k-mer keys
- `include/bruijn/kmer_counter.hpp`
  - Multi-threaded k-mer counting pipeline
- `include/bruijn/debruijn_graph.hpp`
  - Graph construction, cleaning, and contig traversal
- `include/bruijn/assembler.hpp`
  - High-level orchestration
- `src/main.cpp`
  - CLI entry point

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

## Run

```bash
./build/bruijn_assembler \
  --input data/example_reads.fasta \
  --output data/example_contigs.fasta \
  --k 15 \
  --min-count 2 \
  --canonical-kmers \
  --threads 2 \
  --tip-length 4 \
  --bubble-depth 32
```

Supported `k` values are compile-time dispatched: `15,17,19,21,23,25,27,29,31`.

## Notes on Memory

- Reads are stored in one contiguous byte buffer (`std::vector<uint8_t>`) with POD spans.
- The Bloom filter reduces peak memory for exact hash-map counting by filtering likely singleton/noise k-mers.
- K-mers are 2-bit encoded in `uint64_t` for compact representation and fast shifts.

## Tests

```bash
ctest --test-dir build --output-on-failure
```

The test target checks:

- k-mer extraction behavior and canonical reverse-complement collapsing
- tip clipping on dead-end sequencing artifacts
- bounded-BFS bubble popping on a multi-branch reconvergent structure

## Input Formats

- FASTA
- FASTQ (single-line sequence supported; quality can span multiple lines)

## Output

FASTA file with records like:

```text
>contig_1 len=1234 cov=42
ACGT...
```

where `cov` is mean node coverage across the contig path.
