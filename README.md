> [!NOTE]
> **Archive Notice**: This project is an educational and reference implementation of a memory-aware De Bruijn graph assembler in modern C++20. It is finished and archived; no active maintenance or feature development is planned.

# De Bruijn Graph Assembler

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

The bundled `data/example_reads.fasta` contains overlapping reads derived from the
*Escherichia coli* K-12 MG1655 reference genome (`U00096.3`, positions 21-65),
including one single-base mismatch read (`*_err1`) to mimic sequencing noise.

Supported `k` values are compile-time dispatched: `15,17,19,21,23,25,27,29,31`.

## Toy Example Walkthrough

A clean, minimal toy dataset is provided in [`data/toy_reads.fasta`](file:///home/andreassag/repo/bruijn/data/toy_reads.fasta) for demonstration and testing.

Target sequence (22 bp):
```text
TACGATCGATCGATCGAGCTAG
```

The dataset consists of 17-bp reads sliding across the target with single-base offsets (16-bp overlap):
1. `TACGATCGATCGATCGA`
2. `ACGATCGATCGATCGAG`
3. `CGATCGATCGATCGAGC`
4. `GATCGATCGATCGAGCT`
5. `ATCGATCGATCGAGCTA`
6. `TCGATCGATCGAGCTAG`
7. `TACGATCGATCGATCGA` (duplicate terminal read)
8. `TCGATCGATCGAGCTAG` (duplicate terminal read)

### Step-by-Step Graph Assembly with $k = 15$

1. **K-mer Extraction & Counting ($k=15$, `--min-count 2`)**:
   - Each 17-bp read yields $(17 - 15 + 1) = 3$ consecutive 15-mers.
   - Sliding through all 8 reads produces 8 distinct 15-mers spanning the 22-bp target:
     - $k_1$: `TACGATCGATCGATC` (count = 2, from reads 1, 7)
     - $k_2$: `ACGATCGATCGATCG` (count = 3, from reads 1, 2, 7)
     - $k_3$: `CGATCGATCGATCGA` (count = 4, from reads 1, 2, 3, 7)
     - $k_4$: `GATCGATCGATCGAG` (count = 3, from reads 2, 3, 4)
     - $k_5$: `ATCGATCGATCGAGC` (count = 3, from reads 3, 4, 5)
     - $k_6$: `TCGATCGATCGAGCT` (count = 4, from reads 4, 5, 6, 8)
     - $k_7$: `CGATCGATCGAGCTA` (count = 3, from reads 5, 6, 8)
     - $k_8$: `GATCGATCGAGCTAG` (count = 2, from reads 6, 8)
   - All 8 k-mers have count $\ge 2$, so all 8 are retained. Total coverage sum across all 8 nodes is $2 + 3 + 4 + 3 + 3 + 4 + 3 + 2 = 24$.

2. **Graph Construction**:
   - 8 nodes are created (one per k-mer).
   - Overlapping consecutive k-mers form directed edges via $(k-1)$-mer overlaps (14 bp):
     $k_1 \xrightarrow{\text{G}} k_2 \xrightarrow{\text{A}} k_3 \xrightarrow{\text{G}} k_4 \xrightarrow{\text{C}} k_5 \xrightarrow{\text{T}} k_6 \xrightarrow{\text{A}} k_7 \xrightarrow{\text{G}} k_8$.
   - Yields a linear chain of 8 nodes and 7 edges with no branches.

3. **Graph Cleaning**:
   - Since the graph is a single unbranched path of length 8 ($> 4$), `--tip-length 4` retains the path intact.
   - No bubbles are present.

4. **Contig Extraction**:
   - The linear path is walked from source $k_1$ to sink $k_8$, appending each transition base.
   - Result: 1 contig of length $15 + 7 = 22\text{ bp}$ (`TACGATCGATCGATCGAGCTAG`) with average coverage $\text{cov} = 3$.

Execute the toy assembly:
```bash
./build/bruijn_assembler \
  --input data/toy_reads.fasta \
  --output data/toy_contigs.fasta \
  --k 15 \
  --min-count 2 \
  --threads 1 \
  --tip-length 4
```

## K-mer 2-Bit Representation

K-mers are encoded into compact 64-bit unsigned integers (`uint64_t`), supporting values of $k$ from 1 to 31 ($2k \le 62\text{ bits}$):

### Base Encoding & Bit Operations
Nucleotides are mapped to 2-bit values:
- `A` / `a` $\rightarrow 00_2$ (`0`)
- `C` / `c` $\rightarrow 01_2$ (`1`)
- `G` / `g` $\rightarrow 10_2$ (`2`)
- `T` / `t` $\rightarrow 11_2$ (`3`)
- Non-ACGT characters (`N`, etc.) invalidate the rolling window.

Key bitwise operations in [`include/bruijn/dna.hpp`](file:///home/andreassag/repo/bruijn/include/bruijn/dna.hpp):
- **Masking**: `kmer_mask<K>() = (1ULL << (2 * K)) - 1ULL` isolates the $2k$ active bits.
- **Base complement**: `complement_bits(b) = 0x3U - (b & 0x3U)` (e.g. $00_2 (\text{A}) \leftrightarrow 11_2 (\text{T})$, $01_2 (\text{C}) \leftrightarrow 10_2 (\text{G})$).
- **Shift-left append**: `((kmer << 2) & mask) | new_base` shifts left to add a successor base at the 3' end.
- **Shift-right prepend**: `(kmer >> 2) | (new_base << (2 * (K - 1)))` shifts right to add a predecessor base at the 5' end.

### Rolling Window Extraction
In `bruijn::dna::for_each_kmer<K>(sequence, fn)`:
- Scans `sequence` linearly in $O(N)$ time.
- Maintains `rolling = ((rolling << 2) | bits) & mask`.
- Consecutive valid bases increment a counter `valid_run`. When `valid_run >= K`, the current rolling integer is emitted without string allocations. If an invalid base is encountered, `valid_run` resets to 0.

### Canonicalization
DNA sequencing reads can originate from either forward or reverse strands:
- `reverse_complement<K>(kmer)` iterates through the $K$ 2-bit bases in reverse order and inverts each base via `complement_bits`.
- `canonical_kmer<K>(kmer) = std::min(kmer, reverse_complement<K>(kmer))`.
- When `--canonical-kmers` is enabled, abundance filtering aggregates forward and reverse-complement counts together, while directional graph nodes are preserved during assembly.

## Memory Footprint & Data Structures

The assembler is engineered for low memory overhead across all processing phases:

### 1. Read Storage (`ReadBuffer`)
- **Contiguous Byte Buffer**: `std::vector<uint8_t> bases_` stores all raw nucleotide bytes consecutively without per-string heap metadata or null terminators.
- **POD Read Spans**: `std::vector<ReadSpan> spans_` stores lightweight POD descriptors:
  ```cpp
  struct ReadSpan {
      uint32_t offset;
      uint32_t length;
  };
  ```
- **Memory**: Exactly 8 bytes per read in `spans_` plus 1 byte per base in `bases_`. For example, 10 million 150-bp reads require $\approx 80\text{ MB} + 1.5\text{ GB} \approx 1.58\text{ GB}$.

### 2. Counting Bloom Filter Prefilter (`CountingBloomFilter`)
- Used in Pass 1 of k-mer counting to filter out singleton/noise k-mers before inserting into hash tables.
- Implemented as an array of 8-bit saturating counters (`std::vector<uint8_t> counters_`).
- Default sizing: $2^{24}$ counters (16,777,216 bytes = 16 MB) with 4 hash functions generated via `splitmix64`.
- Multi-threading: Each thread allocates a local Bloom filter; thread-local filters are merged via saturating addition into a shared global filter.

### 3. Compact Open-Addressed Hash Map (`CompactKmerMap`)
- A power-of-two open-addressed hash table with linear probing and a 70% max load factor.
- Avoids pointer indirection, cache misses, and per-node overhead inherent in `std::unordered_map`.
- Slot layout:
  ```cpp
  struct Slot {
      uint64_t key = 0;   // 8 bytes (2-bit encoded k-mer)
      TValue value{};     // e.g. uint32_t count or node index (4 bytes)
      uint8_t used = 0;   // 1 byte occupancy flag (+ 3 bytes padding)
  }; // total: 16 bytes per slot
  ```
- **Memory**: At 70% load factor, an index of $N$ unique k-mers uses $\approx \frac{16\text{ bytes}}{0.7} \approx 22.8\text{ bytes}$ per unique k-mer.

### 4. Graph Representation (`DeBruijnGraph`)
- **Node Storage**: Nodes are packed into a flat `std::vector<Node>`:
  ```cpp
  struct Node {
      uint64_t kmer = 0;  // 8 bytes (k-mer sequence)
      uint32_t count = 0; // 4 bytes (abundance coverage)
      uint8_t in_mask = 0; // 4-bit bitmask of incoming edge bases (A=1, C=2, G=4, T=8)
      uint8_t out_mask = 0;// 4-bit bitmask of outgoing edge bases (A=1, C=2, G=4, T=8)
      uint8_t alive = 1;  // 1 byte deletion marker
  }; // total: 16 bytes per node (with alignment padding)
  ```
- **4-bit Edge Bitmasks**: Edges are not stored as pointer lists or separate 8-byte edge objects. Instead, incoming and outgoing edges for A, C, G, T are encoded directly in 4 bits each (`in_mask` and `out_mask`), allowing $O(1)$ in/out degree calculation via `std::popcount` and $O(1)$ edge existence checks.
- **K-mer to Node Index**: A `CompactKmerMap<uint32_t> index_` maps each $k$-mer key to its `uint32_t` index in `nodes_` (16 bytes per slot).

### Peak Memory Across Pipeline Phases

| Phase | Active Data Structures | Dominant Memory Components |
|---|---|---|
| **Phase 1: Read Loading** | `ReadBuffer` | $\text{total\_bases} \times 1\text{ B} + \text{read\_count} \times 8\text{ B}$ |
| **Phase 2: Bloom Pre-counting** | `ReadBuffer` + thread-local Bloom filters | ReadBuffer + $\text{threads} \times 16\text{ MB}$ |
| **Phase 3: Exact Counting** | `ReadBuffer` + merged Bloom filter + `CompactKmerMap` | ReadBuffer + $16\text{ MB} + \text{candidate\_kmers} \times 22.8\text{ B}$ |
| **Phase 4: Graph Construction & Cleaning** | `ReadBuffer` + `nodes_` (`std::vector<Node>`) + `index_` (`CompactKmerMap<uint32_t>`) | ReadBuffer + $|V| \times 16\text{ B} + (|V| / 0.7) \times 16\text{ B} \approx \text{ReadBuffer} + |V| \times 38.8\text{ B}$ |
| **Phase 5: Contig Extraction** | `nodes_` + `visited_out` (`std::vector<uint8_t>`) + output contigs | Negligible additional memory ($|V| \times 1\text{ B}$) |

*Note*: Because Bloom filters drop low-frequency k-mers before Phase 3, the exact map only expands for candidate k-mers, dramatically curbing peak memory on large, noisy sequencing runs.

## Tests

Run the test suite via `ctest` or directly via the compiled binary:

```bash
ctest --test-dir build --output-on-failure
# or directly:
./build/bruijn_tests
```

The test target checks:

- k-mer extraction behavior and canonical reverse-complement collapsing
- tip clipping on dead-end sequencing artifacts
- bounded-BFS bubble popping on a multi-branch reconvergent structure
- end-to-end toy dataset assembly verifying exact contig reconstruction and graph invariants

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
