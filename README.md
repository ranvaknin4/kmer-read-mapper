# K-mer Read Mapper

A command-line tool for mapping DNA sequencing reads to reference genomes using k-mers.

The tool builds a k-mer index from a set of reference genomes (FASTA), then maps sequencing reads (FASTQ) to those genomes and classifies each read as **uniquely mapped**, **ambiguously mapped** or **unmapped**.

Developed as the final project of the Introduction to Computer Science course at the Hebrew University of Jerusalem (grade: 94).

## How it works

1. **Indexing:** each reference genome is split into overlapping k-mers. A dictionary maps every k-mer to the genomes and positions where it appears.
2. **Mapping:** each read is split into k-mers, which are looked up in the index. K-mers that appear in only one genome are *specific*; the rest are *non-specific*.
3. **Classification:**
   - A read with no specific k-mers is **unmapped**.
   - A read whose specific k-mers point to one genome, or clearly favor one genome (by at least the unique threshold), is **uniquely mapped**.
   - Otherwise it is **ambiguously mapped**.
   - An optional second step re-checks uniquely mapped reads against all their k-mers, using the ambiguous threshold.

## Requirements

- Python 3.10+
- NumPy

```bash
pip install numpy
```

## Usage

The tool has four tasks, selected with `-t`:

| Task | Description |
|------|-------------|
| `reference` | Build a k-mer index from a FASTA file and save it to a compressed file |
| `dumpref` | Print a summary of the index as JSON |
| `align` | Map reads from a FASTQ file and save the results to a compressed file |
| `dumpalign` | Print the mapping statistics as JSON |

### Examples

Build a reference index with k-mers of length 31:
```bash
python3 main.py -t reference -g genomes.fa -k 31 -r reference.kdb
```

Print the index summary:
```bash
python3 main.py -t dumpref -r reference.kdb
```

Map reads to the reference:
```bash
python3 main.py -t align -r reference.kdb --reads reads.fq -a results.aln
```

Print the mapping statistics:
```bash
python3 main.py -t dumpalign -a results.aln
```

`dumpref` and `dumpalign` can also run directly from the input files (a FASTA file with `-k`, and a FASTQ file with `--reads`) without saving an intermediate file.

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-g`, `--genomefile` | Reference genomes FASTA file | |
| `-k`, `--kmer-size` | K-mer length | |
| `-r`, `--referencefile` | Reference index file (input or output) | |
| `--reads` | Reads FASTQ file | |
| `-a`, `--alignfile` | Mapping results file (input or output) | |
| `-m`, `--unique-threshold` | Minimum difference in specific k-mers for a unique mapping | 1 |
| `-p`, `--ambiguous-threhold` | Threshold for the second classification step | 1 |

### Optional filters

| Option | Description |
|--------|-------------|
| `--min-read-quality` | Skip reads whose mean quality score is below this value |
| `--min-kmer-quality` | Skip k-mers whose mean quality score is below this value |
| `--max-genomes` | Skip k-mers that appear in more than this number of genomes |
| `--filter-similar` | Remove highly similar genomes from the reference |
| `--similarity-threshold` | Similarity score (0–1) above which a genome is removed (default: 0.95) |

## Project structure

```
main.py                 Command-line interface and argument validation
commands_implement.py   Implementation of the four tasks
helper.py               FASTA/FASTQ parsing, indexing, filtering and read classification
classes.py              Genome, Read and KmerReference classes
```
