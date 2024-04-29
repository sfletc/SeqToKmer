## SeqToKmer: a fast nucleotide sequence to kmer tool.        

Processes a FASTA or fastq.gz file containing off-target sequences for to a kmer binary for use by [dsRNAmax](https://github.com/sfletc/dsRNAmax)


## Key Features

- Configurable K-mer Length, with the selected kmer length written to file for decoding upon load
- Filtering: Apply a minimum kmer count filter to reduce retention of erroneous kmers due to sequencing errors
- Efficient Binary Output: Store k-mer counts compactly in a binary file format for later retrieval and analysis.
- Sequence Reconstruction: Convert stored k-mer numerical representations back to their original nucelotide sequences.

## Usage

```
seqToKmer -h
Usage of seqToKmer:
  -input string
    	Path to the input FASTA/FASTQ file
  -k int
    	K-mer length (default 21)
  -load
    	Load k-mers from binary file and print to stdout
  -minKmerCount int
    	Minimum k-mer count for filtering
  -output string
    	Path to the output binary file
```

## Coming

- optionally write interim maps to disk to save RAM at the cost of speed
