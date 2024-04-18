## SeqToKmer: a fast nucleotide sequence to kmer tool.        

Developed to process FASTA or fastq.gz files containing off-target sequences for use by [dsRNAmax](https://github.com/sfletc/dsRNAmax)


### Key Features

- Parallel Processing: Utilizes Go's concurrency primitives (goroutines and channels) to optimize performance by distributing sequence processing and k-mer counting across multiple workers.
- Configurable K-mer Length: Adapt the 'k' parameter to analyze kmers of varying lengths.
- Filtering: Apply a minimum kmer count filter to reduce retention of erroneous kmers due to sequencing errors
- Efficient Binary Output: Store k-mer counts compactly in a binary file format for later retrieval and analysis.
- Sequence Reconstruction: Convert stored k-mer numerical representations back to their original DNA sequences.