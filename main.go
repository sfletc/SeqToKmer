// MIT License

// Copyright (c) 2023 Stephen Fletcher

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
)

// fastaLineReader reads FASTA formatted data from an input reader and sends each sequence
// to a channel. It skips the sequence identifier lines (which start with '>') and
// concatenates all subsequent lines until the next identifier is reached. The sequences
// are sent to the 'sequences' channel without their FASTA headers. The reader buffers data
// up to 4MB to accommodate large lines.
//
// Parameters:
// input: an io.Reader from which the FASTA data is read
// sequences: a channel of strings to which sequences (not including headers) are sent.
// done: a channel of bool that signals when the reading process is complete.
//
// Usage:
// - The function should be run as a goroutine to allow concurrent processing of the input.
// - It assumes that input data is correctly formatted according to the FASTA standard.
// - If a line exceeds the 4MB buffer limit, it logs an error and returns early.
// - It is the caller's responsibility to close the 'done' channel once all sequences have been processed.
//
// It is essential to handle the 'sequences' and 'done' channels correctly to avoid deadlocks or lost data.
func fastaLineReader(input io.Reader, sequences chan<- string, done chan<- bool) {
	scanner := bufio.NewScanner(input)
	capacity := 4000 * 1000 // 4MB buffer
	scanner.Buffer(make([]byte, capacity), capacity)

	var builder strings.Builder // For efficient string concatenation

	for scanner.Scan() {
		line := scanner.Text()
		if len(line) == 0 {
			continue // Skip empty lines
		}
		if line[0] == '>' {
			// Send the current sequence to the channel if it's not empty
			if builder.Len() > 0 {
				sequences <- strings.ToUpper(builder.String()) // Convert to uppercase before sending
				builder.Reset()                                // Reset the builder for the next sequence
			}
		} else {
			// Append the line to the current sequence
			builder.WriteString(line)
		}
	}

	// Handle scanner errors
	if err := scanner.Err(); err != nil {
		if err == bufio.ErrTooLong {
			fmt.Println("Error: Encountered a line too long to fit in buffer")
		} else {
			fmt.Printf("Error while scanning: %v\n", err)
		}
		return // Exit the function on error
	}

	// Send any remaining sequence to the channel, making sure it is in uppercase
	if builder.Len() > 0 {
		sequences <- strings.ToUpper(builder.String()) // Convert to uppercase before sending
	}

	close(sequences) // Close the sequences channel when done
	done <- true     // Signal that the process is done
}

// fastqGzSeqLineReader reads sequences from a gzipped FASTQ formatted file provided as an input reader.
// It extracts the sequence lines, skipping the identifier, '+', and quality lines. The extracted sequences
// are then sent to the 'sequences' channel. This function is specifically tailored to parse gzipped FASTQ files,
// where every fourth line starting from the second line (1-based) is the sequence line.
//
// Parameters:
// - input: an io.Reader from which the gzipped FASTQ data is read, such as a file or stdin.
// - sequences: a channel of strings to which the extracted sequences are sent.
// - done: a channel of bool that signals when the reading and extraction process is complete.
//
// The function panics if it fails to initialize a gzip reader with the provided input or encounters
// an error during scanning the input.
//
// Usage:
//   - The function is intended to be run as a goroutine alongside other concurrent processes that
//     consume the 'sequences' channel output.
//   - The calling code is responsible for ensuring that the 'done' channel is checked after the 'sequences'
//     channel is closed to confirm that the reading process has completed successfully.
//
// Note:
//   - It is crucial to provide a correctly formatted gzipped FASTQ file to avoid unexpected behavior.
//   - The function panics on errors, so it should be used within a context where panics are recovered,
//     or the input is guaranteed to be error-free.
func fastqGzSeqLineReader(input io.Reader, sequences chan<- string, done chan<- bool) {
	gz, err := gzip.NewReader(input)
	if err != nil {
		panic(err)
	}
	defer gz.Close()

	scanner := bufio.NewScanner(gz)
	lineCount := 0

	for scanner.Scan() {
		lineCount++
		text := scanner.Text()

		if lineCount%4 == 2 {
			sequences <- strings.TrimSpace(text)
		}
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}
	close(sequences)
	done <- true
}

// splitNonACGT scans through sequences from a channel and segments them into substrings
// excluding any runs of characters that are not 'A', 'C', 'G', or 'T'. The segments
// that meet or exceed the specified minimum length are then sent to an output channel.
//
// Parameters:
// - minLength: the minimum length that a sequence segment must have to be considered valid.
// - sequences: a channel from which input sequences are received for segmentation.
// - splitSequences: a channel to which the resulting segments are sent if they meet the length criteria.
// - done: a channel used to signal the completion of the segmentation process.
//
// The function reads sequences from the 'sequences' channel, iterating over each rune (character)
// in the sequence. When it encounters a rune that is not 'A', 'C', 'G', or 'T', it checks if the
// segment prior to this rune meets the minimum length requirement. If it does, the segment is sent
// to 'splitSequences'. This continues until the function has processed all characters in the sequence.
//
// Usage:
//   - This function is designed to operate concurrently as a goroutine, processing sequences in parallel.
//   - Once all sequences have been read and appropriate segments have been sent to 'splitSequences',
//     the function signals its completion on the 'done' channel.
//
// Note:
//   - Sequences are handled in a streaming manner, suitable for processing large datasets without
//     requiring all data to be loaded into memory at once.
//   - Segments containing any non-AGCT characters or shorter than 'minLength' are excluded
//     from the output.
func splitNonACGT(minLength int, sequences <-chan string, splitSequences chan<- string, done chan<- bool) {
	for sequence, ok := <-sequences; ok; sequence, ok = <-sequences {
		start := 0
		for i, r := range sequence {
			if r != 'A' && r != 'C' && r != 'G' && r != 'T' {
				if start <= i-minLength {
					splitSequences <- sequence[start:i]
				}
				start = i + 1
			}
		}

		if start <= len(sequence)-minLength {
			splitSequences <- sequence[start:]
		}
	}

	done <- true

}

// makeKmers computes k-mers and their reverse complements from a given set of DNA sequences.
// Each k-mer is a substring of length 'k', and its reverse complement is also calculated.
// The smaller of each k-mer and its reverse complement is determined using bitwise operations
// to be efficiently stored and processed.
//
// The function takes four parameters:
// k: The size of the k-mers to generate.
// sequences: A channel that provides DNA sequences from which to generate k-mers.
// kmersOut: An array of channels through which the generated k-mers are sent to subscribers.
// done: A channel that signals the completion of the k-mer generation process.
//
// The function processes each sequence from the 'sequences' channel in turn.
// For each base in a sequence, it computes a numerical value that represents
// the k-mer ending at that base position (as well as the reverse complement).
// These values are stored in a slice of uint64s, which is sent through each output channel in 'kmersOut'.
// Once all sequences have been processed, the function signals completion by sending 'true' to the 'done' channel.
//
// Note: This function is designed for use in a concurrent environment, where sequences are supplied
// and k-mers are consumed by separate goroutines. Proper synchronization via channels is assumed.
func makeKmers(k int, sequences <-chan string, kmersOut []chan []uint64, done chan<- bool) {
	// toShift calculates the number of bit positions to shift for reverse complement calculation.
	toShift := uint((k - 1) * 2)
	// mask is used to isolate the bits corresponding to the k-mer.
	mask := (^uint64(0)) >> uint(64-k*2)

	// Iterate over all sequences provided via the channel.
	for seq := range sequences {
		// Initialize a slice to store k-mers for the current sequence.
		kmers := make([]uint64, len(seq)-k+1)
		var next, nextRC uint64 // next holds the current k-mer, nextRC holds its reverse complement.

		// Iterate over each base in the sequence.
		for i, b := range []byte(seq) {
			val := ((b >> 1) ^ ((b & 4) >> 2)) & 3                    // Convert base to 2-bit representation.
			next = ((next << 2) | uint64(val)) & mask                 // Calculate the next k-mer.
			nextRC = ((nextRC >> 2) | ((^uint64(val))<<toShift)&mask) // Calculate the reverse complement.

			// Once we have seen at least 'k' bases, start storing k-mers.
			if i >= k-1 {
				// Store the lexicographically smaller of the k-mer and its reverse complement.
				kmers[i-k+1] = min(next, nextRC)
			}
		}
		// Send the slice of k-mers to each channel in 'kmersOut'.
		for _, kout := range kmersOut {
			kout <- append(make([]uint64, 0, len(kmers)), kmers...)
		}
	}
	// Signal that k-mer generation is complete.
	done <- true
}

// min returns the smaller of two uint64 values.
// It's a utility function commonly used to compare two numerical values
// such as integers representing encoded DNA k-mers and their reverse complements.
//
// Parameters:
// a: The first uint64 value to compare.
// b: The second uint64 value to compare.
//
// Returns:
// The smaller of the two values 'a' and 'b'.
func min(a, b uint64) uint64 {
	if a < b {
		return a
	}
	return b
}

// countKmers counts the occurrences of k-mers in a channel of k-mer sequences
// applying a filtering mask to select only specific k-mers of interest.
//
// The function iterates over a channel of k-mer sequences, each represented as a slice of uint64 values,
// and counts the number of times each unique k-mer appears, subject to the specified filter criteria.
// It updates a map with the count of each k-mer.
//
// Parameters:
// - filter: The specific binary pattern that a k-mer must match to be counted.
// - filterMask: A bitmask to isolate bits of interest in a k-mer for comparison with the filter.
// - counts: A map from k-mer (uint64) to its count (int32) that stores the occurrence count of each k-mer.
// - kmerSequences: A channel that provides sequences of k-mers to be counted.
// - done: A channel used to signal the completion of the counting process.
// - minKmerCount: The minimum threshold for a k-mer count to be incremented. Once a k-mer reaches this count,
//                 it will no longer be incremented to prevent integer overflow and optimize performance.
//
// The function ignores k-mers that do not match the filter when the filterMask is applied.
// If a k-mer has previously been counted and its count is less than minKmerCount, its count is incremented.
// If it's a new k-mer, it's added to the map with an initial count of 1.
//
// Note that the function operates on a channel of k-mer sequences, allowing it to work efficiently
// in a concurrent context with multiple producers of k-mer sequences and a single consumer counting them.

func countKmers(filter, filterMask uint64, counts map[uint64]int32, kmerSequences <-chan []uint64, done chan<- bool, minKmerCount int) {
	for kmers := range kmerSequences {
		for _, kmer := range kmers {
			if (kmer & filterMask) != filter {
				continue
			}
			if val, ok := counts[kmer]; ok {
				if val < int32(minKmerCount) {
					counts[kmer]++
				}
			} else {
				counts[kmer] = 1
			}
		}
	}
	done <- true
}

// splitInputSequences distributes the task of splitting input strings based on non-ACGT characters
// across multiple worker goroutines. It sets up a pipeline where each worker receives input strings,
// processes them to split at any character that is not A, C, G, or T, and sends the resulting sequences
// through a channel to be further processed or consumed.
//
// Parameters:
// - k: The minimum length of a sequence segment to be considered valid and thus retained.
// - inputStrings: A channel of strings representing the input sequences that need to be processed.
// - sequences: A channel where processed sequences are sent after splitting.
// - numWorkers: The number of concurrent worker goroutines to be spawned for processing.
//
// Returns:
//   - doneSplits: A channel of boolean values used to signal the completion of processing by each worker.
//     The caller can listen on this channel to know when all workers have finished processing.
func splitInputSequences(k int, inputStrings chan string, sequences chan string, numWorkers int) chan bool {
	doneSplits := make(chan bool, numWorkers)
	for i := 0; i < numWorkers; i++ {
		go splitNonACGT(k, inputStrings, sequences, doneSplits)
	}
	return doneSplits
}

// generateKmers orchestrates the concurrent generation of k-mers from a stream of input sequences.
// It launches a specified number of worker goroutines that each execute the makeKmers function to
// convert sequences into their constituent k-mers.
//
// Parameters:
// - k: The size of the k-mers to generate.
// - sequences: A channel that provides the input sequences from which to generate k-mers.
// - kmerSequences: An array of channels that workers use to output the generated k-mers.
// - kmerWorkers: The number of worker goroutines to spawn for k-mer generation.
//
// Returns:
//   - doneKmers: A channel of boolean values used to signal the completion of k-mer generation by each worker.
//     Consumers of this function can listen on this channel to determine when all workers have
//     finished their k-mer generation tasks.
func generateKmers(k int, sequences chan string, kmerSequences []chan []uint64, kmerWorkers int) chan bool {
	doneKmers := make(chan bool, kmerWorkers)
	for i := 0; i < kmerWorkers; i++ {
		go makeKmers(k, sequences, kmerSequences, doneKmers)
	}
	return doneKmers
}

// countAllKmers launches multiple goroutines to count the occurrences of k-mers from multiple channels of k-mer sequences.
// Each goroutine applies a specific filter to only count relevant k-mers and populates a corresponding map with the count
// of each k-mer that passes the filter.
//
// Parameters:
// - numMaps: The number of k-mer count maps and corresponding goroutines to be used.
// - filterMask: A bitmask to be applied to each k-mer to filter out irrelevant ones.
// - kmerSequences: A slice of channels that supply the k-mer sequences to be counted.
// - allCounts: A slice of maps, where each map will hold the count of k-mers as determined by its corresponding goroutine.
// - minKmerCount: The minimum count threshold for a k-mer to be included in the final count map.
//
// Returns:
// - doneCounting: A channel that signals the completion of the counting process for each goroutine.
func countAllKmers(numMaps int, filterMask uint64, kmerSequenceChans []chan []uint64, allCounts []map[uint64]int32, minKmerCount int) chan bool {
	doneCounting := make(chan bool, numMaps)
	for i, seqs := range kmerSequenceChans {
		go countKmers(uint64(i), filterMask, allCounts[i], seqs, doneCounting, minKmerCount)
	}
	return doneCounting
}

// writeCountsToFile writes the counts of k-mers to a file, ensuring that only k-mers with a count
// equal to or above a specified threshold are written. The function is designed to serialize k-mer
// counts into a binary format, which can be later used for further analysis or for generating reports.
//
// Parameters:
// - allCounts: A slice of maps where each map contains k-mer counts.
// - finalFilePath: The path to the file where the k-mer counts should be written.
// - minKmerCount: The minimum count threshold for a k-mer to be included in the output file.
func writeCountsToFile(allCounts []map[uint64]int32, finalFilePath string, minKmerCount int) {
	finalFile, err := os.Create(finalFilePath)
	if err != nil {
		log.Fatal(err)
	}
	defer finalFile.Close()
	bufferedWriter := bufio.NewWriter(finalFile)

	var buffer bytes.Buffer
	chunkSize := 1024
	for _, counts := range allCounts {
		for kmer, count := range counts {
			if count >= int32(minKmerCount) {
				binary.Write(&buffer, binary.LittleEndian, kmer)
				if buffer.Len() >= chunkSize {
					if _, err := bufferedWriter.Write(buffer.Bytes()); err != nil {
						log.Fatal(err)
					}
					buffer.Reset()
				}
			}
		}
	}
	if _, err := bufferedWriter.Write(buffer.Bytes()); err != nil {
		log.Fatal(err)
	}
	bufferedWriter.Flush()
}

// kmerSearch performs a complete k-mer counting workflow, from reading input sequences
// to writing the final k-mer counts to a file. The function is designed to be efficient
// by using concurrent processing and channel communication.
//
// Parameters:
// - k: The length of the k-mers to be searched.
// - inputSequences: A channel that receives the input sequences to be processed.
// - doneInputSequences: A channel that signals the completion of input sequence processing.
// - numMaps: The number of concurrent maps to use for k-mer counting.
// - finalFilePath: The path to the file where the final k-mer counts should be written.
// - minKmerCount: The minimum count threshold for a k-mer to be included in the output file.
//
// The function sets up a multi-stage pipeline:
// 1. Splitting the input sequences into subsequences that do not contain non-AGCT bases.
// 2. Generating k-mers from the split subsequences.
// 3. Counting the occurrences of each unique k-mer.
// 4. Writing the counts to a file, filtering out k-mers with counts below minKmerCount.
func kmerSearch(k int, inputSequences chan string, doneInputSequences chan bool, numMaps int, finalFilePath string, minKmerCount int) {
	mask := uint64(numMaps - 1)
	splitWorkers := 4
	splitSequences := make(chan string, 3)
	doneSplitSequences := splitInputSequences(k, inputSequences, splitSequences, splitWorkers)
	kmerSearchWorkers := 8
	kmerSequences := make([]chan []uint64, numMaps)
	for i := 0; i < numMaps; i++ {
		kmerSequences[i] = make(chan []uint64, 3)
	}
	doneKmerID := generateKmers(k, splitSequences, kmerSequences, kmerSearchWorkers)
	allCountMaps := make([]map[uint64]int32, numMaps)
	for i := range allCountMaps {
		allCountMaps[i] = make(map[uint64]int32)
	}
	doneKmerCounting := countAllKmers(numMaps, mask, kmerSequences, allCountMaps, minKmerCount)

	<-doneInputSequences

	for i := 0; i < splitWorkers; i++ {
		<-doneSplitSequences
	}
	close(splitSequences)

	for i := 0; i < kmerSearchWorkers; i++ {
		<-doneKmerID
	}
	for _, s := range kmerSequences {
		close(s)
	}

	for i := 0; i < numMaps; i++ {
		<-doneKmerCounting
	}

	writeCountsToFile(allCountMaps, finalFilePath, minKmerCount)
}

// kmerToSequence converts a numeric representation of a k-mer back into its
// sequence of nucleotides as a string. Each k-mer is represented as a 64-bit
// unsigned integer where every 2 bits correspond to a nucleotide in the order
// of ACGT (i.e., A=00, C=01, G=10, T=11).
//
// Parameters:
// - kmer: The 64-bit unsigned integer representing the k-mer to be converted.
// - k: The length of the k-mer sequence to be reconstructed.
//
// Returns:
// - A string representing the sequence of nucleotides for the given k-mer.
func kmerToSequence(kmer uint64, k int) string {
	bases := "ACGT"
	sequence := make([]byte, k)

	for i := k - 1; i >= 0; i-- {
		sequence[i] = bases[kmer&3]
		kmer >>= 2
	}
	return string(sequence)
}

// readKmersFromDisk reads k-mers encoded as unsigned 64-bit integers from a binary file
// and converts each k-mer into its respective DNA sequence string. This function assumes
// that the binary file contains a sequence of 64-bit integers where each integer
// represents a k-mer. It returns a slice of strings representing the sequences and an error, if any.
//
// Parameters:
// - filename: The name of the binary file containing the k-mers.
// - k: The length of the k-mers that are stored in the file.
//
// Returns:
// - A slice of strings representing the sequences of nucleotides for the k-mers.
// - An error if any occurs during file operation or reading.
func readKmersFromDisk(filename string, k int) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	var sequences []string
	for {
		var kmer uint64
		err := binary.Read(file, binary.LittleEndian, &kmer)
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, fmt.Errorf("error reading binary data: %v", err)
		}

		sequence := kmerToSequence(kmer, k)
		sequences = append(sequences, sequence)
	}

	return sequences, nil
}

// main orchestrates the execution of a k-mer counting pipeline for DNA sequences.
// It processes command-line arguments to determine the mode of operation, input and output file paths,
// k-mer length, and the minimum count of k-mers for filtering.
//
// If the loadAndPrint flag is set, it will read and print k-mers from the specified binary output file.
// Otherwise, it will initiate the k-mer counting process, which involves:
//   - Opening and validating the input file.
//   - Reading sequences from the input file in a format-specific manner.
//   - Counting k-mers across the sequences using a parallelized approach.
//   - Saving the counted k-mers to a binary file.
//
// Input and output file paths are mandatory unless the loadAndPrint flag is set, in which case
// only the output file path is required. The function provides error messages for unsupported file formats
// and other input validation failures.
func main() {

	numMaps := 32

	inputFilePath := flag.String("input", "", "Path to the input FASTA/FASTQ file")
	kmerLength := flag.Int("k", 21, "K-mer length")
	minKmerCount := flag.Int("minKmerCount", 0, "Minimum k-mer count for filtering")
	outputFilePath := flag.String("output", "", "Path to the output binary file")
	loadAndPrint := flag.Bool("load", false, "Load k-mers from binary file and print to stdout")
	flag.Parse()

	if *loadAndPrint {
		if *outputFilePath == "" {
			fmt.Println("Please provide an output binary file path to load k-mers from.")
			return
		}
	} else {
		if *inputFilePath == "" || *outputFilePath == "" {
			fmt.Println("Please provide both an input file path and an output file path.")
			return
		}
	}

	var inputStrings chan string = make(chan string)
	var doneLines chan bool = make(chan bool)
	file, err := os.Open(*inputFilePath)
	if err != nil {
		fmt.Printf("Error opening file: %v\n", err)
		return
	}
	defer file.Close()

	if strings.HasSuffix(*inputFilePath, ".fa") || strings.HasSuffix(*inputFilePath, ".fna") || strings.HasSuffix(*inputFilePath, ".fasta") {
		go fastaLineReader(file, inputStrings, doneLines)
	} else if strings.HasSuffix(*inputFilePath, ".fq.gz") || strings.HasSuffix(*inputFilePath, ".fastq.gz") {
		go fastqGzSeqLineReader(file, inputStrings, doneLines)
	} else {
		fmt.Println("Error: Unsupported file format.")
		return
	}

	kmerSearch(*kmerLength, inputStrings, doneLines, numMaps, *outputFilePath, *minKmerCount)
	fmt.Println("Done counting and saving to disk")

	// Optionally load and print if the flag is set
	if *loadAndPrint {
		fmt.Printf("Loading output file: %v\n", *outputFilePath)
		readKmersFromDisk(*outputFilePath, *kmerLength)
	}
}
