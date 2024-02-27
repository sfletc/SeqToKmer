package main

import (
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"reflect"
	"strings"
	"sync"
	"testing"
	"time"

	"os"
)

func TestSeqLineReader(t *testing.T) {
	// Simulating the content of test.fa
	input := strings.NewReader(">sequence_1\nATCG\n>sequence_2\nGGTT\n")
	lines := make(chan string)
	done := make(chan bool)

	go fastaLineReader(input, lines, done)

	received := make([]string, 0)
	loopDone := false
	for !loopDone {
		select {
		case line, ok := <-lines:
			if ok {
				received = append(received, line)
			}
		case success := <-done:
			if !success {
				t.Errorf("An error occurred during scanning")
			}
			loopDone = true // End the loop once the done channel has signaled
		}
	}

	// Check the number of received sequences
	if len(received) != 2 {
		t.Errorf("Expected 2 sequences, got %d", len(received))
	}

	// Check the actual sequences
	if received[0] != "ATCG" || received[1] != "GGTT" {
		t.Errorf("Received sequences are incorrect, got %v", received)
	}
}
func TestFastqGzSeqLineReader(t *testing.T) {
	// Create a buffer containing the gzipped fastq content
	var buf bytes.Buffer
	writer := gzip.NewWriter(&buf)
	_, err := writer.Write([]byte("@seq1\nACTG\n+\n!!@@\n@seq2\nTGAC\n+\n@@!!\n"))
	if err != nil {
		t.Fatalf("Failed to write to buffer: %v", err)
	}
	if err := writer.Close(); err != nil {
		t.Fatalf("Failed to close writer: %v", err)
	}

	// Create channels
	lines := make(chan string)
	done := make(chan bool)

	// Start the reader function
	go fastqGzSeqLineReader(&buf, lines, done)

	// Collect lines
	var receivedLines []string
	timeout := time.After(2 * time.Second)
	for {
		select {
		case line, ok := <-lines:
			if ok {
				t.Logf("Received line: %s", line) // Debugging line
				receivedLines = append(receivedLines, line)
			}
		case <-done:
			t.Logf("Received done signal") // Debugging line
			goto doneReading               // Using goto to jump out of select and for loop
		case <-timeout:
			t.Fatal("Test timed out")
			return
		}
	}

doneReading:

	// Check if the lines received are as expected
	expectedLines := []string{"ACTG", "TGAC"}
	if !reflect.DeepEqual(receivedLines, expectedLines) {
		t.Errorf("Expected lines %v, got %v", expectedLines, receivedLines)
	}
}

func TestSplitNs(t *testing.T) {
	// Set up channels and data
	lines := make(chan string, 3)
	splitLines := make(chan string, 3)
	done := make(chan bool, 1)

	// Fill the lines channel with example data
	lines <- "AACTGA"
	lines <- "TACG"
	lines <- "TTNNGGT"
	close(lines)
	var wg sync.WaitGroup
	wg.Add(1) // Add a job for the wait group
	go func() {
		splitNonACGT(3, lines, splitLines, done)
		wg.Done() // Signal job completion
	}()

	go func() {
		wg.Wait()
		close(splitLines)
	}()

	output := []string{}
	for line := range splitLines {
		output = append(output, line)
	}

	<-done // Wait for the done signal

	// Expected output
	expected := []string{"AACTGA", "TACG", "GGT"}

	// Validate
	if !reflect.DeepEqual(output, expected) {
		t.Errorf("Expected %v, but got %v", expected, output)
	}
}

func TestMakeKmers(t *testing.T) {
	// Setup
	numMaps := 4
	kmerWorkers := 8
	sequences := make(chan string, 3)
	k := 2

	// Sample input sequence and expected kmers
	inputData := "ACGT"
	expectedKmers := []uint64{1, 6, 1} // TODO: Update with correct expected kmers

	sequences <- inputData

	kmerSequences := make([]chan []uint64, numMaps)
	for i := 0; i < numMaps; i++ {
		kmerSequences[i] = make(chan []uint64, 3)
	}

	doneKmers := make(chan bool, kmerWorkers)
	for i := 0; i < kmerWorkers; i++ {
		go makeKmers(k, sequences, kmerSequences, doneKmers)
	}

	close(sequences)

	// Wait for all kmer workers to finish
	for i := 0; i < kmerWorkers; i++ {
		<-doneKmers
	}
	for _, s := range kmerSequences {
		close(s)
	}

	// Collect and validate results
	collectedKmers := make([]uint64, 0, len(expectedKmers))
	for i := 0; i < numMaps; i++ {
		for kmer := range kmerSequences[i] {
			collectedKmers = append(collectedKmers, kmer...)
			if len(collectedKmers) >= len(expectedKmers) {
				break
			}
		}
		if len(collectedKmers) >= len(expectedKmers) {
			break
		}
	}

	if !reflect.DeepEqual(expectedKmers, collectedKmers) {
		t.Errorf("Expected kmers: %v, Got: %v", expectedKmers, collectedKmers)
	}
}

func TestCountKmers(t *testing.T) {
	// Setup test cases
	tests := []struct {
		filter        uint64
		filterMask    uint64
		kmerSequences [][]uint64
		minKmerCount  int
		expected      map[uint64]int32
	}{
		{
			filter:        1,
			filterMask:    1,
			kmerSequences: [][]uint64{{1, 5, 5, 3, 3, 3, 7}},
			minKmerCount:  2,
			expected:      map[uint64]int32{1: 1, 3: 2, 5: 2, 7: 1},
		},
	}

	for _, tt := range tests {
		counts := make(map[uint64]int32)
		kmerSequencesChan := make(chan []uint64, len(tt.kmerSequences))

		// Populate the kmerSequences channel.
		for _, seq := range tt.kmerSequences {
			kmerSequencesChan <- seq
		}
		// Close the channel to signal the end of data.
		close(kmerSequencesChan)

		// Now, start the countKmers function.
		done := make(chan bool)
		go countKmers(tt.filter, tt.filterMask, counts, kmerSequencesChan, done, tt.minKmerCount)

		// Wait for countKmers to complete.
		<-done

		// Validate the results.
		if !reflect.DeepEqual(counts, tt.expected) {
			t.Errorf("Expected: %+v, Got: %+v", tt.expected, counts)
		}

		// Teardown (clean-up).
		close(done)
	}
}

func TestSplitInputSequences(t *testing.T) {
	// Setup
	numWorkers := 4
	k := 2
	inputStrings := make(chan string, 3)
	sequences := make(chan string, 3)

	// Sample input data
	inputData := []string{"ACGT", "TGCA", "ANCGTT"}

	// Send input data to inputStrings channel
	for _, data := range inputData {
		inputStrings <- data
	}
	close(inputStrings)

	// Call the function to be tested
	doneSplits := splitInputSequences(k, inputStrings, sequences, numWorkers)

	// Wait for all splits to be done
	for i := 0; i < numWorkers; i++ {
		<-doneSplits
	}

	close(sequences)

	// Collect results from sequences channel
	collectedSequences := make([]string, 0, len(inputData)*2)
	for seq := range sequences {
		collectedSequences = append(collectedSequences, seq)
	}

	// Validate results
	expectedSequences := []string{"ACGT", "TGCA", "CGTT"}
	if !reflect.DeepEqual(expectedSequences, collectedSequences) {
		t.Errorf("Expected sequences: %v, Got: %v", expectedSequences, collectedSequences)
	}
}

func TestGenerateKmers(t *testing.T) {
	// Setup
	numMaps := 1
	kmerWorkers := 8
	sequences := make(chan string, 3)
	k := 2

	// Sample input sequence and expected kmers
	inputData := "ACGT"
	expectedKmers := []uint64{1, 6, 1}

	sequences <- inputData
	close(sequences)

	kmerSequences := make([]chan []uint64, numMaps)
	for i := 0; i < numMaps; i++ {
		kmerSequences[i] = make(chan []uint64, 3)
	}

	// Call the function to be tested
	doneKmers := generateKmers(k, sequences, kmerSequences, kmerWorkers)

	// Wait for all kmer generation to be done
	for i := 0; i < kmerWorkers; i++ {
		<-doneKmers
	}
	for _, s := range kmerSequences {
		close(s)
	}

	// Collect and validate results
	collectedKmers := make([]uint64, 0, len(expectedKmers)*numMaps)
	for i := 0; i < numMaps; i++ {
		for kmer := range kmerSequences[i] {
			collectedKmers = append(collectedKmers, kmer...)
		}
	}

	if !reflect.DeepEqual(expectedKmers, collectedKmers) {
		t.Errorf("Expected kmers: %v, Got: %v", expectedKmers, collectedKmers)
	}
}

func TestCountAllKmers(t *testing.T) {
	// Setup
	numMaps := 4
	filterMask := uint64(3)
	kmerSequences := make([]chan []uint64, numMaps)
	allCounts := make([]map[uint64]int32, numMaps)
	minKmerCount := 2

	for i := 0; i < numMaps; i++ {
		kmerSequences[i] = make(chan []uint64, 3)
		allCounts[i] = make(map[uint64]int32)
	}

	// Sample input kmers
	inputData := []uint64{1, 2, 3, 4, 5, 6, 7, 8}
	expectedCounts := []map[uint64]int32{
		{4: 1, 8: 1},
		{1: 1, 5: 1},
		{2: 1, 6: 1},
		{3: 1, 7: 1},
	}

	for i := 0; i < numMaps; i++ {
		go func(i int) {
			kmerSequences[i] <- inputData
			close(kmerSequences[i])
		}(i)
	}

	// Call the function to be tested
	doneCounting := countAllKmers(numMaps, filterMask, kmerSequences, allCounts, minKmerCount)

	// Wait for all counting to be done
	for i := 0; i < numMaps; i++ {
		<-doneCounting
	}

	// Validate results
	for i := 0; i < numMaps; i++ {
		if !reflect.DeepEqual(expectedCounts[i], allCounts[i]) {
			t.Errorf("Expected counts: %v, Got: %v", expectedCounts[i], allCounts[i])
		}
	}
}

func TestWriteCountsToFile(t *testing.T) {
	// Setup
	minKmerCount := 2
	allCounts := []map[uint64]int32{
		{0: 2, 1: 1},
		{2: 3, 3: 4},
	}
	tmpFile, err := os.CreateTemp("", "test_write_counts_*.bin")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpFile.Name())

	// Call function under test
	writeCountsToFile(allCounts, tmpFile.Name(), minKmerCount)

	// Validate results
	content, err := os.ReadFile(tmpFile.Name())
	if err != nil {
		t.Fatal(err)
	}
	expected := map[uint64]struct{}{
		0: {},
		2: {},
		3: {},
	}
	reader := bytes.NewReader(content)
	var kmer uint64
	readKmers := map[uint64]struct{}{}
	for reader.Len() > 0 {
		if err := binary.Read(reader, binary.LittleEndian, &kmer); err != nil {
			t.Fatal(err)
		}
		readKmers[kmer] = struct{}{}
	}
	if !reflect.DeepEqual(expected, readKmers) {
		t.Errorf("Expected kmers %v, got %v", expected, readKmers)
	}
}
