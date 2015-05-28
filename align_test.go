package gobioinfo

import (
	"fmt"
	"testing"
)

func TestSGAlign(t *testing.T) {

	// create FASTA reads

	// create predicted alignments

	// test if the alignments match up

	subject := "GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTGCTTG"
	query := "GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTGTTCTACTTGAA"

	result := Align(subject, query)
	fmt.Println("results:")
	fmt.Println(result.SubjectStart)
	fmt.Println(result.QueryStart)
	fmt.Println(result.SubjectAlignLen)
	fmt.Println(result.QueryAlignLen)
	fmt.Println(result.Subject)
	fmt.Println(result.GappedSubject)
	fmt.Println(result.AlignmentRepresentation)
	fmt.Println(result.GappedQuery)
	fmt.Println(result.Query)
}
