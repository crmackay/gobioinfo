package gobioinfo

import (
	"fmt"
	"testing"
)

func TestSGAlign(t *testing.T) {

	// create FASTQ reads

	queryRead := NewFASTQRead(
		"test Query",
		[]rune("GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAA"),
		"",
		[]rune("@@@FFFFFHHFFFFFHGHJ@FH?BFHF<HIGGIJIGJJGG=CCGGGHIC@=DDECHHED3>@CDCDCACC>>@A:9>99@)<>?@>@5)8<@CC:A>A<A"),
	)
	subjectRead := FASTARead{
		DNASequence: DNASequence{
			Sequence: NucleotideSequence([]rune("GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTGCTTG")),
		},
		ID: "test Subject",
	}

	// create predicted alignments

	// GTGTCAGTCACTTCCAGCGGTCGTATGCCGTC-T-TGCTTG
	// |||| ||||||||||||||||||||||||||| | ||||||
	// GTGT-AGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG
	//

	// test if the alignments match up

	//subject := "GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTGCTTG"
	//query := "GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTGTTCTACTTGAA"

	result := queryRead.Sequence.Align(subjectRead.Sequence)

	if result.SubjectStart != 0 {
		t.Error("SubjectStart should be 0, got ", result.SubjectStart)
	}
	if result.QueryStart != 58 {
		t.Error("SubjectStart should be 57, got ", result.QueryStart)
	}
	if result.GappedSubject != "GTGTCAGTCACTTCCAGCGGTCGTATGCCGTC-T-TGCTTG" {
		t.Error("Gapped Subject not matching predicted string")
	}
	if result.GappedQuery != "GTGT-AGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG" {
		t.Error("Gapped Query is not matching predicted string")
	}
	if result.AlignmentRepresentation != "|||| ||||||||||||||||||||||||||| | ||||||" {
		t.Error("Alignment Representation is not matching predicted value")
	}

	fmt.Println("results:")
	fmt.Println("Subject Start\t", result.SubjectStart)
	fmt.Println("Query Start\t", result.QueryStart)
	fmt.Println("SubjectAlignLen\t", result.SubjectAlignLen)
	fmt.Println("QueryAlignLen\t", result.QueryAlignLen)
	fmt.Println("Subject\t", string(result.Subject))
	fmt.Println("GappedSubject\t", result.GappedSubject)
	fmt.Println("AlignRepr\t", result.AlignmentRepresentation)
	fmt.Println("GappedQuery\t", result.GappedQuery)
	fmt.Println("Query\t", string(result.Query))
}

// func alignmentRepr(alignment PairWiseAlignment) PairWiseAlignment {}
