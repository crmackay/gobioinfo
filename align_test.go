package gobioinfo

import (
	"bufio"
	"bytes"
	"fmt"
	"strings"
	"testing"
)

func TestSGAlign(t *testing.T) {

	fmt.Println("testing SGAlign...")

	numRead := 6

	rawTestData := bytes.NewBufferString("@HWI-ST560:155:C574EACXX:3:1101:2214:1998 1:N:0:\nAGCAGGGAGGACGATGCGGAACTGATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCGTGTCAGTCACTTCCAGCGGGCGTATGCCGT\n+\nCCCFFFFFHHHHHJJJJJJJJJJJJJJJJIIGHEIIJIGIGGIGFFCCDDECCDDDDDDDDDDDDCDCC-9?28A@>4@AAC3>4439B###########\n@HWI-ST560:155:C574EACXX:3:1101:2300:1939 1:N:0:\nCACAGGGAGGACGATGCGGGAGTGAGACCGTCTTGCTTACTTGTCCGATGAAATGAATGAAATAGAAAGTGGGAAAATAATGTGTCAGTCACTTCCAGCG\n+\nBBCFFFFFHHHHHJJJJJJJIJHIGIHHIJJJJJIGIIEHIIGHHHHFFFCEDDEDCCCDDDACDCDDD@CDDACBCCCDDCDEDDCDCCCCCDCC:@99\n@HWI-ST560:155:C574EACXX:3:1101:2409:1942 1:N:0:\nCACAGGGAGGACGATGCGGAAAAGAATGTGAATCATGGTGTTCTTGTGGTTGGCTATGGGACTCTTGATGGCAAAGATTACTGGCTTGTGAAAAAAGGGT\n+\nBBBFFFFFHHHHHFIJIJJJJJJJJIHIFIJJJJJJJJ=FGIJJJJHIJJJJHHHHFFFFFDEEEEDDDDDDDDDDDDDDDDDDDDD5?CCDCDD#####\n@HWI-ST560:155:C574EACXX:3:1101:2433:1960 1:N:0:\nAGCAGGGAGGACGATGCGGACAAGTCCCTGAGGAGCCCTTTGAGCCTGGTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAA\n+\nBCCFFFFFHHHHHJJJJJJJIJJJGJIJJJIJJJJJJJJJJGJIJIJHHHGFFFFFFEEEEEEEDDDBDDDDDDDCDDDDDDDDDDDDDDDDDDD@B@DB\n@HWI-ST560:155:C574EACXX:3:1101:2381:1976 1:N:0:\nAGCAGGGAGGACGATGCGGTGATGTTCACAGTGGCTAAGTTCCGCGGTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAA\n+\nBCCFFFFFHHHHHIJJJJJJJJJIIHIIJJJFHJJJJJJJJJJJJJJHHFFFFFFEEEEEEEDDDDDDDDDDCDDBBBDDDDDDDDCDDDDDDD9BB>BD\n@HWI-ST560:155:C574EACXX:3:1101:2403:1977 1:N:0:\nGCTAGGGAGGACGATGCGGCTAAGTGGTTGGAACCCGATTGCCTCTCTGGAGCGTGTCAGTCACTTCCAGCGGGTGTCAGTCACTTCCAGCGGTCGTATG\n+\n@@@FFFFFHHGHHJJJJJJGIEFHFHGDHGIEGGHIIJIICHHIJHEFHGDDDCDD@BCCDDDDDDA@CDDDDD@><ACDDCCCCCC>CC?>B9@B>833\n")

	//rawTestData := bytes.NewBufferString("@HWI-ST560:155:C574EACXX:3:1101:2403:1977 1:N:0:\nGCTAGGGAGGACGATGCGGCTAAGTGGTTGGAACCCGATTGCCTCTCTGGAGCGTGTCAGTCACTTCCAGCGGGTGTCAGTCACTTCCAGCGGTCGTATG\n+\n@@@FFFFFHHGHHJJJJJJGIEFHFHGDHGIEGGHIIJIICHHIJHEFHGDDDCDD@BCCDDDDDDA@CDDDDD@><ACDDCCCCCC>CC?>B9@B>833\n")
	//																									                                                             GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG
	scanner := FASTQScanner{Scanner: bufio.NewScanner(rawTestData), File: nil}

	testReads := make([]FASTQRead, numRead)

	for i := 0; i < numRead; i++ {
		newRead := scanner.NextRead()
		//fmt.Println(newRead.DNASequence.Sequence)
		testReads[i] = newRead
	}

	subject := NewDNASequence("GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG")

	type testGroup struct {
		alignment          PairWiseAlignment
		read               FASTQRead
		expectedCIGAR      []string
		expectedQueryStart int
	}

	testSuite := make([]testGroup, numRead)

	rawCIGARS := []string{
		"mmmmmmmmmmmmmmmmmmmmxmmmmmmmmmm",
		"mmmmmmmmmmmmmmmmmmm",
		"mxmm",
		"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm",
		"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm",
		"mmmmmmmmmmmmmmmmmmmmmmmmmmm",
	}

	expectedCIGARS := make([][]string, numRead)

	for i, s := range rawCIGARS {
		expectedCIGARS[i] = strings.Split(s, "")
	}

	for i, elem := range testReads {
		testSuite[i] = testGroup{
			alignment:          elem.Sequence.Align(subject.Sequence),
			read:               elem,
			expectedCIGAR:      expectedCIGARS[i],
			expectedQueryStart: 69,
		}
	}

	for _, elem := range testSuite {
		if strings.Join(elem.expectedCIGAR, "") != strings.Join(elem.alignment.ExpandedCIGAR, "") {
			t.Error(
				"expected: ", elem.expectedCIGAR, "\n",
				"but got: ", elem.alignment.ExpandedCIGAR, "\n",
			)
		}
	}

	// create FASTQ reads
	/*
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
	*/
}

// func alignmentRepr(alignment PairWiseAlignment) PairWiseAlignment {}
