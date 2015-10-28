package gobioinfo

import (
	"bytes"
	"fmt"
	"testing"
)

func TestSGAlign(t *testing.T) {

	fmt.Println("testing SGAlign...")

	rawTestData := bytes.NewBufferString(`@HWI-ST560:155:C574EACXX:3:1101:2214:1998 1:N:0:
AGCAGGGAGGACGATGCGGAACTGATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCGTGTCAGTCACTTCCAGCGGGCGTATGCCGT
+
CCCFFFFFHHHHHJJJJJJJJJJJJJJJJIIGHEIIJIGIGGIGFFCCDDECCDDDDDDDDDDDDCDCC-9?28A@>4@AAC3>4439B###########
@HWI-ST560:155:C574EACXX:3:1101:2300:1939 1:N:0:
CACAGGGAGGACGATGCGGGAGTGAGACCGTCTTGCTTACTTGTCCGATGAAATGAATGAAATAGAAAGTGGGAAAATAATGTGTCAGTCACTTCCAGCG
+
BBCFFFFFHHHHHJJJJJJJIJHIGIHHIJJJJJIGIIEHIIGHHHHFFFCEDDEDCCCDDDACDCDDD@CDDACBCCCDDCDEDDCDCCCCCDCC:@99
@HWI-ST560:155:C574EACXX:3:1101:2409:1942 1:N:0:
CACAGGGAGGACGATGCGGAAAAGAATGTGAATCATGGTGTTCTTGTGGTTGGCTATGGGACTCTTGATGGCAAAGATTACTGGCTTGTGAAAAAAGGGT
+
BBBFFFFFHHHHHFIJIJJJJJJJJIHIFIJJJJJJJJ=FGIJJJJHIJJJJHHHHFFFFFDEEEEDDDDDDDDDDDDDDDDDDDDD5?CCDCDD#####
@HWI-ST560:155:C574EACXX:3:1101:2433:1960 1:N:0:
AGCAGGGAGGACGATGCGGACAAGTCCCTGAGGAGCCCTTTGAGCCTGGTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAA
+
BCCFFFFFHHHHHJJJJJJJIJJJGJIJJJIJJJJJJJJJJGJIJIJHHHGFFFFFFEEEEEEEDDDBDDDDDDDCDDDDDDDDDDDDDDDDDDD@B@DB
@HWI-ST560:155:C574EACXX:3:1101:2381:1976 1:N:0:
AGCAGGGAGGACGATGCGGTGATGTTCACAGTGGCTAAGTTCCGCGGTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAA
+
BCCFFFFFHHHHHIJJJJJJJJJIIHIIJJJFHJJJJJJJJJJJJJJHHFFFFFFEEEEEEEDDDDDDDDDDCDDBBBDDDDDDDDCDDDDDDD9BB>BD
@HWI-ST560:155:C574EACXX:3:1101:2403:1977 1:N:0:
GCTAGGGAGGACGATGCGGCTAAGTGGTTGGAACCCGATTGCCTCTCTGGAGCGTGTCAGTCACTTCCAGCGGGTGTCAGTCACTTCCAGCGGTCGTATG
+
@@@FFFFFHHGHHJJJJJJGIEFHFHGDHGIEGGHIIJIICHHIJHEFHGDDDCDD@BCCDDDDDDA@CDDDDD@><ACDDCCCCCC>CC?>B9@B>833
@HWI-ST560:155:C574EACXX:3:1101:4452:1943 1:N:0:
CACAGGGAGGACGATGCGGAGGAGAAGACCACATATGTGAAGGCCCCTGGTTGACTGGTTGTGGGCTCAGCTGACCAGCTGGGCTTGCCTGCTGCAGGCG
+
BBBFFFFFHHHHHJJJJJJJIJGGIJFBGH@GGHIGGDHIICHIJJHHHGFFFCDEEDCED@@@DDDDDDCCACDD@BBCCDDBCDDDDDDDDCA99<9@
@another test
AGCAGGGAGGACGATGCGGTTGTGATATAATACAACCTGCTAA
+
B@CDFFFFHHDDFHHIJJJIIIDGHBFGFHIGHGHG@EGHJHH
`)

	//rawTestData := bytes.NewBufferString("@HWI-ST560:155:C574EACXX:3:1101:2403:1977 1:N:0:\nGCTAGGGAGGACGATGCGGCTAAGTGGTTGGAACCCGATTGCCTCTCTGGAGCGTGTCAGTCACTTCCAGCGGGTGTCAGTCACTTCCAGCGGTCGTATG\n+\n@@@FFFFFHHGHHJJJJJJGIEFHFHGDHGIEGGHIIJIICHHIJHEFHGDDDCDD@BCCDDDDDDA@CDDDDD@><ACDDCCCCCC>CC?>B9@B>833\n")
	//																									                                                             GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG
	scanner := NewFASTQScanner(rawTestData)

	var testReads []FASTQRead
	for {
		newRead, err := scanner.NextRead()
		if err != nil {
			if err.Error() == "EOF" {
				break
			}
		}
		//	fmt.Println(newRead.DNASequence.Sequence)
		testReads = append(testReads, newRead)
	}
	subject := NewDNASequence("GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG")

	type testGroup struct {
		alignment          PairWiseAlignment
		read               FASTQRead
		expectedCIGAR      string
		expectedQueryStart int
	}

	var testSuite []testGroup

	expectedCIGARS := []string{
		"mmmmmmmmmmmmmmmmmmmmxmmmmmmmmmm",
		"mmmmmmmmmmmmmmmmmmm",
		"mm",
		"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm",
		"mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm",
		"mmmmmmmmmmmmmmmmmmmmmmmmmmm",
		"immmxxmjmmmjmmjxmmmmmjmmxmimimmimiimmxmmmmxxm",
		"",
	}

	for i, elem := range testReads {
		testSuite = append(testSuite, testGroup{
			alignment:          elem.Sequence.SG3pAlign(subject.Sequence),
			read:               elem,
			expectedCIGAR:      expectedCIGARS[i],
			expectedQueryStart: 69,
		})
	}

	for _, elem := range testSuite {
		fmt.Print("HERE")
		if elem.expectedCIGAR != elem.alignment.ExpandedCIGAR {
			t.Error(
				"expected: \t",
				elem.expectedCIGAR,
				"\nbut got: \t\t",
				elem.alignment.ExpandedCIGAR,
				"\n",
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
