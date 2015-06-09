package gobioinfo

import (
//"fmt"
)

// NucleotideSequence is
type NucleotideSequence []rune

// DNASequence is
type DNASequence struct {
	Sequence NucleotideSequence
}

// PHRED is
type PHRED struct {
	Encoded  []rune
	Decoded  []uint8
	Encoding string
}

// Decode takes the
func (p *PHRED) Decode() {
	decodedPHRED := make([]uint8, len(p.Encoded))

	for i, rune := range p.Encoded {
		decodedPHRED[i] = PHREDEncodings[p.Encoding][string(rune)]
	}

	p.Decoded = decodedPHRED
}

// DecodePHRED takes a
func DecodePHRED(encoded []rune, encoding string) (decoded []uint8) {
	decoded = make([]uint8, len(encoded))
	for i, rune := range encoded {
		decoded[i] = PHREDEncodings[encoding][string(rune)]
	}
	return decoded
}

// FASTQRead is
type FASTQRead struct {
	DNASequence
	PHRED
	ID   string
	Misc string
}

// FASTARead is
type FASTARead struct {
	ID string
	DNASequence
}

// PairWiseAlignment creates an pairwise alignment structure
type PairWiseAlignment struct {
	Subject                 NucleotideSequence
	Query                   NucleotideSequence
	ExpandedCIGAR           []string
	SubjectStart            int
	QueryStart              int
	SubjectAlignLen         int
	QueryAlignLen           int
	GappedSubject           string
	GappedQuery             string
	AlignmentRepresentation string
}

// PHREDEncodings
var PHREDEncodings = map[string]map[string]uint8{
	"illumina_1.8": {
		"!":  0,
		"\"": 1,
		"#":  2,
		"$":  3,
		"%":  4,
		"&":  5,
		"'":  6,
		"(":  7,
		")":  8,
		"*":  9,
		"+":  10,
		",":  11,
		"-":  12,
		".":  13,
		"/":  14,
		"0":  15,
		"1":  16,
		"2":  17,
		"3":  18,
		"4":  19,
		"5":  20,
		"6":  21,
		"7":  22,
		"8":  23,
		"9":  24,
		":":  25,
		";":  26,
		"<":  27,
		"=":  28,
		">":  29,
		"?":  30,
		"@":  31,
		"A":  32,
		"B":  33,
		"C":  34,
		"D":  35,
		"E":  36,
		"F":  37,
		"G":  38,
		"H":  39,
		"I":  40,
		"J":  41,
	},
}
