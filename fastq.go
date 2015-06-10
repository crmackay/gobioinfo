package gobioinfo

import (
//	"fmt"
)

// FASTQRead is
type FASTQRead struct {
	DNASequence
	PHRED
	ID   string
	Misc string
}

// PHRED is
type PHRED struct {
	Encoded  []rune
	Decoded  []uint8
	Encoding string
}

// NewFASTQRead ...
func NewFASTQRead(ln1 string, ln2 []rune, ln3 string, ln4 []rune) (newRead FASTQRead) {

	phredEncoding := "illumina_1.8"

	decodedQuality := DecodePHRED(ln4, phredEncoding)

	newSequence := NucleotideSequence(ln2)

	newRead = FASTQRead{
		ID: ln1,
		DNASequence: DNASequence{
			Sequence: newSequence,
		},
		Misc: ln3,
		PHRED: PHRED{
			Encoded:  ln4,
			Decoded:  decodedQuality,
			Encoding: phredEncoding,
		},
	}

	return (newRead)
}

// DecodePHRED takes a
func DecodePHRED(encoded []rune, encoding string) (decoded []uint8) {
	decoded = make([]uint8, len(encoded))
	for i, rune := range encoded {
		decoded[i] = PHREDEncodings[encoding][string(rune)]
	}
	return decoded
}

// Decode turns the Encoded part of a PHRED struct and in-place decodes it
// and stores in the the Decoded element of the PHRED
func (p *PHRED) Decode() {
	decodedPHRED := make([]uint8, len(p.Encoded))

	for i, rune := range p.Encoded {
		decodedPHRED[i] = PHREDEncodings[p.Encoding][string(rune)]
	}

	p.Decoded = decodedPHRED
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
