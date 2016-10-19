package gobioinfo

//	"fmt"

// FASTQRead is structure holding all of the elements of a FASTQ read, which includes a sequences,
// a quality (PHRED) sequence, an ID sequence, and a miscallaneous string (line 3 of each FASTQRead).
type FASTQRead struct {
	DNASequence
	PHRED
	ID   string
	Misc string
}

// TODO: refactoring -
// type Phred int
// type Qual struct {
// 		Encoded  []Phred
//		Encoding string
// }

// PHRED is
type PHRED struct {
	Encoded []rune
	Offset  int
}

// NewFASTQRead ...
func NewFASTQRead(ln1 string, ln2 string, ln3 string, ln4 string) (newRead FASTQRead) {

	// "illumina_1.8", this is the default until we make this user configurable
	phredOffset := 33

	newSequence := NucSeq(ln2)

	newRead = FASTQRead{
		ID: ln1,
		DNASequence: DNASequence{
			Sequence: newSequence,
		},
		Misc: ln3,
		PHRED: PHRED{
			Encoded: []rune(ln4),
			Offset:  phredOffset,
		},
	}

	return (newRead)
}

// DecodePHRED takes a
func DecodePHRED(encoded []rune, phredOffset int) (decoded []uint8) {
	decoded = make([]uint8, len(encoded))
	for i, rune := range encoded {
		decoded[i] = PHREDSymToNum(rune, phredOffset)
	}
	return decoded
}

// PHREDSymToNums is is a simple function that takes the symbolic representation of a PHRED score (ala
// an encoded PHRED score) and its encoding offset, converts the symbol to its integer value, and
// subtracts the offset, arriving at the decoded PHRED score. The default offset if 33, which is used
// by Illumina in version 1.8+
func PHREDSymToNum(value rune, offset int) uint8 {
	return uint8(int(value) - offset)
}

func PHREDNumToSym(value uint8, offset int) rune {
	return rune(int(value) + offset)
}
