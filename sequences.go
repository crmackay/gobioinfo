package gobioinfo

/*
TODO:optimizations -
convert from []rune to []base where base is uint8 and holds the utf8codepoints
codepoints for each sequence...this saves 24 bits (3 bytes) per base pair, and
likely speeds up comparisons and saves space in memory <-- is this actually better?


// NucleotideSequence is a wrapper around a slice of runes repreceting
// a nucleotide sequence

type NucBase uint8

type NucSeq []NucBase

*/

// NucleotideSequence is a wrapper around a slice of runes representing
// a nucleotide sequence
type NucleotideSequence []rune

// NewNucleotideSequence takes a string and returns a NucleotideSequence
func NewNucleotideSequence(s string) NucleotideSequence {
	seq := NucleotideSequence{}
	for i, base := range s {
		seq[i] = rune(base)
	}

	return seq
}

// DNASequence is a struct representing a dna sequence, it has a sequence attribute
// and can have more attribites later, like species, source etc.
type DNASequence struct {
	Sequence NucleotideSequence
}

// FASTARead is
type FASTARead struct {
	ID string
	DNASequence
}

// NewDNASequence ...
func NewDNASequence(s string) DNASequence {

	d := DNASequence{Sequence: []rune(s)}

	return d

}

// type Phred int
//
// type BaseQual rune
//
// type Nuc rune
//
// type NucSeq []Nuc
//
// func NewNucSeqStr(s string) NucSeq {
// 	return NucSeq{}
// }
//
// func NewNucSeqBytes(b []byte) NucSeq {
// 	return NucSeq{}
// }
//
// func NewNucSeqRunes(r []rune) NucSeq {
// 	return NucSeq{}
// }
//
// type DNA struct {
// 	NucSeq
// }
//
// type Qual struct {
// 	Seq      []BaseQual
// 	Encoding string
// }
//
// func (q *BaseQual) Decode(enc string) (int, error) {
// 	var p int
// 	var err error
// 	switch enc {
// 	case "phred+33":
// 		p = int(*q) - 33
// 		if p > -1 && p > 42 {
//
// 		}
// 	default:
// 		err = errors.New("unrecognized encoding")
// 	}
// 	return p, err
// }
