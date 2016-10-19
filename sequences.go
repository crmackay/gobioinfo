package gobioinfo

/*
TODO:optimizations -
convert from []rune to []base where base is uint8 and holds the utf8codepoints
codepoints for each sequence...this saves 24 bits (3 bytes) per base pair, and
likely speeds up comparisons and saves space in memory <-- is this actually better?


// NucSeq is a wrapper around a slice of runes repreceting
// a nucleotide sequence

type NucBase uint8

type NucSeq []NucBase

*/

// NucSeq is a wrapper around a slice of runes representing
// a nucleotide sequence
type NucSeq []uint8

// NewNucSeq takes a string and returns a NucSeq
func NewNucSeq(s string) NucSeq {
	seq := NucSeq{}
	for i, base := range s {
		seq[i] = uint8(base)
	}

	return seq
}

// DNASequence is a struct representing a dna sequence, it has a sequence attribute
// and can have more attribites later, like species, source etc.
type DNASequence struct {
	Sequence NucSeq
}

// FASTARead is
type FASTARead struct {
	ID string
	DNASequence
}

// NewDNASequence ...
func NewDNASequence(s string) DNASequence {

	d := DNASequence{Sequence: []uint8(s)}

	return d

}
