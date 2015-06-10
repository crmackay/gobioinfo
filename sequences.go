package gobioinfo

import (
//"fmt"
)

// NucleotideSequence is a wrapper around a slice of runes repreceting
// a nucleotide sequence
type NucleotideSequence []rune

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
