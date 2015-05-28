package gobioinfo

import (
//"fmt"
)

type nucleotideSequence struct {
	Sequence []byte
	Name     string
	Misc     string
}

type DNASequence struct {
	nucleotideSequence
}

type RNASequence struct {
	nucleotideSequence
}

type QSequence struct {
	Quality  []byte
	PHRED    []byte
	Encoding string
}
