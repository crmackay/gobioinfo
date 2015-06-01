package gobioinfo

import (
//"fmt"
)

type NucleotideSequence struct {
	Sequence []rune
}

type QSequence struct {
	QualByte []rune
	PHRED    []uint8
	Encoding string
}
