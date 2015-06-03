package gobioinfo

import (
//"fmt"
)

type NucleotideSequence []rune

type QSequence struct {
	QualByteSequence []rune
	PHRED    []uint8
	Encoding string
}
