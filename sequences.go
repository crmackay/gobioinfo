package gobioinfo

import (
//"fmt"
)

type NucleotideSequence []rune

type QSequence struct {
	QualByte []rune
	PHRED    []uint8
	Encoding string
}
