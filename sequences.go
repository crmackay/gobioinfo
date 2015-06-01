package gobioinfo

import (
//"fmt"
)

type NucleotideSequence struct {
	Sequence []byte
}

type QSequence struct {
	QualByte []byte
	PHRED    []uint8
	Encoding string
}
