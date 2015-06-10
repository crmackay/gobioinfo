package gobioinfo

import (
	"fmt"
	"testing"
)

// func NewDNASequence(s string) DNASequence {}
func TestNewDNASequence(t *testing.T) {
	fmt.Println("testing NewDNASequence()")
	type testPair struct {
		input  string
		output []rune
	}

	testSuite := []testPair{
		{"AT", []rune{65, 84}},
	}

	for _, elem := range testSuite {
		newDNA := NewDNASequence(elem.input)
		for i, base := range newDNA.Sequence {
			if base != elem.output[i] {
				t.Error("NewDNASequence is not working, got ", elem.output, "but was expecting ", NewDNASequence(elem.input))
			}
		}
	}
}
