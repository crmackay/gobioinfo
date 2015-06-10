package gobioinfo

import (
	"fmt"
	//"os/user"
	"testing"
)

// NewFASTQRead(ln1 string, ln2 []rune, ln3 string, ln4 []rune) (newRead FASTQRead)
func TestNewFASTQRead(t *testing.T) {
	fmt.Println("testing NewFASTQRead")
}

//func DecodePHRED(encoded []rune, encoding string) (decoded []uint8)
func TestDecodePHRED(t *testing.T) {
	fmt.Println("testing DecodePHRED")
}

// func (p *PHRED) Decode()
func TestDecode(t *testing.T) {
	fmt.Println("testing *PHRED.Decode()")
}

/*
func TestDecodeQualByteSequence(t *testing.T) {

	testSequence := []rune("@@@FFFFFHHFFFFFHGHJ@FH?")

	result := DecodePHRED(testSequence, "illumina_1.8")
	//                         @    @   @   F   F   F   F   F   H   H   F   F   F   F   F   H   G   H   J   @   F   H   ?
	predictedResult := []uint8{31, 31, 31, 37, 37, 37, 37, 37, 39, 39, 37, 37, 37, 37, 37, 39, 38, 39, 41, 31, 37, 39, 30}

	for i, value := range result {
		if value != predictedResult[i] {
			t.Error("PHRED decoding not as predicted: got", value, " expeceted ", predictedResult[i])
		}
	}
}
*/
