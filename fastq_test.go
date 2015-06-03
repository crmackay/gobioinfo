package gobioinfo

import (
	"fmt"
	"os/user"
	"testing"
)

func TestFASTQReader(t *testing.T) {

	// read FASTQ file

	// test FASTQ read values

	// t.Error("this is a test error")

	usr, _ := user.Current()
	homeDir := usr.HomeDir
	testPath := homeDir + "/Desktop/coding/golang/src/github.com/crmackay/gobioinfo/testData/sample_50.fastq"

	fastqscanner := NewFASTQScanner(testPath)

	defer fastqscanner.Close()

	for {
		myRead := fastqscanner.NextRead()
		//fmt.Println(myRead.Id)
		//fmt.Println(myRead.Sequence)
		if myRead.Id == "" {
			break
		}
	}
}

func TestFASTQWriter(t *testing.T) {

	// create FASTQ objects

	// write reads to a file

	// is file the same?

	// delete file

	var err error

	usr, _ := user.Current()

	homeDir := usr.HomeDir

	testPath := homeDir + "/Desktop/coding/golang/src/testing/test.fastq"

	newRead := FASTQRead{
		Id: "this is my read name",
		DNASequence: DNASequence{
			Sequence: NucleotideSequence([]rune("AATCGATCGATGAGATAGTC")),
		},
		Misc:      "+",
		QSequence: QSequence{QualByteSequence: []rune("*(&*^%^%$^%#Q(*&0(&(*&^&^%^%$"), PHRED: []uint8{10, 20, 18, 10, 16, 25, 35, 35, 40, 35, 35, 36, 27, 32, 34, 23, 34, 23, 23, 34, 4, 5, 45, 45, 5, 45, 45, 5, 45}},
	}

	w := NewFASTQWriter(testPath)

	err = w.Write(newRead)
	err = w.Write(newRead)

	if err != nil {
		fmt.Println(err)
	}

	w.Close()
}

func TestFASTAReader(t *testing.T) {
}

func TestFASTAWriter(t *testing.T) {
}

func TestDecodeQualByteSequence(t *testing.T) {

	testSequence := []rune("@@@FFFFFHHFFFFFHGHJ@FH?")

	result := DecodeQualByteSequence(testSequence, "illumina 1.8")
	//                         @    @   @   F   F   F   F   F   H   H   F   F   F   F   F   H   G   H   J   @   F   H   ?
	predictedResult := []uint8{31, 31, 31, 37, 37, 37, 37, 37, 39, 39, 37, 37, 37, 37, 37, 39, 38, 39, 41, 31, 37, 39, 30}

	for i, value := range result {
		if value != predictedResult[i] {
			t.Error("PHRED decoding not as predicted: got", value, " expeceted ", predictedResult[i])
		}
	}
}
