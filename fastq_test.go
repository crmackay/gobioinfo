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
		fmt.Println(myRead.Id)
		fmt.Println(myRead.Sequence.Sequence)
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
		Id:       "this is my read name",
		Sequence: NucleotideSequence{Sequence: []byte("AATCGATCGATGAGATAGTC")},
		Misc:     "+",
		Quality:  QSequence{QualByte: []byte("*(&*^%^%$^%#Q(*&0(&(*&^&^%^%$"), PHRED: []uint8{10, 20, 18, 10, 16, 25, 35, 35, 40, 35, 35, 36, 27, 32, 34, 23, 34, 23, 23, 34, 4, 5, 45, 45, 5, 45, 45, 5, 45}},
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
