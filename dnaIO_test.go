package gobioinfo

import (
	"fmt"
	"os/user"
	"testing"
)

// func NewFASTQScanner(filePath string) FASTQScanner
func TestNewFASTQScanner(t *testing.T) {

}

//func (s *FASTQScanner) NextRead() FASTQRead {}
func TestFASTQScannerNextRead(t *testing.T) {

}

// func NewFASTQWriter(filePath string) FASTQWriter {}

// func (w *FASTQWriter) Write(r FASTQRead) error {}

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
		if myRead.ID == "" {
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
		ID: "this is my read name",
		DNASequence: DNASequence{
			Sequence: NucleotideSequence([]rune("AATCGATCGATGAGATAGTC")),
		},
		Misc:  "+",
		PHRED: PHRED{Encoded: []rune("*(&*^%^%$^%#Q(*&0(&(*&^&^%^%$"), Decoded: []uint8{10, 20, 18, 10, 16, 25, 35, 35, 40, 35, 35, 36, 27, 32, 34, 23, 34, 23, 23, 34, 4, 5, 45, 45, 5, 45, 45, 5, 45}},
	}

	w := NewFASTQWriter(testPath)

	err = w.Write(newRead)
	err = w.Write(newRead)

	if err != nil {
		fmt.Println(err)
	}

	w.Close()
}

// func NewFASTAWriter(filePath string) FASTAWriter {}
// func (w *FASTAWriter) Write(r FASTQRead) error {}
// func (w *FASTAWriter) Close() {}
