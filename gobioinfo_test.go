package gobioinfo

import (
	"fmt"
	"os/user"
	"testing"
)

func TestFASTQReader(t *testing.T) {

	usr, _ := user.Current()
	homeDir := usr.HomeDir
	testPath := homeDir + "/Desktop/coding/golang/src/testing/sample_50_2.fastq"

	fastqscanner := NewFASTQScanner(testPath)

	defer fastqscanner.Close()

	for {
		myRead := fastqscanner.NextRead()
		fmt.Print(myRead.Name)
		if myRead.Name == "" {
			break
		}
	}
}


func TestFASTQWriter(t *testing.T) {
	var err error
	usr, _ := user.Current()
	homeDir := usr.HomeDir
	testPath := homeDir + "/Desktop/coding/golang/src/testing/test.fastq"
	newRead := FASTQRead{
		Name:     "this is my read name",
		Sequence: "AATCGATCGATGAGATAGTC",
		Misc:     "+",
		QualityString:  "*(&*^%^%$^%#Q(*&0(&(*&^&^%^%$",
		QualityPHRED: []uint8{10, 20, 18, 10, 16, 25, 35, 35, 40 ,35, 35, 36, 27, 32, 34, 23, 34, 23, 23, 34, 4, 5,45, 45, 5, 45, 45, 5, 45},
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

func TestSGAlign(t *testing.T) {

	subject := "GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTGCTTG"
	query := "GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTGTTCTACTTGAA"

	result := Align(subject, query)
	fmt.Println("results:")
	fmt.Println(result.SubjectStart)
	fmt.Println(result.QueryStart)
	fmt.Println(result.SubjectAlignLen)
	fmt.Println(result.QueryAlignLen)
	fmt.Println(result.Subject)
	fmt.Println(result.GappedSubject)
	fmt.Println(result.AlignmentRepresentation)
	fmt.Println(result.GappedQuery)
	fmt.Println(result.Query)
}
