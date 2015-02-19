package gobioinfo

import (
	"bufio"
	"fmt"
	"os"
)

type FASTQRead struct {
	Sequence string
	Quality  string
	Name     string
	Misc     string
}

type FASTQScanner struct {
	*os.File
	*bufio.Scanner
}

func NewFASTQScanner(filePath string) FASTQScanner {

	file, err := os.Open(filePath)

	if err != nil {
		fmt.Println("error opening file= ", err)
		os.Exit(1)
	}
	fastqscanner := FASTQScanner{Scanner: bufio.NewScanner(file), File: file}
	return (fastqscanner)
}

func (s *FASTQScanner) Close() {
	s.File.Close()
	fmt.Println("closing fastqscanner")
}

func (s *FASTQScanner) NextRead() FASTQRead {

	var ln1, ln2, ln3, ln4 string
	var newRead FASTQRead

	if s.Scan() {
		ln1 = s.Text()
		fmt.Print(ln1)
	} else {
		return FASTQRead{}
	}
	if s.Scan() {
		ln2 = s.Text()
	} else {
		return FASTQRead{}
	}
	if s.Scan() {
		ln3 = s.Text()
	} else {
		return FASTQRead{}
	}
	if s.Scan() {
		ln4 = s.Text()
	} else {
		return FASTQRead{}
	}

	newRead = FASTQRead{Name: ln1, Sequence: ln2, Misc: ln3, Quality: ln4}

	return newRead

}

func Testfastq() {

	testPath := "/Users/christophermackay/Desktop/deepseq_data/pir1/hits-clip/working_data/sample_data/fastq/sample_50_2.fastq"

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

//TODO unit testing?? error handling???
