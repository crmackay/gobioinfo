package gobioinfo

import (
	"bufio"
	"fmt"
	"os"
	"os/user"
	"strings"
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

	if s.Scanner.Scan() {
		ln1 = s.Scanner.Text()
		fmt.Print(ln1)
	} else {
		return FASTQRead{}
	}
	if s.Scanner.Scan() {
		ln2 = s.Scanner.Text()
	} else {
		return FASTQRead{}
	}
	if s.Scanner.Scan() {
		ln3 = s.Scanner.Text()
	} else {
		return FASTQRead{}
	}
	if s.Scanner.Scan() {
		ln4 = s.Scanner.Text()
	} else {
		return FASTQRead{}
	}

	newRead = FASTQRead{Name: ln1, Sequence: ln2, Misc: ln3, Quality: ln4}

	return newRead

}

type FASTQWriter struct {
	*os.File
	*bufio.Writer
}

func NewFASTQWriter(filePath string) FASTQWriter {

	file, err := os.Create(filePath)

	if err != nil {
		fmt.Println("error opening file= ", err)
		os.Exit(1)
	}
	fastqwriter := FASTQWriter{Writer: bufio.NewWriter(file), File: file}
	return (fastqwriter)
}

func (w *FASTQWriter) Write(r FASTQRead) error {
	//compose FASTQRead struct into the proper format

	for_writing := strings.Join([]string{"@" + r.Name, r.Sequence, r.Misc, r.Quality, ""}, "\n")
	fmt.Println(for_writing)
	//check if there are enough bytes lefts in the writer buffer
	var err error
	var i int
	if w.Writer.Available() < len(for_writing) {
		w.Writer.Flush()

		//flush the write buffer, then write
		_, err = w.Writer.WriteString(for_writing)
	} else {
		fmt.Println("here")
		//write to the write buffer
		i, err = w.Writer.WriteString(for_writing)
		fmt.Println("return from writestring: ", i)
		fmt.Println("here")
	}
	return (err)
}

func (w *FASTQWriter) Close() {
	w.Writer.Flush()
	w.File.Close()
	fmt.Println("closing the writer")
}

/*

NOTE...there is a lot of repetitive code here (between FASTA and FASTQ).
TODO: look into a better way of doing this a bit...they both have the same
functions, they only differ in the exact string that is passed to the io buffer

*/

type FASTAWriter struct {
	*os.File
	*bufio.Writer
}

func NewFASTAWriter(filePath string) FASTAWriter {
	file, err := os.Create(filePath)

	if err != nil {
		fmt.Println("error opening file= ", err)
		os.Exit(1)
	}
	fastawriter := FASTAWriter{Writer: bufio.NewWriter(file), File: file}
	return (fastawriter)

}

func (w *FASTAWriter) Write(r FASTQRead) error {
	//put the read name and sequence in the proper format(adding a carrot) and
	//splitting by newlines (\n)

	//TODO: add a method of splitting up the sequence into 80 character long
	//lines as is the *spec* for fasta
	for_writing := strings.Join([]string{">" + r.Name, r.Sequence}, "\n")
	var err error
	if w.Writer.Available() < len(for_writing) {
		w.Writer.Flush()

		_, err = w.Writer.WriteString(for_writing)
	} else {
		_, err = w.Writer.WriteString(for_writing)
	}
	return (err)

}

func (w *FASTAWriter) Close() {
	w.Writer.Flush()
	w.File.Close()
	fmt.Println("closing FASTA writer")

}

func Testfastq() {

	usr, _ := user.Current()
	homeDir := usr.HomeDir
	testPath := homeDir + "/Desktop/coding/golang/src/testing/sample_50_2.fastq"

	fastqscanner := NewFASTQScanner(testPath)

	defer fastqscanner.Close()

	for {
		myRead := fastqscanner.NextRead()
		//fmt.Print(myRead.Name)
		if myRead.Name == "" {
			break
		}
	}
}

func Testfastqwriter() {

	var err error
	usr, _ := user.Current()
	homeDir := usr.HomeDir
	testPath := homeDir + "/Desktop/coding/golang/src/testing/test.fastq"
	newRead := FASTQRead{
		Name:     "this is my read name",
		Sequence: "AATCGATCGATGAGATAGTC",
		Misc:     "+",
		Quality:  "*(&*^%^%$^%#Q(*&0(&(*&^&^%^%$",
	}
	w := NewFASTQWriter(testPath)

	err = w.Write(newRead)
	err = w.Write(newRead)
	if err != nil {
		fmt.Println(err)
	}

	w.Close()
}

//TODO unit testing?? error handling???
