package gobioinfo

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
)

// FASTQScanner is a wrapper around bufio.Scanner, which allows for easiliy accessing each
// read in a FASTQ file in a iterative manner via the Next() function
type FASTQScanner struct {
	*os.File
	*bufio.Scanner
}

// NewFASTQScanner take a string for a path of a fastq file, and returns a FASTQScanner
func NewFASTQScanner(filePath string) FASTQScanner {

	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("error opening file= ", err)
		os.Exit(1)
	}

	var reader io.Reader

	if filepath.Ext(filePath) == ".gz" {
		reader, err = gzip.NewReader(file)
	} else {
		reader = file
	}

	fastqscanner := FASTQScanner{Scanner: bufio.NewScanner(reader), File: file}
	return (fastqscanner)
}

// Close closes the FASTQScanner and the underlying file. Idiomatically this should
// be deferred everytime a new FASTQScanner is created
func (s *FASTQScanner) Close() {
	s.File.Close()
	fmt.Println("closing fastqscanner")
}

// NextRead returns the next read from a FASTQScanner
func (s *FASTQScanner) NextRead() (FASTQRead, error) {

	var ln1, ln3 string
	var ln2, ln4 []rune
	var newRead FASTQRead

	if s.Scanner.Scan() {
		ln1 = s.Scanner.Text()
	} else {
		return FASTQRead{}, errors.New("EOF")
	}
	if s.Scanner.Scan() {
		ln2 = []rune(s.Scanner.Text())
	} else {
		return FASTQRead{}, errors.New("EOF")
	}
	if s.Scanner.Scan() {
		ln3 = s.Scanner.Text()
	} else {
		return FASTQRead{}, errors.New("EOF")
	}
	if s.Scanner.Scan() {
		ln4 = []rune(s.Scanner.Text())
	} else {
		return FASTQRead{}, errors.New("EOF")
	}

	newRead = NewFASTQRead(ln1, ln2, ln3, ln4)

	return newRead, nil

}

// FASTQWriter defines the FASTQWriter structure, which contains a pointer to the file to
// which FASTQReads will be written and to a bufio.Writer instance
type FASTQWriter struct {
	*os.File
	*bufio.Writer
}

//NewFASTQWriter creates a new FASTQWriter instance, which contains the file and writer objects
//from a file path string
func NewFASTQWriter(filePath string) FASTQWriter {

	file, err := os.Create(filePath)

	if err != nil {
		fmt.Println("error opening file= ", err)
		os.Exit(1)
	}
	fastqwriter := FASTQWriter{Writer: bufio.NewWriter(file), File: file}
	return (fastqwriter)
}

//writes a FASTQRead to file, returns an error
func (w *FASTQWriter) Write(r FASTQRead) error {

	//compose FASTQRead struct into the proper format
	forWriting := strings.Join([]string{"@" + r.ID, string(r.Sequence), r.Misc,
		string(r.PHRED.Encoded), ""}, "\n")
	//fmt.Println(for_writing)

	var err error

	//check if there are enough bytes lefts in the writer buffer
	if w.Writer.Available() < len(forWriting) {

		//flush the write buffer
		w.Writer.Flush()

		//then write the current read to the buffer
		_, err = w.Writer.WriteString(forWriting)
	} else {
		//write to the write buffer
		_, err = w.Writer.WriteString(forWriting)
		// fmt.Println("return from writestring: ", i)
	}
	//returns any write error that might have occurred
	return (err)
}

// Close flushes the FASTQWriter buffer and closes it
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

// FASTAWriter is a wrapper around bufio.Writer which allows for FASTA
// formatting writing of FASTAReads
type FASTAWriter struct {
	*os.File
	*bufio.Writer
}

// NewFASTAWriter creates a new FASTAWriter, which is a wrapper around the
// bufio.Writer it takes a string file path as input
func NewFASTAWriter(filePath string) FASTAWriter {
	file, err := os.Create(filePath)

	if err != nil {
		fmt.Println("error opening file= ", err)
		os.Exit(1)
	}
	fastawriter := FASTAWriter{Writer: bufio.NewWriter(file), File: file}
	return (fastawriter)

}

// Writes a FASTA read to file
func (w *FASTAWriter) Write(r FASTQRead) error {

	//put the read name and sequence in the proper format(adding a carrot) and
	//splitting by newlines (\n)

	//TODO: add a method of splitting up the sequence into 80 character long
	//lines as is the *spec* for fasta
	forWriting := strings.Join([]string{">" + r.ID, string(r.Sequence)}, "\n")
	var err error
	if w.Writer.Available() < len(forWriting) {
		w.Writer.Flush()

		_, err = w.Writer.WriteString(forWriting)
	} else {
		_, err = w.Writer.WriteString(forWriting)
	}
	return (err)
}

// Close flushes the bufio.Writer buffer to file and then closes the file
func (w *FASTAWriter) Close() {
	w.Writer.Flush()
	w.File.Close()
	fmt.Println("closing FASTA writer")
}

/*
TODO: add proper error handling
TODO: change the filepath parameters to io.Reader, or io.Writer objects?
*/
