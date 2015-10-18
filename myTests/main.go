package main

import (
	"fmt"
	bio "github.com/crmackay/gobioinfo"
)

func main() {
	reader := bio.NewFASTQScanner("/Users/christophermackay/Desktop/deepseq_data/pir1/hits-clip/working_data/sample_data/fastq/sblade-test/sample_25000_2.fastq")

	read, _ := reader.NextRead()

	fmt.Println(string(read.Sequence))

	for i := 0; i < 10; i++ {
		read, _ := reader.NextRead()

		fmt.Println(len(string(read.Sequence)))
	}
}
