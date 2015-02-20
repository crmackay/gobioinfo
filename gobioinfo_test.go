package gobioinfo

import "testing"

func TestFASTQReader(t *testing.T) {

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

func TestFASTQWriter(t *testing.T) {
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

func TestFASTQReader(t *testing.T) {
}

func TestFASTQWriter(t *testing.T) {
}

funf TestSGAlign(t *testing.T) {

	subject := "GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTGCTTG"
	query := "GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTGTTCTACTTGAA"

	result := align(subject, query)
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