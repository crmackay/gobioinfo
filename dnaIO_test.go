package gobioinfo

import (
	"bufio"
	"bytes"
	"compress/gzip"
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

	// read gzipped file
	testGZip := bytes.NewBuffer(gzipRawData)
	reader, _ := gzip.NewReader(testGZip)

	fastqscanner := FASTQScanner{Scanner: bufio.NewScanner(reader), File: nil}

	defer fastqscanner.Close()

	for {
		myRead, err := fastqscanner.NextRead()
		if err != nil {
			if err.Error() == "EOF" {
				break
			}
		}
		fmt.Println(myRead.ID)
		fmt.Println(myRead.Sequence)
	}

	// read unzipped file

	// compare both

	// if any errors fail test
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

// this is a 10 read FASTQ file that has been gzipped and then represented as a slice of bytes
var gzipRawData = []byte{
	0x1f, 0x8b, 0x08, 0x08, 0xe2, 0xf4, 0xec, 0x53, 0x00, 0x03, 0x73, 0x61, 0x6d, 0x70, 0x6c, 0x65,
	0x5f, 0x31, 0x30, 0x2e, 0x66, 0x61, 0x73, 0x74, 0x71, 0x00, 0x9d, 0x55, 0xcb, 0x72, 0xa4, 0x38,
	0x10, 0xbc, 0xfb, 0x2b, 0x36, 0xc2, 0x87, 0x0d, 0xcf, 0x8e, 0x23, 0xba, 0x4d, 0xe3, 0x06, 0x9a,
	0xd6, 0xab, 0x24, 0x4a, 0x70, 0xe8, 0xcb, 0x12, 0x31, 0xf3, 0xff, 0x5f, 0x32, 0x99, 0x85, 0xdf,
	0xde, 0xb1, 0xbd, 0x23, 0x40, 0x12, 0x25, 0xa0, 0x3b, 0x2b, 0xb3, 0x52, 0xa1, 0xfe, 0x98, 0x6f,
	0xff, 0x5d, 0xdb, 0xfb, 0xdd, 0xb0, 0x6f, 0xdb, 0x41, 0xda, 0xe3, 0xa1, 0x44, 0xf9, 0xf9, 0x73,
	0x68, 0x86, 0xfd, 0x7e, 0xb7, 0x47, 0xd7, 0xf6, 0xc3, 0xbe, 0x6f, 0x8e, 0x7f, 0xed, 0x87, 0xcb,
	0xb0, 0x1b, 0xae, 0x54, 0xd6, 0xa8, 0xaa, 0xb8, 0xa2, 0x68, 0x5c, 0x55, 0x54, 0x57, 0x9e, 0x9c,
	0xe2, 0x94, 0x28, 0x71, 0xc5, 0x85, 0xa7, 0xa2, 0xe2, 0xd4, 0x55, 0x10, 0xe6, 0xb3, 0x82, 0x87,
	0xf0, 0xea, 0xca, 0x35, 0x04, 0xa3, 0xbd, 0x89, 0x20, 0x5f, 0xc4, 0xc0, 0x20, 0x3f, 0x81, 0x2f,
	0xc5, 0xab, 0x7f, 0xae, 0x42, 0x08, 0x13, 0x5b, 0xad, 0xdb, 0xa0, 0x75, 0x09, 0x53, 0xf5, 0x69,
	0xaa, 0xd3, 0x58, 0x67, 0xd5, 0x79, 0x99, 0x75, 0x59, 0x54, 0xcf, 0x78, 0x57, 0xb5, 0xce, 0x12,
	0xce, 0x39, 0x17, 0xa9, 0xb5, 0xe4, 0xc6, 0x05, 0xc9, 0x38, 0xa2, 0x88, 0x73, 0x21, 0x0e, 0xbd,
	0xeb, 0xfb, 0x70, 0x33, 0x3a, 0x1f, 0x5c, 0x68, 0x6f, 0xba, 0x31, 0x88, 0x0c, 0xd1, 0xc5, 0x31,
	0x5e, 0x85, 0x4f, 0xc1, 0xef, 0x3b, 0x80, 0x3f, 0x74, 0x4f, 0xe0, 0x57, 0x7d, 0x0b, 0x1e, 0x18,
	0x81, 0x8b, 0xa8, 0x38, 0x51, 0xc3, 0xbd, 0x41, 0xb9, 0xac, 0x04, 0x7b, 0x41, 0x7b, 0x07, 0xf8,
	0xb2, 0x35, 0x59, 0x2f, 0x06, 0xf8, 0x12, 0xb7, 0x26, 0x11, 0x8f, 0x03, 0xbc, 0x0f, 0x3e, 0xb3,
	0x4d, 0x75, 0x90, 0x29, 0x88, 0x4e, 0x33, 0x5a, 0x1a, 0x70, 0x7c, 0xf3, 0x59, 0x81, 0xbe, 0xa2,
	0x2b, 0x85, 0x79, 0x38, 0x87, 0xeb, 0xdb, 0x36, 0x14, 0xb9, 0x46, 0xbb, 0xfd, 0x7e, 0x4a, 0x24,
	0xc0, 0xfb, 0xef, 0xcd, 0xe9, 0x1c, 0x2d, 0xf6, 0xdb, 0xf6, 0x39, 0xf8, 0xbb, 0x3d, 0x66, 0xfd,
	0x7d, 0xf7, 0x01, 0xf3, 0x4a, 0xd4, 0xe4, 0x1b, 0x24, 0x73, 0x66, 0xd9, 0xe1, 0xad, 0xb1, 0x4f,
	0x51, 0x98, 0x2c, 0x78, 0x87, 0x44, 0x59, 0x1a, 0x56, 0xe6, 0x30, 0xf2, 0x9d, 0x95, 0xd0, 0x4d,
	0x15, 0x4c, 0x9d, 0xd8, 0xea, 0x0b, 0xe6, 0x55, 0x73, 0xc5, 0x6f, 0xd4, 0x0a, 0xec, 0x0b, 0x31,
	0xcf, 0xb3, 0x4d, 0x17, 0x64, 0xa0, 0xc6, 0x50, 0x6e, 0x8f, 0xa1, 0xd4, 0x32, 0x95, 0x52, 0x22,
	0x52, 0x25, 0x29, 0x63, 0x2a, 0x05, 0xf8, 0x93, 0xf0, 0x3e, 0x07, 0x27, 0x60, 0x9a, 0xf3, 0x82,
	0xbb, 0xd4, 0x8d, 0x58, 0x10, 0x87, 0x48, 0x90, 0xfe, 0x0b, 0xcc, 0x1f, 0x29, 0xfb, 0x63, 0xff,
	0x08, 0x1e, 0xec, 0xbd, 0x03, 0xbf, 0xfd, 0xef, 0x67, 0x39, 0xaf, 0xef, 0x22, 0xbf, 0x13, 0xf8,
	0x53, 0x13, 0x66, 0x62, 0x65, 0x87, 0xa9, 0x81, 0x17, 0x03, 0x9f, 0x6b, 0xae, 0xd3, 0x0c, 0x96,
	0xeb, 0x3c, 0xf9, 0x5c, 0xab, 0x9f, 0xa8, 0xf8, 0xad, 0x9d, 0xa6, 0x58, 0x54, 0x02, 0xf8, 0xaf,
	0xda, 0x97, 0x8c, 0x87, 0xd3, 0xdf, 0xfd, 0xa9, 0x23, 0x09, 0x31, 0x9d, 0x00, 0x74, 0x70, 0xcd,
	0xa1, 0x71, 0x87, 0xae, 0x19, 0xe2, 0x9f, 0x31, 0xbf, 0x3f, 0x10, 0x7c, 0xff, 0x04, 0xde, 0xc4,
	0xfd, 0xae, 0xe6, 0xc1, 0xb0, 0xac, 0x26, 0x00, 0x13, 0x06, 0x01, 0x83, 0x52, 0x12, 0x19, 0x31,
	0x44, 0x2a, 0x51, 0x0c, 0x9b, 0x01, 0x16, 0x73, 0x09, 0x7a, 0x04, 0x57, 0x8d, 0xfc, 0x57, 0xf9,
	0x7a, 0x66, 0x1e, 0xb2, 0x47, 0xd5, 0xd7, 0x0a, 0xbc, 0x4b, 0x61, 0x0a, 0xa6, 0x3a, 0x63, 0x04,
	0xfb, 0xcb, 0xb2, 0xcc, 0x79, 0x56, 0x81, 0xfc, 0x21, 0x82, 0x19, 0x23, 0x0e, 0xa4, 0xa1, 0x92,
	0xe1, 0x09, 0x95, 0x52, 0x0a, 0x3c, 0x40, 0xda, 0x73, 0x8a, 0x29, 0xa7, 0x0c, 0x35, 0xb0, 0x8a,
	0x82, 0xf8, 0x81, 0xc9, 0x09, 0x07, 0x27, 0xee, 0x0b, 0xb2, 0xef, 0x76, 0xaf, 0x0c, 0xef, 0xf2,
	0x5f, 0xe0, 0xa3, 0xa1, 0x26, 0x68, 0x53, 0x30, 0x83, 0x42, 0x90, 0x26, 0xef, 0x75, 0x73, 0x44,
	0x31, 0x9c, 0x22, 0x9b, 0x41, 0xc4, 0xcd, 0x1e, 0xc5, 0x34, 0x60, 0xba, 0x78, 0x81, 0x1e, 0xe0,
	0xaf, 0xf7, 0x30, 0x30, 0x12, 0xcf, 0x66, 0x3a, 0x27, 0xf9, 0xe8, 0xe6, 0x92, 0x51, 0x04, 0x45,
	0x6b, 0xd2, 0xda, 0xc5, 0x0a, 0x64, 0xe1, 0xe4, 0x5d, 0xe4, 0x71, 0x1e, 0xc3, 0x18, 0x08, 0x0d,
	0x3f, 0x82, 0xd1, 0xc1, 0xf4, 0x5c, 0x18, 0x06, 0x27, 0x3e, 0xe5, 0x31, 0x65, 0xaa, 0x2b, 0x25,
	0x77, 0xd7, 0x3d, 0x5a, 0xc1, 0xe7, 0xe0, 0x0f, 0x07, 0x82, 0x6f, 0x9b, 0x0f, 0x98, 0x37, 0xcb,
	0x23, 0xa9, 0xf8, 0x55, 0x35, 0xf7, 0xdd, 0x98, 0xb7, 0xa9, 0xe5, 0x20, 0xd2, 0x02, 0x8d, 0x61,
	0xd5, 0x87, 0xb8, 0x7d, 0x83, 0x95, 0x80, 0xc0, 0xdb, 0x4a, 0x01, 0xf8, 0xf4, 0x20, 0x7b, 0x03,
	0xcf, 0x4a, 0x87, 0xb5, 0x5b, 0x06, 0x90, 0x00, 0x92, 0x0e, 0x3d, 0xe4, 0xe7, 0x96, 0x78, 0x09,
	0xfb, 0x8c, 0xfa, 0x2e, 0x1c, 0x7d, 0x7a, 0x5a, 0x0d, 0x63, 0xef, 0xa8, 0xbd, 0x90, 0xd3, 0x88,
	0x83, 0x12, 0x88, 0x03, 0xf6, 0x81, 0xcf, 0xc1, 0x37, 0xcd, 0x1d, 0x0d, 0x6f, 0xf7, 0x91, 0xe1,
	0x09, 0xb1, 0x19, 0xde, 0xb8, 0x6e, 0xd2, 0xe7, 0xca, 0x6b, 0x27, 0xfb, 0x1f, 0xf5, 0xff, 0xb2,
	0xe6, 0x89, 0x5e, 0x6d, 0x6b, 0xa3, 0xf2, 0x97, 0xe9, 0x5b, 0x0f, 0x12, 0x27, 0x85, 0xcd, 0xcf,
	0xf8, 0x0f, 0x13, 0xb6, 0xba, 0x10, 0x54, 0x72, 0x2c, 0xe5, 0x0c, 0xd4, 0x1e, 0x17, 0x04, 0x1e,
	0x71, 0x0d, 0x21, 0xa5, 0xd4, 0xf9, 0x1e, 0x85, 0x1f, 0x6e, 0x5a, 0xd8, 0x5c, 0xe3, 0xe8, 0x7b,
	0x5d, 0x60, 0xef, 0x6e, 0xda, 0x3e, 0xe5, 0x2f, 0x80, 0xbf, 0xe7, 0x56, 0x77, 0x7f, 0x78, 0x32,
	0xbc, 0xc8, 0x4d, 0x4d, 0x49, 0xb1, 0xa1, 0x31, 0x9c, 0xd1, 0x0c, 0xdf, 0xe0, 0x18, 0x8e, 0x47,
	0xe4, 0x6a, 0x8b, 0x51, 0x1e, 0xac, 0x8c, 0x59, 0xda, 0xe2, 0xb2, 0xcd, 0xad, 0xfe, 0x59, 0x2b,
	0xcc, 0xa5, 0x39, 0xa9, 0x59, 0x02, 0xc0, 0x77, 0x3e, 0x18, 0x7f, 0x40, 0x3f, 0x70, 0x77, 0xd7,
	0xba, 0x9d, 0x61, 0xd2, 0x04, 0xcf, 0x87, 0xf7, 0x49, 0xd6, 0x54, 0xf1, 0xa5, 0x89, 0x1e, 0x90,
	0x10, 0x9c, 0x03, 0x2a, 0x3b, 0x96, 0x3a, 0xa1, 0xcb, 0xf9, 0x94, 0x52, 0x94, 0x94, 0xbc, 0xef,
	0x3d, 0x12, 0x09, 0xfd, 0x5b, 0x41, 0xf8, 0xbb, 0x31, 0x74, 0xe3, 0x38, 0xfa, 0x31, 0x7e, 0x49,
	0xf6, 0x47, 0x03, 0xdf, 0x7e, 0xe0, 0xf6, 0x9c, 0xe8, 0x6b, 0x8a, 0xff, 0xd0, 0xf0, 0xb7, 0xfe,
	0x41, 0xf6, 0x79, 0x23, 0x9e, 0x7b, 0xfd, 0xe6, 0x6f, 0x33, 0xac, 0xad, 0x82, 0x76, 0x98, 0x1e,
	0xdd, 0x4d, 0x03, 0x2a, 0x7f, 0x61, 0x46, 0x8a, 0x79, 0x23, 0x76, 0x39, 0x07, 0x0b, 0x70, 0xe7,
	0x33, 0x24, 0xef, 0xb9, 0xaf, 0x71, 0x00, 0x66, 0x94, 0x42, 0xb2, 0xdd, 0x2e, 0xbc, 0x71, 0xfb,
	0x5f, 0x32, 0xab, 0x05, 0x91, 0xe2, 0x09, 0x00, 0x00,
}

var rawUnzipped = `@HWI-ST560:155:C574EACXX:3:1101:1159:1937 1:N:0:
GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAA
+
@@@FFFFFHHFFFFFHGHJ@FH?BFHF<HIGGIJIGJJGG=CCGGGHIC@=DDECHHED3>@CDCDCACC>>@A:9>99@)<>?@>@5)8<@CC:A>A<A
@HWI-ST560:155:C574EACXX:3:1101:1118:1948 1:N:0:
GTGAGGGAGGACGATGCGGAAGTGTTCACAGTGGCTAAGTTCTGCNTGTCANNNNCTTCCAGCGGTCGTATGNNNNNNNCTNCTTGANAAAAAAACAANN
+
?@?DDDDDFH:CF@CGFIIIIB:B:B*?DGIGGHGIGEEHIGG=@#-5@EC####-,;BCCAC??,3;=AC#############################
@HWI-ST560:155:C574EACXX:3:1101:1213:1968 1:N:0:
GCTAGGGAGGACGATGCGGGACAGATACCCTCAGAGAGGGTACCTAAGATGATATGCTTAAGGAACTTCCTGTGAAATACTTTGACGTGTCAGTCCTTCC
+
@@@FFFFFGGDHGCGHHIIIJIGGHIIHIIIJJDGIHA@E-7@EHEFEEEADDDCBDFEECECCABCDDDCD@>C@CCCDDDEDCDB8<ABCC>CCC@C9
@HWI-ST560:155:C574EACXX:3:1101:1179:1979 1:N:0:
AGCAGGGAGGACGATGCGGGTCAGTCACTTCCAGCGTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAACAATATAATAACA
+
@@CFFFFFDHDHFIGIGHIF?DHH?FGIJIIIIIII;FAEGC@EEHHG9EDFFDB'9;8ACCCAB;B8<:>343>483:A####################
@HWI-ST560:155:C574EACXX:3:1101:1149:1999 1:N:0:
CACAGGGAGGACGATGCGGTGGAGACTACAGACTAGGCAGCTTTTCCTATTTACCACCATAATGAAAACTGGTGGCTGATTTCTTTGTGTCAGTCACTTC
+
@@@FFFFDFHHHFHHJIIJEGIGHFHIEGIIHIJJJIDIGCHGIDGIICHGHGHEHHHEDCDFF@CEEDECC5=BABDBDCECDDDD@C?:ACCC@4>C>
@HWI-ST560:155:C574EACXX:3:1101:1280:1937 1:N:0:
NACAGGGAGGACGATGCGGTGAAGACTGGCAACGTGGCGGCAAACCTGTGTGTGGTGCTTTCTCCAGTGTTAGCTGCCCCTGCTGTCGTGTCAGTCACTT
+
#1=DDFFDHHHHHIIIJJIGHIJIGIEDGHHEGHBGH8AHABD@;?>A>A>A=<@<@ACCCCCA@AC>CDC>@::>C?BD<BDACAABB>28AC######
@HWI-ST560:155:C574EACXX:3:1101:1440:1953 1:N:0:
CACAGGGAGGACGATGCGGGGAAGATGACCCGCCGGGCAGCTTCCGGGAAACCAAAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCGTGTCAGTCACTTC
+
B@CFFFFFHHHHHGHIIJJGIJJIGIJIJJJIHHFDDDDDDDDDDDBDDDBCBDDDDCCDEDDDD?BDDDDDDDDD@<9>CCAC@DB<B<BCDDDA:>@A
@HWI-ST560:155:C574EACXX:3:1101:1332:1960 1:N:0:
GCTAGGGAGGACGATGCGGGCAAGTAGCTTATCAGACTGATGTTGACGTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAA
+
@@CFFFFFFHHHGIGGIJJIIJJF*9D<BFGEHIIGAGFG=C@@GCDAEE=CDE?CDECECACEC:@BBB8?9<:>@)5ABC3>CCDD8@CCDD>)59BD
@HWI-ST560:155:C574EACXX:3:1101:1368:1964 1:N:0:
AATTCAGCGCCTTCCATGATGAGACAGGCCGTTTGAATGTTGACGGGATGAACATAATAAGCAATGACGGCAGCAATAAACTCAACAGGAGCAGGAAAAC
+
8?@DDDD?FFH:?BFHGHHGHHGH@FGBGGH?DHCDGBHGGCFEGIIBBGGI@DCEAEHFEAEDD;BBACBB??9?@@CC>@@ACCC?2<@8<<<?<A##
@HWI-ST560:155:C574EACXX:3:1101:1478:1965 1:N:0:
AGCAGGGAGGACGATGCGGACGAGGTGTCAGTCACTTCGTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAACAAAAACAAC
+
B@CDFFFFHHDDFHHIJJJIIIDGHBFGFHIGHGHG@EGHJHHGHEFFFFDDFE>ABD>==?BD?CDB8BD?C@ACDCCBCCDDD@##############
`
