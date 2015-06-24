/* semi-global (aka *glocal* or global/local) alignment algorithm with affine
gaps penalties

semi-global alignment does not add gap penalties for gaps at the end of each
strand being aligned, and begins alignment matrix traceback from an edge.
Together these two changes to a standard *local* alignment algorithm
(aka Smith-Waterman) make it ideal for aligning sequences of differing
lengths, or sequences which might only partially overlap.

*/

package gobioinfo

import (
	"fmt"
)

/*

# Example:

linker: GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG
fastq read:

ln1: @HWI-ST560:155:C574EACXX:3:1101:1159:1937 1:N:0:
ln2: GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAA
ln3: +
ln4: @@@FFFFFHHFFFFFHGHJ@FH?BFHF<HIGGIJIGJJGG=CCGGGHIC@=DDECHHED3>@CDCDCACC>>@A:9>99@)<>?@>@5)8<@CC:A>A<A


result:

linker: start-GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG-end
              |||| ||||||||||||||||||||||||||||||||||||
read     [...]GTGT-AGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG-end


*/

// PairWiseAlignment creates an pairwise alignment structure
type PairWiseAlignment struct {
	Subject                 NucleotideSequence
	Query                   NucleotideSequence
	ExpandedCIGAR           []string
	SubjectStart            int
	QueryStart              int
	SubjectAlignLen         int
	QueryAlignLen           int
	GappedSubject           string
	GappedQuery             string
	AlignmentRepresentation string
}

type matrixMovement struct {
	Score  int
	Origin string
}

type matrixPosition struct {
	i int
	j int
}

func max(list []matrixMovement) matrixMovement {

	max := list[0]

	for i := 1; i < len(list); i++ {
		if list[i].Score > max.Score {
			max = list[i]
		}
	}
	return (max)
}

func maxInt(list []int) int {

	max := list[0]
	for i := 1; i < len(list); i++ {
		if list[i] > max {
			max = list[i]
		}
	}
	//	fmt.Println(list)
	//	fmt.Println(max)
	return (max)
}

func alignmentRepr(alignment PairWiseAlignment) PairWiseAlignment {

	subject := alignment.Subject
	query := alignment.Query
	CIGAR := alignment.ExpandedCIGAR
	subjectStart := alignment.SubjectStart
	queryStart := alignment.QueryStart

	//parse CIGAR

	ins := 0
	dels := 0

	var subjectRepr string
	var queryRepr string
	var alignmentRepr string

	for i := 0; i < len(CIGAR); i++ {
		subjectPosition := subjectStart + i - dels
		queryPosition := queryStart + i - ins
		switch {
		case string(CIGAR[i]) == "m":
			subjectRepr = subjectRepr + string(subject[subjectPosition])
			queryRepr = queryRepr + string(query[queryPosition])
			alignmentRepr = alignmentRepr + "|"
		case string(CIGAR[i]) == "x":
			subjectRepr = subjectRepr + string(subject[subjectPosition])
			queryRepr = queryRepr + string(query[queryPosition])
			alignmentRepr = alignmentRepr + " "
		case string(CIGAR[i]) == "n":
			subjectRepr = subjectRepr + string(subject[subjectPosition])
			queryRepr = queryRepr + string(query[queryPosition])
			alignmentRepr = alignmentRepr + " "
		case string(CIGAR[i]) == "i":
			subjectRepr = subjectRepr + string(subject[subjectPosition])
			queryRepr = queryRepr + "-"
			alignmentRepr = alignmentRepr + " "
			ins = ins + 1
		case string(CIGAR[i]) == "j":
			subjectRepr = subjectRepr + "-"
			queryRepr = queryRepr + string(query[queryPosition])
			alignmentRepr = alignmentRepr + " "
			dels = dels + 1
		}

	}
	alignment.GappedSubject = subjectRepr
	alignment.GappedQuery = queryRepr
	alignment.AlignmentRepresentation = alignmentRepr

	alignment.SubjectAlignLen = len(CIGAR) - dels
	alignment.QueryAlignLen = len(CIGAR) - ins

	return (alignment)

}

// alignment algorithm

// Align applies a semi-global alignment algorithm to the query and subject sequences
func (query NucleotideSequence) Align(subject NucleotideSequence) PairWiseAlignment {

	// get the length of the input strings
	lenSubject := len(subject)

	lenI := lenSubject + 1

	lenQuery := len(query)

	lenJ := lenQuery + 1

	const (
		h             = 30 //gap opening penalty
		g             = 5  //gap extension penalty
		matchScore    = 15
		mismatchScore = -13
	)

	/* TODO: add visualization of the alignment matrix

	   subjectString


	*/

	// create matrices (j dimension = position along query)

	H := make([][]int, lenJ)
	I := make([][]int, lenJ)
	J := make([][]int, lenJ)
	D := make([][]string, lenJ)

	// create the inner arrays (i dimension = position along subject)

	for index := range H {

		H[index] = make([]int, lenI)
		J[index] = make([]int, lenI)
		I[index] = make([]int, lenI)
		D[index] = make([]string, lenI)

	}

	match := func(i int, j int) matrixMovement {

		var returnValue matrixMovement
		switch {
		case string(subject[i-1]) == string(query[j-1]):
			/*if the position is a match*/
			returnValue.Score = H[j-1][i-1] + matchScore
			returnValue.Origin = "m"

		case string(query[j-1]) == "N":
			/*if the base is undefined treat is a neutral*/
			returnValue.Score = H[j-1][i-1]
			returnValue.Origin = "n"

		default:
			/*otherwise it is a mismatch*/
			returnValue.Score = H[j-1][i-1] + mismatchScore
			returnValue.Origin = "x"
		}
		return (returnValue)
	}

	traceback := func(start matrixPosition) matrixPosition {
		currentVector := D[start.j][start.i]
		var nextI int
		var nextJ int
		switch {
		case string(currentVector) == "m":
			nextI = start.i - 1
			nextJ = start.j - 1
		case string(currentVector) == "x":
			nextI = start.i - 1
			nextJ = start.j - 1
		case string(currentVector) == "n":
			nextI = start.i - 1
			nextJ = start.j - 1
		case string(currentVector) == "i":
			nextI = start.i - 1
			nextJ = start.j
		case string(currentVector) == "j":
			nextI = start.i
			nextJ = start.j - 1
		}

		nextPosition := matrixPosition{i: nextI, j: nextJ}

		return (nextPosition)
	}

	// fill matrices

	for j := range H {
		for i := range H[j] {
			if i != 0 && j != 0 {

				//fill in I

				I[j][i] = maxInt([]int{
					H[j][i-1] - h,
					I[j][i-1] - g,
				})

				// fill in J
				J[j][i] = maxInt([]int{
					H[j-1][i] - h,
					J[j-1][i] - g,
				})

				// fill in H and D
				bestMove := max([]matrixMovement{
					match(i, j),
					matrixMovement{I[j][i], "i"},
					matrixMovement{J[j][i], "j"},
				})
				H[j][i] = bestMove.Score
				D[j][i] = bestMove.Origin
				//fmt.Println(bestMove.Score, bestMove.Origin)

			} else {

				H[j][i] = 0
				D[j][i] = "-"

			}
		}

	}

	//fmt.Println(H)
	//fmt.Println(I)
	//fmt.Println(J)
	//fmt.Println(D)

	// create traceback

	//find max score in the last row or column
	var maxScore int
	var maxPosition matrixPosition

	for j := 0; j < lenJ; j++ {
		if H[j][lenI-1] > maxScore {
			maxScore = H[j][lenI-1]
			maxPosition.i = lenI - 1
			maxPosition.j = j
		}
	}

	for i := 0; i < lenI; i++ {

		if H[lenJ-1][i] > maxScore {
			maxScore = H[lenJ-1][i]
			maxPosition.i = i
			maxPosition.j = lenJ - 1
		}
		//	fmt.Println("i: ", i)
		//	fmt.Println("H[lenJ-1][i]: ", H[lenJ-1][i])
	}
	//	fmt.Println("max score: ", maxScore)
	//	fmt.Println("maxPosition.i: ", maxPosition.i)
	//	fmt.Println("maxPosition.j: ", maxPosition.j)

	//fmt.Println(string("max position"),maxPosition.i, maxPosition.j)

	//build reverse cigar string

	var currentPosition matrixPosition

	// TODO: change CIGAR, revCIGAR, and the movement matrix to []rune, and also rename to revCigar
	var revCIGAR []string

	completedTraceback := false

	currentPosition = maxPosition

	for completedTraceback == false {

		currentVector := D[currentPosition.j][currentPosition.i]

		revCIGAR = append(revCIGAR, currentVector)

		nextPosition := traceback(currentPosition)

		if nextPosition.i == 0 || nextPosition.j == 0 {
			completedTraceback = true
		} else {
			currentPosition = nextPosition
		}
	}
	//fmt.Println("current position", currentPosition)
	//fmt.Println(revCIGAR)

	// create an forward cigar

	var CIGAR []string

	for i := 0; i < len(revCIGAR); i++ {
		nextLetter := revCIGAR[len(revCIGAR)-1-i]
		CIGAR = append(CIGAR, nextLetter)
	}
	fmt.Println(CIGAR)

	// TODO: create new alignment object

	// cigar start = currentPosition
	// cigar end = maxPosition

	newAlignment := PairWiseAlignment{
		Subject:       subject,
		Query:         query,
		ExpandedCIGAR: CIGAR,
		SubjectStart:  currentPosition.i - 1,
		QueryStart:    currentPosition.j - 1,
	}
	newAlignment = alignmentRepr(newAlignment)

	fmt.Println(newAlignment.GappedQuery)
	fmt.Println(newAlignment.AlignmentRepresentation)
	fmt.Println(newAlignment.GappedSubject)

	// TODO: create print method for alignment object

	// TODO: resolve same scores on traceback

	// return the new alignment object
	return (newAlignment)
}
