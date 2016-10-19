package gobioinfo

/*
semi-global (aka *glocal* or global/local) alignment algorithm with affine
gaps penalties

semi-global alignment does not add gap penalties for gaps at the end of each
strand being aligned, and begins alignment matrix traceback from an edge.
Together these two changes to a standard *local* alignment algorithm
(aka Smith-Waterman) make it ideal for aligning sequences of differing
lengths, or sequences which might only partially overlap.

Aside from a generic semi-global alignment algorithm, there are also two usage specific
implementations that allow gap penelties at either the 5'-end or 3'end of the query string,
which is useful when looking for a substring that you know should be at one end or another of the
subject string.

	example:
		3' Linker : GTGTCAGCACA
		Subject: CACATACACTAAGAAGGTCCTGGACGCGTGTAGTCAC

		the 3'-specific adjustment function `SGAlign3p()` would favor this alignment

					                           5'-GTGTCAGCACA-3'
					                              ||||
					5'-CACATACACTAAGAAGGTCCTGGACGCGTGT-3'

		over this alignment:

					5'-GTGTCAGCACA-3'
					          ||||
					       5'-CACATACACTAAGAAGGTCCTGGACGCGTGT-3'

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

import (
	"fmt"
	"strings"
)

const (
	ntA = 65 //
	ntC = 67 //
	ntG = 71 //
	ntT = 84 //
	ntN = 78 // nucleotide N as represented by the unicode/ascii code point for "N"
)

// matrix directions
const (
	match      = iota // m 00000001
	mismatch          // x 00000010
	neutral           // n 00000100
	openInsI          // i	etc...
	extendInsI        //
	openInsJ          // open a gap in the j direction
	extendInsJ        // extend a gap in the j direction
	gap               // "-"
	edge              // "$"
	fromI             // movement of a value from the I matrix to the H matrix
	fromJ             // movement of a value from the J matrix to the H matrix
)

// the ids for three DP matrixes
const (
	matrixH = iota
	matrixI
	matrixJ
)

// PWAlignment creates an pairwise alignment structure
type PWAlignment struct {
	Subject       NucSeq
	Query         NucSeq
	ExpandedCIGAR string // can be used to get represtantations of the gappedsubject, gappedquery, and the alignment itself
	// SubjectStart  int
	// QueryStart    int
	// QueryAlignLen   int
	// SubjectAlignLen int
	// GappedSubject   string --> these can be reconstructed from the cigar if needed
	// GappedQuery     string
	//AlignmentRepresentation string --> same as above
}

type matrixMove struct {
	Score  int
	Origin int
}

type derivedFrom int

type matrixPos struct {
	matrix int
	i      int
	j      int
}

type tracebackBranch struct {
	dir      derivedFrom
	branches []tracebackBranch
}

// the components of the alignment matrices
type matrixElem struct {
	score   int
	numDirs int
	dirs    [3]derivedFrom // each matrix position can be derived from a max of three starting positions
	// this is only possible if they have the same score...this will result in a branch during traceback
	matrix int
}

// a DP alignment matric
type alignMatrix [][]matrixElem

// String creates a printable representation of the alignment matrix
// TODO: implement this
func (m alignMatrix) String() string {
	return ("not sure how to do this yet")
}

// take the potential movements for each matrix element and return the proper elements for
// that matrix cells, including recording the several directions if the highest score is shared
// by several potential movements
func calcElem(list ...matrixElem) matrixElem {
	idxOfBestElems := []int{0}
	idxOfBestElems[0] = 0
	highScore := list[0].score

	for i := 1; i < len(list); i++ {
		if list[i].score > highScore {
			idxOfBestElems = []int{i}
		} else if list[i].score == highScore {
			idxOfBestElems = append(idxOfBestElems, i)
		}
	}
	bestDirs := [3]derivedFrom{}
	for i := 0; i < len(idxOfBestElems); i++ {
		bestDirs[i] = list[idxOfBestElems[i]].dirs[0]
	}

	finalElem := matrixElem{
		score:   highScore,
		numDirs: len(idxOfBestElems),
		dirs:    bestDirs,
	}
	return finalElem
}

// Align applies a semi-global alignment algorithm to the query and subject sequences
// and returns a slide of pairwise alignments (this slice is larger than one, or several alignments
// result in the exact same alignment score)
func (q NucSeq) sGAlign(s NucSeq) []PWAlignment {
	//func (q NucSeq) sGAlign(s NucSeq, end string) []PWAlignment {

	// get the length of the input strings
	lenS := len(s)

	lenI := lenS + 1

	lenQ := len(q)

	lenJ := lenQ + 1

	const (
		h             = 6 // gap opening penalty
		g             = 3 // gap extension penalty
		matchScore    = 3
		mismatchScore = -4
	)

	type expCIGAR []int // expanded cigar string

	// TODO: add visualization of the alignment matrix

	// set up matrices: H is the match matrix, I is the matrix of gaps in the *i* dimentions (eg
	// along the subject), J is the matric for gaps in the *j* dimension (eg along the query)

	// TODO: optimization...used a fixed spot in memory instead of creating and destroying
	// pointers for every iteration of this functions (accept an array as an argument?)

	// create the outer arrays (j dimension = position along query)
	H := make([][]matrixElem, lenJ)
	I := make([][]matrixElem, lenJ)
	J := make([][]matrixElem, lenJ)

	// create the inner arrays (i dimension = position along subject)
	for j := range H {
		H[j] = make([]matrixElem, lenI)
		J[j] = make([]matrixElem, lenI)
		I[j] = make([]matrixElem, lenI)
	}

	// matcher determines the approriate matrix movement value from i-1, j-1 (match, mismatch, nuetral)
	// and returns it
	matcher := func(i int, j int) matrixElem {

		var returnValue matrixElem
		switch {

		case s[i-1] == q[j-1]:
			// if the position is a match
			returnValue.score = H[j-1][i-1].score + matchScore
			returnValue.dirs[0] = match

		case q[j-1] == ntN || s[i-1] == ntN:
			// if the base is undefined treat it is a neutral
			returnValue.score = H[j-1][i-1].score
			returnValue.dirs[0] = neutral

		default:
			// otherwise it is a mismatch
			returnValue.score = H[j-1][i-1].score + mismatchScore
			returnValue.dirs[0] = mismatch
		}
		return (returnValue)
	}

	// fill in the matrices
	for j := range H {
		for i := range H[j] {
			switch {

			// if we are not in the first column or first row
			case i != 0 && j != 0:

				//fill in I
				I[j][i] = calcElem(
					matrixElem{score: H[j][i-1].score - h, dirs: [3]derivedFrom{openInsI}},
					matrixElem{score: I[j][i-1].score - g, dirs: [3]derivedFrom{extendInsI}},
				)

				// fill in J
				J[j][i] = calcElem(
					matrixElem{score: H[j-1][i].score - h, dirs: [3]derivedFrom{openInsJ}},
					matrixElem{score: J[j-1][i].score - g, dirs: [3]derivedFrom{extendInsJ}},
				)

				// fill in H and D
				H[j][i] = calcElem(
					matcher(i, j),
					matrixElem{score: I[j][i].score, dirs: [3]derivedFrom{fromI}},
					matrixElem{score: J[j][i].score, dirs: [3]derivedFrom{fromJ}},
				)
				//fmt.Println(bestMove.Score, bestMove.Origin)

			// if we are in the first position (0,0):
			case i == 0 && j == 0:

				H[j][i] = matrixElem{score: 0, dirs: [3]derivedFrom{fromI}}
				I[j][i] = matrixElem{score: 0, dirs: [3]derivedFrom{fromI}}
				J[i][j] = matrixElem{score: 0, dirs: [3]derivedFrom{fromI}}

			// if we are in the first column
			case i == 0 && j != 0:

				// TODO: fill this out
				// H[j][i] =
				// I[j][i] =
				// J[j][i] =

			// if we are in the first row
			case i != 0 && j == 0:

			}
		}

	}

	fmt.Println(H)
	fmt.Println("---")

	// carries out the traceback operation through the three alignment matrices
	// TODO: this needs to actually traverse the other two matrices properly
	// TODO: this needs to recursivley traceback all "equivalent" matches
	// traceback takes a starting position in a matrix and returns a list of directions
	// to achieve the best alignment
	// this function is recursive, so it calls itself and progressively builds the direction list
	var traceback func(int, int, alignMatrix) []expCIGAR
	traceback = func(startI, startJ int, matrix alignMatrix) []expCIGAR {

		// <<<<<<<<<<<<<<<<<<<<<HERE>>>>>>>>>>>>>>>>>>>>>>

		// outline:
		// for each dir in current elemt
		// - find the next elemement
		// - run traceback
		// - append the right cigar element to the array of results
		// - return the array of results

		var revCIGAR = new(expCIGAR) // the reverse, expanded CIGAR string

		currentElem := matrix[startJ][startI]

		var results []expCIGAR

		for i := 0; i < currentElem.numDirs; i++ {
			var tracebackI, tracebackJ, matrix int

			switch currentElem.dirs[i] {
			case openInsI:
				tracebackI = startI - 1
				tracebackJ = startJ
				matrix = matrixH
			case extendInsI:
				tracebackI = startI - 1
				tracebackJ = startJ
				matrix = matrixI
			case openInsJ:
				tracebackI = startI
				tracebackJ = startJ - 1
				matrix = matrixH
			case extendInsJ:
				tracebackI = startI
				tracebackJ = startJ - 1
				matrix = matrixJ
			case fromI:
				tracebackI = startI
				tracebackJ = startJ
				matrix = matrixI
			case fromJ:
				tracebackI = startI
				tracebackJ = startJ
				matrix = matrixI

			case match:

			case mismatch:

			case neutral:

			}
			prelimResults := traceback(tracebackI, tracebackJ, matrix)
			for i, _ := range prelimResults {
				prelimResults[i] = append(prelimResults[i], currentElem.dirs[i])
			}
			results = append(results, prelimResults...)
		}

		// the end condition is if the starting coordinates are at any edge
		if startI == 0 || startJ == 0 {
			for i, _ := range results {
				results[i] = append(results[i], edge)
			}

		}
		return results
	}

	tbackToPWAlign := func(r expCIGAR, traceStart matrixPos) PWAlignment {
		var cIGAR string

		tracePos := traceStart
		// this check protects against the situation where no alignment takes place
		if r != nil {
			for i := range r {
				nextVector := r[len(r)-1-i]
				var nextLetter string
				switch {
				case nextVector == match:
					nextLetter = "m"
					tracePos.i--
					tracePos.j--
				case nextVector == mismatch:
					nextLetter = "x"
					tracePos.i--
					tracePos.j--
				case nextVector == neutral:
					nextLetter = "n"
					tracePos.i--
					tracePos.j--
				case nextVector == openInsI:
					nextLetter = "i"
					tracePos.i--
				case nextVector == extendInsI:
					nextLetter = "i"
					tracePos.i--
				case nextVector == openInsJ:
					nextLetter = "j"
					tracePos.j--
				case nextVector == extendInsJ:
					nextLetter = "j"
					tracePos.i--
				case nextVector == edge:
					nextLetter = ""
				}
				cIGAR += nextLetter
			}
		} else {
			cIGAR = ""
		}

		// TODO: create new alignment object

		// cigar start = currentPosition
		// cigar end = maxPosition
		var newAlignment PWAlignment

		subjectStart := tracePos.i - 1

		queryStart := tracePos.j - 1

		if cIGAR != "" {
			switch {
			case string(cIGAR[0]) == "i":
				subjectStart = tracePos.i - 1
				queryStart = tracePos.j
			case string(cIGAR[0]) == "j":
				subjectStart = tracePos.i
				queryStart = tracePos.j - 1

			}
		}

		// fmt.Println(subjectStart)
		// fmt.Println(queryStart)

		if tracePos.i != 0 || tracePos.j != 0 {
			newAlignment = PWAlignment{
				Subject:       s,
				Query:         q,
				ExpandedCIGAR: cIGAR,
				SubjectStart:  subjectStart,
				QueryStart:    queryStart,
			}
		} else {
			newAlignment = PWAlignment{
				Subject:       s,
				Query:         q,
				ExpandedCIGAR: cIGAR,
			}
		}
		return newAlignment
	}

	// create traceback

	//find max score in the last row or column
	maxScore := H[lenJ-1][lenI-1]
	maxPos := matrixPos{
		i: lenI - 1,
		j: lenJ - 1,
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
	for j := 0; j < lenJ; j++ {
		if H[j][lenI-1] > maxScore {
			maxScore = H[j][lenI-1]
			maxPosition.i = lenI - 1
			maxPosition.j = j
		}
	}

	completedTraceback := false

	currentPosition = maxPosition
	var nextPosition matrixPosition

	// check to make sure you don't start a traceback at a i=0 or j=0 position
	// (which arises when there essentially is no alignment and the max alignment score is zero)
	if currentPosition.i != 0 && currentPosition.j != 0 {
		for completedTraceback == false {

			currentVector := D[currentPosition.j][currentPosition.i]
			revCIGAR = append(revCIGAR, currentVector)
			nextPosition = traceback(currentPosition)

			// check to see if we have reach the top or left of the matrix
			if nextPosition.i == 0 || nextPosition.j == 0 {
				completedTraceback = true
			} else {
				currentPosition = nextPosition
			}
		}
	} else {
		revCIGAR = nil
	}

	// TODO update this for more alignments if they are "equivalent"
	return tbackToPWAlign(traceback(maxPos))

}

// a convenience function to print create an easily printable representation of a
// PairWiseAlignment struct
func (a PWAlignment) String() string {

	subject := a.Subject
	query := a.Query
	CIGAR := a.ExpandedCIGAR
	subjectStart := a.SubjectStart
	queryStart := a.QueryStart

	//parse CIGAR

	ins := 0
	dels := 0

	var subjectString string
	var queryString string
	var alignString string
	if CIGAR == "" {
		return ("")
	}
	for i, base := range CIGAR {
		subjectPos := subjectStart + i - dels
		queryPos := queryStart + i - ins
		switch string(base) {
		case "m":
			subjectString += string(subject[subjectPos])
			queryString += string(query[queryPos])
			alignString += "|"
		case "x":
			subjectString += string(subject[subjectPos])
			queryString += string(query[queryPos])
			alignString += " "
		case "n":
			subjectString += string(subject[subjectPos])
			queryString += string(query[queryPos])
			alignString += " "
		case "i":
			subjectString += string(subject[subjectPos])
			queryString += "-"
			alignString += " "
			ins++
		case "j":
			subjectString += "-"
			queryString += string(query[queryPos])
			alignString += " "
			dels++
		}

	}

	// alignment.SubjectAlignLen = len(CIGAR) - dels
	// alignment.QueryAlignLen = len(CIGAR) - ins
	// return alignmentRepr
	repr := strings.Join([]string{subjectString, alignString, queryString}, "\n")

	return (repr)
}
