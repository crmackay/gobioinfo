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
	match      = 1 << iota // m 00000001
	mismatch               // x 00000010
	neutral                // n 00000100
	openInsI               // i	etc...
	extendInsI             //
	openInsJ               // j
	extendInsJ             //
	gap                    // "-"
	edge                   // "$"
)

const (
	matrixH = 1 << iota
	matrixI
	matrixJ
)

type alignMatrix [][]matrixElem

// String creates a printable representation of the alignment matrix
func (m alignMatrix) String() string {
	return ("not sure how to do this yet")
}

// PWAlignment creates an pairwise alignment structure
type PWAlignment struct {
	Subject         NucSeq
	Query           NucSeq
	ExpandedCIGAR   string
	SubjectStart    int
	QueryStart      int
	QueryAlignLen   int
	SubjectAlignLen int
	// GappedSubject   string --> these can be reconstructed from the cigar if needed
	// GappedQuery     string
	//AlignmentRepresentation string --> same as above
}

type matrixMovement struct {
	Score  int
	Origin int
}

type matrixDir int

type matrixPos struct {
	matrix int
	i      int
	j      int
}

type matrixElem struct {
	score   int
	numDirs int
	dirs    [3]matrixDir
}

// finds the maxmimum move from a list of possible matrix moves (when filling in the
// alignment matrix)
func maxMove(list ...matrixMovement) matrixMovement {

	max := list[0]

	for i := 1; i < len(list); i++ {
		if list[i].Score > max.Score {
			max = list[i]
		}
	}
	return (max)
}

func maxInt(list ...int) int {

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

func bestMove(list ...matrixElem) matrixElem {
	bestScore := list[0].score
	bestMove := list[0]

	for i := 1; i < len(list); i++ {
		if list[i].score > bestScore {
			bestMove = list[i]
		}
	}
	return bestMove
}

// Align applies a semi-global alignment algorithm to the query and subject sequences
func (q NucSeq) sGAlign(s NucSeq) []PWAlignment {
	//func (q NucSeq) sGAlign(s NucSeq, end string) []PWAlignment {

	// get the length of the input strings
	lenS := len(s)

	lenI := lenS + 1

	lenQ := len(q)

	lenJ := lenQ + 1

	const (
		h             = 6 //gap opening penalty
		g             = 3 //gap extension penalty
		matchScore    = 3
		mismatchScore = -4
	)
	type revCIGAR []matrixDir

	// TODO: add visualization of the alignment matrix

	// set up matrices: H is the match matrix, I is the matrix of gaps in the *i* dimentions (eg
	// along the subject), J is the matric for gaps in the *j* dimension (eg along the query), and
	// D is the matrix of which direction (vector) was chosen to fill the H matrix.
	// --> this approach with the "D" matrix doesnt work...need a direction for every value in all
	// matrices...and potentially a slice of directions, if there is branching

	// matrix position: score,[]directions{insertion extension, new insertion, match, mismatch}

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

	// matcher determines the approriate matrix movement valie from i-1, j-1 (match, mismatch, nuetral)
	// and returns it
	matcher := func(i int, j int) matrixMovement {

		var returnValue matrixMovement
		switch {

		case s[i-1] == q[j-1]:
			// if the position is a match
			returnValue.Score = H[j-1][i-1].score + matchScore
			returnValue.Origin = match

		case q[j-1] == ntN || s[i-1] == ntN:
			// if the base is undefined treat it is a neutral
			returnValue.Score = H[j-1][i-1].score
			returnValue.Origin = neutral

		default:
			// otherwise it is a mismatch
			returnValue.Score = H[j-1][i-1].score + mismatchScore
			returnValue.Origin = mismatch
		}
		return (returnValue)
	}

	// carries out the traceback operation through the three alignment matrices
	// TODO: this needs to actualy traverse the other two matrices properly

	// TODO: this needs to recursivley traceback all "equivalent" matches

	// traceback takes a starting position in a matrix and returns a list of directions
	// to achieve the best alignment
	// this function is recursive, so it calls itself and progressively builds the direction list
	traceback := func(startI, startJ, startMarix int) revCIGAR {

		var results revCIGAR

		// TODO: fix this to not just take the first one...
		currentVector := H[startJ][startI].dirs[0]

		// the end condition is if the starting coordinates are at an edge
		if startI == 0 || startJ == 0 {
			results = append(results, edge)
			return results
		}

		var nextI int
		var nextJ int
		var nextMatrix int

		switch currentVector {
		case match:
			nextI = startI - 1
			nextJ = startJ - 1
			nextMatrix = matrixH
		case mismatch:
			nextI = startI - 1
			nextJ = startJ - 1
			nextMatrix = matrixH
		case neutral:
			nextI = startI - 1
			nextJ = startJ - 1
			nextMatrix = matrixH
		case openInsI:
			nextI = startI - 1
			nextJ = startJ
			nextMatrix = matrixH
		case extendInsI:
			nextI = startI - 1
			nextJ = startJ
			nextMatrix = matrixI
		case openInsJ:
			nextI = startI
			nextJ = startJ - 1
			nextMatrix = matrixH
		case extendInsJ:
			nextI = startI
			nextJ = startJ - 1
			nextMatrix = matrixJ
		}

		results[0] = currentVector
		nextPos := traceback(nextI, nextJ, nextMatrix)
		results = append(results, nextPos)
		return results
	}

	tbackToPWAlign := func(r revCIGAR, traceStart matrixPos) PWAlignment {
		var CIGAR string

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
				CIGAR += nextLetter
			}
		} else {
			CIGAR = ""
		}

		// TODO: create new alignment object

		// cigar start = currentPosition
		// cigar end = maxPosition
		var newAlignment PWAlignment

		subjectStart := tracePos.i - 1

		queryStart := tracePos.j - 1

		if CIGAR != "" {
			switch {
			case string(CIGAR[0]) == "i":
				subjectStart = tracePos.i - 1
				queryStart = tracePos.j
			case string(CIGAR[0]) == "j":
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
				ExpandedCIGAR: CIGAR,
				SubjectStart:  subjectStart,
				QueryStart:    queryStart,
			}
		} else {
			newAlignment = PWAlignment{
				Subject:       s,
				Query:         q,
				ExpandedCIGAR: CIGAR,
			}
		}
		return newAlignment
	}

	// fill in the matrices
	for j := range H {
		for i := range H[j] {
			switch {

			// if we are not in the first column or first row
			case i != 0 && j != 0:

				//fill in I
				I[j][i] = bestMove(
					matrixElem{score: H[j][i-1].score - h, dirs: [3]matrixDir{openInsI}},
					matrixElem{score: I[j][i-1].score - g, dirs: [3]matrixDir{extendInsI}},
				)

				// fill in J
				J[j][i] = bestMove(
					matrixElem{score: H[j-1][i].score - h, dirs: [3]matrixDir{openInsJ}},
					matrixElem{score: J[j-1][i].score - g, dirs: [3]matrixDir{extendInsJ}},
				)

				// fill in H and D
				bestMove := maxMove(
					matcher(i, j),
					matrixMovement{I[j][i], insI},
					matrixMovement{J[j][i], insJ},
				)
				H[j][i] = bestMove.Score
				D[j][i] = bestMove.Origin
				//fmt.Println(bestMove.Score, bestMove.Origin)

			// if we are in the first position (0,0):
			case i == 0 && j == 0:

				H[j][i] = 0
				I[j][i] = 0
				J[i][j] = 0
				D[i][j] = gap

			// if we are in the first column
			case i == 0 && j != 0:

			// if we are in the first row
			case i != 0 && j == 0:

			}
		}

	}

	fmt.Println(H)
	fmt.Println("---")

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

	// TODO update this for more alignments if the are "equivalent"
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
