/* semi-global (aka *glocal* or global/local) alignment algorithm with affine
gaps penalties

semi-global alignment does not add gap penalties for gaps at the end of each
strand being aligned, and begins alignment matrix traceback from an edge.
Together these two changes to a standard *local* alignment algorithm
(aka Smith-Waterman) make it ideal for aligning sequences of differing
lengths, or sequences which might only partially overlap.

Aside from a generic semi-global alignment algorithm, there are also two usage specific
implementations that allow gap penelties at the 5' end of either the query string or the subject,
which is useful when looking for a substring that you know should be at one end or another of the
subject string.

	example:
		3' Linker : GTGTCAGCACA
		Subject: CACATACACTAAGAAGGTCCTGGACGCGTGTAGTCAC

		the 3'-specific adjustment would favor this alignment:

									   GTGTCAGCACA
									   ||||
			CACATACACTAAGAAGGTCCTGGACGCGTGT

		over this alignment:

			GTGTCAGCACA
		           ||||
		   		   CACATACACTAAGAAGGTCCTGGACGCGTGT



*/

package gobioinfo

//"fmt"

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

const (
	ntA = 65 //
	ntC = 67 //
	ntG = 71 //
	ntT = 84 //
	ntN = 78 // nucleotide N as represented by the unicode/ascii code point for "N"
)

const (
	match    = 1 << iota // m
	mismatch             // x
	neutral              // n
	insI                 // i
	insJ                 // j
	gap                  // "-"
)

// PairWiseAlignment creates an pairwise alignment structure
type PairWiseAlignment struct {
	Subject                 NucleotideSequence
	Query                   NucleotideSequence
	ExpandedCIGAR           string
	SubjectStart            int
	QueryStart              int
	QueryAlignLen           int
	SubjectAlignLen         int
	GappedSubject           string
	GappedQuery             string
	AlignmentRepresentation string
}

// PairWiseRepresentation is a convenience struct for printing and displaying a pairwise alignment
type PairWiseRepresentation struct {
	Subject string
	Query   string
	Matches string
}

type matrixMovement struct {
	Score  int
	Origin int
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
	if CIGAR == "" {
		return (alignment)
	}
	for i := 0; i < len(CIGAR); i++ {
		subjectPosition := subjectStart + i - dels
		queryPosition := queryStart + i - ins
		switch string(CIGAR[i]) {
		case "m":
			subjectRepr = subjectRepr + string(subject[subjectPosition])
			queryRepr = queryRepr + string(query[queryPosition])
			alignmentRepr = alignmentRepr + "|"
		case "x":
			subjectRepr = subjectRepr + string(subject[subjectPosition])
			queryRepr = queryRepr + string(query[queryPosition])
			alignmentRepr = alignmentRepr + " "
		case "n":
			subjectRepr = subjectRepr + string(subject[subjectPosition])
			queryRepr = queryRepr + string(query[queryPosition])
			alignmentRepr = alignmentRepr + " "
		case "i":
			subjectRepr = subjectRepr + string(subject[subjectPosition])
			queryRepr = queryRepr + "-"
			alignmentRepr = alignmentRepr + " "
			ins = ins + 1
		case "j":
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

// SG5pAlign aligns...
func (q NucleotideSequence) SG5pAlign(s NucleotideSequence) PairWiseAlignment {
	return (q.sGAlign(s, "five"))
}

// SG3pAlign aligns...
func (q NucleotideSequence) SG3pAlign(s NucleotideSequence) PairWiseAlignment {
	return (q.sGAlign(s, "three"))
}

// Align applies a semi-global alignment algorithm to the query and subject sequences
func (q NucleotideSequence) sGAlign(s NucleotideSequence, end string) PairWiseAlignment {

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

	// TODO: add visualization of the alignment matrix

	//   subjectString

	//TODO: optimization - use a matrix (github.com/gonum/matrix/mat64) instead of a [][]int
	// 	this removes half off all the matric/array lookups

	// set up matrices: H is the match matrix, I is the matrix of gaps in the *i* dimentions (eg
	// along the subject), J is the matric for gaps in the *j* dimension (eg along the query), and
	// D is the matrix of which direction (vector) was chosen to fill the H matrix.
	//
	// create the outer arrays (j dimension = position along query)
	H := make([][]int, lenJ)
	I := make([][]int, lenJ)
	J := make([][]int, lenJ)
	D := make([][]int, lenJ)

	// create the inner arrays (i dimension = position along subject)
	for index := range H {
		H[index] = make([]int, lenI)
		J[index] = make([]int, lenI)
		I[index] = make([]int, lenI)
		D[index] = make([]int, lenI)
	}

	//TODO:optimization , remove all of the casting here, and do comparisons against constants
	matcher := func(i int, j int) matrixMovement {

		var returnValue matrixMovement
		switch {

		case s[i-1] == q[j-1]:
			// if the position is a match
			returnValue.Score = H[j-1][i-1] + matchScore
			returnValue.Origin = match

		case q[j-1] == ntN:
			// if the base is undefined treat it is a neutral
			returnValue.Score = H[j-1][i-1]
			returnValue.Origin = neutral

		default:
			// otherwise it is a mismatch
			returnValue.Score = H[j-1][i-1] + mismatchScore
			returnValue.Origin = mismatch
		}
		return (returnValue)
	}

	traceback := func(start matrixPosition) matrixPosition {
		currentVector := D[start.j][start.i]
		var nextI int
		var nextJ int
		switch currentVector {
		case match:
			nextI = start.i - 1
			nextJ = start.j - 1
		case mismatch:
			nextI = start.i - 1
			nextJ = start.j - 1
		case neutral:
			nextI = start.i - 1
			nextJ = start.j - 1
		case insI:
			nextI = start.i - 1
			nextJ = start.j
		case insJ:
			nextI = start.i
			nextJ = start.j - 1
		}

		nextPos := matrixPosition{i: nextI, j: nextJ}

		return (nextPos)
	}

	// fill matrices

	for j := range H {
		for i := range H[j] {
			switch {
			case i != 0 && j != 0:
				// if we are not in the first column or first row

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
					matcher(i, j),
					matrixMovement{I[j][i], insI},
					matrixMovement{J[j][i], insJ},
				})
				H[j][i] = bestMove.Score
				D[j][i] = bestMove.Origin
				//fmt.Println(bestMove.Score, bestMove.Origin)

			case i == 0 && j == 0:
				// if we are in the first position (0,0)
				H[j][i] = 0
				I[j][i] = 0
				J[i][j] = 0
				D[i][j] = gap

			case i == 0 && j != 0:
				// if we are in the first column
				switch {
				case end == "three":
					// if we are conducting a three-prime bias'ed semi-global alignment we will
					// put all zeros
					H[j][i] = 0
					I[j][i] = 0
					J[j][i] = 0
					D[j][i] = gap
				case end == "five":
					// if we are conducting a five-prim bias'ed semi-global alignment we will
					// fill in the insertions along J matrix, and the H and D matrices accordingly
					I[j][i] = 0

					// fill in J
					J[j][i] = maxInt([]int{
						H[j-1][i] - h,
						J[j-1][i] - g,
					})

					// fill in H and D
					bestMove := max([]matrixMovement{
						// matrixMovement{I[j][i], "i"},
						matrixMovement{J[j][i], insJ},
					})
					H[j][i] = bestMove.Score
					D[j][i] = bestMove.Origin
				}

			case i != 0 && j == 0:
				switch {
				case end == "three":

					I[j][i] = maxInt([]int{
						H[j][i-1] - h,
						I[j][i-1] - g,
					})

					// fill in J
					J[j][i] = 0

					// fill in H and D
					bestMove := max([]matrixMovement{
						matrixMovement{I[j][i], insI},
						// matrixMovement{J[j][i], "j"},
					})
					H[j][i] = bestMove.Score
					D[j][i] = bestMove.Origin

				case end == "five":
					H[j][i] = 0
					I[j][i] = 0
					J[j][i] = 0
					D[j][i] = gap
				}

			}
		}

	}

	// create traceback

	//find max score in the last row or column

	maxScore := H[lenJ-1][lenI-1]
	maxPosition := matrixPosition{
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

	// fmt.Println("max score: ", maxScore)
	// fmt.Println("maxPosition.i: ", maxPosition.i)
	// fmt.Println("maxPosition.j: ", maxPosition.j)

	//fmt.Println(string("max position"), maxPosition.i, maxPosition.j)

	//build reverse cigar string

	var currentPosition matrixPosition

	var revCIGAR []int

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
	// fmt.Println("current position", currentPosition)
	// fmt.Println(revCIGAR)

	// create an forward cigar

	var CIGAR string

	// this check protects against the situation where not alignment takes place
	if revCIGAR != nil {
		for i := range revCIGAR {
			nextVector := revCIGAR[len(revCIGAR)-1-i]
			var nextLetter string
			switch {
			case nextVector == match:
				nextLetter = "m"
			case nextVector == mismatch:
				nextLetter = "x"
			case nextVector == neutral:
				nextLetter = "n"
			case nextVector == insI:
				nextLetter = "i"
			case nextVector == insJ:
				nextLetter = "j"
			case nextVector == gap:
				nextLetter = "-"
			}
			CIGAR += nextLetter
		}
	} else {
		CIGAR = ""
		// fmt.Println("here")
	}

	// fmt.Println(CIGAR)

	// TODO: create new alignment object

	// cigar start = currentPosition
	// cigar end = maxPosition
	var newAlignment PairWiseAlignment

	subjectStart := currentPosition.i - 1

	queryStart := currentPosition.j - 1

	if CIGAR != "" {
		switch {
		case string(CIGAR[0]) == "i":
			subjectStart = currentPosition.i - 1
			queryStart = currentPosition.j
		case string(CIGAR[0]) == "j":
			subjectStart = currentPosition.i
			queryStart = currentPosition.j - 1

		}
	}

	// fmt.Println(subjectStart)
	// fmt.Println(queryStart)

	if currentPosition.i != 0 || currentPosition.j != 0 {
		newAlignment = PairWiseAlignment{
			Subject:       s,
			Query:         q,
			ExpandedCIGAR: CIGAR,
			SubjectStart:  subjectStart,
			QueryStart:    queryStart,
		}
	} else {
		newAlignment = PairWiseAlignment{
			Subject:       s,
			Query:         q,
			ExpandedCIGAR: CIGAR,
		}
	}

	newAlignment = alignmentRepr(newAlignment)

	// TODO: create print method for alignment object

	// TODO: resolve same scores on traceback

	// return the new alignment object
	return (newAlignment)
}

// TODO add a string method for alignment

// func(a *PairWiseAlignment) String {
//
// }

// func (a *PairWiseAlignment) GetSubjLen() int {}
// func (a *PairWiseAlignment) GetSubjLen() int {}
