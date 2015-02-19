/* semi-global (aka *glocal* or global/local) alignment algorithm with affine
gaps penalties

semi-global alignment does not add gap penalties for gaps at the end of each
strand being aligned, and begins alignment matrix traceback from an edge.
Together these two changes to a standard *local* alignment algorithm
(aka Smith-Waterman) make it ideal for aligning sequences of differing
lengths, or sequences which might only partially overlap.

*/

package gobioinfo

//package RNASeq

import (
	"fmt"
	//"strings"
)

/* semi-local alignment
    - allow local


test linker: GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG
test read:

@HWI-ST560:155:C574EACXX:3:1101:1159:1937 1:N:0:
GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAA
+
@@@FFFFFHHFFFFFHGHJ@FH?BFHF<HIGGIJIGJJGG=CCGGGHIC@=DDECHHED3>@CDCDCACC>>@A:9>99@)<>?@>@5)8<@CC:A>A<A

*/

type MatrixMovement struct {
	Score  int
	Origin string
}

type MatrixPosition struct {
	i int
	j int
}

// create an pairwise alignment structure

type PairWiseAlignment struct {
	Subject                 string
	Query                   string
	ExpandedCIGAR           []string
	SubjectStart            int
	QueryStart              int
	SubjectAlignLen         int
	QueryAlignLen           int
	GappedSubject           string
	GappedQuery             string
	AlignmentRepresentation string
}

const (
	h = 2
	g = 4
)

func max(list []MatrixMovement) MatrixMovement {

	var max MatrixMovement
	for i := 0; i < len(list); i++ {
		if i == 0 {
			max = list[i]
		} else if list[i].Score > max.Score {
			max = list[i]
		}
	}
	return (max)
}

func maxInt(list []int) int {

	max := 0
	for i := 0; i < len(list); i++ {
		if list[i] > max {
			max = list[i]
		}
	}
	return (max)
}

func AlignmentRepr(alignment PairWiseAlignment) PairWiseAlignment {

	subject := alignment.Subject
	query := alignment.Query
	CIGAR := alignment.ExpandedCIGAR
	subject_start := alignment.SubjectStart
	query_start := alignment.QueryStart

	//parse CIGAR

	ins := 0
	dels := 0

	var subject_repr string
	var query_repr string
	var alignment_repr string

	for i := 0; i < len(CIGAR); i++ {
		subject_position := subject_start + i - dels
		query_position := query_start + i - ins
		switch {
		case string(CIGAR[i]) == "m":
			subject_repr = subject_repr + string(subject[subject_position])
			query_repr = query_repr + string(query[query_position])
			alignment_repr = alignment_repr + "|"
		case string(CIGAR[i]) == "x":
			subject_repr = subject_repr + string(subject[subject_position])
			query_repr = query_repr + string(query[query_position])
			alignment_repr = alignment_repr + " "
		case string(CIGAR[i]) == "n":
			subject_repr = subject_repr + string(subject[subject_position])
			query_repr = query_repr + string(query[query_position])
			alignment_repr = alignment_repr + " "
		case string(CIGAR[i]) == "i":
			subject_repr = subject_repr + string(subject[subject_position])
			query_repr = query_repr + "-"
			alignment_repr = alignment_repr + " "
			ins = ins + 1
		case string(CIGAR[i]) == "j":
			subject_repr = subject_repr + "-"
			query_repr = query_repr + string(query[query_position])
			alignment_repr = alignment_repr + " "
			dels = dels + 1
		}

	}
	alignment.GappedSubject = subject_repr
	alignment.GappedQuery = query_repr
	alignment.AlignmentRepresentation = alignment_repr

	alignment.SubjectAlignLen = len(CIGAR) - dels
	alignment.QueryAlignLen = len(CIGAR) - ins

	return (alignment)

}

// alignment algorithm

func align(subject string, query string) PairWiseAlignment {

	// get the length of the input strings
	len_subject := len(subject)

	len_i := len_subject + 1

	len_query := len(query)

	len_j := len_query + 1

	const (
		h = 5
		g = 2
	)

	/* visualization of the alignment matrix

	   subject_string


	*/

	// create matrices (i dimension = position along subject)

	H := make([][]int, len_j)
	I := make([][]int, len_j)
	J := make([][]int, len_j)
	D := make([][]string, len_j)

	// create the inner arrays (j dimension = position along query)

	for index := range H {

		H[index] = make([]int, len_i)
		J[index] = make([]int, len_i)
		I[index] = make([]int, len_i)
		D[index] = make([]string, len_i)

	}

	match := func(i int, j int) MatrixMovement {

		var return_value MatrixMovement

		if string(subject[i-1]) == string(query[j-1]) {
			/*if the position is a match*/
			return_value.Score = H[j-1][i-1] + 4
			return_value.Origin = "m"

		} else if string(query[j-1]) == "N" {
			/*if the base is undefined treat is a neutral*/
			return_value.Score = H[j-1][i-1]
			return_value.Origin = "n"

		} else {
			/*otherwise it is a mismatch*/
			return_value.Score = H[j-1][i-1] - 5
			return_value.Origin = "x"
		}
		return (return_value)
	}

	traceback := func(start MatrixPosition) MatrixPosition {
		current_vector := D[start.j][start.i]
		var next_i int
		var next_j int
		switch {
		case string(current_vector) == "m":
			next_i = start.i - 1
			next_j = start.j - 1
		case string(current_vector) == "x":
			next_i = start.i - 1
			next_j = start.j - 1
		case string(current_vector) == "n":
			next_i = start.i - 1
			next_j = start.j - 1
		case string(current_vector) == "i":
			next_i = start.i - 1
			next_j = start.j
		case string(current_vector) == "j":
			next_i = start.i
			next_j = start.j - 1
		}

		next_position := MatrixPosition{i: next_i, j: next_j}

		return (next_position)
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
				best_move := max([]MatrixMovement{
					match(i, j),
					MatrixMovement{I[j][i], "i"},
					MatrixMovement{J[j][i], "j"},
				})
				H[j][i] = best_move.Score
				D[j][i] = best_move.Origin
				//fmt.Println(best_move.Score, best_move.Origin)

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
	var max_score int
	var max_position MatrixPosition

	for j := 0; j < len_j; j++ {
		if H[j][len_i-1] > max_score {
			max_score = H[j][len_i-1]
			max_position.i = len_i - 1
			max_position.j = j
		}
	}

	for i := 0; i < len_i; i++ {
		if H[len_j-1][i] > max_score {
			max_score = H[len_j-1][i]
			max_position.i = i
			max_position.j = len_j - 1
		}
	}

	//fmt.Println(string("max position"),max_position.i, max_position.j)

	//build reverse cigar string

	var current_position MatrixPosition

	var rev_CIGAR []string

	completed_traceback := false

	current_position = max_position

	for completed_traceback == false {

		current_vector := D[current_position.j][current_position.i]

		rev_CIGAR = append(rev_CIGAR, current_vector)

		next_position := traceback(current_position)

		if next_position.i == 0 || next_position.j == 0 {
			completed_traceback = true
		} else {
			current_position = next_position
		}
	}
	//fmt.Println("current position", current_position)
	//fmt.Println(rev_CIGAR)

	// create an forward cigar

	var CIGAR []string

	for i := 0; i < len(rev_CIGAR); i++ {
		next_letter := rev_CIGAR[len(rev_CIGAR)-1-i]
		CIGAR = append(CIGAR, next_letter)
	}
	fmt.Println(CIGAR)

	// TODO: create new alignment object

	// cigar start = current_position
	// cigar end = max_position

	newAlignment := PairWiseAlignment{
		Subject:       subject,
		Query:         query,
		ExpandedCIGAR: CIGAR,
		SubjectStart:  current_position.i - 1,
		QueryStart:    current_position.j - 1,
	}
	newAlignment = AlignmentRepr(newAlignment)

	// TODO: create print method for alignment object

	// TODO: resolve same scores on traceback

	// return the new alignment object
	return (newAlignment)
}

func Testalign() {

	subject := "GTGTCAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTG"
	query := "GCTAGGGAGGACGATGCGGTGGTGATGCTGCCACATACACTAAGAAGGTCCTGGACGCGTGTAGTCACTTCCAGCGGTCGTATGCCGTCTTCTGCTTGAA"

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
