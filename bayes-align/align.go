package main

import (
	"fmt"
	"math"
)

var probPCRError = 0.00001490711984999862 + 0.00001038461538461538*30

func phredToProbIncorrect(p byte) float64 {
	phred := int(p) + 33
	prob := math.Pow(10, (-float64(phred) / 10))
	return prob
}

func max(d, i, j float64) (float64, int) {
	max := d
	dir := 1

	if i > max {
		max = i
		dir = 2
	}
	if j > max {
		max = j
		dir = 3
	}
	// fmt.Println(m)
	// fmt.Println(max)
	return max, dir
}

func match(s, q, qual byte) float64 {
	probMisscall := phredToProbIncorrect(qual)

	probCorrcall := 1 - probMisscall

	var value float64
	if s == q {
		// fmt.Println(probPCRError * probMisscall)
		value = (probPCRError * probMisscall) + (probCorrcall * (1 - probPCRError))
	} else {
		value = (float64(1) / 3 * (1 - probMisscall)) + (float64(1) / 3 * probCorrcall * probPCRError) + (float64(2) / 9 * probMisscall * probPCRError)
	}
	return value
}

func chance(s, q, qual byte) float64 {
	if s == q {
		return float64(1) / 4
	} else if s != q && q != 'N' {
		return float64(3) / 4
	} else {
		return 1
	}
}

var subject = []byte("GTGTCAGTCACTTCC")

var query = []byte("NTGTGCANNNNCTTCC")

// var query = []byte("GTGTCAGTCACTTCC")

var quals = []byte("@#-%@EC####-,;BC")

func align(subj, query, quals []byte) [][]float64 {

	matrix := make([][]float64, len(subj)+1)

	matrixSL := make([][]float64, len(subj)+1)

	matrixSLP := make([][]float64, len(subj)+1)

	matrixDir := make([][]int, len(subj)+1)

	for i := range matrix {
		matrix[i] = make([]float64, len(query)+1)
		matrixSL[i] = make([]float64, len(query)+1)
		matrixSLP[i] = make([]float64, len(query)+1)
		matrixDir[i] = make([]int, len(query)+1)
	}

	for i := range matrix {
		for j := range matrix[i] {
			if i == 0 || j == 0 {
				matrix[i][j] = 1

				// matrixSL[i][j] =
				// matrixSLP[i][j] =
				// matrixDir[i][j] =

			} else if query[j-1] == byte('N') {
				matrix[i][j] = matrix[i-1][j-1]

				// matrixSL[i][j] =
				// matrixSLP[i][j] =
				// matrixDir[i][j] =

			} else {
				probSC, dir := max(
					matrix[i-1][j-1]*match(query[j-1], subj[i-1], quals[j-1]),
					matrix[i-1][j]*probPCRError,
					matrix[i][j-1]*probPCRError,
				)
				if dir == 1 {

				}
				probSCP := chance(query[j-1], subj[i-1], quals[j-1])

				// TODO: this is wrong: we need to store each of these values individually... and
				// recalculated at each one... the alignment matrix could be probSC/probSCP
				matrix[i][j] = probSC / probSCP

				// matrixSL[i][j] =
				// matrixSLP[i][j] =
				// matrixDir[i][j] =

			}
		}
	}

	// find highest SL/SLP ratio

	// traceback

	// return alignment representation, and P(S|L), and P(S|L')

	return matrix

}

//
// func isAlignmentLinker(a alignment) bool {

// take alignment, calculate P(L|S) and P(L'|S) and compare, return true/false
// }

func printMatrix(m [][]float64) {
	for i := range m {
		for j := range m[i] {
			fmt.Printf("%2.3f\t", m[i][j])
		}
		fmt.Printf("\n")
	}
	fmt.Printf("\n")

}

func main() {
	m := align(subject, query, quals)
	printMatrix(m)
}
