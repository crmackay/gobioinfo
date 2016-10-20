package main

import "fmt"

type seq []byte

type read struct {
	seq  seq
	qual []byte
}

type matrixElem struct {
	fromDir int
	probSL  float64
	probSLP float64
}

type matrix [][]matrixElem

func (m *matrix) print() {
	fmt.Println(m)
}

func calcI() float64 {
	return 0
}

func calcJ() float64 {
	return 0
}

func calcH() float64 {
	return 0
}

func align(s, q read) matrix {

	lenS := len(s.seq)
	lenQ := len(q.seq)
	lenI := lenS + 1
	lenJ := lenQ + 1

	var h matrix
	h = make([][]matrixElem, lenI)
	for i, _ := range h {
		h[i] = make([]matrixElem, lenJ)
	}

	h.print()

	// fill in matrix

	for i := 0; i < lenI; i++ {
		for j := 0; j < lenJ; j++ {
			h[i][j] = matrixElem{
				fromDir: 10,
				probSL:  0.025,
				probSLP: 0.0001,
			}
		}
	}

	return h
}

// hasLinker take a matrix, and checks it to see if there is an alignment that exists that has a probSL
// which is greater than a probSL. If yes, the traceback is performed on the matrix, and True, as well
// as the index in the sunject from which linker begins is returned
func findLinker(h matrix) (bool, int) {
	return true, 1
}

func processRead(r read, l seq) read {
	newRead := r
	return newRead
}

func main() {

	// get sequences

	// process sequences

	// save new sequences

}
