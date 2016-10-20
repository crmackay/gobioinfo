#include <stdio.h>
#include <math.h>
#include <string.h>

// RT enzyme error plus the DNA polyermase enzyme error (per base) time the number of cycles of PCR
const double probPCRError = 0.00001490711984999862 + (0.00001038461538461538 * 30);

// takes a PHRED+33 encoded quality score from the sequencer and returns the probability of a
// sequencer miscall, as a double
double phredToProbIncorrect(char p) {
	int phred = (int)p + 33;
	double prob = pow(10, (-((double)phred) / 10));
	return prob;
}

double max_match(double a, double b, double c) {
	double max = a;
	if (b > max) {
			max = b;
	}
    if (c > max){
        max = c;
    }

    return max;
}

double match(char s, char q, char qual) {
	double probMisscall = phredToProbIncorrect(qual);

	double probCorrcall = 1 - probMisscall;

	double value = 0;

	if (s == q) {
		value = (probPCRError * probMisscall) + (probCorrcall * (1 - probPCRError));
	} else {
		value = ((double)1 / 3 * (1 - probMisscall)) + ((double)1 / 3 * probCorrcall * probPCRError) + ((double)2 / 9 * probMisscall * probPCRError);
	}
	return value;
}

double chance(char s, char q, char qual) {
	if (s == q) {
		return (double)1 / 4;
	} else if (s != q && q != 'N') {
		return (double)3 / 4;
	} else {
		return (double)1;
	}
}


void printMatrix(int rows, int cols, double m[rows][cols]) {
	for (int i = 0; i<rows; i++){
		for (int j = 0; j<cols; j++){
			printf("%2.3f\t", m[i][j]);
		}
		printf("\n");
	}
	printf("\n");

}


char* subject = "GTGTCAGTCACTTCC";

char* query = "NTGTGCANNNNCTTCC";

// var query = []byte("GTGTCAGTCACTTCC")

char* quals = "@#-%@EC####-,;BC";

void align(char* subj, char* query, char* quals ) {

    int max_i = strlen(subj)+1;
    int max_j = strlen(query)+1;

    double matrix[max_i][max_j];

    //fill the matrix

	for (int i = 0; i < max_i; i++) {
		for (int j = 0; j < max_j; j++){
			if (i == 0 || j == 0) {
				matrix[i][j] = (double) 1;
			} else if (query[j-1] == 'N') {
				matrix[i][j] = matrix[i-1][j-1];
			} else {
				double probSC = max_match(
					matrix[i-1][j-1]*match(query[j-1], subj[i-1], quals[j-1]),
					matrix[i-1][j]*probPCRError,
					matrix[i][j-1]*probPCRError
				);
				double probSCP = chance(query[j-1], subj[i-1], quals[j-1]);
				matrix[i][j] = probSC / probSCP;
			}
		}
	}

	printMatrix(max_i, max_j, matrix);

}

int main(void) {
	align(subject, query, quals);
}
