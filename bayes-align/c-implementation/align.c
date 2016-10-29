#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <sys/time.h>
#include <sys/resource.h>


typedef enum {
	DIAG_MATCH,			// diagonal move because of a MATCH
	DIAG_MISMATCH,		// diagonal move because of a mismatch
	DIAG_N,				// diagonal move because of the query string has an "N"
	INS_I,				// insertion in the "i" direction
	INS_J,				// insertion in the "j" direction
	EDGE				// the edge of the matrix
} matrix_move;

// RT enzyme error plus the DNA polyermase enzyme error (per base) time the number of cycles of PCR
const double probPCRError = 0.00001490711984999862 + (0.00001038461538461538 * 30.0);

// takes a PHRED+33 encoded quality score from the sequencer and returns the probability of a
// sequencer miscall, as a double
double phredToProbIncorrect(char p)
{
	int phred = (int)p - 33;
	// printf("pred int for %c: %i\n", p, phred);
	double prob = pow(10, (-((double)phred) / 10));
	// printf("pred prob of miscall for %c: %f\n", p, prob);
	return prob;
}

int max_match(double values[3])
{

	// // printf("0: %f\n", values[0]);
	// // printf("1: %f\n", values[1]);
	// // printf("2: %f\n", values[2]);

	int maxIndex = 0;

    for (int i = 1; i < 3; i++)
    {
        if (values[i] > values[maxIndex])
        {
            maxIndex = i;
        }
    }

	// // printf("max: %i\n", maxIndex);

    return maxIndex;
}


matrix_move match(char s, char q, char q_qual)
{
	matrix_move new_move;
	if (q == s)
	{
		new_move = DIAG_MATCH;
	}
	else if (q =='N')
	{
		new_move = DIAG_N;
	}
	else
	{
		new_move = DIAG_MISMATCH;
	}

	return(new_move);
}

double probSLMatch(char s, char q, char q_qual)
{
	double probMisscall = phredToProbIncorrect(q_qual);

	double probCorrcall = 1.0 - probMisscall;

	double value = (probPCRError * probMisscall) + (probCorrcall * (1.0 - probPCRError));

	// // printf("prob SL for %c, %c, %c, is: %2.10f\n", s, q, q_qual, value);
	return(value);
}


double probSLMismatch(char s, char q, char q_qual)
{
	double probMisscall = phredToProbIncorrect(q_qual);

	double probCorrcall = 1.0 - probMisscall;

	double value = (1.0 / 3.0 * probMisscall * (1.0 - probPCRError)) + (1.0 / 3.0 * probCorrcall * probPCRError) + (2.0 / 9.0 * probMisscall * probPCRError);

	// // printf("prob SLM for %c, %c, %c, is: %2.10f\n", s, q, q_qual, value);

	return value;
}

void printMatrix(int rows, int cols, double m[rows][cols])
{
	for (int i = 0; i<rows; i++)
	{
		for (int j = 0; j<cols; j++)
		{
			printf("%2.2e\t", m[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void printMatrixInt(int rows, int cols, int m[rows][cols])
{
	for (int i = 0; i<rows; i++)
	{
		for (int j = 0; j<cols; j++)
		{
			printf("%i\t", m[i][j]);
		}
		printf("\n");
	}
	printf("\n");

}
// refactor code a bit:
// make alignment matrices, traceback, test/trim, all seperate...testable things


int align(char* subj, char* query, char* query_quals)
{

    int max_i = strlen(subj)+1;
    int max_j = strlen(query)+1;

    double H[max_i][max_j];
	int D[max_i][max_j];

    double SL[max_i][max_j];
    double SLP[max_i][max_j];

	for (int i = 0; i< max_i; i++)
	{
		for (int j = 0; j < max_j; j++)
		{
			if (i == 0 || j == 0)
			{
				H[i][j] = 0;
				SL[i][j] = 1;
				SLP[i][j] = 1;
				D[i][j] = EDGE;
			}
			else
			{
				matrix_move diag_move = match(subj[i-1], query[j-1], query_quals[j-1]);

				double diag_value = 0;

				double diag_PSL = 0;

				double diag_PSLP = 0;
				double next_factor;
				switch (diag_move)
				{
					case DIAG_MATCH:
						diag_PSL = probSLMatch(subj[i-1], query[j-1], query_quals[j-1]);
						diag_PSLP = 0.25;

						next_factor = log10(diag_PSL / diag_PSLP);

						diag_value = H[i-1][j-1] + next_factor;
						break;

				 	 case DIAG_MISMATCH:
					 	diag_PSL = probSLMismatch(subj[i-1], query[j-1], query_quals[j-1]);
						diag_PSLP = 0.75;

						next_factor = log10(diag_PSL / diag_PSLP);

						diag_value = H[i-1][j-1] + next_factor;
						break;

					case DIAG_N:
						diag_PSL = 1;
						diag_PSLP = 1;
						diag_value = 1;
						break;

					default:
						break;
				}

				// // printf("diag_value: %f\n", diag_value);

				// // printf("diag_PSL: %f\n", diag_PSL);

				// // printf("diag_PSLP: %f\n", diag_PSLP);

				double gap_i = H[i-1][j] + log10(probPCRError);

				double gap_j = H[i][j-1] + log10(probPCRError);

				int best_move = max_match( (double[]){diag_value, gap_i, gap_j});
				switch (best_move)
				{
					case 0:
						H[i][j] = diag_value;
						D[i][j] = diag_move;
						SL[i][j] = SL[i-1][j-1] * diag_PSL;
						SLP[i][j] = SLP[i-1][j-1] * diag_PSLP;

						break;

					case 1:
						H[i][j] = gap_i;
						D[i][j] = INS_I;
						SL[i][j] = SL[i-1][j] * probPCRError;
						SLP[i][j] = SLP[i-1][j];

						break;

					case 2:
						H[i][j] = gap_j;
						D[i][j] = INS_J;
						SL[i][j] = SL[i][j-1] * probPCRError;
						SLP[i][j] = SLP[i][j-1];
						break;
				}
			}
		}
	}

	printMatrix(max_i, max_j, H);

	printMatrix(max_i, max_j, SL);

	printMatrix(max_i, max_j, SLP);

    printf("dir matrix:\n");



	printMatrixInt(max_i, max_j, D);

	// find the highest value in the last column or last row of matrix H
	int highest_score_i = max_i-1;
	int highest_score_j = max_j-1;
	double highest_score = H[highest_score_i][highest_score_j];
	for (int i = 1; i < max_i; i++)
	{
		if (H[i][max_j-1] > highest_score)
		{
			highest_score = H[i][max_j-1];
			highest_score_i = i;
			highest_score_j = max_j-1;
		}
	}
	for (int j = 1; j < max_j-1; j++)
	{
		if (H[max_i-1][j] > highest_score)
		{
			highest_score = H[max_i-1][j];
			highest_score_i = max_i-1;
			highest_score_j = j;
		}
	}
	// // printf("here\n");

	int tback_pos_i = highest_score_i;
	int tback_pos_j = highest_score_j;
	matrix_move curr_dir = D[highest_score_i][highest_score_j];
	matrix_move rev_cigar[max_i + max_j];
	int rev_cigar_pos = 0;
	rev_cigar[rev_cigar_pos] = curr_dir;
	// while both traceback cursors are greater than zero (ie we are not at an edge yet)
	while (curr_dir != EDGE)
	{
		// // printf("tback_pos_i: %i \ntback_pos_j: %i\n\n", tback_pos_i, tback_pos_j);
		switch(curr_dir)
		{
			case DIAG_MATCH:			// diagonal move because of a MATCH
				tback_pos_i = tback_pos_i - 1;
				tback_pos_j = tback_pos_j - 1;
				break;

			case DIAG_MISMATCH:			// diagonal move because of a mismatch
				tback_pos_i = tback_pos_i - 1;
				tback_pos_j = tback_pos_j - 1;
				break;

			case DIAG_N:				// diagonal move because of the query string has an "N"
				tback_pos_i = tback_pos_i - 1;
				tback_pos_j = tback_pos_j - 1;
				break;

			case INS_I:					// insertion in the "i" direction
				tback_pos_i = tback_pos_i - 1;
				break;

			case INS_J:					// insertion in the "j" direction
				tback_pos_j = tback_pos_j - 1;
				break;

			case EDGE:
				// this should never run
				break;
		}
		curr_dir = D[tback_pos_i][tback_pos_j];
		rev_cigar_pos ++;
		rev_cigar[rev_cigar_pos] = curr_dir;

	}

	matrix_move cigar[rev_cigar_pos+1];

	// printf("rev-cigar_pos: %i\n", rev_cigar_pos);

	char* cigar_str = malloc((rev_cigar_pos+1)*sizeof(char));
	// error check on the malloc call:
	if (cigar_str == NULL){printf("there was an error allocating memory\n"); return(1);}

	// printf("here\n");
	// invert the reverse cigar into a proper cigar string
	for (int i = 0; i < rev_cigar_pos; i++)
	{
		cigar[i] = rev_cigar[rev_cigar_pos-1-i];
		// printf("rev cigar val no: %i is %i\n", i, rev_cigar[rev_cigar_pos-1-i]);
		cigar_str[i] = (char) rev_cigar[rev_cigar_pos - 1 - i]+65;
	}
	cigar_str[rev_cigar_pos] = '\0';

	// printf("%s\n", cigar_str);

	// TODO: alignment representation




	// TODO: bayesian prob test

	// calc P(L|S) and P(L'|S), for the alignment, print out true or false


	// TODO: make more robust so that a read and linker can be inputted, they are processed, and
	// then the clean linker is returned...

	// TODO: make a verbose version that prints to stdout data re: cutting

	// TODO: deal with conflicts on matrix creation? two dirs with same value?



	return(0);
}


double get_time(void)
{
	struct timeval t;
	struct timezone tzp;
	gettimeofday(&t, &tzp);
	return t.tv_sec + t.tv_usec*1e-6;
}


int main(int argc, char *argv[])
{



	double start = get_time();
	// printf("start:\t%1.5f\n", start);
	char *subject = argv[1];
	//char* subject = "GTGTCAG";
	char* query = "GTGTCAGT";
	char* quals = "DDDDDFFF";

	int err  = align(subject, query, quals);
	if (err != 0){return(1);}

	double end = get_time();
	printf("end:\t%1.5f\n", end-start);

	return(0);

}
