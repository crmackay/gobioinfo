#include <stdio.h>
#include <math.h>
#include <string.h>

#define MATCH 0
#define DIR_I 1
#define DIR_J 2
#define END   3

// RT enzyme error plus the DNA polyermase enzyme error (per base) time the number of cycles of PCR
const double probPCRError = 0.00001490711984999862 + (0.00001038461538461538 * 30.0);

// takes a PHRED+33 encoded quality score from the sequencer and returns the probability of a
// sequencer miscall, as a double
double phredToProbIncorrect(char p)
{
	int phred = (int)p + 33;
	double prob = pow(10, (-((double)phred) / 10));
	return prob;
}

int max_match(double values[3])
{
	int maxIndex = 0;
    for (int i = 1; i < 3; i++)
    {
        if (values[i] > values[maxIndex])
        {
            maxIndex = i;
        }
    }

    return maxIndex;
}

double match(char s, char q, char qual)
{
	double probMisscall = phredToProbIncorrect(qual);

	double probCorrcall = 1.0 - probMisscall;

	double value = 0;

    if (q == 'N')
    {
        value = 1;
    }
	else if (s == q)
	{
		value = (probPCRError * probMisscall) + (probCorrcall * (1.0 - probPCRError));
	}
	else
	{
		value = ((double)1 / 3 * (1.0 - probMisscall)) + ((double)1 / 3 * probCorrcall * probPCRError) + ((double)2 / 9 * probMisscall * probPCRError);
	}
    printf("value: %f\n", value);
	return value;
}

double chance(char s, char q, char qual)
{
	if (s == q)
	{
		return (double)1 / 4;
	}
	else if (s != q && q != 'N')
	{
		return (double)3 / 4;
	}
	else
	{
		return (double)1;
	}
}


void printMatrix(int rows, int cols, double m[rows][cols])
{
	for (int i = 0; i<rows; i++)
	{
		for (int j = 0; j<cols; j++)
		{
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

void align(char* subj, char* query, char* quals )
{

    int max_i = strlen(subj)+1;
    int max_j = strlen(query)+1;

    double matrixRatio[max_i][max_j];

    double matrixSL[max_i][max_j];
    double matrixSLP[max_i][max_j];
    double matrixDir[max_i][max_j];

    //fill the matrix

	for (int i = 0; i < max_i; i++)
	{
		for (int j = 0; j < max_j; j++)
		{
			if (i == 0 || j == 0)
			{
				matrixRatio[i][j] = 1.0;

                matrixSL[i][j] = 1.0;
                matrixSLP[i][j] = 1.0;
                matrixDir[i][j] = END;

			}
			else
			{
                double potentialMoves[3] =
                {
                    matrixSL[i-1][j-1]*match(query[j-1], subj[i-1], quals[j-1]),
                    matrixSL[i-1][j]*probPCRError,
                    matrixSL[i][j-1]*probPCRError
                };
                printf("%f, %f, %f\n", potentialMoves[0], potentialMoves[1], potentialMoves[2]);

				int bestMoveIndex = max_match(potentialMoves);

                int dir = 0;
                int lastDir_i = 0;
                int lastDir_j = 0;
                switch(bestMoveIndex)
                {
                    case 0:
                        dir = MATCH;
                        lastDir_i = i-1;
                        lastDir_j = j-1;
                        break;
                    case 1:
                        dir = DIR_I;
                        lastDir_i = i-1;
                        lastDir_j = j;
                        break;
                    case 2:
                        dir = DIR_J;
                        lastDir_i = i;
                        lastDir_j = j-1;
                        break;
                }

                double probSL = potentialMoves[bestMoveIndex];

				double probSLP = chance(query[lastDir_j], subj[lastDir_i], quals[lastDir_j])*matrixSLP[lastDir_i][lastDir_j];

				matrixRatio[i][j] = probSL / probSLP;

                matrixSL[i][j] = probSL;
                matrixSLP[i][j] = probSLP;
                matrixDir[i][j] = dir;
			}
		}
	}

	printMatrix(max_i, max_j, matrixRatio);
    printf("dir matrix:\n");
    printMatrix(max_i, max_j, matrixDir);

}

int main(void)
{
	align(subject, query, quals);
}
