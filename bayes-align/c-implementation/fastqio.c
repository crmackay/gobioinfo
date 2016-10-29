#include <string.h>
#include <stdio.h>

// fastq io

// open file

// create fastq reads

typedef struct sb_fastq_read
{
    char seq[100];
    int seq_length;
    char quals[100];
    int quals_length;
    char id[100];
    int id_length;
    char misc[100];
    int misc_length;
    int cleaned_start;
    int cleaned_end;
} sb_read;

sb_read new_sb_read(char* sequence, char* qualities, char* id, char* misc)
{
    sb_read new_read;

    int seq_length = strlen(sequence);

    int quals_length = strlen(qualities);

    for (int i = 0; i < seq_length ; i++)
    {
        new_read.seq[i] = sequence[i];
    }

    for (int i = 0; i < quals_length ; i++)
    {
        new_read.quals[i] = qualities[i];
    }

    new_read.seq_length = seq_length;
    new_read.quals_length = quals_length;
    new_read.cleaned_start = 0;
    new_read.cleaned_end = seq_length;

    return(new_read);
}

// next read returns error if there's an error
int next_read_from_file(FILE *fp, sb_read *next_read)
{
    int error = ferror(fp);
    int eof = feof(fp);
    char myChar = fgetc(fp);
    printf("3: %c\n", myChar);
    return(0);
}

// close file
int close (FILE *fp)
{
    fclose(fp);
    return(0);
}

int main(void)
{

    FILE *my_fp;

    my_fp = fopen("sample.fastq", "r");

    sb_read *my_read;

    int err = next_read_from_file(my_fp, my_read);

    // printf("%s", my_read->seq);
    //
    // err = close(my_fp);


}
