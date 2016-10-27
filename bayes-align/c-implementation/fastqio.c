#include <string.h>

// fastq io

// open file

// create fastq reads

typedef struct sb_fastq_read
{
    char seq[];
    int seq_length;
    char quals[];
    int quals_length;
    char id[];
    int id_length;
    char misc[];
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
        new.seq[i] = sequence[i];
    }

    for (int i = 0; i < qual_length ; i++)
    {
        new.quals[i] = qualities[i];
    }

    new_read.seq_length = seq_length;
    new_read.quals_length = quals_length;
    new_read.cleaned_start = 0;
    new_read.cleaned_end = seq_length;

    return(new_read);
}

// open file
int open(char* path, FILE *fp)
{
    FILE *fp;
    char buff[255];
    fp = fopen(path, "r");


}

// next read returns error if there's an error
int next_read_from_file(FIPE fp, sb_read *next_read)
{
    int ferror(fp);
    int feof(fp);
    fgetc(fp);
    printf("3: %s\n", buff );

}

// close file
int close (FILE *fp)
{
    fclose(fp);
}
