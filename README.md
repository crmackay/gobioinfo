# gobioinfo
my golang bioinformatics package

version 0.0.-1

**this is all experimental code** which is not all totally *idiomatic* go. I am learning as I go along, so take this as you will.

I would recommend the great [**Biogo** package](http://github.com/biogo) if you are looking for a more mature bioinformatics package (although they don't have a semiglobal alignment algorithm), which is why I struck out on my own.

## Components so far:

- semiglobal alignment algorithm (and pairwise alignment struct)
- a FASTQ scanner structure for scanning a FASTQ file read by read

## To Be Added

- FASTQ Reader and Writer that work with unix pipes (eg are general io.Reader and io.Writers)
- better sequence structs
