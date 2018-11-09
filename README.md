# Homework_4

Due December 7, 12pm.
Homework must be uploaded to a new GitHub repository you create specifically for Homework 4. Submission of the link to the repository suffices for submitting the homework. Please don't submit real paper. Please use either markdown a file or a mix of a markdown file and scripts to answer your homework. One markdown file should serve as the "master file" describing your answers and which files have the code that demonstrate your answers (if necessary). Remember, markdown itself can display code snippets. The timestamp of your files will serve to determine whether your assignment is on time.

## Summarize partitions of a genome assembly
We will be revisiting the Drosophila melanogaster genome. As with Homework 3, start at flybase.org. Go to the most current download genomes section and download the gzipped fasta file for all chromosomes.
Hint: The partitioning can be accomplished in many different ways. In my opinion, the easiest way is by using bioawk and faSize. The bioawk tool can be found in the module jje/jjeutils and the fa* utilities can be found in the module jje/kent.

```sh
qrsh
cd /pub/jje/ee282/${USER}
mkdir hmwk4
cd hmwk4
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz                         
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/gff/dmel-all-r6.24.gff.gz

# Unzip
ls 
gunzip *.gz

# download bioawk and faSize
module load jje/jjeutils/0.1a
module load jje/kent/2014.02.19
```


### Calculate the following for all sequences ≤ 100kb and all sequences > 100kb:
1. Total number of nucleotides
2. Total number of Ns
3. Total number of sequences
Because the calculations will be for two genome partitions, there will be **6 total responses**.

```sh

```

Plots of the following for the whole genome, for all sequences ≤ 100kb, and all sequences > 100kb:
Hint: bioawk has a function called gc(). Don't forget about the CDF plotting utility we used in class.
    Sequence length distribution
    Sequence GC% distribution
    Cumulative genome size sorted from largest to smallest sequences
Because the calculations will be for the whole genome and two genome partitions, there will be **9 total plots**.

```sh

```
## Genome Assembly
