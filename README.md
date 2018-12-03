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

1. Total number of nucleotides
```sh
ls

# Need to first partition the genomes to ≤ 100kb and > 100kb
# First for ≤ 100kb

faFilter -maxSize=100000 dmel-all-chromosome-r6.24.fasta \
> faSize dmel_less.fasta

# For > 100kb
faFilter -minSize=100001 dmel-all-chromosome-r6.24.fasta \
> faSize dmel_more.fasta
ls


# Trying another way for partition (will only use faSize partition version for remaining homework):
bioawk -c fastx 'length($seq) <= 100000{ print ">"$name; print $seq }'  *.fasta \
> sample_less 
bioawk -c fastx 'length($seq) > 100000{ print ">"$name; print $seq }'  *.fasta \
> sample_more

# To calculate total nucleotides, Ns, and sequences
faSize dmel_less.fasta
faSize dmel_more.fasta

# Another way to solve for total number of nucleotides
bioawk -c fastx '{ print $name, length($seq) }' dmel_less.fasta
bioawk -c fastx '{ print $name, length($seq) }' dmel_more.fasta

# Another way to solve for total number of Ns
bioawk -t -c fastx 'END {print NR}' dmel_less.fasta
bioawk -t -c fastx 'END {print NR}' dmel_more.fasta

```
**For ≤ 100 kb**
1. Total nucleotides: 6178042
2. Total Ns: 662593
3. Total sequences: 1863

**For > 100 kb**
1. Total nucleotides: 137547960
2. Total Ns: 490385
3. Total sequences: 7


## Plots of the following for the whole genome, for all sequences ≤ 100kb, and all sequences > 100kb:
Hint: bioawk has a function called gc(). Don't forget about the CDF plotting utility we used in class.
    Sequence length distribution
    Sequence GC% distribution
    Cumulative genome size sorted from largest to smallest sequences

Because the calculations will be for the whole genome and two genome partitions, there will be **9 total plots**.

```sh
# sequence length distribution ## maybe??
bioawk -c fastx '{ print $name, length($seq) }' input.fasta

# GC% Distribution
bioawk -c fastx '{ print $name, gc($seq) }' dmel_less.fasta
bioawk -c fastx '{ print $name, gc($seq) }' dmel_more.fasta
bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.24.fasta

# Cumulative genome size largest to smallest -rn?

```
## Genome Assembly
**Note:** This part of homework 4 is still being arranged. When this note is gone, it should be ready.

# Assemble a genome from MinION reads
**Hint:** Read up on miniasm here. We're using one of the simplest assembly approaches possible. This assembly can literally be accomplished with three lines of code. This will literally take only 3 command lines.
Download the reads from here
Use minimap to overlap reads
Use miniasm to construct an assembly

```sh
pwd
mkdir assembly
cd assembly

### THIS PART TAKES A RIDICULOUSLY LONG TIME
module load jje/jjeutils
wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz 
# Trying script from class
minimap=$(which minimap)
miniasm=$(which miniasm)
basedir=/pub/jje/ee282/$USER
projname=nanopore_assembly1
basedir=$basedir/$projname
raw=$basedir/$projname/data/raw
processed=$basedir/$projname/data/processed
figures=$basedir/$projname/output/figures
reports=$basedir/$projname/output/reports

createProject $projname $basedir
ln -sf /bio/share/solarese/hw4/rawdata/iso1_onp_a2_1kb.fastq $raw/reads.fq

$minimap -t 32 -Sw5 -L100 -m0 $raw/reads.fq{,} \
| gzip -1 \
> $processed/onp.paf.gz

$miniasm -f $raw/reads.fq $processed/onp.paf.gz \
> $processed/reads.gfa

```

## Assembly assessment
**Hint:** For MUMmer, you should run nucmer, delta-filter, and mummerplot.
1. Calculate the N50 of your assembly (this can be done with only faSize+awk+sort or with bioawk+awk+sort) and compare it to the Drosophila community reference's contig N50 (here)

```sh
# For calculating n50
n50 () {
  bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
  | sort -rn \
  | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
}

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' $processed/reads.gfa \
| tee >(n50 /dev/stdin > $reports/n50.txt) \
| fold -w 60 \
> $processed/unitigs.fa

```
2. Compare your assembly to the contig assembly (not the scaffold assembly!) from Drosophila melanogaster on FlyBase using a dotplot constructed with MUMmer (Hint: use faSplitByN as demonstrated in class)

```sh
# Need to first make a contig assembly

module load jje/jjeutils perl
faSplitByN dmel-all-chromosome-r6.24.fasta dmel-all-chromosome-cntg-r6.24.fasta 10

# To make the dotplot
source /pub/jje/ee282/bin/.qmbashrc
module load gnuplot

## Sample
mummer -mum -b -c ref.fasta qry.fasta > ref_qry.mums
mummerplot --postscript --prefix=ref_qry ref_qry.mums
gnuplot ref_qry.gp

## Trying; Need to be in the right working directory (/pub/jje/ee282/$USER/hmwk4
mummer -mum -b -c dmel-all-chromosome-cntg-r6.24.fasta  /pub/jje/ee282/bsorouri/nanopore_assembly1/nanopore_assembly1/data/processedreads.gfa > ref_qry.mums

### Trying with edwin
REF="dmel-all-chromosome-r6.24.fasta"
PREFIX="flybase"
SGE_TASK_ID=1
QRY=$(ls pub/jje/ee282/bsorouri/nanopore_assembly1/nanopore_assembly1/data/processed/unitigs.gfa | head -n $SGE_TASK_ID | tail -n 1)
PREFIX=${PREFIX}_$(basename ${QRY} .fa)

```

3. Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot (Hint: use plotCDF2 as demonstrated in class and see this example)

```sh

```

4. Calculate BUSCO scores of both assemblies and compare them

```sh
# First load BUSCO
module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc

INPUTTYPE="geno"
MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
##OPTIONS="${OPTIONS} -sp 4577"
QRY="unitigs.fasta"
###Please change this based on your qry file. I.e. .fasta or .fa or .gfa
MYEXT=".fasta" 

#my busco run
#you can change the value after -c to tell busco how many cores to run on. Here we are using only 1 core.
BUSCO.py -c 1 -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}


```
