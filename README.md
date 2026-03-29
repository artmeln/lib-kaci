# Description

This is a library for estimating decoding probabilities from a **query nucleotide sequence** based on a set of **reference amino acid kmers**.

The library is written in C and is meant to be used via the provided Python wrapper. The full list of functions and their descriptions can be found in **kaci.h**

# Installation

The library can be built on Linux (or WSL).

1. Install the dependencies for the library:
```
sudo apt update
sudo apt install build-essential
sudo apt install zlib1g-dev
```

2. Clone this repo.
```
git clone https://github.com/artmeln/lib-kaci.git
```

3. Install python dependencies (these are used only by the tests and examples, not by the library itself):
```
pip install -r requirements.txt
```

4. Build the library:
```
cd src
make libkaci
```

# Getting started

The best way to start is to run the provided examples (both examples will take a few seconds to execute). First build a reference:

```
cd examples
python3 build_reference.py
```

To process queries run:

```
mkdir results
wget -P genomes/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
python3 process_query.py
```

This will produce typical KACI output files for several genomes located in **test** directory.

You should look at these examples to familarize yourself with the inputs that control KACI. 

## Inputs

The following inputs should be given to the library:

* amino acid kmer length

* translation table - a string of single letter amino acid decodings for the 64 possible codons. The order of decodings in the table is based on the ascending ordering of codons where T<C<A<G. The translation table for the standard genetic code is `FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG`

* maximum query count - how many times an amino acid kmer that is found repeatedly in the query and matches a reference kmer will be used in estimating decoding probabilities

* the name of the file containing query nucleotide sequences in fasta format (*.fna or *.fna.gz) OR a list of such file names when working in the batch mode

* the name of the file containing the reference kmers (the format of this file is described below)

If the reference file is not available, it should be constructed first. The library settings that influence this process are

* maximum allowed count for any single amino acid in a family reference kmer

* minimum allowed number of proteins in a family

* minimum allowed link number - the number of kmers in one family that are linked by a single amino acid substitution at a specific position

* name of the file containing representative sequences of protein domains for all available protein families in fasta format. One example of such file is **Pfam-A.fasta.gz** which accompanies releases of Pfam domain database.

## Outputs

* a json file containing the decoding probabilities for every codon and the number of times each codon contributed to the decoding probability calculation

* optionally a csv file listing all query kmers that contributed to the calculation

## Tests

If you make modifications to the library, you may want to check that everything is still working as expected. To run the tests:

1. Download genomes from NCBI:
```
wget -P test/genomes/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
wget -P test/genomes/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/661/245/GCA_001661245.1_Pacta1_2/GCA_001661245.1_Pacta1_2_genomic.fna.gz
wget -P test/genomes/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/965/GCF_000182965.3_ASM18296v3/GCF_000182965.3_ASM18296v3_genomic.fna.gz
```

2. Run the tests:
```
python3 tests.py
```
**Note:** the last test needs [seqkit](https://bioinf.shenwei.me/seqkit/)


# Acknowledgments
This project was made possible by a [very fast hash table C library](https://github.com/attractivechaos/khashl)

In addition, it utilizes a C implementation of [logaddexp](https://github.com/horta/logaddexp)


# Citation

Artem V Melnykov. "New genetic codes in bacteria and archaea identified with a fast k-mer based algorithm".

If you want to reproduce results from this paper, follow these steps:

1. Make sure you are running on a machine with at least 16 GB of available memory.

2. Build the library as described above.

3. Download the [reference](https://zenodo.org/records/19318166), unzip it and place into `examples/ref/`:
```
cd examples
wget -P ref/ https://zenodo.org/records/19318166/files/ref_k11_link20.tar.gz?download=1
tar -xvzf ref/ref_k11_link20.tar.gz ref/ref_k11_link20/
```

4. Download genomes of interest and place them into `examples/genomes/` (you don't have to unzip them).

5. To process all files located in `examples/genomes/` use `1_batch_process_genomes.py`. You may want to edit this file and change the number of threads to be used.
```
python3 1_batch_process_genomes.py
```

6. To extract the inferrences:
```
python3 2_evaluate_results.py
```

