sega
====

Self-orgnizing approach for metagenome

# Installtion
	cd src; make sega
# Usage
## Data Reprocessing
`sega` has a command `seq2qu` to transform the input sequences in FastA or FastQ format to their quanternary representation and concatenates them into a single sequence.
#### Examples:
	bin/sega seq2qu test/ecoli_apec_o1.fasta test/ecoli_apec_o1.fasta.qu
## Training
#### Models
See [reference]
#### Examples
If we want to know if there are two different codon usage patterns in E.coli, the option -c should be specified to be 2. We can run 10 iterations to see whether it converges. The training sequence is the first 50000*99 bps.
1. initilize phase labelling
	bin/sega train test/ecoli_apec_o1.fasta.qu tmp_dir -c 2 -k 10 n 50000 -w 99
2. initilize frequency tables
	bin/sega train test/ecoli_apec_o1.fasta.qu tmp_dir -p test/init_cut -c 2 -k 10 n 50000 -w 99

[reference]: http://benthamscience.com/open/tobioij/articles/V006/28TOBIOIJ.pdf
