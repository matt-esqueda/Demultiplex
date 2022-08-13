Demultiplexing and Index Swapping Assignment the first

	Part 1: Quality Score Distribution-per-nucleotide

	1. initial data exploration
	
	i.	bash commands used:
			zcat 1294_S1_L008_R1_001.fastq.gz | head
				- used on all four files (R1, R2, R3, R4)
				- determined that R1 = read 1, R2 = read 2
				- determined that R3 = index 1, R4 = index 2
			

	ii. 	zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
				- used on all four files 
				- determined that the read length for the read files is 101 nt
				- determined that the read length for the index files is 8 


	iii.	zcat 1294_S1_L008_R2_001.fastq.gz | head -1000000 | tail -20 > /projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/test_index.fq.gz 
				- used to create test index fastq files
			
			zcat 1294_S1_L008_R1_001.fastq.gz | head -1000000 | tail -20 > /projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/test_read.fq.gz 
				- used to create test read fastq files
				

	2. Used code from PS4 part 1 to create part_1.py python code to determine the average quality scores at each postion for all reads and indexes. 
	Generated a per nucleotide mean distribution for each of the four files:

		import argparse
		import gzip
		import bioinfo 

		python script:
		part_1.py = /projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/part_1.py

		slurm script:
		part_1.sh = /projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/part_1.sh

		1294_S1_L008_R1_001.fastq.gz
		1294_S1_L008_R2_001.fastq.gz
		1294_S1_L008_R3_001.fastq.gz
		1294_S1_L008_R4_001.fastq.gz
		
		- all located in /projects/bgmp/shared/2017_sequencing/

		ii. Good quality scores for indexes and reads...

		Used these commands to find how many indexes have undetermined (N) base calls:

		zcat 1294_S1_L008_R2_001.fastq.gz | grep -A1 "^@" | grep -v "^--" | grep -v "^@" | grep "N" | wc -l	
		3976613

		zcat 1294_S1_L008_R3_001.fastq.gz | grep -A1 "^@" | grep -v "^--" | grep -v "^@" | grep "N" | wc -l
		3328051	


	Part 2: Develop an algorithm to de-multiplex the samples

		- Wrote psuedocode to develop algorithm for de-multiplexing files and reporting the amount of index hopping. 
			/projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-first/Part_2/psuedocode.txt



Demultiplexing and Index Swapping Assignment the Second
	
	Reviewed three peers psuedocode and provided feedbach through github.

	Gained some insight regarding how to make my own code more efficient. 


Demultiplexing and Index Swapping Assignment the Third

	Used psuedocode from part 1 to write python script to demultiplex the samples. 

	*Look up itertools, permutations, for creating mismatched indexes dictionary
		https://docs.python.org/3/library/itertools.html (very useful)

	Biggest challenges were opening all four read files simultaneously and writing to all the approiate fastq files.
		-used gzip to open zipped files, good to know for the future

	Added reverse_complement() to bioinfo.py

	Figured out the general logic for separating the indexes appropiately, and they code seems to be working for the test files
		-hopped and unknown put into their own folders

	bash command for test:
	./demultiplex.py -i test_indexes  -t1 R1_test_input.fq  -t2 R2_test_input.fq  -t3 R3_test_input.fq  -t4 R4_test_input.fq -ic 30
		-will hardcode the actual files after finished testing, rather than using argparse. Command becoming way too long

	bash command for interactive node:
	srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=2:00:00 --cpus-per-task=1 --pty bash
	

	python ./demultiplex.py -ic 30
		*set index cutoff to 30

	Ran with the read files, failed the first time failed because index_mismatch_dict was not set up properly
		-fixed dictionary, run again
	
	Run failed because incorrect variable used in readfile(), need to be careful to not mix up variables, better naming...

	Run failed becuase tsv variables got mixed up, 'pair' is not defnined. Easy fix, run again..

	Run worked, all files created properly and correct counts with appropiate tables and tsv created

	Used pigz to zip all fastq files (works fast!)
		pigz *fq 
	
	Tried creating graphs for stats output. The graphs were either uniformative or too many variables to be readable. Tried to create heatmap but was too difficult.
		-For the future, look to pandas and seaborn for creating heatmap