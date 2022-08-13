# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | +33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | +33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | +33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | +33 |


2. Per-base NT distribution
    1.  Use markdown to insert your 4 histograms here.
        ![read_1](https://github.com/matt-esqueda/Demultiplex/blob/bfa12ce98800595b094e8121482323b3c91375b4/Assignment-the-first/Part_1/read1.png)
        ![read_2](https://github.com/matt-esqueda/Demultiplex/blob/bfa12ce98800595b094e8121482323b3c91375b4/Assignment-the-first/Part_1/read2.png)
        ![index_1](https://github.com/matt-esqueda/Demultiplex/blob/bfa12ce98800595b094e8121482323b3c91375b4/Assignment-the-first/Part_1/index1.png)
        ![index_2](https://github.com/matt-esqueda/Demultiplex/blob/bfa12ce98800595b094e8121482323b3c91375b4/Assignment-the-first/Part_1/index2.png)
    
    2.  Index quality score cutoff: 35 
        Based off the graph this score will exclude the low quality indexs
        Read quality score cutoff: 30 
        Based off the graph this score will exclude the low quality reads
    
    3.  zcat 1294_S1_L008_R2_001.fastq.gz | grep -A1 "^@" | grep -v "^--" | grep -v "^@" | grep "N" | wc -l
            3976613 N base calls in index 1

        zcat 1294_S1_L008_R3_001.fastq.gz | grep -A1 "^@" | grep -v "^--" | grep -v "^@" | grep "N" | wc -l
            3328051 N base calls in index 2
    
## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
