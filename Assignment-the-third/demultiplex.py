#!/usr/bin/env python

import bioinfo
import argparse 
import matplotlib.pyplot as plt
import gzip
import itertools

def get_args():
    parser=argparse.ArgumentParser(description='Demultplex the input fastq files')
    parser.add_argument('-ic', help='index_cutoff_score')
    return parser.parse_args()

args=get_args()

f = args.f
ic = int(args.ic)


### Files and Folders

i = '/projects/bgmp/shared/2017_sequencing/indexes.txt'
r1 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz'
r2 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz'
r3 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz'
r4 = '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz'
u1 = '/projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-third/unknown_indexes/unknown.R1.fq'
u2 = '/projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-third/unknown_indexes/unknown.R2.fq'
h1 = '/projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-third/hopped_indexes/hopped.R1.fq'
h2 = '/projects/bgmp/mesqueda/bioinfo/Bi622/Demultiplex/Assignment-the-third/hopped_indexes/hopped.R2.fq'


### List of known indexes placed into index set

index_set = set()

with open(i, 'r') as index:
    for line in index:
        if line[0] == 's':
            continue      
        index_set.add(line.strip('\n').split('\t')[4])



### Dictionaries

# Dictionary of index matched to self (R2 = R3)
index_self_match_dict = {}

# Dictionary of counts for each group of matched, index hopped, and unknown indexes
index_type_dict = {'dual_matched' : 0, 'index_hopped' : 0, 'unknown' : 0}

# Dictionary for how many times barcdodes got matched to incorrect barcode 
index_mismatch_dict = {}

# Dictionary for output files
output_file_handles = {}

# Dictionary for reverse complent function
reversed_comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}


### This section will set up the index mismatch dictionary to keep counts of all mismatched pairs, including self matched indexes

# Uses itertools to create a list of all index mismatched permutations
index_mismatch_list = list(itertools.permutations(index_set, 2))

# Populates the mismatch dictionary with each mismatched index pair and an inital count = 0
for i in index_mismatch_list:
    index_mismatch_dict[i] = 0
 
# populates the mismatch dictionary with self matched index pairs and an initial count = 0
for i in index_set:
    index_mismatch_dict[(i,i)] = 0 


### Functions

def reverse_complement(index):
    '''Will take an index as an argument and return the reverse complement of that index'''

    reversed_strand = ""
    length = len(index)
    for i in range(length):
        c = index[length - 1 - i]
        reversed_strand += reversed_comp_dict.get(c, c)
    
    return(reversed_strand)


### Main

# opens the write files for unknown and hopped indexes
U1 = open(u1, mode = 'w') 
U2 = open(u2, mode = 'w')
H1 = open(h1, mode = 'w')
H2 = open(h2, mode = 'w')

# opens all write files for the 48 properly matched indexes
for i in index_set:
    output_file_handles[i] = (open(i + '.R1.fq', mode = 'w'), open(i + '.R2.fq', mode = 'w'))

# opens the read files (R1, R2, R3, R4)
with gzip.open(r1, mode = 'rt') as read_1, gzip.open(r2, mode= 'rt') as read_2, gzip.open(r3, mode= 'rt') as read_3, gzip.open(r4, mode= 'rt') as read_4:
     
    while True:

        # creates tuples to temporarily store each read line for R1-R4 records while looping through all lines in the file
        header = read_1.readline().strip()
        if header == '': 
            # EOF
            break
        seq = read_1.readline().strip()
        plus = read_1.readline().strip()
        q_score = read_1.readline().strip()
        record_R1 = (header, seq, plus, q_score)

        header = read_2.readline().strip()
        seq = read_2.readline().strip()
        plus = read_2.readline().strip()
        q_score = read_2.readline().strip()
        record_index1 = (header, seq, plus, q_score)

        header = read_3.readline().strip()
        seq = read_3.readline().strip()
        plus = read_3.readline().strip()
        q_score = read_3.readline().strip()
        record_index2 = (header, seq, plus, q_score)

        header = read_4.readline().strip()
        seq = read_4.readline().strip()
        plus = read_4.readline().strip()
        q_score = read_4.readline().strip()
        record_read2 = (header, seq, plus, q_score)
        

        # Variables
        header1 = record_R1[0]
        index1 = record_index1[1]
        index2 = record_index2[1]
        header2 = record_read2[0]
        index2_rev = reverse_complement(index2)
        index1_qs = bioinfo.qual_score(record_index1[3])
        index2_qs = bioinfo.qual_score(record_index2[3])


        ### Unknown Indexes

        # checks for N's in index
        if 'N' in index1:
            header1 += '_' + index1 + '-' + index2_rev
            header2 += '_' + index1 + '-' + index2_rev
            
            # add R1 and R2 to unknown
            U1.write(header1 + '\n' + record_R1[1] + '\n' + record_R1[2] + '\n' + record_R1[3] + '\n')
            U2.write(header2 + '\n' + record_read2[1] + '\n' + record_read2[2] + '\n' + record_read2[3] + '\n')
            
            # count number of reads in unknown
            index_type_dict['unknown'] += 1

        #  checks for low quality score
        elif (index1_qs < ic) or (index2_qs < ic):
            header1 += '_' + index1 + '-' + index2_rev
            header2 += '_' + index1 + '-' + index2_rev
            
            # add R1 and R2 to unknown
            U1.write(header1 + '\n' + record_R1[1] + '\n' + record_R1[2] + '\n' + record_R1[3] + '\n')
            U2.write(header2 + '\n' + record_read2[1] + '\n' + record_read2[2] + '\n' + record_read2[3] + '\n')
            
            # count number of reads in unknown
            index_type_dict['unknown'] += 1

        # checks if indexes are not in known indexes set
        elif (index1 not in index_set) or (index2_rev not in index_set):
            header1 += '_' + index1 + '-' + index2_rev
            header2 += '_' + index1 + '-' + index2_rev
            
            # add R1 and R2 to unknown
            U1.write(header1 + '\n' + record_R1[1] + '\n' + record_R1[2] + '\n' + record_R1[3] + '\n')
            U2.write(header2 + '\n' + record_read2[1] + '\n' + record_read2[2] + '\n' + record_read2[3] + '\n')
            
            # count number of reads in unknown
            index_type_dict['unknown'] += 1 

        
        ### Index Hopped

        else:

            if index1 != index2_rev:
                header1 += '_' + index1 + '-' + index2_rev
                header2 += '_' + index1 + '-' + index2_rev
              
                # add R1 and R2 to hopped
                H1.write(header1 + '\n' + record_R1[1] + '\n' + record_R1[2] + '\n' + record_R1[3] + '\n')
                H2.write(header2 + '\n' + record_read2[1] + '\n' + record_read2[2] + '\n' + record_read2[3] + '\n')
                
                # count number of reads in index hopped 
                index_type_dict['index_hopped'] += 1

                # index_mismatch will store every unique index pair with count
                if (index1, index2_rev) in index_mismatch_dict:
                    index_mismatch_dict[(index1, index2_rev)] += 1
                
                 # counts of mismatched pairs where index is matched to self and adds count to index_mismatch
                elif (index1 == index2):
                    index_mismatch_dict[(index1, index2)] += 1
            
            
            ### Dual Matched

            else:

                header1 += '_' + index1 + '-' + index2_rev
                header2 += '_' + index1 + '-' + index2_rev
                
                # add R1 and R2 to dual matched
                output_file_handles[index1][0].write(header1 + '\n' + record_R1[1] + '\n' + record_R1[2] + '\n' + record_R1[3] + '\n')
                output_file_handles[index1][1].write(header2 + '\n' + record_read2[1] + '\n' + record_read2[2] + '\n' + record_read2[3] + '\n')
                    
                # count number of reads in dual matched
                index_type_dict['dual_matched'] += 1

                # index_mismatch will store every unique index pair with count
                if (index1, index2_rev) in index_mismatch_dict:
                    index_mismatch_dict[(index1, index2_rev)] += 1


### Close all write files

U1.close()
U2.close()
H1.close()
H2.close()
for i in index_set:
    output_file_handles[i][0].close()
    output_file_handles[i][1].close()
           

### Stats output

# creates tsv for index pairs and counts
mismatch_list = sorted(index_mismatch_dict(), key=lambda x:x[1], reverse=True)
sorted_mismatch_dict = dict(mismatch_list)

tot = sum(sorted_mismatch_dict.values())

print('Index Pair','\t','Count','\t','Percentage','\n')
for key in sorted_mismatch_dict:
    percent = (sorted_mismatch_dict[key] / tot) * 100
    print(key, sorted_mismatch_dict[key], percent,sep='\t')
print('\n')
print('\n')


# creates tsv for index categories and counts
total = sum(index_type_dict.values())

print('Index Category','\t','Count','\t','Percentage','\n')
for key in index_type_dict:
    percent = (index_type_dict[key] / total) * 100
    print(key, index_type_dict[key], percent,sep='\t')

print('total', total, sep='\t')





                    
                    




       










    
        
    

