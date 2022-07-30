#!/usr/bin/env python

import bioinfo
import argparse 
import gzip


def get_args():
    parser=argparse.ArgumentParser(description='Caclulates the Q score Distribution-per nucleotide')
    parser.add_argument('-f', help='file_name')
    parser.add_argument('-l', help='read_length')
    parser.add_argument('-o', help='output_file')
    return parser.parse_args()

args=get_args()

f = args.f
l = int(args.l)
o = args.o

# print(f)


def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    
    for element in range(l):
        lst.append(value)
    return lst


my_list: list = []
my_list = init_list(my_list)


    
with gzip.open(f,'rt') as fq:
    line_count=0
    for line in fq:
        line_count+=1
        line=line.strip('\n')
        
        if line_count%4 == 0:   # gives the phred score line
            phred_score = line    
        
        
            for count, p_count in enumerate(phred_score):   # calculates the phred score and populates q_list
                my_list[count] += bioinfo.convert_phred(p_count)

# print(my_list)  

for position, count in enumerate(my_list):
    mean_score = count/(line_count/4)
    my_list[position] = mean_score

# print(my_list)

import matplotlib.pyplot as plt

x = range(len(my_list))
y = (my_list)

plt.bar(x,y)
plt.xlabel("Base Pair Position")
plt.ylabel("Mean Quality Score")
plt.title("Quality Score Mean Distribution")
plt.savefig(o)



  
    
    
    
    
    
    
    
    
   