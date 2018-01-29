### CSC 480 Assignment 1
### Author: Nick Frogley

import os

# load in file contents
infile = open("sequences.txt")
seqs = infile.read()
infile.close()

lines = seqs.split("\n")

i = 1;

# look for sequence in each line of file and print info about it
for the_line in lines:
    colon_position = the_line.find(":")
    position = the_line.find("GGGAAC")
    if position != -1:
        print ("Sample: " + str(i) + "\tPosition:"+ str(position - colon_position - 1))
        print ("Sequence after: " + the_line[position + 6::])
    i += 1
    
