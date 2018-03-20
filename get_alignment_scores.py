# CSC 480 Assignment 2
# Nick Frogley
# 2/3/18
#
# Gets information about direct matches for a target sequence in a provided
# FASTA file of sequences and writes it to a CSV file

from Bio import SeqIO # to parse sequence data
from Bio.Seq import Seq

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

import sys


# this class stores info about an alignment
class AlignmentRecord :

    id_str = ""         # the id from the FASTA record
    score = 0           # alignment score
    strand = -1         # strand; 0 = original, 1 = reverse complement
    
    # set values
    def __init__(self, new_id, new_score, new_strand) :    
        self.id_str = new_id
        self.score = new_score
        self.strand = new_strand

    # print values
    def printValues(self) :
        print("ID: " + self.id_str)
        print("SCORE: " + str(self.score))
        if (self.strand == 0) :
            print("FOUND IN ORIGINAL")
        else :
            print("FOUND IN REVERSE COMPLEMENT")

    # return values as CSV row
    def getCSV(self) :
        csv_str =   self.id_str + ","
        csv_str +=  str(self.score) + ","
        if (self.strand == 0) :
            csv_str += "ORIGINAL"
        else :
            csv_str += "REVERSE COMPLEMENT"
        return csv_str

              

# Get pairwise alignments with a target sequence for each sequence in a FASTA file 
def getAlignments(target_seq, file) :

    rc_target_seq = Seq(target_seq).reverse_complement()
    
    # parse sequences in 'fasta' format; this returns an iterator, which stores
    # a sequence of elements
    sequences = SeqIO.parse(file, 'fasta')

    records = list()    # stores our AlignmentRecord objects (alignment scores)

    target_seq_len = len(target_seq)

    # loop through sequence records in the file
    counter = 1
    for s in sequences:
    
        if len(s) < target_seq_len:
            print("length skip")
            continue

        # get alignment for the original version of this sequence
        if counter % 20 == 0 :
            print ("GETTING ALIGNMENT FOR SEQ " + str(counter) + "...")

        alignments = pairwise2.align.globalms(target_seq, s.seq, 2,-1, -2, -1, 
                                      penalize_end_gaps = (False, True),
                                      score_only = True)

        the_record = AlignmentRecord(s.id, alignments, 0)
        #the_record.printValues()
        records.append(the_record)

        # get alignment for the reverse complement of this sequence
        alignments = pairwise2.align.globalms(rc_target_seq, s.seq, 2,-1, -2, -1, 
                                      penalize_end_gaps = (False, True),
                                      score_only = True)

        the_record = AlignmentRecord(s.id, alignments, 1)
        #the_record.printValues()
        records.append(the_record)
        
        counter += 1
        
    return records                


def outputRecords(records, file = None, alignments = False) :
    if file == None :
        for r in records :
            r.printValues()        
    else :
        # output matches to CSV file
        outfile = open(file,"w")        
        outfile.write("ID,SCORE,STRAND\n")
        for r in records :
            outfile.write(r.getCSV() + "\n")
        outfile.close();

seq_dict = {}

seq_dict['el1'] = "ccgaggtgagtccggaaatgggctcaaaactgcggtgaaacc".upper()
seq_dict['el2'] = "actgacatccggacagcgttgcgacagtggcgcttttagcgcagcccgggggtttttacaggatacc".upper()        
seq_dict['el3'] = "gtggcgcttttagcgcagcccgggggtttttacaggatacca".upper()
seq_dict['el312'] = "AATTGAGGTGGATCGGTGGATCGGTGGATCAGTTCATTTCGGAACTGAAATGAGCCGTGTCCGAGGTGAGTCCGGAAATGGGCTCAAAACTGCGGTGAAACCACTGACATCCGGACAGCGTTGCGACAGTGGCGCTTTTAGCGCAGCCCGGGGGTTTTTACAGGATACC"


if len(sys.argv) == 4 :
    print("Running with provided args")
    alignments = getAlignments(seq_dict[sys.argv[2]], sys.argv[1])        
    outputRecords(alignments, sys.argv[3], True)
else:
    print("Missing args (FASTA file, contig name, output file)")
                   
