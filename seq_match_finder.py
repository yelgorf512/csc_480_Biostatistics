# CSC 480 Assignment 2
# Nick Frogley
# 2/3/18
#
# Gets information about direct matches for a target sequence in a provided
# FASTA file of sequences and writes it to a CSV file

from Bio import SeqIO # to parse sequence data
from Bio.Seq import Seq
import sys
import os

import owcheck      # custom module

# this class stores info about a match we find and allows for it to be
# printed or returned as a CSV row string
class MatchRecord :

    id_str = ""         # the id from the FASTA record
    match_pos = 0       # the index of the match in the sequence (starts at 0)
    strand = -1         # strand; 0 = original, 1 = reverse complement
    downstream_str = ""   # the sequence 200 bp downstream from the match (if 200 bp downstream exists)

    # set values
    def __init__(self, new_id, new_pos, new_strand, new_downstream = "") :    
        self.id_str = new_id
        self.match_pos = new_pos
        self.strand = new_strand
        self.downstream_str = new_downstream
        self.downstream_len = len(self.downstream_str)
        
    # print values
    def printValues(self) :
        print("ID: " + self.id_str)
        print("POSITION: " + str(self.match_pos + 1))
        if (self.strand == 0) :
            print("FOUND IN ORIGINAL")
        else :
            print("FOUND IN REVERSE COMPLEMENT")
        print("DWNSTR LENGTH: " + self.downstream_len)
        print("DOWNSTREAM (first 100bp): " + self.downstream_str[:100:])
        

    # return values as CSV row
    def getCSV(self) :
        csv_str =   self.id_str + ","
        csv_str +=  str(self.match_pos + 1) + ","
        if (self.strand == 0) :
            csv_str += "ORIGINAL" + ","
        else :
            csv_str += "REVERSE COMPLEMENT" + ","
        csv_str += str(self.downstream_len) + ","
        csv_str += self.downstream_str
        return csv_str

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
            csv_str += "ORIGINAL" + ","
        else :
            csv_str += "REVERSE COMPLEMENT" + ","
        return csv_str

# get matches between target sequences and sequences from a file that are
# in fasta format
def getHits(target_seq, file) :

    # parse sequences in 'fasta' format; this returns an iterator, which stores
    # a sequence of elements
    sequences = SeqIO.parse(file, 'fasta')

    records = list()    # stores our MatchRecord objects (sequence matches)

    target_seq_len = len(target_seq)

    # loop through sequence records in the file
    for s in sequences:
    
        if len(s) < target_seq_len:    
            continue

        # look for a match in the original version of this sequence
        hit_pos = s.seq.find(target_seq)
        if (hit_pos != -1) :    # match found, record it
            downstream = ""
            if hit_pos + target_seq_len + 200 < len(s.seq):
                downstream = s.seq[hit_pos + target_seq_len + 200::]
            records.append(MatchRecord(s.id, hit_pos, 0, str(downstream)))

        # look for a match in the reverse complement of this sequence
        rev_comp = s.seq.reverse_complement()
        hit_pos = rev_comp.find(target_seq)
        if (hit_pos != -1) :    
            downstream = ""
            if hit_pos + target_seq_len + 200 < len(s.seq):
                downstream = s.seq[hit_pos + target_seq_len + 200::]
            records.append(MatchRecord(s.id, hit_pos, 1, str(downstream)))
    print("# matches found:", len(records))
    return records                



def outputRecords(records, file = None, alignments = False) :
    if file == None :
        for r in records :
            r.printValues()        
    else :
        # output matches to CSV file
        outfile = open(file,"w")
        if (alignments == False) :
            outfile.write("ID,POSITION,STRAND,DWNSTR LENGTH,DOWNSTREAM\n")
        else :
            outfile.write("ID,SCORE,STRAND\n")
        for r in records :
            outfile.write(r.getCSV() + "\n")
        outfile.close();

seq_dict = {}

seq_dict['el1'] = "ccgaggtgagtccggaaatgggctcaaaactgcggtgaaacc".upper()
seq_dict['el2'] = "actgacatccggacagcgttgcgacagtggcgcttttagcgcagcccgggggtttttacaggatacc".upper()        
seq_dict['el3'] = "gtggcgcttttagcgcagcccgggggtttttacaggatacca".upper()
seq_dict['el312'] = "AATTGAGGTGGATCGGTGGATCGGTGGATCAGTTCATTTCGGAACTGAAATGAGCCGTGTCCGAGGTGAGTCCGGAAATGGGCTCAAAACTGCGGTGAAACCACTGACATCCGGACAGCGTTGCGACAGTGGCGCTTTTAGCGCAGCCCGGGGGTTTTTACAGGATACC"

if len(sys.argv) == 3 :
    print("Running with provided args")
    input_filename = os.path.split(sys.argv[1])[1]
    output_filename = "matches-" + input_filename + "-" + sys.argv[2] + ".csv"
    owcheck.overwriteFile(output_filename)
    target_seq = str(list(SeqIO.parse(sys.argv[2], 'fasta'))[0].seq)
    hits = getHits(target_seq, sys.argv[1])
    outputRecords(hits, output_filename, False)
    print("Your results have been saved to the following file: " + output_filename)

else:
    print("Missing args (sequences FASTA file, contig FASTA file), running with hardcoded values")
    hits = getHits(seq_dict['el312'], "seqs.fa")
    outputRecords(hits, "direct_matches.csv", False)
