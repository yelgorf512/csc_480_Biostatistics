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

import io # to convert string to 'handle'
import sys

# this class stores info about a match we find and allows for it to be
# printed or returned as a CSV row string
class MatchRecord :

    id_str = ""         # the id from the FASTA record
    match_pos = 0       # the index of the match in the sequence (starts at 0)
    strand = -1         # strand; 0 = original, 1 = reverse complement
    upstream_str = ""   # the sequence 200 bp upstream from the match (if 200 bp upstream exists)

    # set values
    def __init__(self, new_id, new_pos, new_strand, new_upstream = "") :    
        self.id_str = new_id
        self.match_pos = new_pos
        self.strand = new_strand
        self.upstream_str = new_upstream
        self.upstream_len = len(self.upstream_str)
        
    # print values
    def printValues(self) :
        print("ID: " + self.id_str)
        print("POSITION: " + str(self.match_pos + 1))
        if (self.strand == 0) :
            print("FOUND IN ORIGINAL")
        else :
            print("FOUND IN REVERSE COMPLEMENT")
        print("UPST LENGTH: " + self.upstream_len)
        print("UPSTREAM (first 100bp): " + self.upstream_str[:100:])
        

    # return values as CSV row
    def getCSV(self) :
        csv_str =   self.id_str + ","
        csv_str +=  str(self.match_pos + 1) + ","
        if (self.strand == 0) :
            csv_str += "ORIGINAL" + ","
        else :
            csv_str += "REVERSE COMPLEMENT" + ","
        csv_str += str(self.upstream_len) + ","
        csv_str += self.upstream_str
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
            upstream = ""
            if hit_pos + target_seq_len + 200 < len(s.seq):
                upstream = s.seq[hit_pos + target_seq_len + 200::]
            records.append(MatchRecord(s.id, hit_pos, 0, str(upstream)))

        # look for a match in the reverse complement of this sequence
        rev_comp = s.seq.reverse_complement()
        hit_pos = rev_comp.find(target_seq)
        if (hit_pos != -1) :    
            upstream = ""
            if hit_pos + target_seq_len + 200 < len(s.seq):
                upstream = s.seq[hit_pos + target_seq_len + 200::]
            records.append(MatchRecord(s.id, hit_pos, 1, str(upstream)))
    print("# matches found:", len(records))
    return records                

# Get pairwise alignments with a target sequence for each sequence in a FASTA file 
def getAlignments(target_seq, rc_target_seq, file) :

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
        if (alignments == False) :
            outfile.write("ID,POSITION,STRAND,UPST LENGTH,UPSTREAM\n")
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

seq_dict['rc_el1'] = Seq(seq_dict['el1']).reverse_complement()
seq_dict['rc_el2'] = Seq(seq_dict['el2']).reverse_complement()
seq_dict['rc_el3'] = Seq(seq_dict['el3']).reverse_complement()
seq_dict['rc_el312'] = Seq(seq_dict['el312']).reverse_complement()

#res1 = getHits(seq_dict['el1'], "seqs.fa")
#outputRecords(res1, "el1_results.csv")

#res2 = getHits(seq_dict['el2'], "seqs.fa")
#outputRecords(res2, "el2_results.csv")

#res3 = getHits(seq_dict['el3'], "seqs.fa")
#outputRecords(res3, "el3_results.csv")

#res321 = getHits(seq_dict['el312'], "seqs.fa")
#outputRecords(res312, "el312_results.csv")

if len(sys.argv) == 4 :
    print("Running with provided args")
    #alignments = getAlignments(seq_dict[sys.argv[2]], seq_dict['rc_' + sys.argv[2]], sys.argv[1])
    #outputRecords(alignments, sys.argv[3], True)
    hits = getHits(seq_dict[sys.argv[2]], sys.argv[1])
    outputRecords(hits, sys.argv[3], True)
else:
    print("Missing args (FASTA file, contig name, output file), running with hardcoded values")
    #alignments = getAlignments(seq_dict['el312'], seq_dict['rc_el312'], "seqs.fa")
    #outputRecords(alignments, "el312_alignment_results.csv", True)
    hits = getHits(seq_dict['el312'], "seqs_4.fa")
    outputRecords(hits, "direct_matches.csv", False)
                   
#alignments = pairwise2.align.globalms(the_target, contig_1810, 2,-1, -2, -1, 
#                                      penalize_end_gaps = (False, True))

#print("Optimal alignment has score of: " + str(alignments[0][2]) + "\n")
#print("Optimal alignments are below: ")
#for a in alignments:
#    print(format_alignment(*a))


    
