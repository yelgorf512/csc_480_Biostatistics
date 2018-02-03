# CSC 480 Assignment 2
# Nick Frogley
# 2/3/18
#
# Gets information about direct matches for a target sequence in a provided
# FASTA file of sequences and writes it to a CSV file

from Bio import SeqIO # to parse sequence data

import io # to convert string to 'handle'

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

    # print values
    def printValues(self) :
        print("ID: " + self.id_str)
        print("POSITION: " + str(self.match_pos + 1))
        if (self.strand == 0) :
            print("FOUND IN ORIGINAL")
        else :
            print("FOUND IN REVERSE COMPLEMENT")
        print("UPSTREAM (first 100bp): " + self.upstream_str[:100:])

    # return values as CSV row
    def getCSV(self) :
        csv_str =   self.id_str + ","
        csv_str +=  str(self.match_pos + 1) + ","
        if (self.strand == 0) :
            csv_str += "ORIGINAL" + ","
        else :
            csv_str += "REVERSE COMPLEMENT" + ","
        csv_str += self.upstream_str
        return csv_str

# open and initialize sequence file
handle = open("seqs.fa", "r")
handle = handle.read()
handle = io.StringIO(handle)

# parse sequences in 'fasta' format; this returns an iterator, which stores
# a sequence of elements
sequences = SeqIO.parse(handle, "fasta")

# each sequence is stored as a SeqRecord object
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:seq_features

# EL 3, 1, 2 SEQ: AATTGAGGTGGATCGGTGGATCGGTGGATCAGTTCATTTCGGAACTGAAATGAGCCGTGTCCGAGGTGAGTCCGGAAATGGGCTCAAAACTGCGGTGAAACCACTGACATCCGGACAGCGTTGCGACAGTGGCGCTTTTAGCGCAGCCCGGGGGTTTTTACAGGATACC

num = 1
target_seq = "AATTGAGGTGGATCGGTGGATCGGTGGATCAGTTCATTTCGGAACTGAAATGAGCCGTGTCCGAGGTGAGTCCGGAAATGGGCTCAAAACTGCGGTGAAACCACTGACATCCGGACAGCGTTGCGACAGTGGCGCTTTTAGCGCAGCCCGGGGGTTTTTACAGGATACC"
target_seq_len = len(target_seq)

records = list()    # stores our MatchRecord objects (sequence matches)

# loop through sequence records in the file
for s in sequences:
    
    if len(s) < target_seq_len :
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
        
    #num = num + 1
    #if num == 4 :
        #break

for r in records :
    r.printValues()
    
print("TOTAL MATCHES: " + str(len(records)))

print("WRITING TO CSV FILE...")

# output matches to CSV file
outfile = open("direct_matches.csv","w")
outfile.write("ID,POSITION,STRAND,UPSTREAM\n")
for r in records :
    outfile.write(r.getCSV() + "\n")
outfile.close();

print("DONE")
