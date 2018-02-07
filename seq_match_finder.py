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


def outputRecords(records, file = None) :
    if file == None :
        for r in records :
            r.printValues()        
    else :
        # output matches to CSV file
        outfile = open(file,"w")
        outfile.write("ID,POSITION,STRAND,UPSTREAM\n")
        for r in records :
            outfile.write(r.getCSV() + "\n")
        outfile.close();

el1 = "ccgaggtgagtccggaaatgggctcaaaactgcggtgaaacc".upper()
el2 = "actgacatccggacagcgttgcgacagtggcgcttttagcgcagcccgggggtttttacaggatacc".upper()        
el3 = "gtggcgcttttagcgcagcccgggggtttttacaggatacca".upper()
el321 = "AATTGAGGTGGATCGGTGGATCGGTGGATCAGTTCATTTCGGAACTGAAATGAGCCGTGTCCGAGGTGAGTCCGGAAATGGGCTCAAAACTGCGGTGAAACCACTGACATCCGGACAGCGTTGCGACAGTGGCGCTTTTAGCGCAGCCCGGGGGTTTTTACAGGATACC"

res1 = getHits(el1, "seqs.fa")
outputRecords(res1, "el1_results.csv")

res2 = getHits(el2, "seqs.fa")
outputRecords(res2, "el2_results.csv")

res3 = getHits(el3, "seqs.fa")
outputRecords(res3, "el3_results.csv")

res321 = getHits(el321, "seqs.fa")
outputRecords(res321, "el321_results.csv")




    
