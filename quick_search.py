from Bio import SeqIO # to parse sequence data
from Bio.Seq import Seq
import sys

if len(sys.argv) != 2 :
    raise Exception("Usage: python quick_search contigs.fa")

file = sys.argv[1]

print("Argument is:", sys.argv[1], file = sys.stderr)

sequences = SeqIO.parse(file, 'fasta')

target  = Seq("AATTGAGGTGGATCGGTGGATCGGTGGATCAGTTCATTTCGGAACTGAAATGAGCCGTGTCCGAGGTGAGTCCGGAAATGGGCTCAAAACTGCGGTGAAACCACTGACATCCGGACAGCGTTGCGACAGTGGCGCTTTTAGCGCAGCCCGGGGGTTTTTACAGGATACC")

rc = target.reverse_complement()

for s in sequences :

    if ("18999-AE09B-Vi40b4011-GCCAAT_S3_L001_R1_001_(paired)_contig_8533" in s.id) :
        print(">", s.id, sep = "")
        print(s.seq)
    #if target in s.seq :
    #    print(s.id)
       # print(s.seq)
       # exit(-1)
    #if rc in s.seq :
        #print(s.id, "(RC)")
    
