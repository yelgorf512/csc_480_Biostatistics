from Bio import SeqIO # to parse sequence data
from Bio.Seq import Seq

sequences = SeqIO.parse("Vi41b.fa", 'fasta')
target = Seq("AATTGAGGTGGATCGGTGGATCGGTGGATCAGTTCATTTCGGAACTGAAATGAGCCGTGTCCGAGGTGAGTCCGGAAATGGGCTCAAAACTGCGGTGAAACCACTGACATCCGGACAGCGTTGCGACAGTGGCGCTTTTAGCGCAGCCCGGGGGTTTTTACAGGATACC")
rc = target.reverse_complement()

for s in sequences :
    if target in s.seq :
        print(s.id)
    if rc in s.seq :
        print(s.id, "(RC)")
    