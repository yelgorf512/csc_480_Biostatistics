from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML

import sys

# Assignment #2:

# write a script that reads in a results file from the
# match script, and for each sequence, finds the top BLAST hit.

# Create a new file that copies the original input file and adds the
# following columns: the definition (.hit_def).

# Then for the first high-scoring pair (from topHit.hsps), get the alignment
# length and the percent identity (the number of identities / the alignment
# length).

# This script should take two arguments – the input file name and
# output file name. An example BLAST script is below:

def get_csv_lines(data_file):
    infile = open(data_file)
    contents_raw = infile.read()
    infile.close()
    contents_split_nl = contents_raw.split("\n")
    return contents_split_nl

def output_csv_lines(csv_lines, file = None) :
    if file == None :
        for csv_line in csv_lines :
            print(csv_line)        
    else :
        # output matches to CSV file
        outfile = open(file,"w")
        for csv_line in csv_lines :
            outfile.write(csv_line + "\n")
        outfile.close();

seq_len_cutoff = -1 # if this gets set to something else later sequence will cutoff at this index

if len(sys.argv) == 4 :
    print("Running with provided args")
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    if sys.argv[3] != "ALL":
        seq_len_cutoff = int(sys.argv[3])

elif len(sys.argv) == 3 :
    print("Running with provided args")
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
else:
    print("Missing args (input file, output file, [seq len cutoff]), running with hardcoded values")
    input_file = "direct_matches.csv"
    output_file = "blast_results_2_redo.csv"
    #seq_len_cutoff = 15000

seq_lines = get_csv_lines(input_file)
print("GOT SEQ LINES")

seq_list = list()

skip = 1
counter = 0
limit = -1  

for the_line in seq_lines:
    comma_split_line = the_line.split(",")
    if counter >= skip and len(comma_split_line) == 5 :
        if seq_len_cutoff != -1 and len(the_line) > seq_len_cutoff:
            print('Adding sequence (cutoff)' + str(counter))
            seq_list.append(comma_split_line[4][0:seq_len_cutoff:])
        else: 
            print('Adding sequence ' + str(counter))
            seq_list.append(comma_split_line[4])
    else :
        seq_list.append("")
    counter += 1
    if counter == limit :
        break

new_csv_lines = list()
new_csv_lines.append("ID, POSITION, STRAND, UPST LENGTH, UPSTREAM, TOP DEFINITION, ALIGN LENGTH, PERCENT IDENTITY, Q START")

counter = 0

for seq in seq_list :
    if len(seq) > 0 :
        print("BLASTING SEQUENCE " + str(counter))
        result_handle = NCBIWWW.qblast("blastn", "nt", seq)
        print("BLAST DONE")

        # save results to file
        file = open("blast_result_ " + str(counter) + ".xml", "w")
        file.write(result_handle.read())
        file.close()

        # read in results and parse them,
        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:parsing-blast
        result_handle = open("blast_result_ " + str(counter) + ".xml")
        #result_handle = open("my_blast.xml")
        blast_record = NCBIXML.read(result_handle)

        topHit = blast_record.alignments[0]

        # look at first high scoring pair (hsp)
        hsp = topHit.hsps[0]

        new_csv_line = ""

        comma_split_line = seq_lines[counter].split(",")
        cell_counter = 0
        for cell in comma_split_line:
            if (cell_counter == 4 and seq_len_cutoff != -1):
                new_csv_line += comma_split_line[cell_counter][0:seq_len_cutoff:] + ","
            else:        
                new_csv_line += comma_split_line[cell_counter] + ","
            cell_counter += 1
            
        new_csv_line += topHit.hit_def.replace(",", "|") + ","
        new_csv_line += str(hsp.align_length) + ","
        new_csv_line += str(hsp.identities/hsp.align_length * 100) + ","
        new_csv_line += str(hsp.query_start)
        new_csv_lines.append(new_csv_line)

    elif (counter >= skip) :
        new_csv_line = seq_lines[counter]
        new_csv_lines.append(new_csv_line)

    counter += 1


output_csv_lines(new_csv_lines, output_file)
