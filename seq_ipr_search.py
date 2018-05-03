# example usage - python seq_ipr_search.py direct_matches.csv nickfro@gmail.com 15000

import os
import sys
import xml.etree.ElementTree as ET

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import owcheck 		# custom module

seq_len_cutoff = -1

def get_csv_lines(data_file):
    infile = open(data_file)
    contents_raw = infile.read()
    infile.close()
    contents_split_nl = contents_raw.split("\n")
    return contents_split_nl

if len(sys.argv) >= 3 :
    print("Running with provided args")
    input_file_path = sys.argv[1]
    input_filename = os.path.split(sys.argv[1])[1]
    email = sys.argv[2]
    output_filename = "ipr_results_" + input_filename + ".csv"
    if (len(sys.argv) == 4) :
        seq_len_cutoff = int(sys.argv[3])

else:
    print("Missing args (input file, your email, [cutoff point; optional]), exiting...")
    exit(-1)

owcheck.overwriteFile(output_filename)

seq_lines = get_csv_lines(input_file_path)

print("GOT CONTIGS FROM INPUT FILE")

seq_list = list()
contig_list = list()
protein_list = list()

skip = 1
counter = 0
limit = -1

for the_line in seq_lines:
    comma_split_line = the_line.split(",")
    if counter >= skip :
        if len(comma_split_line) == 5 :
            contig_list.append(comma_split_line[0])     # add the contig name at the same index as the sequence
            if seq_len_cutoff != -1 and len(the_line) > seq_len_cutoff:
                print('Adding sequence (cutoff) ' + str(counter))
                seq_list.append(comma_split_line[4][0:seq_len_cutoff:])
            else: 
                print('Adding sequence ' + str(counter))
                seq_list.append(comma_split_line[4])
    
    counter += 1
    if counter == limit :
        break

for dna_seq in seq_list :
    the_seq = Seq(dna_seq)
    the_rc_seq = the_seq.reverse_complement()
    
    frame_list = list()
    
    for frame_offset in range(0, 3) :       # get first 3 reading frames (original seq with offset 0, 1, 2)
        the_offset_seq = str(the_seq)[frame_offset::]   # offset seq for a different reading frame
        len_mod_three = len(the_offset_seq) % 3    # see if the seq length is divisble by 3 - if it's not, we need to fill in N's to change that
        
        if (len_mod_three > 0) :
            for i in range(0, (3 - len_mod_three)) :
                the_offset_seq += "N"

        the_protein_seq = str(Seq(the_offset_seq).translate())
        the_protein_seq = the_protein_seq.replace("*", "") # Interpro won't take stop codons, which are *s in the protein sequence, so we get rid of them here
        frame_list.append(the_protein_seq)

    for frame_offset in range(0, 3) :       # get last 3 reading frames (reverse complement seq with offset 0, 1, 2)
        the_offset_rc_seq = str(the_rc_seq)[frame_offset::]   # offset seq for a different reading frame
        len_mod_three = len(the_offset_rc_seq) % 3    # see if the seq length is divisble by 3 - if it's not, we need to fill in N's to change that
        
        if (len_mod_three > 0) :
            for i in range(0, (3 - len_mod_three)) :
                the_offset_rc_seq += "N"

        the_protein_seq = str(Seq(the_offset_rc_seq).translate())
        the_protein_seq = the_protein_seq.replace("*", "") # Interpro won't take stop codons, which are *s in the protein sequence, so we get rid of them here
        frame_list.append(the_protein_seq)      

    protein_list.append(frame_list)


counter = 0
str_csv = ""    # this is the string that the CSV output file will be built in

for protein_set in protein_list:
    if (len(protein_set[0]) >= 11) :     # 11 is minimum length protein allowed by Interpro Scan
        print("HANDLING CONTIG " + str(counter + 1) + "/" + str(len(protein_list)))
        str_csv += (contig_list[counter] + ",FRAME,AC,DESC,NAME\n")
        frame_counter = 1
        for protein_seq in protein_set :
            print("CALLING IPR FOR FRAME " + str(frame_counter) + "/6...")
            
            # put the protein sequence in a temporary FASTA file for IPR to read in
            str_protein_fasta = "> Protein sequence for IPR search\n" + protein_seq
            outfile = open("ipr_input_protein.fa","w")
            outfile.write(str_protein_fasta)
            outfile.close()

            str_ipr_call = "python iprscan5_urllib3.py --goterms --pathways --appl=PfamA,CDD --outfile=ipr_output --outformat=xml"
            str_ipr_call += " --email=" + email + " ipr_input_protein.fa"
            
            os.system(str_ipr_call)

            print("PARSING IPR RESULTS")

            tree = ET.parse('ipr_output.xml.xml')      # ending up with .xml twice for some reason
            the_iter = tree.getiterator()

            for elem in the_iter :
                tailing_tag = elem.tag.split("}")[1]    # IPR XML results end up with tags of form {http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5}tag
                if (tailing_tag == "signature") :
                    str_csv += ("," + str(frame_counter) + "," + elem.get("ac") + "," + elem.get("desc") + "," + elem.get("name") + "\n")

            frame_counter += 1
    elif (len(protein_set[0]) > 0) :
        print("PROTEIN TOO SHORT FOR ANALYSIS FOR CONTIG" + str(counter + 1) + ", skipping...")
        str_csv += (contig_list[counter] + ",PROTEIN TOO SHORT FOR ANALYSIS\n")
    else :
        print("NO DOWNSTREAM SEQUENCE FOR CONTIG " + str(counter + 1) + ", skipping...")
        str_csv += (contig_list[counter] + ",NO DOWNSTREAM SEQUENCE\n")
    counter += 1    

outfile = open(output_filename,"w")
outfile.write(str_csv)
outfile.close()
print("Your results have been saved to the following file: " + output_filename)
