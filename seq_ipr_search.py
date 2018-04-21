import os
import sys
import xml.etree.ElementTree as ET

import owcheck 		# custom module

if len(sys.argv) == 3 :
    print("Running with provided args")
    input_filename = os.path.split(sys.argv[1])[1]
    email = sys.argv[2]
    output_filename = "ipr_results_" + input_filename + ".csv"
else:
    print("Missing args (input file, your email), exiting...")
    exit(-1)

owcheck.overwriteFile(output_filename)

str_ipr_call = "python iprscan5_urllib3.py --goterms --pathways --appl=PfamA,CDD --outfile=ipr_output --outformat=xml"
str_ipr_call += " --email=" + email + " " + input_filename

print("CALLING IPR...")
os.system(str_ipr_call)

print("PARSING IPR RESULTS")

outfile = open(output_filename,"w")
outfile.write("AC,DESC,NAME\n")

tree = ET.parse('ipr_output.xml.xml')      # ending up with .xml twice for some reason
the_iter = tree.getiterator()

for elem in the_iter :
	tailing_tag = elem.tag.split("}")[1]	# IPR XML results end up with tags of form {http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5}tag
	if (tailing_tag == "signature") :
		outfile.write(elem.get("ac") + "," + elem.get("desc") + "," + elem.get("name") + "\n")

outfile.close()

print("Your results have been saved to the following file: " + output_filename)