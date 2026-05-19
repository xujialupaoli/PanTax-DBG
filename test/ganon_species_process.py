

import sys

report_file = sys.argv[1]

tax_profile_dict = {}
with open(report_file, "r") as f_in:
    for line in f_in:
        if line.strip().startswith("species"):
            tokens = line.strip().split("\t")
            species_taxid = tokens[2].split("|")[-1]
            abundance = float(tokens[8]) / 100
            tax_profile_dict[species_taxid] = abundance

sorted_tax_profile_dict = dict(sorted(tax_profile_dict.items(), key=lambda item: item[1], reverse=True))
with open("species_abundance.txt", "w") as f_out:          
    f_out.write("species_taxid\tpredicted_abundance\n")    
    for k,v in sorted_tax_profile_dict.items():
        f_out.write(f"{k}\t{v}\n")
        