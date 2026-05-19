import argparse
import pandas as pd

usage = "Convert seqid2tax and result from sylph to species abundance profile"

def process(sylph_result_file, seqid2tax_file, profile_level="species"):
    sylph_result = pd.read_csv(sylph_result_file, sep="\t")
    sylph_result["seq_id"] = sylph_result["Contig_name"].str.split(" ").str[0]
    if profile_level == "species":
        seqid2tax = pd.read_csv(seqid2tax_file, sep="\t", header=None, dtype=object)
    elif profile_level == "strain":
        seqid2tax = pd.read_csv(seqid2tax_file, sep="\t", header=None, dtype=object, usecols=[0,3])
        cols = seqid2tax.columns[[1,0]]
        seqid2tax = seqid2tax[cols]
    seqid2tax.columns = ["seq_id", "taxonomy"]
    merge_df = pd.merge(sylph_result, seqid2tax, on="seq_id", how="left")
    selected_df = merge_df[["taxonomy", "Taxonomic_abundance"]]
    new_df = selected_df.groupby("taxonomy", as_index=False)["Taxonomic_abundance"].sum()
    new_df["Taxonomic_abundance"] = new_df["Taxonomic_abundance"]/100
    new_df.columns = ["taxonomy", "abundance"]
    new_df.to_csv("sylph_abundance.txt", index=False, sep="\t")
    
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="sylph_convert.py", description=usage)
    parser.add_argument("sylph_result", type=str, help="sylph result file")
    parser.add_argument("seqid2tax", type=str, help="seqid2tax(from centrifuge)")
    parser.add_argument("profile_level", type=str, help="profile_level")
    args = parser.parse_args()
    sylph_result_file = args.sylph_result
    seqid2tax_file = args.seqid2tax
    profile_level = args.profile_level
    process(sylph_result_file, seqid2tax_file, profile_level)