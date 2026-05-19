

import sys
import pandas as pd

kmcp_profile_file = sys.argv[1]
if len(sys.argv) == 3:
    genomes_info_file = sys.argv[2]
else:
    genomes_info_file = None
if genomes_info_file and genomes_info_file != "-":
    genomesID = pd.read_csv(genomes_info_file, sep="\t", usecols=[0])
    genomesID.columns = ["strain_taxid"]
else:
    genomesID = None
kmcp_profile = pd.read_csv(kmcp_profile_file, sep="\t", usecols=[0,1])
kmcp_profile.columns = ["strain_taxid", "abundance"]
kmcp_profile["abundance"] = kmcp_profile["abundance"] / 100
if isinstance(genomesID, pd.DataFrame) and not genomesID.empty:
    kmcp_profile = pd.merge(kmcp_profile, genomesID)
    kmcp_profile["abundance"] = kmcp_profile["abundance"] / kmcp_profile["abundance"].sum()
kmcp_profile_sorted = kmcp_profile.sort_values(by="abundance", ascending=False)
kmcp_profile_sorted.to_csv("strain_abundance.txt", sep="\t", index=False)


