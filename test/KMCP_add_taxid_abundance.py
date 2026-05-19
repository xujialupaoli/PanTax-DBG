import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description='Process strain abundance and genome info files.')
    parser.add_argument('--abundance', required=True, help='Path to strain_abundance.txt file')
    parser.add_argument('--genome', required=True, help='Path to genome info file')
    parser.add_argument('--output', required=True, help='Path to output file')
    args = parser.parse_args()

    
    abundance_df = pd.read_csv(args.abundance, sep='\t')

    #
    genome_df = pd.read_csv(
        args.genome, sep='\t', header=0,
        dtype={'genome_ID': str, 'strain_taxid': str, 'species_taxid': str,
               'organism_name': str, 'id': str}
    )

    # 
    merged_df = pd.merge(abundance_df, genome_df,
                         left_on='strain_taxid', right_on='genome_ID', how='left')

    # 
    species_abundance = (
        merged_df.groupby('species_taxid')['abundance']
        .sum().reset_index()
        .rename(columns={'abundance': 'species_abundance'})
    )

    #
    final_df = pd.merge(merged_df, species_abundance, on='species_taxid', how='left')

    # 
    final_df = final_df[['strain_taxid_x', 'abundance', 'strain_taxid_y',
                         'species_taxid', 'species_abundance']]
    final_df.columns = ['genome_ID', 'strain_abundance',
                        'strain_taxid', 'species_taxid', 'species_abundance']

    # 
    final_df['species_taxid'] = pd.to_numeric(final_df['species_taxid'], errors='coerce')
    final_df = final_df.sort_values(by='species_taxid', na_position='last')

    # 
    final_df.to_csv(args.output, sep='\t', index=False)
    print(f"✅ Results saved to {args.output}")

if __name__ == '__main__':
    main()
