import pandas as pd

# Load the bin data
bins_df = pd.read_csv('updated_merged_chromosomes_bins.csv')

# Load the gene annotation data
genes_df = pd.read_csv('annotated_genes_plot.txt', sep='\t')

# Create a key for exact match comparison
bins_df['key'] = bins_df['chrom'].astype(str) + '_' + bins_df['start'].astype(str) + '_' + bins_df['end'].astype(str)
genes_df['key'] = genes_df['chrom'].astype(str) + '_' + genes_df['start'].astype(str) + '_' + genes_df['end'].astype(str)

# Group gene info by key (in case multiple genes match the same bin)
grouped_genes = genes_df.groupby('key').agg({
    'gene_name': lambda x: ','.join(x),
    'arabidopsis_gene': lambda x: ','.join(x)
}).reset_index()

# Merge the gene data into the bin data
merged_df = bins_df.merge(grouped_genes, on='key', how='left')

# Replace NaNs with empty strings if needed
merged_df['gene_name'] = merged_df['gene_name'].fillna('')
merged_df['arabidopsis_gene'] = merged_df['arabidopsis_gene'].fillna('')

# Output to new CSV
merged_df.drop(columns=['key']).to_csv('bins_with_gene_annotations.csv', index=False)

# Report genes that didn't match any bin
matched_keys = set(bins_df['key'])
unmatched_genes = genes_df[~genes_df['key'].isin(matched_keys)]

if not unmatched_genes.empty:
    print("Genes that did not match any bin:")
    print(unmatched_genes[['chrom', 'start', 'end', 'gene_name', 'arabidopsis_gene']])
else:
    print("All genes matched successfully.")

