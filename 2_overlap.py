import pandas as pd

# Load data
bins_df = pd.read_csv('updated_merged_chromosomes_bins.csv', low_memory=False)
genes_df = pd.read_csv('annotated_genes_plot.txt', sep='\t')

# Strip whitespace and remove rows with missing or invalid start/end
for col in ['start', 'end']:
    bins_df[col] = bins_df[col].astype(str).str.strip()
    bins_df = bins_df[bins_df[col].str.isnumeric()]  # Keep only numeric rows
    bins_df[col] = bins_df[col].astype(int)

    genes_df[col] = genes_df[col].astype(str).str.strip()
    genes_df = genes_df[genes_df[col].str.isnumeric()]
    genes_df[col] = genes_df[col].astype(int)

# Now safe to continue
bins_df['key'] = bins_df['chrom'].astype(str) + '_' + bins_df['start'].astype(str) + '_' + bins_df['end'].astype(str)
genes_df['key'] = genes_df['chrom'].astype(str) + '_' + genes_df['start'].astype(str) + '_' + genes_df['end'].astype(str)

# (Rest of the script remains unchanged)

# Exact match first
grouped_genes = genes_df.groupby('key').agg({
    'gene_name': lambda x: ','.join(x),
    'arabidopsis_gene': lambda x: ','.join(x)
}).reset_index()

merged_df = bins_df.merge(grouped_genes, on='key', how='left')
merged_df['gene_name'] = merged_df['gene_name'].fillna('')
merged_df['arabidopsis_gene'] = merged_df['arabidopsis_gene'].fillna('')

# Now handle unmatched genes
matched_keys = set(merged_df['key'])
unmatched_genes = genes_df[~genes_df['key'].isin(matched_keys)]

# Try to match within ±1000 bp range
def find_nearby_bin(gene_row, bins):
    chrom_match = bins[bins['chrom'] == gene_row['chrom']]
    condition = (
        (chrom_match['start'] >= gene_row['start'] - 10000) &
        (chrom_match['end'] <= gene_row['end'] + 10000)
    )
    return chrom_match[condition]

# Track unmatched after relaxation
still_unmatched = []

# Add found nearby matches
for _, gene in unmatched_genes.iterrows():
    nearby_bins = find_nearby_bin(gene, bins_df)
    if not nearby_bins.empty:
        for idx in nearby_bins.index:
            if merged_df.at[idx, 'gene_name']:
                merged_df.at[idx, 'gene_name'] += ',' + gene['gene_name']
                merged_df.at[idx, 'arabidopsis_gene'] += ',' + gene['arabidopsis_gene']
            else:
                merged_df.at[idx, 'gene_name'] = gene['gene_name']
                merged_df.at[idx, 'arabidopsis_gene'] = gene['arabidopsis_gene']
    else:
        still_unmatched.append(gene)

# Output final file
merged_df.drop(columns=['key']).to_csv('bins_with_gene_annotations_rescued.csv', index=False)

# Report final unmatched genes
if still_unmatched:
    print("Genes that could not be matched even within ±1000 bp:")
    print(pd.DataFrame(still_unmatched)[['chrom', 'start', 'end', 'gene_name', 'arabidopsis_gene']])
else:
    print("All genes were successfully matched either exactly or within ±1000 bp.")

