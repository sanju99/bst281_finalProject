import numpy as np
import pandas as pd

# get the database of interactions
virus_human = pd.read_csv("Processed/virus_human.csv")

# all human proteins are in the B columns now
virus_human_proteins = virus_human["Protein B"].unique()

# save to a text file
pd.DataFrame(virus_human_proteins).to_csv("Data/human_proteins.txt", index=False, header=None, sep="\t")

# converted UniProt to UniProt KB ID by searching online to update names. This is a full dataframe, not a text file
human_fixed_df = pd.read_csv("Data/human_proteins_fixed.tab", sep="\t")

print(f"{len(human_fixed_df)} human proteins interact with viral proteins")

# indicates that the second column is what we want
assert (sum(human_fixed_df.iloc[:, 0].values != human_fixed_df.iloc[:, 1].values)) != 0

## Get gene names from the UniProt KB IDS

# downloaded from UniProt. This is the human proteome
human_proteins_all = pd.read_csv("Data/human_proteins_all.gz", sep="\t", compression="gzip")

print(f"{len(human_proteins_all)} proteins in the human interactome")

### Separate the gene names from the columns in the dataframe 

def get_genes(df, col):
    
    gene_names = []

    for _, row in df.iterrows():

        if not pd.isnull(row["Gene names"]):
            gene_names.append(row["Gene names"].split(" ")[0])
        else:
            gene_names.append(row[col].split("_")[0])
            
    return gene_names

# full database
gene_names_all = get_genes(human_proteins_all, "Entry name")

# makes sense, ~20,000 genes in the human genome
print(f"{len(np.unique(gene_names_all))} genes in the human genome")

# genes from the virus database
gene_names = get_genes(human_fixed_df, "Entry name")

# there are a few in gene_names that are not in gene_names_all. It's because of orf names or deprecated/mislabeled names
gene_names = list(set(gene_names).intersection(gene_names_all))
print(f"{len(gene_names)} genes involved in viral infection")

# save both to files to search using GOrilla
pd.DataFrame(np.unique(gene_names_all)).to_csv("Data/human_genes_all.txt", sep="\t", header=None, index=False)
pd.DataFrame(gene_names).to_csv("Data/human_genes.txt", sep="\t", header=None, index=False)
