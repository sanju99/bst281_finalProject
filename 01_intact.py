import numpy as np
import pandas as pd

virus_db = pd.read_csv("Data/virus_intact.txt", 
                             usecols=['#ID(s) interactor A', 'ID(s) interactor B', 'Interaction detection method(s)', 'Publication Identifier(s)', 
                                      'Taxid interactor A', 'Taxid interactor B', 'Interaction type(s)', 'Source database(s)', 
                                      'Interaction identifier(s)','Confidence value(s)'],
                             sep="\t")

print(f"Original database shape: {virus_db.shape}")

# make a copy to work with so that we don't chane the original
virus_intact = virus_db.copy()

# rename columns to make it easier to work with
virus_intact = virus_intact.rename(columns=dict(zip(list(virus_intact.columns),
    ["Protein A", "Protein B", "Method", "Publication", "Taxid A", "Taxid B", "Interaction type", "Source DB", "Interaction identifier", "Confidence"])))

# missing values have a dash
virus_intact = virus_intact.replace("-", np.nan)

# remove any rows where a protein is missing. I don't know what that means honestly. Because there are rows where Protein A = Protein B, so that's
# a self interaction
virus_intact = virus_intact.dropna(subset=["Protein A", "Protein B"]).reset_index(drop=True)
print(f"Shape after dropping protein NaNs: {virus_intact.shape}")

## Separate protein IDs so that the source and ID aren't together in a string.

def fix_protein_ids(col1, col2):
    
    split_df = virus_intact[col1].str.split(":", expand=True)
    
    types = []
    ids = []

    for i, row in split_df.iterrows():

        types.append(row[0])

        if pd.isnull(row[2]):
            ids.append(row[1])
        else:
            ids.append(row[2])

    virus_intact[col2] = types
    virus_intact[col1] = ids
    
    
# run the function on the two columns
fix_protein_ids("Protein A", "Protein A DB")
fix_protein_ids("Protein B", "Protein B DB")


## Do the same thing with the taxids to separate the IDs and organism names

temp_A = np.unique([val.split(":")[1].split("|")[0] for val in virus_intact["Taxid A"].values])
temp_2_A = pd.DataFrame([val.rstrip(")").split("(") for val in temp_A])

non_taxa_A = temp_2_A.loc[temp_2_A[0].astype(int) < 0].iloc[:, :2]
non_taxa_mapping = dict(zip(non_taxa_A[0], non_taxa_A[1]))

temp_B = np.unique([val.split(":")[1].split("|")[0] for val in virus_intact["Taxid B"].values])
temp_2_B = pd.DataFrame([val.rstrip(")").split("(") for val in temp_B])

non_taxa_B = temp_2_B.loc[temp_2_B[0].astype(int) < 0].iloc[:, :2]
non_taxa_mapping = {**non_taxa_mapping, **dict(zip(non_taxa_B[0], non_taxa_B[1]))}
non_taxa_df = pd.DataFrame(non_taxa_mapping, index=[0]).T.reset_index()
non_taxa_df.columns = ["taxid", "taxname"]

taxids = list(set(temp_2_B.loc[temp_2_B[0].astype(int) > 0].iloc[:, 0].values.astype(int)).union(temp_2_A.loc[temp_2_A[0].astype(int) > 0].iloc[:, 0].values.astype(int)))

# save the list of taxids to a dataframe
pd.DataFrame(pd.Series(taxids)).to_csv("Data/taxids.txt", sep="\t", index=False, header=None)

# searched in the NCBI taxonomy browser and got the following file
taxids_with_names = pd.read_csv("Data/taxids_with_names.txt", usecols=["taxid", "taxname"], sep="\t")
taxids_with_names = pd.concat([non_taxa_df, taxids_with_names]).reset_index(drop=True)

### Separate the taxids for the two columns. The code above gets only the unique names for searching. But here, do it for every row and save it.

taxids_A = np.array([val.split(":")[1].split("|")[0].split("(")[0] for val in virus_intact["Taxid A"].values]).astype(int)
taxids_B = np.array([val.split(":")[1].split("|")[0].split("(")[0] for val in virus_intact["Taxid B"].values]).astype(int)

# replace the columns with only the taxids
virus_intact['Taxid A'] = taxids_A
virus_intact['Taxid B'] = taxids_B

# create a dictionary to map between taxid and taxname
dict_for_taxa = dict(zip(taxids_with_names["taxid"], taxids_with_names["taxname"]))

# make the columns for the taxonomies
virus_intact['Taxname A'] = virus_intact['Taxid A'].map(dict_for_taxa)
virus_intact['Taxname B'] = virus_intact['Taxid B'].map(dict_for_taxa)

# get interactions with human proteins
virus_human = virus_intact.loc[(virus_intact["Taxname A"] == 'Homo sapiens') | (virus_intact["Taxname B"] == 'Homo sapiens')].reset_index(drop=True)
print(f"{len(virus_human) / len(virus_intact)} proportion involves human proteins")

# not sure why, but there are some other species in this dataset: human proteins interacting with proteins from other organisms
interact_taxa = pd.Series(list(set(virus_human["Taxname A"]).union(virus_human["Taxname B"]))).dropna()
interact_taxa = list(interact_taxa)

keep_taxa = [name for name in interact_taxa if "virus" in name or "HIV" in name or "SARS" in name or "Homo sapiens" in name]
drop_taxa = [name for name in interact_taxa if name not in keep_taxa]

virus_human = virus_human.loc[(virus_human["Taxname A"].isin(keep_taxa)) & (virus_human["Taxname B"].isin(keep_taxa))].reset_index(drop=True)

## Some human-human interactions were found to be related to viral infection, which is why they're in this database. Keep them separately for now.

human_human = virus_human.loc[(virus_human["Taxname A"] == 'Homo sapiens') & (virus_human["Taxname B"] == 'Homo sapiens')]

# separate virus-host and host-host that's important for viruses
contains_viruses_indices = list(set(virus_human.index) - set(human_human.index))
virus_human = virus_human.iloc[contains_viruses_indices, :]

# save both files to CSVs
virus_human.to_csv("Processed/virus_human.csv", index=False)
human_human.to_csv("Processed/human_important_for_virus.csv", index=False)
