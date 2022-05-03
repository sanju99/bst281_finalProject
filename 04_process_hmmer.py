import glob
import os
import sys


_, input_dir, fasta_file, output_file = sys.argv


def combine_hmmer_results(input_dir, fasta_file, output_file):
    
    fNames = glob.glob(os.path.join(input_dir, "*.tsv"))
    print(f"{len(fNames)} proteins in the dataset")
    
    nums = [int(os.path.basename(fName).split("_")[1].split(".")[0]) for fName in fNames]
    nums_files_dict = dict(zip(nums, fNames))
    nums_files_dict = dict(sorted(nums_files_dict.items()))
    
    protein_names = []
    with open(fasta_file, "r") as f:

        for line in f:
            if ">" in line:
                protein_names.append(line.lstrip(">").rstrip("\n"))

    # protein_names is 0-indexed, nums_files_dict is 1-indexed
    present_proteins = [protein_names[i-1] for i in list(nums_files_dict.keys())]

    assert len(present_proteins) == len(fNames)

    dfs_lst = []

    for i in range(len(present_proteins)):

        file = pd.read_csv(list(nums_files_dict.values())[i], usecols=["Family id", "Family Accession"], sep="\t")
        file["UniProt"] = [present_proteins[i]]*len(file)
        dfs_lst.append(file)

    df = pd.concat(dfs_lst).reset_index(drop=True)

    res = pd.DataFrame(df.groupby(["UniProt", "Family id"]).count()).reset_index()
    res.columns = ["UniProt", "PFAM", "Count"]
    
    # in case another file extension was provided
    output_name = os.path.splitext(output_file)[0] + ".csv"
    
    res.to_csv(output_name, index=False)