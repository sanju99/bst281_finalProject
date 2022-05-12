# Final Project for BST 281, Spring 2022

## Directories in this repository:

<b>Data</b>: Raw data downloaded from various databases

<b>Processed</b>: Tidy data files for analysis

<b>Analysis</b>: iPython notebooks used for analyses

## Data files by folder:

### Data

<ul>
  <li><code>virus_intact.txt</code>: database of 42,154 pairwise interactions between virus and host cell proteins. This computaionally-maintained database was downloaded from <a href="https://www.ebi.ac.uk/intact/download/datasets#computationally" target="_blank">IntAct</a> via this <a href="https://www.ebi.ac.uk/intact/query/annot:%22dataset:virus%22?conversationContext=7" target="_blank">search query</a> generated by the website to point to the virus database and then selecting <b>MI-TAB 2.5</b> format for downloading.</li>
  <li><code>virus_genome_sizes.csv</code>: File of viral genome sizes downloaded from the <a href="https://www.ncbi.nlm.nih.gov/genome/browse/#!/viruses/" target="_blank">NCBI</a>. </li>
  <li><code>taxids.txt</code>: List of taxonomic IDs from <code>virus_intact.txt</code> to search to get names for. <code>virus_intact.txt</code> only contains taxonomic IDs, and it is helpful to filter by names.</li>
  <li><code>taxids_with_names.txt</code>: Table of taxonomic IDs with associated names, obtained by searching the <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi" target="_blank">NCBI taxonomy browser</a>.</li>
  <li><code>human_proteins.txt</code>: List of all human proteins from the virus IntAct database</li>
  <li><code>human_proteins_fixed.tab</code>: Table of updated human proteins. It was obtained by searching <code>human_proteins.txt</code> on the UniProt website to update and clean deprecated names.</li>
  <li><code>human_genes.txt</code>: Table of human genes from the host:virus database. Gene names were extracted from <code>human_proteins_fixed.tab</code>.</li>
  <li><code>human_proteins_all.gz</code>: Table of all ~79,000 proteins in the human proteome downloaded from <a href="https://www.uniprot.org/uniprot/?query=proteome:UP000005640" target="_blank">UniProt</a>.</li>
  <li><code>human_genes_all.txt</code>: List of all human genes to use as the background list for GOrilla enrichment searching. Gene names were extracted from <code>human_proteins_all.gz</code>.</li>
  <li><code>GO_function.tsv</code>: GOrilla function results returned from searching <code>human_genes.txt</code> and <code>human_genes_all.txt</code> using the "Two unranked lists of genes" feature</li>
  <li><code>GO_process.tsv</code>: GOrilla process results returned from searching <code>human_genes.txt</code> and <code>human_genes_all.txt</code> using the "Two unranked lists of genes" feature</li>
</ul>

### Processed:

<ul>
  <li><code>virus_human.csv</code>: Protein-protein interactions between viruses and humans. This file includes columns for taxonomies that were added in later steps. All viral proteins are in columns labeled <code>A</code>, and all human proteins are in columns labeled <code>B</code>. </li>
  <li><code>human_important_for_virus.csv</code>: Pairwise interactions between human proteins that are involved in viral infection. This file is not used in this project, but these interactions were present in <code>Data/virus_intact.txt</code> downloaded from IntAct.</li>
  <li><code>GO_function</code: Pickle file that is very similar to <code>Data/GO_function.tsv</code>. The only difference is that the Genes column has been cleaned up to remove gene descriptions and contains lists of genes instead of strings.</li>
  <li><code>GO_process</code: Pickle file that is very similar to <code>Data/GO_process.tsv</code>. The only difference is that the Genes column has been cleaned up to remove gene descriptions and contains lists of genes instead of strings.</li>
  <li><code>enriched_GO_proteins.csv</code>: Tidy dataframe of proteins by Gene name and UniProt ID with their associated gene ontology annotations. This file includes ONLY genes associated with GO functions with enrichment scores greater than the mean, obtained from <code>Processed/GO_function</code>.</li>
  <li><code>protein_domains.csv</code>: Tidy dataframe of proteins in <code>enriched_GO_proteins.csv</code> by their UniProt IDs and the numbers and types of PFAM domains they contain.</li>
  <li><code>GO_distances.csv.gz</code>: Dataframe of the number of GO functions shared by every pair of proteins in <code>enriched_GO_proteins.csv</code></li>
  <li><code>GO_PFAM_pairwise_share.csv.gz</code>Same as <code>GO_distances.csv.gz</code>, but contains another column for the number of PFAM domains shared by every pair of proteins.</li>
  <li><code>viral_taxids.txt</code>: List of taxonomic IDs of the viruses that predominantly infect humans extracted from the "Taxid A" column of <code>Processed/virus_human.csv</code>.</li>
  <li><code>viral_taxonomies.csv</code>: Full taxonomic classifications for the viral taxids in <code>Processed/viral_taxids.txt</code>, obtained using the <code>Taxonomizr</code> package in R.</li>
  <li><code>model_predictors.csv</code>: Dataframe of predictors (genome size, number of publications, and number of unique interactions with human proteins) for a linear model.</li>
  <li><code>model_predictors_human.csv</code>: Same as <code>model_predictors.csv</code>, but excluding viruses that don't normally infect humans. Animal viruses are included in the dataset because many that are closely related to human viruses have been studied in laboratory settings.</li>
</ul>

## Python scripts:

<ul>
  <li><code>01_intact.py</code>: Removes unnecessary columns from <code>Data/virus_intact.txt</code>, renames columns to make them easier to work with, keeps only virus:human interactions, adds taxonomy names from NCBI taxonomy IDs, and switches the columns so that the human protein always appears second (this makes later analyses easier). This script creates <code>Processed/virus_human.csv</code>.</li>
  <li><code>02_GO_files.py</code>: Creates <code>Data/human_genes.txt</code> <code>Data/human_genes_all.txt</code> to search on GOrilla for gene enrichment.</li>
  <li><code>03_download_hmmer.py</code>: Downloads PFAM domain information from <a href="https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan" target="_blank">HMMER</a> using <code>wget</code>. </li>
  <li><code>04_process_hmmer.py</code>: Parses the tsv files downloaded by <code>03_download_hmmer.py</code> and combines them into a single file for each batch of query proteins. All the batches were combined later in <code>PFAM.ipynb</code> to generate <code>Processed/protein_domains.csv</code>.</li>
</ul>

## iPython Notebooks:

<ul>
  <li><code>01_protein_enrichment.ipynb</code>: Creates some summary plots of the GOrilla results in <code>Processed/GO_function</code> and <code>Processed/GO_process</code> and computes pairwise overlap in GO function between genes associated with the most highly enriched GO functions. </li>
  <li><code>02_PFAM.ipynb</code>: Combines PFAM domains from the protein batches into <code>protein_domains.csv</code>. Performs (a very uninformative) principal component analysis and correlates the number of GO functions and PFAM domains shared between every pair of proteins in the filtered dataset.</li>
  <li><code>03_virus_groups.ipynb</code>: Creates <code>model_predictors.csv</code> and <code>model_predictors_human.csv</code> and performs a logit ordinal regression to predict viral severity group (an ordinal variable) from the predictors in <code>model_predictors_human.csv</code>. </li>
</ul>
