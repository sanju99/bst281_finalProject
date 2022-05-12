# Final Project for BST 281, Spring 2022

## Directories contained in this repository:

<b>Data</b>: Raw data downloaded from various databases

<b>Processed</b>: Tidy data files for analysis

<b>Analysis</b>: iPython notebooks used for analyses

## Files contained in this repository

### Data

<ul>
  <li><code>virus_intact.txt</code>: database of 42,154 pairwise interactions between virus and host cell proteins. This computaionally-maintained database was downloaded from <a href="https://www.ebi.ac.uk/intact/download/datasets#computationally" target="_blank">IntAct</a> via <a href="https://www.ebi.ac.uk/intact/query/annot:%22dataset:virus%22?conversationContext=7" target="_blank">this</a> search query generated by the website to point to the virus database and then selecting <b>MI-TAB 2.5</b> format for downloading</li>
  <li><code>virus_genome_sizes.csv</code>: File of viral genome sizes downloaded from the <a href="https://www.ncbi.nlm.nih.gov/genome/browse/#!/viruses/" target="_blank">NCBI</a>. </li>
  <li><code>taxids.txt</code>: List of taxonomic IDs from <code>virus_intact.txt</code> to search to get names for. <code>virus_intact.txt</code> only contains taxonomic IDs, and it is helpful to filter by names.</li>
  <li><code>taxids_with_names.txt</code>: Table of taxonomic IDs with associated names, obtained by searching the <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi" target="_blank">NCBI taxonomy browser</a>.</li>
  <li><code>human_proteins.txt</code>: List of all human proteins from the virus IntAct database</li>
  <li><code>human_proteins_fixed.tab</code>: Table of updated human proteins. It was obtained by searching <code>human_proteins.txt</code> on the UniProt website to update and clean deprecated names.</li>
  <li><code>human_genes.txt</code>: Table of human genes from the host:virus database. Gene names were extracted from <code>human_proteins_fixed.tab</code>.<li>
  <li><code>human_proteins_all.gz</code>: Table of all ~79,000 proteins in the human proteome downloaded from <a href="https://www.uniprot.org/uniprot/?query=proteome:UP000005640" target="_blank">UniProt</a>.</li>
  <li><code>human_genes_all.txt</code>: List of all human genes to use as the background list for GOrilla enrichment searching. Gene names were extracted from <code>human_proteins_all.gz</code>.</li>
  <li><code>GO_function.tsv</code>: GOrilla function results returned from searching <code>human_genes.txt</code> and <code>human_genes_all.txt</code> using the "Two unranked lists of genes" feature</li>
  <li><code>GO_process.tsv</code>: GOrilla process results returned from searching <code>human_genes.txt</code> and <code>human_genes_all.txt</code> using the "Two unranked lists of genes" feature</li>
</ul>

### Processed:

<ul>
  <li><code>virus_human.csv</code>This file includes columns for taxonomies that were added in later steps
