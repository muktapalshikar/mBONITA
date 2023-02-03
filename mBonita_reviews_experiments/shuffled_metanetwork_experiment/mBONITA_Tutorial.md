![GitHub last commit](https://img.shields.io/github/last-commit/mgp13/moBONITA?style=for-the-badge)

# mBONITA: multi-omics Boolean Omics Network Invariant-Time Analysis

## Authors:
Mukta G. Palshikar, Xiaojun Min, Alexander Crystal, Jiayue Meng, Shannon P. Hilchey, Martin Zand, Juilee Thakar

## Abstract:

Multi-omics profiling provides a holistic picture of a condition being examined and capture the complexity of signaling events, beginning from the original cause (environmental or genetic), to downstream functional changes at multiple molecular layers. Pathway enrichment analysis has been used with multi-omics datasets to characterize signaling mechanisms. However, technical and biological variability between these layered datasets are challenges for integrative computational analyses. We present a Boolean network-based method, multi-omics Boolean Omics Network Invariant-Time Analysis (mBONITA) to integrate omics datasets that quantify multiple molecular layers. mBONITA utilizes prior knowledge networks to perform topology-based pathway analysis. In addition, mBONITA identifies genes that are consistently modulated across molecular measurements by combining observed fold-changes and variance with a measure of node (i.e., gene or protein) influence over signaling, and a measure of the strength of evidence for that gene across datasets. We used mBONITA to integrate multi-omics datasets from RAMOS B cells treated with the immunosuppressant drug cyclosporine A under varying oxygen tensions to identify pathways involved in hypoxia-mediated chemotaxis. We compare mBONITA's performance with 6 other pathway analysis methods designed for multi-omics data and show that mBONITA identifies a set of pathways with evidence of modulation across all omics layers.

![Graphical abstract - Light mode](https://github.com/mgp13/moBONITA/blob/main/Picture1.png?raw=true#gh-light-mode-only)
![Graphical abstract - Light mode](https://github.com/mgp13/moBONITA/blob/main/Picture2.png?raw=true#gh-dark-mode-only)

## mBONITA tutorial
multi-omics Boolean Omics Network Invariant-Time Analysis (mBONITA) is a method (and corresponding software package) that builds off our previous work BONITA  


### Requirements

The *mBONITA* tool is written in Python3 and C. I strongly recommend that mBONITA be run on a computing cluster such as the University of Rochester's BlueHive, and that jobs are submitted using a scheduler such as SLURM. Dependencies are listed in the conda environment file (SPECIFY FILENAME HERE).

**Minor caveat** - *mBONITA* is not a Python package like numpy or scipy, which allow users to import individual functions and (re)use them in custom code. mBONITA is an all-in-one pipeline that doesn't allow function import or much customization beyond the pre-specified parameters. I welcome advanced users to modify code and submit pull requests, but this is beyond what most users will need. 

*mBONITA* requires the following inputs (Step 0):

- A pre-preprocessed multi-omics dataset from matched samples, prepared in a combined matrix format as in (link to Python notebook here)
- A conditions file in matrix format, which specfies the experimental conditions for each sample in the training dataset above
- A contrast file that specifies which pairs of experimental conditions are to be compared during pathway analysis

to perform the following tasks:

- Download and prepare KEGG pathways for pathway analysis (Step 1)
- Infer Boolean regulatory/signaling rules for KEGG pathways using the combined multi-omics dataset (Step 2)
- Perform topology-informed pathway analysis for user-specified pairs of experimental conditions (Step 3)

This tutorial will go through the mBONITA pipeline using a multi-omics dataset of transcriptomics, proteomics, and phosphoproteomics from RAMOS B cells, as described in the mBONITA publication.

### Step 0: Process multi-omics data and generate conditions and contrast files

I expect that most users will begin with 2 or more processed datasets from separate multi-omics datasets. These datasets will usually be log2-normalized. The R notebook (INSERT LINK) outlines how to combine log2-normalized proteomics, phosphoproteomics and transcriptomics datasets as in the mBONITA publication and prepare them in a matrix format for mBONITA.

mBONITA also requires a condition and contrast file for pathway analysis. An example of how to prepare these files is in (INSERT LINK).

Briefly, if your dataset looks something like this:

| Genes | Condition1_replicate1_proteomics  | Condition1_replicate2_proteomics | Condition2_replicate1_proteomics  | Condition2_replicate2_proteomics | Condition1_replicate1_phosphoproteomics | Condition1_replicate2_phosphoproteomics | Condition2_replicate1_phosphoproteomics | Condition2_replicate2_phosphoproteomics |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| Gene1 | - | - | - | - | - | - | - | - |
| Gene2  | - | - | - | - | - | - | - | - |
| Gene3  | - | - | - | - | - | - | - | - |
| Gene4  | - | - | - | - | - | - | - | - |

Then your condition file will look like this:

| Sample |  Condition1 | Condition2  | 
| ------------- | ------------- | ------------- | 
| Condition1_replicate1_proteomics | 1  | 0  | 
| Condition1_replicate2_proteomics  | 1  | 0  |
| Condition2_replicate1_proteomics | 0  | 1  | 
| Condition2_replicate2_proteomics  | 0  | 1  |
| Condition1_replicate1_phosphoproteomics | 1  | 0  | 
| Condition1_replicate2_phosphoproteomics  | 1  | 0  |
| Condition2_replicate1_phosphoproteomics | 0  | 1  | 
| Condition2_replicate2_phosphoproteomics  | 0  | 1  |


And your contrast file will look like this:

|  Condition1 | Condition2  |

## Step 1: Download and prepare KEGG pathways for pathway analysis

## Step 2: Infer Boolean regulatory/signaling rules for KEGG pathways using the combined multi-omics dataset

Simply run the script find_rules_pathway_analysis.sh which will automatically submit appropriate jobs to a SLURM queue:

```bash find_rules_pathway_analysis.sh```

Please note that these scripts are written for SLURM. find_rules_pathway_analysis.sh loops over all networks to execute the script calcNodeImportance.sh, which in turn executes the Python script INSERT NAME HERE. I'm open to writing these scripts for other job scheduling managers. The Python script can also be run by itself on a desktop, but I advise doing this only for small networks/training datasets.

## Step 3: Perform topology-informed pathway analysis for user-specified pairs of experimental conditions

Run the Python script pathway_analysis_score_pathways_mBonita.py with the following parameters. An example is listed below. 

- path to training dataset file (concatenated)
- conditions file
- contrast file

For file formats, please refer to Step 0.

Here is an example command:

```python3 pathway_analysis_score_pathways_mBonita.py concatenated_datasets.csv concatenated_conditions.csv contrasts.csv -sep ,```
