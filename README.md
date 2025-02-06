 # About Us

📌 We are a drug discovery team with an interest in the development of publicly available open-source customizable cheminformatics tools
to be used in computer-assisted drug discovery. We belong to the Laboratory of Bioactive Research and Development (LIDeB) of the
National University of La Plata (UNLP), Argentina. Our research group is focused on computer-guided drug repurposing and rational discovery
of new drug candidates to treat epilepsy and neglected tropical diseases.

💻 [LIDeB Web Site](https://lideb.biol.unlp.edu.ar)


-------------------------------------------------------------------------------------------------

## LUDe WebApp

LUDe (LIDEB’s Useful Decoys) is a WebApp that generates, from a set of active compounds, decoys (putative inactive compounds) which
can be used to retrospectively validate virtual screening tools/protocols. Decoys are molecules that have not been tested against
a molecular target of interest but due to their structural features are presumably not prone to bind the target with high affinity.
LUDe finds decoys in a curated ChEMBL30 database; these decoys are paired with the known active compounds in relation to certain
general physicochemical properties (e.g., molecular weight, log P, and others) but are topologically different from the query compounds.
LUDe is conceptually similar to the Directory of Useful Decoys enhanced, but additional filters have been serially implemented
to assure the topological dissimilarity between the decoys and the active compounds used as input.

In this Web App, decoys are obtained through four sequential steps:

1) Searching molecules with similar physicochemical properties of the input active molecules in a curated ChEMBL database.
2) Filtering the selected molecules by dissimilarity against each individual input molecule.
3) Randomly selecting a desired number of decoys for each individual input molecule.
4) Filtering the selected molecules by the dissimilarity against all the input molecules (set of active compounds used as query). The decoys will have low Tanimoto similarity with the input compounds, and also low Maximum Common Substructure (MCS) ratio and distinctive molecular frameworks (Murcko scaffold). All in all, these three serial filters assure that the decoys will have a different topology in relation to the active compound used as a query. Furthermore, dissimilarity to the remaining active compounds is also checked. Finally, you can download a file with your decoys.

![LUDe_worflow](https://github.com/LIDeB/LUDe.v1.0/blob/main/workflow_lude.png?raw=true)


## Installation

Download the LUDe.py file, the Lude-env.yml file and the Base_ChEMBL_final folder (which contains the ChEMBL databases).

The Lude-env.yml file includes all the required dependencies for easily installation. You can set up the environment using Conda:

```
conda env create -f Lude-env.yml
conda activate Lude-env
```


## Running LUDe

Before running LUDe.py, you need to modify specific lines in the script related to folder paths (in the Queries and Folder Paths section). Open the file in a text editor and update the following variables (located between lines 50 and 57) with your file names and paths:

```
smiles_dataset: Name of the .txt file containing the queries (one SMILES per line).
directory: Path to the directory containing the SMILES dataset.
directory_chembl: Path to the directory containing the ChEMBL datasets.
```

If you need to modify any default program settings (e.g., physicochemical feature limits or dissimilarity conditions), you can adjust the parameters in the Customizable Options section (lines 64 to 87) or leave them as they are.

To run LUDe, after making all the modifications to the script, open a command prompt, navigate to the script's folder, and type:

```
python LUDe.py
```
## Contact

📧 If you are looking to contact us, please send a mail to lideb@biol.unlp.edu.ar or contact us by [Twitter](https://twitter.com/LIDeB_UNLP)
