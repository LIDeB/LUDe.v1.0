📌About Us

We are a drug discovery team with an interest in the development of publicly available open-source customizable cheminformatics tools
to be used in computer-assisted drug discovery. We belong to the Laboratory of Bioactive Research and Development (LIDeB) of the
National University of La Plata (UNLP), Argentina. Our research group is focused on computer-guided drug repurposing and rational discovery
of new drug candidates to treat epilepsy and neglected tropical diseases.

💻Web Site https://lideb.biol.unlp.edu.ar



-------------------------------------------------------------------------------------------------

**LUDe WebApp**

LUDe (LIDEB’s Useful Decoys) is a WebApp that generates, from a set of active compounds, decoys (putative inactive compounds) which
can be used to retrospectively validate virtual screening tools/protocols. Decoys are molecules that have not been tested against
a molecular target of interest but due to their structural features are presumably not prone to bind the target with high affinity.
LUDe finds decoys in a curated ChEMBL27 database; these decoys are paired with the known active compounds in relation to certain
general physicochemical properties (e.g., molecular weight, log P, and others) but are topologically different from the query compounds.
LUDe is conceptually similar to the Directory of Useful Decoys enhanced, but additional filters have been serially implemented
to assure the topological dissimilarity between the decoys and the active compounds used as input.

In this Web App, decoys are obtained through four sequential steps:

In this WebApp, decoys are obtained by four sequential steps:

1) Searching molecules with similar physicochemical properties of the input active molecules in a curated ChEMBL database.
2) Filtering the selected molecules by dissimilarity against each individual input molecule.
3) Randomly selecting a desired number of decoys for each individual input molecule.
4) Filtering the selected molecules by the dissimilarity against all the input molecules (set of active compounds used as query). The decoys will have low Tanimoto similarity with the input compounds, and also low Maximum Common Substructure (MCS) ratio and distinctive molecular frameworks (Murcko scaffold). All in all, these three serial filters assure that the decoys will have a different topology in relation to the active compound used as a query. Furthermore, dissimilarity to the remaining active compounds is also checked. Finally, you can download a file with your decoys.

If you are looking to contact us, please send a mail to lideb@biol.unlp.edu.ar or contact us by Twitter (https://twitter.com/LIDeB_UNLP)
