# acid-amine-enumeration-2
Code for "Exploring the Combinatorial Explosion of Amine–Acid Reaction Space via Graph Editing".
Tested on the following package versions: ipython 7.16.1, jupyterlab 3.1.4, matplotlib 3.3.4, numpy 1.19.2, pandas 1.1.3, RDKit 2019.09.3, seaborn 0.11.1, tqdm 4.56.0, umap-learn 0.5.1, circos 0.69-9.

Circos downloaded from http://circos.ca/software/download/circos/

Drugbank structures to be downloaded from https://go.drugbank.com/releases and saved as structures.sdf. Version 5.1.8 is used in the manuscript.

To reproduce results and figures in the manuscript, download the repository and run all notebooks in ascending order of numeric prefixes. While it may be possible to run them out of order, certain code blocks will expect certain files produced by preceding notebooks, and fail if the file is absent, and a substitute is not provided.

Most intended output file sizes exceed Github's support, and are provided at https://drive.google.com/drive/u/1/folders/1ecw-hD_zBuTfuVpWSQYYc_10UhOi3Ul8.



The notebooks in the ./umap/ and ./late_stage_div/ folders will expect input data produced by notebooks in the parent directory, so running those first is recommended. Notebooks within these directories should also be run in ascending order of numeric prefixes.


A copy of the code that produces a significantly smaller dataset, and only requires one core, is present in the ./Demo/ folder. 
