# blood-tcr-pipeline
This pipeline takes either 10X cellranger outputs for two tissues (always needed to generate metadata) or an existing Seurat object stored in an .RDS file, finds matching
clones based on TCR chain sequences between these samples, either creates or
loads a Seurat object, creates a custom metadata table, then write this modified
Seurat object out to an .rds file.
Below are some example commands to run the pipeline.

---

## Unified Metadata Format
After running this pipeline, the modified metadata in the outputted Seurat v3
object accessible at `seurat_object_name@meta.data` will have the following
format:

| Cell Barcode | orig.ident | nCount_RNA | nFeature_RNA | RNA_snn_res.1.2 | seurat_clusters | matching | u1 | u1 | freq | Frequency |


---

## Mount DFCI Kraken Storage To Access Data
Run `mount_smbfs //jluber@shares.dfci.harvard.edu/singerlab ./dfci_share`.

---

## Generate custom TCR metadata

Run `python clone_pipeline.py /Users/jacobluber/dfci_share/jacob/data/K409/K409/K409LNVDJ/filtered_contig_annotations.csv /Users/jacobluber/dfci_share/jacob/data/K409/K409/K409bloodVDJ/filtered_contig_annotations.csv /Users/jacobluber/dfci_share/jacob/data/K409/K409/K409LNGEXouts/outs/analysis/tsne/2_components/projection.csv /Users/jacobluber/dfci_share/jacob/data/K409/K409/K409bloodGEXouts/outs/analysis/tsne/2_components/projection.csv /Users/jacobluber/dfci_share/jacob/data/K409/K409/K409LNVDJ/clonotypes.csv /Users/jacobluber/dfci_share/jacob/data/K409/K409/K409bloodVDJ/clonotypes.csv k409_ln-blood.csv k409_blood-ln.csv`. Note that `clone_pipeline.py` will generate "matching" data for both blood and tumor.

Below is a table explaining the arguments to `clone_pipeline.py`. In the example
the first tissue is tumor and the second tissue is blood. Arguments 3 and 4 are
deprecated (UMAP coordinates generated at later step) but still need to be included currently. All inputs are generated using either `cellranger count` or `cellranger vdj` commands.

| clone_pipeline.py argument # | file details |
|:----------------------------:|:------------:|
|1|first tissue TCR contig annotations|
|2|second tissue TCR contig annotations|
|3|first tissue GEX t-SNE coordinates|
|4|second tissue GEX t-SNE coordinates|
|5|first tissue TCR clonotype info|
|6|second tissue TCR clonotype info|
|7|output file name for TCR metadata file for first tissue|
|8|output file name for TCR metadata file for second tissue|

---

## Append additional TCR metadata to a Seurat object stored in an .rds file
Run `conda activate R` followed by `Rscript append_metadata_to_seurat_object.R /Users/jacobluber/browser/k409.ln-blood.rds k409_ln-blood.csv ~/singer_repos/blood-tcr-pipeline/ k409_test` where the first argument is the existing .rds file, the second argument is the output from generating the TCR metadata above, the third argument is the output path, and the fourth argument is the name of the output file excluding .rds at the end.

---

## Create New Seurat Object and Add TCR metadata
Run `conda activate R` followed by `Rscript make_seurat_object_from_scratch.R /Users/jacobluber/dfci_share/jacob/data/K409/K409/K409LNGEXouts/outs/filtered_gene_bc_matrices/GRCh38/ k409_ln-blood.csv ~/singer_repos/blood-tcr-pipeline/ k409_ln` where the first argument is the full path to the directory containing the cellranger GEX sparse expression matrix, and the remaining arguments are the same as the previous instructions for appending metadata.

This script needs to be run twice: once for blood and once for tumor.

## Load modified Seurat v3 .rds file into R
Run `seurat_object_name <- readRDS("k409_ln.rds")`.
