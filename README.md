# mudRapp-seq

This repository contains the code that accompanies the paper introducing 
"*In-situ sequencing with multiple direct RNA-assisted padlock probing (mudRapp-seq)*":

> Ahmad S, Schaust J, Ambil UB, Fischer SC, Ankenbrand MJ, Smyth RP. *Localization of influenza viral RNA in cells by direct RNA padlock probing and in-situ sequencing.* (in preparation)

Raw data is archived independently in the [BioImage Archive (BIA)](https://www.ebi.ac.uk/bioimage-archive/) with accession number **TODO**.
In order to reproduce our analyses, download the raw data from BIA and put them into the `data/raw` folder.

To re-create the python environment with [`mamba`](https://github.com/mamba-org/mamba) run:

```bash
mamba env create -f environment.yml
mamba activate mudRapp-seq
```