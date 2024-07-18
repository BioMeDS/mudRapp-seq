# mudRapp-seq

This repository contains the code that accompanies the paper introducing 
"*Multiple direct RNA-assisted padlock probing in combination with in-situ sequencing (mudRapp-seq)*":

> Ahmad S, Schaust J, Ambil UB, Fischer SC, Ankenbrand MJ, Smyth RP. *Localization of influenza viral RNA in cells by direct RNA padlock probing and in-situ sequencing.* (in preparation)

## Data

Raw data is archived independently in the [BioImage Archive (BIA)](https://www.ebi.ac.uk/bioimage-archive/) with accession number **TODO**.
In order to reproduce our analyses, download the raw data from BIA and put them into the `data/raw` folder.

All final and some intermediate results are included in the repository to facilitate additional analyses without having to re-process all files from scratch.

## Computational environments

### Python environments

To re-create the python environments with [`mamba`](https://github.com/mamba-org/mamba) run:

```bash
mamba env create -f envs/starfish.yml # mudRapp-seq-starfish
mamba env create -f envs/cellpose.yml # mudRapp-seq-cellpose
```

### R environment

The R environment for this project is managed via [`renv`](https://rstudio.github.io/renv/articles/renv.html). A local environment is automatically created for you, when you run `R` or `Rscript` for the first time in the main project directory. This happens because of the `.Rprofile` file.

### jupyter setting (vscode)

The jupyter files are in sub-folders of `code/` but assume the kernel to run in the project root.
This is necessary to make the R kernel use the local renv, and allows to consistently use paths relative to the project root rather than the specific notebook location.
In VS Code this can be achieved by changing the setting `jupyter.notebookFileRoot` to `${workspaceFolder}`.
For `jupyter lab` there seems to be no simple solution at this moment (see [jupyterlab#11619](https://github.com/jupyterlab/jupyterlab/issues/11619)).

## Analyses

### Data formatting

As [`starfish`](https://github.com/spacetx/starfish) is used, the raw data needs to be restructured in SpaceTx format.

In order to create the formatted data in `data/spacetx` run these steps in the root of the `mudRapp-seq` repo and in the `mudRapp-seq` environment:

```bash
mamba run -n mudRapp-seq-starfish python code/data_formatting/cDNA_vRNA.py
```

### Segmentation

Images were separately segmented for nuclei and cell instances.
Nuclei segmentation is used to separate spots based on their location into nucleus and cytoplasm.
Cell segmentation is used to count spots per cell, filter infected cells and perform single cell analyses.

For nucleus segmentation, a cellpose model (`models/cellpose/nuclei`) was trained and applied to the raw dapi images (in spaceTx format).
The model was trained on a total of X images with human provided sparse labels (list of images for training in file X) **TODO @Joél**.

For cell instance segmentation, two different approaches were used:
1. Watershed of the dapi image with nuclei as seeds
2. A separate cellpose model with manual correction (details below)

The first strategy was used for most data, as it was deemed sufficient for filtering of infected cells and to calculate summary statistics like spots per cell.
However, for single cell analyses the cell borders were not reliable enough.

The separate cellpose model was trained on raw images without ICC (computational clearing by the microscope vendor).
Further, data was preprocessed with intensity scaling and combination of dapi and X channel in a single rgb image (see code in file X) **TODO @Joél**.
In order to maximize the number of correctly detected cells, the following parameters were used: `cell_proba=X`, `flow_threshold=X` and the resulting masks were post-processed, removing small objects and closing small holes and gaps (see code in file X) **TODO @Joél**.
Masks produced this way were manually corrected using label editing tools in `napari`.
Manual correction involved extending cells, shrinking cells, moving cell borders and adding new cells (starting with cell IDs of X) **TODO @Joél**.

#### Strategy 1 (nuclei via cellpose, cells via watershed)

This code performs nuclei segmentation with the cellpose model and watershed for cell segmentation (strategy 1).

```bash
mamba run -n mudRapp-seq-cellpose python code/segmentation/cellpose_nuclei_watershed_cells.py
```

### Spot detection

Spot detection is performed using starfish methods. The following command creates csv and netCDF files in `analysis/spot_detection`

```bash
mamba run -n mudRapp-seq-starfish python code/spot_detection/spot_detection.py
```

This creates separate spots files for each fov, they can be combined to a single `tsv.xz` file for each experiment using

```bash
Rscript code/spot_detection/combine_csvs.R
```

The result of this step is included in the repository:
- `analysis/spot_detection/cDNA_vRNA/all_spots.tsv.xz`

### Spot analysis

#### Sensitivity of cDNA vs direct vRNA probing

The main result is the much higher sensitivity for direct vRNA probing compared to cDNA probing.
For details, see the [analysis notebook](code/spot_analysis/cDNA_vRNA.ipynb).
These figures are included in the manuscript as Supp. Fig. 1b and Fig. 1c, respectively.

![Supplementary Figure 1b](figures/supp-fig1b-cDNA_vRNA-spot_counts-neg_ctrl.svg)
![Figure 1c](figures/fig1c-cDNA_vRNA-spot_counts.svg)

### Segment analysis

**TODO** describe TemporalCorrelation_mRNA_vRNA