# mudRapp-seq

This repository contains the code that accompanies the paper introducing 
"*In-situ sequencing with multiple direct RNA-assisted padlock probing (mudRapp-seq)*":

> Ahmad S, Schaust J, Ambil UB, Fischer SC, Ankenbrand MJ, Smyth RP. *Localization of influenza viral RNA in cells by direct RNA padlock probing and in-situ sequencing.* (in preparation)

## Data

Raw data is archived independently in the [BioImage Archive (BIA)](https://www.ebi.ac.uk/bioimage-archive/) with accession number **TODO**.
In order to reproduce our analyses, download the raw data from BIA and put them into the `data/raw` folder.

## Analyses

### Python environment

To re-create the python environments with [`mamba`](https://github.com/mamba-org/mamba) run:

```bash
mamba env create -f envs/starfish.yml # mudRapp-seq-starfish
mamba env create -f envs/cellpose.yml # mudRapp-seq-cellpose
```

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

