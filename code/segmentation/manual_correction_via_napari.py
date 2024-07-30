#!/usr/bin/env python
# coding: utf-8

import os

from glob import glob
from skimage.io import imread
import napari
from starfish import Experiment
import argparse

parser = argparse.ArgumentParser(description='Open cellpose segmentation results in napari for manual correction')
parser.add_argument('--rep', type=int, default=1, help='replicate number')
parser.add_argument('--moi', type=str, default="0.3MOI", help='magnitude of infection')
parser.add_argument('--hpi', type=int, default=5, help='hours post infection')
parser.add_argument('--fov_index', type=int, default=0, help='index of the field of view to inspect')
args = parser.parse_args()

rep = args.rep
moi = args.moi
hpi = args.hpi
fov_index = args.fov_index
data_folder = f"data/spacetx/seq_2nt/rep{rep}/{moi}/{hpi}hpi"
segmentation_folder = data_folder.replace("data/spacetx", "analysis/segmentation")
prefix = f"{segmentation_folder}/fov_{fov_index}"

# ## Load raw data

exp = Experiment.from_json(os.path.join(data_folder, "experiment.json"))
fov = exp.fovs()[fov_index]
imgs = fov.get_image("primary")
nuclei = fov.get_image("nuclei")
bf = fov.get_image("brightfield")

# ## Load segmentation results

cp_nucleus_mask = imread(f"{prefix}_cp_masks.png")
moi = moi.replace("1.0", "1")
cp_cell_mask = imread(f"{prefix}_cp_cells.png")
# check if manual correction file exists
if os.path.exists(f"{prefix}_cpmc_cells.png"):
    cpmc_cell_mask = imread(f"{prefix}_cpmc_cells.png")
else:
    print(f"File {prefix}_cpmc_cells.png does not exist. Initializing new manual correction from cellpose mask.")
    cpmc_cell_mask = cp_cell_mask.copy()

# ## View results in napari

viewer = napari.Viewer()
for i,b,c in zip(range(3,-1,-1), "AGTC"[::-1], ["magenta", "green", "yellow", "red"][::-1]):
    viewer.add_image(imgs.xarray[:,i], name=b+" (raw)", colormap=c, blending="additive", visible=False, contrast_limits=[0, 0.005])
viewer.add_image(bf.xarray[:,0], name="brightfield", blending="additive", visible=False)
viewer.add_image(nuclei.xarray[:,0], name="dapi", blending="additive", visible=True, contrast_limits=[0, 0.002])
viewer.layers[-1].contrast_limits_range = [0, 0.02]
viewer.add_labels(cp_nucleus_mask, name="nuclei (cp)", opacity=0.4, visible=False)
viewer.add_labels(cp_cell_mask, name="cells (cp)", opacity=0.3, visible=False)
viewer.add_labels(cpmc_cell_mask, name="cells (manually corrected)", opacity=0.4, visible=True)
viewer.dims.set_current_step(0, 0)

@viewer.bind_key('6')
def tg(viewer): viewer.layers[-1].contour = 10 - viewer.layers[-1].contour

napari.run()