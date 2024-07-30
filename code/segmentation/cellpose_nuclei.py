import os
from cellpose import models, io
import skimage as ski
from glob import glob
from tqdm import tqdm
import numpy as np

model = models.CellposeModel(gpu=False, pretrained_model="models/cellpose/nuclei")
files = glob("data/spacetx/seq_2nt/rep*/*MOI/*hpi/nuclei-fov*-c0-r0-z0.tiff")

for file in tqdm(files):
	img = io.imread(file)
	masks, flows, styles = model.eval(img, diameter=None)
	outfile = file.replace("data/spacetx", "analysis/segmentation").replace("nuclei-fov_00","fov_").replace("-c0-r0-z0.tiff", "")
	os.makedirs(os.path.dirname(outfile), exist_ok=True)
	io.save_masks(img[:,:,np.newaxis], masks, flows, outfile, png=True, save_txt=False)