import os
from cellpose import models, io
from skimage import morphology
from glob import glob
from tqdm import tqdm
import numpy as np

model = models.CellposeModel(gpu=False, pretrained_model="models/cellpose/cells")
files = glob("data/intensity_scaled/seq_2nt/rep*/*MOI/*hpi/nuclei-fov*-c0-r0-z0.tiff")

def postprocess_segmentation(masks, min_segment_size=10000):
	processed_image = morphology.remove_small_objects(masks, min_size=min_segment_size, connectivity=2)
	closed_image = morphology.closing(processed_image, footprint=np.ones((5, 5)))
	area_closed_image = morphology.area_closing(closed_image, area_threshold=200)
	return area_closed_image

for file in tqdm(files):
	img = io.imread(file)
	masks, flows, styles = model.eval(img, diameter=None, flow_threshold = 0.7, cellprob_threshold=-4.0)
	pp_masks = postprocess_segmentation(masks)
	outfile = file.replace("data/intensity_scaled", "analysis/segmentation").replace("nuclei-fov_00","fov_").replace("-c0-r0-z0.tiff", "_cp_cells.png")
	print(outfile)
	os.makedirs(os.path.dirname(outfile), exist_ok=True)
	io.imsave(outfile, pp_masks)