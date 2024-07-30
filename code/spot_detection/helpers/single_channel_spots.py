import os
import numpy as np
import skimage.io
import scipy.ndimage as ndi
import starfish.image
from starfish.types import Levels
from starfish.spots import FindSpots
from starfish.core.spots.DecodeSpots.trace_builders import build_spot_traces_exact_match

def spot_detection_segmentation(exp, outdir):
	for fov_index, fov in enumerate(exp.fovs()):
		imgs = fov.get_image("primary")

		# Tophat filtering
		masking_radius = 6
		filt = starfish.image.Filter.WhiteTophat(masking_radius, is_volume=False)
		filtered = filt.run(imgs, verbose=True, in_place=False)

		# Normalization
		cptz= starfish.image.Filter.ClipPercentileToZero(p_min=80, p_max=100, level_method=Levels.SCALE_BY_CHUNK)
		clipped_both_scaled = cptz.run(filtered, in_place=False)

		# Spot detection
		bd = FindSpots.BlobDetector(
			min_sigma=2,
			max_sigma=8,
			num_sigma=7,
			threshold=.05,
			is_volume=True,
			measurement_type='max',
		)

		spots = bd.run(imgs, reference_image=clipped_both_scaled)
		spot_table = build_spot_traces_exact_match(spots)

		feature_df = spot_table.to_features_dataframe()

		# Segmentation
		segdir = outdir.replace("spot_detection", "segmentation")
		nucleus_mask = skimage.io.imread(f"{segdir}/fov_{fov_index}_cp_masks.png")
		cell_mask = skimage.io.imread(f"{segdir}/fov_{fov_index}_cpws_cells.png")

		# Enrich results
		feature_df["nucleus"] = nucleus_mask[feature_df.y, feature_df.x]>0
		distance_to_nucleus = ndi.distance_transform_edt(~(nucleus_mask>0))-ndi.distance_transform_edt(nucleus_mask>0)
		feature_df["nucleus_dist"] = distance_to_nucleus[feature_df.y, feature_df.x]
		feature_df["cell"] = cell_mask[feature_df.y, feature_df.x]
		dist_to_boundary = ndi.distance_transform_edt(~skimage.segmentation.find_boundaries(cell_mask))
		feature_df["boundary_dist"] = dist_to_boundary[feature_df.y, feature_df.x]
		border_cells = np.unique(np.concatenate([cell_mask[:,0], cell_mask[:,-1], cell_mask[0,:], cell_mask[-1,:]]))
		feature_df["border_cell"] = [c in border_cells for c in feature_df["cell"]]

		# Save results
		os.makedirs(outdir, exist_ok=True)
		prefix = f"{outdir}/fov_{fov_index}"
		feature_df.to_csv(f"{prefix}_spots.csv", na_rep=None, index=False)
		clipped_both_scaled.xarray.to_netcdf(f"{prefix}_primary.nc")

