import os
import numpy as np
import pandas as pd
import skimage.io
import scipy.ndimage as ndi
from starfish.types import Axes
from starfish.types import Levels
from starfish.spots import FindSpots
from starfish.core.spots.DecodeSpots.trace_builders import build_spot_traces_exact_match
from starfish.morphology import Binarize, Filter, Merge, Segment
from skimage.filters import threshold_otsu
import starfish.image
from starfish.types import FunctionSource


def spot_detection_segmentation(exp, outdir):
	for fov_index, fov in enumerate(exp.fovs()):
		imgs = fov.get_image("primary")
		nuclei = fov.get_image("nuclei")


		masking_radius = 6
		filt = starfish.image.Filter.WhiteTophat(masking_radius, is_volume=False)
		filtered = filt.run(imgs, verbose=True, in_place=False)

		PR8 = filtered.sel({Axes.CH:2})
		StPt = filtered.sel({Axes.CH:3})

		cptz= starfish.image.Filter.ClipPercentileToZero(p_min=80, p_max=99.9999, level_method=Levels.SCALE_BY_CHUNK)
		clipped_both_scaled_PR8 = cptz.run(PR8, in_place=False)
		clipped_both_scaled_StPt = cptz.run(StPt, in_place=False)

		bd = FindSpots.BlobDetector(
		min_sigma=4,
		max_sigma=8,
		num_sigma=5,
		threshold=.05,
		is_volume=True,
		measurement_type='max',
		)

		spots_PR8 = bd.run(imgs.sel({Axes.CH: 2}), reference_image=clipped_both_scaled_PR8)
		spots_StPt = bd.run(imgs.sel({Axes.CH: 3}), reference_image=clipped_both_scaled_StPt)

		spot_table_PR8 = build_spot_traces_exact_match(spots_PR8)
		spot_table_StPt = build_spot_traces_exact_match(spots_StPt)
		st_PR8 = spot_table_PR8.to_features_dataframe()
		st_PR8['strain'] = 'PR8'
		st_StPt = spot_table_StPt.to_features_dataframe()
		st_StPt['strain'] = 'StPt'
		feature_df = pd.concat([st_PR8, st_StPt])

		smoothed_nuclei = starfish.image.Filter.GaussianLowPass(sigma=5, is_volume=False).run(nuclei)
		dapi_thresh = threshold_otsu(smoothed_nuclei.xarray.values) # e.g. .048, binary mask for cell (nuclear) locations
		stain_thresh = 0#.0001  # binary mask for overall cells // binarization of stain TODO hard coded value - find better way to set this
		min_dist = 50
		hole_size_to_fill = 1000
		min_allowed_size = 5000
		max_allowed_size = 100000


		binarized_nuclei = Binarize.ThresholdBinarize(dapi_thresh).run(smoothed_nuclei)
		filled_nuclei = Filter.Map(FunctionSource.skimage("morphology.remove_small_holes"), hole_size_to_fill, module=None).run(binarized_nuclei)
		labeled_masks = Filter.MinDistanceLabel(min_dist, 1).run(filled_nuclei)
		watershed_markers = Filter.AreaFilter(min_area=min_allowed_size, max_area=max_allowed_size).run(labeled_masks)

		thresholded_stain = Binarize.ThresholdBinarize(stain_thresh).run(nuclei)
		markers_and_stain = Merge.SimpleMerge().run([thresholded_stain, watershed_markers])
		watershed_mask = Filter.Reduce(
		"logical_or",
		lambda shape: np.zeros(shape=shape, dtype=bool)
		).run(markers_and_stain)

		segmenter = Segment.WatershedSegment(connectivity=np.ones((1, 3, 3), dtype=bool))

		masks = segmenter.run(
		smoothed_nuclei,
		watershed_markers,
		watershed_mask,
		)

		nucleus_mask = watershed_markers.to_label_image().xarray.squeeze(Axes.ZPLANE.value).values
		feature_df["nucleus"] = nucleus_mask[feature_df.y, feature_df.x]>0
		distance_to_nucleus = ndi.distance_transform_edt(~(nucleus_mask>0))-ndi.distance_transform_edt(nucleus_mask>0)
		feature_df["nucleus_dist"] = distance_to_nucleus[feature_df.y, feature_df.x]
		cell_mask = masks.to_label_image().xarray.squeeze(Axes.ZPLANE.value).values
		feature_df["cell"] = cell_mask[feature_df.y, feature_df.x]
		dist_to_boundary = ndi.distance_transform_edt(~skimage.segmentation.find_boundaries(cell_mask))
		feature_df["boundary_dist"] = dist_to_boundary[feature_df.y, feature_df.x]
		border_cells = np.unique(np.concatenate([cell_mask[:,0], cell_mask[:,-1], cell_mask[0,:], cell_mask[-1,:]]))
		feature_df["border_cell"] = [c in border_cells for c in feature_df["cell"]]

		os.makedirs(outdir, exist_ok=True)
		prefix = f"{outdir}/fov_{fov_index}"
		skimage.io.imsave(f"{prefix}_nuclei.png", nucleus_mask)
		skimage.io.imsave(f"{prefix}_cells.png", cell_mask)
		feature_df.to_csv(f"{prefix}_spots.csv", na_rep=None, index=False)
		spot_table_PR8.to_netcdf(f"{prefix}_spots_PR8.nc")
		spot_table_StPt.to_netcdf(f"{prefix}_spots_StPt.nc")
		clipped_both_scaled_PR8.xarray.to_netcdf(f"{prefix}_primary_PR8.nc")
		clipped_both_scaled_StPt.xarray.to_netcdf(f"{prefix}_primary_StPt.nc")
