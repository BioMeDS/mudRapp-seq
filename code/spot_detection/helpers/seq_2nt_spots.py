import os
import numpy as np
import skimage.io
import scipy.ndimage as ndi
from starfish import Experiment
from starfish.image import ApplyTransform, LearnTransform
from starfish.types import Axes
from starfish.image import Filter
from starfish.types import Levels
from starfish.spots import FindSpots
from starfish.spots import DecodeSpots

def common_roi(transforms_list):
    # crop the images to the common region
    translations = np.array([x[2].translation for x in transforms_list.transforms])
    lower = -np.floor(translations.min(axis=0)).astype(int)
    upper = -np.ceil(translations.max(axis=0)).astype(int)
    # the `or None` is to handle the case where no cropping is necessary at the upper bound (see: https://stackoverflow.com/a/11337953/4969760)
    return slice(lower[0],upper[0] or None), slice(lower[1],upper[1] or None)#Liste des Prozentualen Verlustes von Spots in der min_sigma range 1.0 bis 2.0:

def analyze_fov(rep, moi, hpi, fov_id):
	exp = Experiment.from_json(os.path.join(f"data/spacetx/seq_2nt/rep{rep}/{moi}MOI/{hpi}hpi", "experiment.json"))
	fov = exp.fovs()[fov_id]
	imgs = fov.get_image("primary")
	nuclei = fov.get_image("nuclei")
	# Spot based registration
	spots = imgs.reduce({Axes.CH}, func="max")
	learn_translation = LearnTransform.Translation(reference_stack=spots.sel({Axes.ROUND: 0}), axes=Axes.ROUND, upsampling=1000)
	warp = ApplyTransform.Warp()
	transforms_list = learn_translation.run(spots)
	registered_imgs = warp.run(imgs, transforms_list=transforms_list)
	transforms_list = learn_translation.run(spots)
	registered_nuclei = warp.run(nuclei, transforms_list=transforms_list)
	transforms_list = learn_translation.run(spots)
	xroi, yroi = common_roi(transforms_list)
	registered_imgs = registered_imgs.sel({Axes.X: xroi, Axes.Y: yroi})
	registered_nuclei = registered_nuclei.sel({Axes.X: xroi, Axes.Y: yroi})
	filtered = Filter.WhiteTophat(masking_radius=6, is_volume=False).run(registered_imgs, in_place=False)
	bleed = np.array(
		[[ 1.  ,  0.  , -0.656,  0.  ],
		[ 0.  ,  1.  ,  0.   ,  0.  ],
		[ 0.  ,  0.  ,  1.   ,  0.  ],
		[ 0.  ,  0.  ,  0.   ,  1.  ]]
	)
	bleed_corrected = Filter.LinearUnmixing(bleed).run(filtered, in_place=False)
	# Normalization
	matched_histograms = Filter.MatchHistograms({Axes.CH})
	scaled_c = matched_histograms.run(bleed_corrected, in_place=False, verbose=False) #n_processed =None)
	cptz= Filter.ClipPercentileToZero(p_min=80, p_max=99.9995, level_method=Levels.SCALE_BY_CHUNK, group_by={Axes.CH})
	clipped_both_scaled = cptz.run(scaled_c, in_place=False)
	# Bandpass filtering
	bandpassed = Filter.Bandpass(lshort=0, llong=17, level_method=Levels.SCALE_BY_CHUNK).run(clipped_both_scaled, in_place=False)
	# further suppress high frequency background noise (in particular in G channel at low hpis)
	median_filtered = bandpassed.apply(skimage.filters.median, in_place=False, selem=skimage.morphology.disk(3))

	# Spot Detection
	dots = median_filtered.reduce([Axes.CH], "max").reduce({Axes.ROUND}, func="mean")
	bd = FindSpots.BlobDetector(
		min_sigma=3,
		max_sigma=10,
		num_sigma=8,
		threshold=.03,
		is_volume=True,
		measurement_type='mean',
	)
	spots = bd.run(median_filtered, reference_image=dots)
	# Spot Decoding
	multi_barcode_decoder = DecodeSpots.MultiBarcodeDecoder(
		codebook=exp.codebook,
		max_distance=.5,
		min_intensity=.1,
		raw_intensity_threshold=.015,
		return_original_intensities=True
	)
	mbd_decoded_spots = multi_barcode_decoder.run(spots=spots)
	# Segmentation
	nucleus_mask = skimage.io.imread(f"analysis/segmentation/seq_2nt/rep{rep}/{moi}MOI/{hpi}hpi/fov_{fov_id}_cp_masks.png")
	nucleus_mask = nucleus_mask[yroi, xroi]
	cell_mask = skimage.io.imread(f"analysis/segmentation/seq_2nt/rep{rep}/{moi}MOI/{hpi}hpi/fov_{fov_id}_cpmc_cells.png")
	cell_mask = cell_mask[yroi, xroi]
	# Enrich Results
	feature_df = mbd_decoded_spots.to_features_dataframe()
	feature_df["nucleus"] = nucleus_mask[feature_df.y, feature_df.x]>0
	distance_to_nucleus = ndi.distance_transform_edt(~(nucleus_mask>0))-ndi.distance_transform_edt(nucleus_mask>0)
	feature_df["nucleus_dist"] = distance_to_nucleus[feature_df.y, feature_df.x]
	feature_df["cell"] = cell_mask[feature_df.y, feature_df.x]
	dist_to_boundary = ndi.distance_transform_edt(~skimage.segmentation.find_boundaries(cell_mask))
	feature_df["boundary_dist"] = dist_to_boundary[feature_df.y, feature_df.x]
	border_cells = np.unique(np.concatenate([cell_mask[:,0], cell_mask[:,-1], cell_mask[0,:], cell_mask[-1,:]]))
	feature_df["border_cell"] = [c in border_cells for c in feature_df["cell"]]
	ambiguous = (mbd_decoded_spots.binarized.sum(dim=["c"]) > 1).sum(dim=["r"]) > 1
	feature_df["ambiguous"] = ambiguous.values
	# Save Results
	outdir = f"analysis/spot_detection/seq_2nt/rep{rep}/{moi}MOI/{hpi}hpi"
	os.makedirs(outdir, exist_ok=True)
	prefix = f"{outdir}/fov_{fov_id}"
	feature_df.to_csv(f"{prefix}_spots.csv", na_rep=None, index=False)
	mbd_decoded_spots.to_netcdf(f"{prefix}_spots.nc")
	median_filtered.xarray.to_netcdf(f"{prefix}_primary.nc")
	transforms_list.to_json(f"{prefix}_registration.json")
