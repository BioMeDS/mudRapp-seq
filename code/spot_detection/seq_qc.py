import os
import numpy as np
import pandas as pd
import starfish
from starfish import Experiment, Codebook
from starfish.types import Axes, Levels
from starfish.image import ApplyTransform, LearnTransform, Filter
from starfish.spots import FindSpots, DecodeSpots

exp = Experiment.from_json(os.path.join("data/spacetx/seq_qc/rep0/A_PB2/", "experiment.json"))

outdir = "analysis/spot_detection/seq_qc/rep0/A_PB2"
os.makedirs(outdir, exist_ok=True)

def common_roi(transforms_list):
    # crop the images to the common region
    translations = np.array([x[2].translation for x in transforms_list.transforms])
    lower = -np.floor(translations.min(axis=0)).astype(int)
    upper = -np.ceil(translations.max(axis=0)).astype(int)
    # the `or None` is to handle the case where no cropping is necessary at the upper bound (see: https://stackoverflow.com/a/11337953/4969760)
    return slice(lower[0],upper[0] or None), slice(lower[1],upper[1] or None)#Liste des Prozentualen Verlustes von Spots in der min_sigma range 1.0 bis 2.0:

def get_decoded_spots_for_fov(fov_id, save=False):
    fov = exp.fovs()[fov_id]
    imgs = fov.get_image("primary")
    spots = imgs.reduce({Axes.CH}, func="max")
    learn_translation = LearnTransform.Translation(reference_stack=spots.sel({Axes.ROUND: 0}), axes=Axes.ROUND, upsampling=1000)
    transforms_list = learn_translation.run(spots)
    warp = ApplyTransform.Warp()
    registered_imgs = warp.run(imgs, transforms_list=transforms_list)
    xroi, yroi = common_roi(transforms_list)
    registered_imgs = registered_imgs.sel({Axes.X: xroi, Axes.Y: yroi})
    masking_radius = 6
    filt = Filter.WhiteTophat(masking_radius, is_volume=False)
    filtered = filt.run(registered_imgs, verbose=False, in_place=False)
    bleed = np.array(
        [[ 1.  ,  0.  , -0.656,  0.  ],
         [ 0.  ,  1.  ,  0.   ,  0.  ],
         [ 0.  ,  0.  ,  1.   ,  0.  ],
         [ 0.  ,  0.  ,  0.   ,  1.  ]]
    )
    lum = starfish.image.Filter.LinearUnmixing(bleed)
    bleed_corrected = lum.run(filtered, in_place=False)
    matched_histograms = starfish.image.Filter.MatchHistograms({Axes.CH})
    scaled_c = matched_histograms.run(bleed_corrected, in_place=False, verbose=False)
    cptz= Filter.ClipPercentileToZero(p_min=80, p_max=99.99995, level_method=Levels.SCALE_BY_CHUNK, group_by={Axes.CH})
    clipped_both_scaled = cptz.run(scaled_c, in_place=False)
    dots = clipped_both_scaled.reduce({Axes.CH}, func="max").reduce({Axes.ROUND}, func="mean")
    bd = FindSpots.BlobDetector(
        min_sigma=1,
        max_sigma=10,
        num_sigma=20,
        threshold=.07,
        is_volume=True,
        measurement_type='max',
    )
    spots = bd.run(clipped_both_scaled, reference_image=dots)
    single_codebook = Codebook.open_json('code/data_formatting/codebooks/codebook_A_PB2.json')
    multi_barcode_decoder = DecodeSpots.MultiBarcodeDecoder(
        codebook=single_codebook,
        max_distance=.5,
        min_intensity=.1,
        return_original_intensities=True,
        raw_intensity_threshold=0.07,
    )
    mbd_decoded_spots = multi_barcode_decoder.run(spots=spots)
    if save:
       prefix = f"{outdir}/fov_{fov_id}"
       feature_df = mbd_decoded_spots.to_features_dataframe()
       feature_df.to_csv(f"{prefix}_spots.csv", na_rep=None, index=False)
       mbd_decoded_spots.to_netcdf(f"{prefix}_spots.nc")
       clipped_both_scaled.xarray.to_netcdf(f"{prefix}_primary.nc")
       transforms_list.to_json(f"{prefix}_registration.json")
    return mbd_decoded_spots

expected_pattern = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0],[0,0,0,1]])
shifted = expected_pattern.copy()
shifted[1:] += shifted[:-1]

def num_correct_by_round(spots, round):
    return np.all(spots["binarized"][:,:(round+1),:] == expected_pattern[:(round+1),:], axis=(1,2)).sum().item()

def num_incorrect_by_round(spots, round):
	ep = expected_pattern[:(round+1),:]
	spb = spots["binarized"][:,:(round+1),:]
	ispb = spb[~np.all(spb == ep, axis=(1,2))]
	additional_only = np.all(np.logical_and(ispb, ep) == ep, axis=(1,2)).sum().item()
	missing_only = np.all(np.logical_or(ispb, ep) == ep, axis=(1,2)).sum().item()
	return additional_only, missing_only, len(ispb)-additional_only-missing_only

def get_metrics_from_spots(spots):
    metrics = {}
    metrics["total"] = len(spots)
    for i in range(6):
       metrics[f"correct_to_{i+1}"] = num_correct_by_round(spots, i)
       metrics[f"additional_to_{i+1}"],metrics[f"missing_to_{i+1}"],metrics[f"additionalAndMissing_to_{i+1}"] = num_incorrect_by_round(spots, i)
    metrics["invalid_consistent"] = len(spots[spots.target == "PB2(invalid)"])
    metrics["invalid_inconsistent"] = len(spots[spots.target == "(invalid)"])
    invalid_spots = spots[[("invalid" in t) for t in spots.target.values]]
    metrics["invalid_shifted"] = np.all(np.logical_or(invalid_spots["binarized"], shifted) == shifted, axis=(1,2)).sum().item()
    missing = len(spots[spots.target == "missing"])
    metrics["missing_consistent"] = np.all(np.logical_or(spots[spots.target=="missing"]["binarized"], expected_pattern) == expected_pattern, axis=(1,2)).values.sum()
    metrics["missing_inconsistent"] = missing - metrics["missing_consistent"]
    return metrics

results = []
for fov_id in range(len(exp.fovs())):
	spots = get_decoded_spots_for_fov(fov_id, save=True)
	metrics = get_metrics_from_spots(spots)
	metrics["fov"] = fov_id
	results.append(metrics)
results_df = pd.DataFrame(results)
results_df.to_csv(f"{outdir}/results.csv", index=False)