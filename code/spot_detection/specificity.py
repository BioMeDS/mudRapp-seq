import os
from starfish import Experiment
from helpers.specificity_spots import spot_detection_segmentation

for sample in ["PR8", "StPt", "PR8_StPt", "NegCtrl"]:
	exp_file = os.path.join(f"data/spacetx/specificity/rep0/{sample}", "experiment.json")
	exp = Experiment.from_json(exp_file)
	dir = f"analysis/spot_detection/specificity/rep0/{sample}"
	spot_detection_segmentation(exp, dir)