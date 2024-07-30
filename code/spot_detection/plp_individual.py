import os
from starfish import Experiment
from tqdm import tqdm
from helpers.single_channel_spots import spot_detection_segmentation

for rep in ["rep1", "rep2"]:
	for segment in ["NA", "HA"]:
		for plp in [c+str(x)+"PLP" for x in range(1,11) for c in ["", "Neg_"]]:
			exp_file = os.path.join(f"data/spacetx/plp_individual/{rep}/{segment}/{plp}", "experiment.json")
			if not os.path.isfile(exp_file):
				print(f"Skipping {exp_file}")
				continue
			exp = Experiment.from_json(exp_file)
			dir = f"analysis/spot_detection/plp_individual/{rep}/{segment}/{plp}"
			spot_detection_segmentation(exp, dir)
