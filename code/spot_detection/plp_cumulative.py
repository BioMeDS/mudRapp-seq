import os
from starfish import Experiment
from helpers.single_channel_spots import spot_detection_segmentation

for rep in ["rep0", "rep1", "rep2", "rep3", "rep4"]:
	for segment in ["NA", "HA", "PB1"]:
		for plp in [str(x)+"PLP" for x in [1,2,3,4,5,6,8,10,"NegCtrl"]]:
			exp_file = os.path.join(f"data/spacetx/plp_cumulative/{rep}/{segment}/{plp}", "experiment.json")
			if not os.path.isfile(exp_file):
				print(f"Skipping {exp_file}")
				continue
			exp = Experiment.from_json(exp_file)
			dir = f"analysis/spot_detection/plp_cumulative/{rep}/{segment}/{plp}"
			spot_detection_segmentation(exp, dir)
