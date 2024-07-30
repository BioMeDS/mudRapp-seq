import os
from starfish import Experiment
from helpers.single_channel_spots import spot_detection_segmentation

for rep in ["rep1", "rep2"]:
	for molecule in ["cDNA", "cDNA_RCAprimer", "vRNA"]:
		for plp in [str(x)+"PLP" for x in [1,5,10,"NegCtrl"]]:
			exp_file = os.path.join(f"data/spacetx/cDNA_vRNA/{rep}/{molecule}/{plp}", "experiment.json")
			if not os.path.isfile(exp_file):
				print(f"Skipping {exp_file}")
				continue
			exp = Experiment.from_json(exp_file)
			dir = f"analysis/spot_detection/cDNA_vRNA/{rep}/{molecule}/{plp}"
			spot_detection_segmentation(exp, dir)
