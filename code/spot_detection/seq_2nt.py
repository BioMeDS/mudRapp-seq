from glob import glob
from tqdm import tqdm
import re
from helpers.seq_2nt_spots import analyze_fov

files = sorted(glob("data/spacetx/seq_2nt/*/*/*/nuclei-fov_*-c0-r0-z0.tiff"))
for file in tqdm(files):
	match = re.search(r"data/spacetx/seq_2nt/rep(\d)/(...)MOI/(\d)hpi/nuclei-fov_00(\d)-c0-r0-z0.tiff", file)
	rep = int(match.group(1))
	moi = match.group(2)
	hpi = int(match.group(3))
	fov_id = int(match.group(4))
	# print(f"analyze_fov({rep}, {moi}, {hpi}, {fov_id})")
	analyze_fov(rep, moi, hpi, fov_id)
