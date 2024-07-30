"""
The 2nt data without ICC is to be used for training and applying the cellpose model.
For this purpose intensity in the dapi channel needs to be rescaled to amplify autofluorescence.
rep0 and the other reps need different scaling factors because of difference in microscope software versions.
"""
import os
from tqdm import tqdm
from skimage import io
import skimage as ski
from helpers.seq_2nt import get_fovs

for rep in tqdm(range(3)):
    dir = "2022.05.09" if rep == 0 else "2023.10.19"
    dir += "_PR8_AllSegments_vRNA_mRNA_8hTimepoints_0.3_1_MOI_noICC"
    intensity_range = (1500, 30000) if rep == 0 else (100, 300)
    for moi in tqdm([0.3, 1], leave=False):
        for hpi in tqdm(range(9), leave=False):
            outputdir = f"data/intensity_scaled/seq_2nt/rep{rep}/{moi:.1f}MOI/{hpi}hpi"
            os.makedirs(outputdir, exist_ok=True)
            fulldir = f"data/raw/{dir}/r{rep}_{moi}MOI"

            dapi_fovs = get_fovs(fulldir, rep, moi, hpi, 0)
            for i,fov in enumerate(dapi_fovs):
                image = io.imread(fov)
                colorscaled = ski.exposure.rescale_intensity(image, in_range=intensity_range)
                new_file_name = f'nuclei-fov_00{i}-c0-r0-z0.tiff'
                output_file_path = os.path.join(outputdir, new_file_name)
                io.imsave(output_file_path, colorscaled)