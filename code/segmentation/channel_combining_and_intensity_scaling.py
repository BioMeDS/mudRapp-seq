import numpy as np 
import skimage as ski
import matplotlib.pyplot as plt
import numpy as np
import os 
from skimage import io, exposure, restoration, filters, morphology

def save_three_channel_image(layer0_path, layer1_path, layer2_path, full_path):
    layer0_data = ski.io.imread(layer0_path)
    layer1_data = ski.io.imread(layer1_path)
    layer2_data = ski.io.imread(layer2_path)
    contrast_limit = 0.05
    adjusted_contrast_0 = exposure.equalize_adapthist(layer0_data, clip_limit=contrast_limit)
    adjusted_contrast_1 = exposure.equalize_adapthist(layer1_data, clip_limit=contrast_limit)
    adjusted_contrast_2 = exposure.equalize_adapthist(layer2_data, clip_limit=contrast_limit)
    img = np.stack([adjusted_contrast_0, adjusted_contrast_1, adjusted_contrast_2], axis=-1)
    ski.io.imsave(full_path, img)

#specify inputs  
layer0_path = (f"path/to/your/image.tif")
layer1_path = (f"path/to/your/image.tif")
layer2_path = (f"path/to/your/image.tif")

filename = f"your_desired_outputname"
directory_path = f"dummy/for/your/save/path"

full_path = f"{directory_path}/{filename}.tiff"

save_three_channel_image(layer0_path, layer1_path, layer2_path, full_path)




## With actual paths for this project:


layer0_path = (f"/storage/biomeds/data/in-situ-seq-influenza/raw/2023.09.15_PR8_Nepal_1MOI4-6hTimepoint_3R_50C_AllSeg_BSPEG_repeat1_2_raw/{Inc} Inc/r{round}/{hpi} hpi/{Inc}Inc_PR8_Nepal_{Moi}MOI_{hpi}hpi_AllSegments_3R_50C__Region {fov}_Processed001_s0{s}_ch00.tif")
layer1_path = (f"/storage/biomeds/data/in-situ-seq-influenza/raw/2023.09.15_PR8_Nepal_1MOI4-6hTimepoint_3R_50C_AllSeg_BSPEG_repeat1_2_raw/{Inc} Inc/r{round}/{hpi} hpi/{Inc}Inc_PR8_Nepal_{Moi}MOI_{hpi}hpi_AllSegments_3R_50C__Region {fov}_Processed001_s0{s}_ch02.tif")
layer2_path = (f"/storage/biomeds/data/in-situ-seq-influenza/raw/2023.09.15_PR8_Nepal_1MOI4-6hTimepoint_3R_50C_AllSeg_BSPEG_repeat1_2_raw/{Inc} Inc/r{round}/{hpi} hpi/{Inc}Inc_PR8_Nepal_{Moi}MOI_{hpi}hpi_AllSegments_3R_50C__Region {fov}_Processed001_s0{s}_ch04.tif")


filename = f"3nt3chan_rep{rep}_{Moi}MOI_{hpi}hpi_fov{fov}_{Inc}Inc_round0{round}_ch0-2-4_s{s}"
directory_path = f"dummy/for/your/save/path"
full_path = f"{directory_path}/{filename}.tiff"
save_three_channel_image(layer0_path, layer1_path, layer2_path, full_path)
