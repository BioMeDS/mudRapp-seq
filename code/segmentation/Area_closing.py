import os
import numpy as np
from skimage import io, morphology, measure
from skimage.morphology import remove_small_objects

input_path = "path/to/your/directory"
output_path = "path/to/your/saving/directory"
file_name = f"name_of_input_image.tif"
input_image_path = os.path.join(input_path, file_name)


segmentation_image = io.imread(input_image_path)
processed_image = morphology.remove_small_objects(segmentation_image, min_size=min_size, connectivity=2)
closed_image = morphology.closing(processed_image, selem=np.ones((5, 5)))
area_closed_image = morphology.area_closing(closed_image, area_threshold=200) #adjust area_threshold to account for bigger holes
output_image_path = os.path.join(output_path, f"RSO_{file_name}")
io.imsave(output_image_path, area_closed_image)




