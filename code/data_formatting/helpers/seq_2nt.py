import os
import re
from glob import glob
from typing import Mapping, Union
from starfish.experiment.builder import FetchedTile, TileFetcher
from starfish.types import Axes
from .spots import metadata_from_xml, metadata_filename, ISSITile

def get_fovs(input_dir: str, rep: int, moi: float, hpi: str, channel: int) -> list:
	fovs = []
	if rep == 0:
		pattern = os.path.join(input_dir, f"1Inc/{hpi}hpi/1Inc_*MOI_{hpi}hpi_PR8_TileScan 1_ICC_Processed001_s*_ch0{channel}.tif")
		fovs = sorted(glob(pattern), reverse=True)

		if hpi == "6" and moi == 0.3:
			glob1 = os.path.join(input_dir, f"1Inc/{hpi}hpi/1Inc_*MOI_{hpi}hpi_PR8_TileScan 1_ICC_Processed001_s*_ch0{channel}.tif")
			glob2 = os.path.join(input_dir, f"1Inc/{hpi}hpi/1Inc_*MOI_{hpi}hpi_PR8_TileScan 2_ICC_Processed001_s*_ch0{channel}.tif")
			fovs = sorted(glob(glob2), reverse=True) + sorted(glob(glob1), reverse=True)[1:]
		
	else:
		fovs = sorted(glob(os.path.join(input_dir, f"1Inc/{hpi}hpi/1Inc_PR8_*MOI_{hpi}hpi_2nt_mRNA_vRNA*_ICC_Processed001_s*_ch0{channel}.tif")))
	return fovs

class PrimaryTileFetcher(TileFetcher):

    def __init__(self, fovs: list) -> None:
        self.fovs = fovs
        self.channel_offset = 1
    
    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        filename = self.fovs[fov_id]
        match = re.search(r"/.*_s(\d+)_ch...tif", filename)
        tile = match.group(1)
        tile_id = int(tile)
        filename = filename.replace("1Inc",f"{round_label+1}Inc")
        filename = filename.replace("_ch01", f"_ch0{ch_label+self.channel_offset}")
        metadata = metadata_from_xml(metadata_filename(filename))
        extras = {"original_path": filename}
        return ISSITile(filename, tile_id, metadata, extras)


class AuxTileFetcher(TileFetcher):

    def __init__(self, fovs: list) -> None:
        self.fovs = fovs

    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        filename = self.fovs[fov_id]
        match = re.search(r"/.*_s(\d+)_ch...tif", filename)
        tile = match.group(1)
        tile_id = int(tile)
        filename = filename.replace("1Inc",f"{round_label+1}Inc")
        metadata = metadata_from_xml(metadata_filename(filename))
        extras = {"original_path": filename}
        return ISSITile(filename, tile_id, metadata, extras)




primary_image_dimensions: Mapping[Union[str, Axes], int] = {
    Axes.ROUND: 2,
    Axes.CH: 4,
    Axes.ZPLANE: 1,
}

aux_images_dimensions: Mapping[str, Mapping[Union[str, Axes], int]] = {
    "nuclei": {
        Axes.ROUND: 2,
        Axes.CH: 1,
        Axes.ZPLANE: 1,
    },
    "brightfield": {
        Axes.ROUND: 2,
        Axes.CH: 1,
        Axes.ZPLANE: 1,
    },
}
