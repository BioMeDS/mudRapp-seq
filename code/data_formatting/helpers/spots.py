import os
import re
import tifffile
import numpy as np
import xml.etree.ElementTree as ET
from glob import glob
from typing import Mapping, Union
from starfish.experiment.builder import FetchedTile, TileFetcher
from starfish.types import Axes, Coordinates

def metadata_from_xml(fname):
    tree = ET.parse(fname)
    root = tree.getroot()
    return {
        "locations": [t.attrib for t in root.find(".//Attachment[@Name='TileScanInfo']")],
        "width": float(root.find(".//DimensionDescription[@DimID='X']").attrib['Length'])/1e6, # assumption: unit=μm, should be checked via .attrib['Unit']
        "height": float(root.find(".//DimensionDescription[@DimID='Y']").attrib['Length'])/1e6, # assumption: unit=μm, should be checked via .attrib['Unit']
    }


def get_position_from_metadata(metadata, tile, axis):
	"""Parse the imagej metadata to get X, Y or Z location for the tile

	Parameters
	----------
	metadata : list<dict>
		Parsed from metadata xml file via `metadata_from_xml`
	tile : int
		The tile number
	axis : str
		The axis to get the position for, either 'X' or 'Y'
	"""
	return float(metadata['locations'][tile][f"Pos{axis}"])


def metadata_filename(image_path):
    elements = image_path.split('/')
    elements.insert(-1, 'MetaData')
    return re.sub(r'_s\d\d_ch\d\d.tif', '_Properties.xml', '/'.join(elements))

# subclass FetchedTile
class ISSITile(FetchedTile):

    def __init__(
            self,
            file_path: str,
            tile_id: int,
            metadata,
            extras,
    ) -> None:
        """Parser for a tile.

        Parameters
        ----------
        file_path : str
            location of the tiff
        tile_id : int
            id of the current tile (to find position in metadata)
        """
        self.file_path = file_path
        self.tile_id = tile_id
        self.metadata = metadata
        self.tfile = tifffile.TiffFile(self.file_path)
        self._extras = extras

        # coordinates must match shape
        # get coordinates from metadata
        x = get_position_from_metadata(self.metadata, self.tile_id, 'X')
        y = get_position_from_metadata(self.metadata, self.tile_id, 'Y')
        z = get_position_from_metadata(self.metadata, self.tile_id, 'Z')
        width = metadata['width']
        height = metadata['height']
        self._coordinates = {
            Coordinates.X: (x, x+width),
            Coordinates.Y: (y, y+height),
            Coordinates.Z: (z, z)
        }

    @property
    def shape(self) -> Mapping[Axes, int]:
        image_shape = self.tfile.asarray().shape
        return {Axes.Y: image_shape[0], Axes.X: image_shape[1]}  # assuming shape to be (row, col) not (x,y)

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def extras(self):
        return self._extras

    def tile_data(self) -> np.ndarray:
        return self.tfile.asarray()


class PrimaryTileFetcher(TileFetcher):

    def __init__(self, input_dir: str) -> None:
        self.input_dir = os.path.join(input_dir)
        self.fovs = sorted(glob(os.path.join(self.input_dir, f"*_ICC_Processed001_s*_ch02.tif")))
    
    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        filename = self.fovs[fov_id]
        match = re.search(r"/(r\d)/.*Region (\d+).*_s(\d+)_ch...tif", filename)
        tile = match.group(3)
        tile_id = int(tile)
        filename = self.fovs[fov_id]
        metadata = metadata_from_xml(metadata_filename(filename))
        extras = {"original_path": filename}
        return ISSITile(filename, tile_id, metadata, extras)


class AuxTileFetcher(TileFetcher):

    def __init__(self, input_dir: str, channel: int) -> None:
        self.input_dir = os.path.join(input_dir)
        self.channel = channel
        self.fovs = sorted(glob(os.path.join(self.input_dir, f"*_ICC_Processed001_s*_ch0{channel}.tif")))

    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        filename = self.fovs[fov_id].replace("1Inc",f"{round_label+1}Inc")
        match = re.search(r"/(r\d)/.*Region (\d+).*_s(\d+)_ch...tif", filename)
        tile = match.group(3)
        tile_id = int(tile)
        filename = self.fovs[fov_id]
        metadata = metadata_from_xml(metadata_filename(filename))
        extras = {"original_path": filename}
        return ISSITile(filename, tile_id, metadata, extras)




primary_image_dimensions: Mapping[Union[str, Axes], int] = {
    Axes.ROUND: 1,
    Axes.CH: 1,
    Axes.ZPLANE: 1,
}


aux_images_dimensions: Mapping[str, Mapping[Union[str, Axes], int]] = {
    "nuclei": {
        Axes.ROUND: 1,
        Axes.CH: 1,
        Axes.ZPLANE: 1,
    },
    "brightfield": {
        Axes.ROUND: 1,
        Axes.CH: 1,
        Axes.ZPLANE: 1,
    },
}
