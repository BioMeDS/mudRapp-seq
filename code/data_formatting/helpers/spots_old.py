import os
import re
import tifffile
import numpy as np
import xml.etree.ElementTree as ET
from glob import glob
from typing import Mapping
from starfish.experiment.builder import FetchedTile, TileFetcher
from starfish.types import Axes, Coordinates

def metadata_from_xml_old(fname):
    tree = ET.parse(fname)
    root = tree.getroot()
    return {
        "x": float(root.find(".//ATLCameraSettingDefinition").attrib['StagePosX'].split(" ")[0])/1000,
        "y": float(root.find(".//ATLCameraSettingDefinition").attrib['StagePosY'].split(" ")[0])/1000,
        "z": float(root.find(".//AdditionalZPosition").attrib['ZPosition']),
        "width": float(root.find(".//DimensionDescription[@DimID='X']").attrib['Length'])/1e6, # assumption: unit=μm, should be checked via .attrib['Unit']
        "height": float(root.find(".//DimensionDescription[@DimID='Y']").attrib['Length'])/1e6, # assumption: unit=μm, should be checked via .attrib['Unit']
    }


def metadata_filename_old(image_path):
    elements = image_path.split('/')
    elements.insert(-1, 'MetaData')
    return re.sub(r'_ch\d\d.tif', '_Properties.xml', '/'.join(elements))

# subclass FetchedTile
class ISSITileOld(FetchedTile):

    def __init__(
            self,
            file_path: str,
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
        self.metadata = metadata
        self.tfile = tifffile.TiffFile(self.file_path)
        self._extras = extras

        # coordinates must match shape
        # get coordinates from metadata
        x = metadata['x']
        y = metadata['y']
        z = metadata['z']
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


class PrimaryTileFetcherOld(TileFetcher):

    def __init__(self, input_dir: str) -> None:
        self.input_dir = os.path.join(input_dir)
        self.fovs = sorted(glob(os.path.join(self.input_dir, f"*_ch01.tif")))
    
    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        filename = self.fovs[fov_id]
        metadata = metadata_from_xml_old(metadata_filename_old(filename))
        extras = {"original_path": filename}
        return ISSITileOld(filename, metadata, extras)


class AuxTileFetcherOld(TileFetcher):

    def __init__(self, input_dir: str, channel: int) -> None:
        self.input_dir = os.path.join(input_dir)
        self.channel = channel
        self.fovs = sorted(glob(os.path.join(self.input_dir, f"*_ch0{channel}.tif")))

    def get_tile(
            self, fov_id: int, round_label: int, ch_label: int, zplane_label: int) -> FetchedTile:
        filename = self.fovs[fov_id].replace("1Inc",f"{round_label+1}Inc")
        metadata = metadata_from_xml_old(metadata_filename_old(filename))
        extras = {"original_path": filename}
        return ISSITileOld(filename, metadata, extras)
