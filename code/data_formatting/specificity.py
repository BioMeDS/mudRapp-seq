import os
import shutil
from glob import glob
from tqdm import tqdm
from slicedimage import ImageFormat
from starfish.experiment.builder import write_experiment_json
from starfish.types import Axes
from helpers.spots import PrimaryTileFetcher, AuxTileFetcher, primary_image_dimensions_4channel, aux_images_dimensions

for dir in tqdm([
    "2023.11.10_PR8_StPt_1MOI_6hpi_HA_Specificity_experiment/NegCtrl",
    "2023.11.10_PR8_StPt_1MOI_6hpi_HA_Specificity_experiment/PR8",
    "2023.11.10_PR8_StPt_1MOI_6hpi_HA_Specificity_experiment/PR8_StPt",
    "2023.11.10_PR8_StPt_1MOI_6hpi_HA_Specificity_experiment/StPt",
    ]):
    _, strain = dir.split("/")
    outputdir = f"data/spacetx/specificity/rep0/{strain}"
    os.makedirs(outputdir, exist_ok=True)
    fulldir = f"data/raw/{dir}"

    primary_tile_fetcher = PrimaryTileFetcher(fulldir, channel_offset=1)
    nuclei_tile_fetcher = AuxTileFetcher(fulldir, 0)
    bf_tile_fetcher = AuxTileFetcher(fulldir, 5)
    write_experiment_json(
        path=outputdir,
        fov_count=len(glob(f"{fulldir}/*_ch02.tif")),
        tile_format=ImageFormat.TIFF,
        primary_image_dimensions=primary_image_dimensions_4channel,
        aux_name_to_dimensions=aux_images_dimensions,
        primary_tile_fetcher=primary_tile_fetcher,
        aux_tile_fetcher={"nuclei": nuclei_tile_fetcher, "brightfield": bf_tile_fetcher},
        dimension_order=(Axes.ROUND, Axes.CH, Axes.ZPLANE)
    )
    shutil.copy("code/data_formatting/codebooks/codebook_specificity.json", outputdir+"/codebook.json")
