import os
import shutil
from glob import glob
from tqdm import tqdm
from slicedimage import ImageFormat
from starfish.experiment.builder import write_experiment_json
from starfish.types import Axes
from helpers.spots import PrimaryTileFetcher, AuxTileFetcher, get_primary_image_dimensions, get_aux_images_dimensions

for dir in tqdm([f"2023_09_15_ATGC_test/{x}" for x in ["A_PB2", "C_PB1", "G_PA", "T_HA"]]):
    _, sample = dir.split("/")
    outputdir = f"data/spacetx/seq_qc/rep0/{sample}"
    os.makedirs(outputdir, exist_ok=True)
    fulldir = f"data/raw/{dir}/1Inc"

    rounds = (6 if sample == "A_PB2" else 2)
    primary_tile_fetcher = PrimaryTileFetcher(fulldir, channel_offset=1)
    nuclei_tile_fetcher = AuxTileFetcher(fulldir, 0)
    bf_tile_fetcher = AuxTileFetcher(fulldir, 5)
    fovs = len(glob(f"{fulldir}/*_ICC_Processed001_s*_ch02.tif"))
    write_experiment_json(
        path=outputdir,
        fov_count=fovs,
        tile_format=ImageFormat.TIFF,
        primary_image_dimensions=get_primary_image_dimensions(rounds),
        aux_name_to_dimensions=get_aux_images_dimensions(rounds),
        primary_tile_fetcher=primary_tile_fetcher,
        aux_tile_fetcher={"nuclei": nuclei_tile_fetcher, "brightfield": bf_tile_fetcher},
        dimension_order=(Axes.ROUND, Axes.CH, Axes.ZPLANE)
    )
    shutil.copy(f"code/data_formatting/codebooks/codebook_{sample}.json", outputdir+"/codebook.json")
