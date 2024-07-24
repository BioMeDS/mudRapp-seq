import os
import shutil
from glob import glob
from tqdm import tqdm
from slicedimage import ImageFormat
from starfish.experiment.builder import write_experiment_json
from starfish.types import Axes
from helpers.spots import PrimaryTileFetcher, AuxTileFetcher, primary_image_dimensions, aux_images_dimensions

for dir in [
    "2023.11.10_PR8_HA_1MOI_6hpi_PLPs1_10_separate/r1",
    "2023.11.10_PR8_HA_1MOI_6hpi_PLPs1_10_separate/r2",
    "2023.11.10_PR8_NA_1MOI_6hpi_PLPs1_10_separate/r1",
    "2023.11.10_PR8_NA_1MOI_6hpi_PLPs1_10_separate/r2",
    ]:
    for i in tqdm([c+"PLP "+str(x) for x in range(1,11) for c in ["", "Neg "]]):
        plp = i.split(" ")[-1] + "PLP"
        fovs = 10
        if i.startswith("Neg"):
            plp = "Neg_" + plp
            fovs = 6
        segment, rep = dir.split("/")
        segment = segment.split("_")[2]
        outputdir = f"data/spacetx/plp_individual/{rep.replace('r','rep')}/{segment}/{plp}"
        os.makedirs(outputdir, exist_ok=True)
        fulldir = f"data/raw/{dir}/{i}"

        primary_tile_fetcher = PrimaryTileFetcher(fulldir)
        nuclei_tile_fetcher = AuxTileFetcher(fulldir, 1)
        bf_tile_fetcher = AuxTileFetcher(fulldir, 0)
        write_experiment_json(
            path=outputdir,
            fov_count=len(glob(f"{fulldir}/*_ch02.tif")),
            tile_format=ImageFormat.TIFF,
            primary_image_dimensions=primary_image_dimensions,
            aux_name_to_dimensions=aux_images_dimensions,
            primary_tile_fetcher=primary_tile_fetcher,
            aux_tile_fetcher={"nuclei": nuclei_tile_fetcher, "brightfield": bf_tile_fetcher},
            dimension_order=(Axes.ROUND, Axes.CH, Axes.ZPLANE)
        )
        shutil.copy("code/data_formatting/codebooks/codebook_spot.json", outputdir+"/codebook.json")
