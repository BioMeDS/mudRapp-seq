import os
import shutil
from glob import glob
from tqdm import tqdm
from slicedimage import ImageFormat
from starfish.experiment.builder import write_experiment_json
from starfish.types import Axes
from helpers.spots import PrimaryTileFetcher, AuxTileFetcher, primary_image_dimensions, aux_images_dimensions
from helpers.spots_old import PrimaryTileFetcherOld, AuxTileFetcherOld

for dir in [
    "2022.02.23_PR8_1MOI_HA_NA_PB1_DiffPLP_newFormat/HA",
    "2022.02.23_PR8_1MOI_HA_NA_PB1_DiffPLP_newFormat/PB1",
    "2022.02.23_PR8_1MOI_HA_NA_PB1_DiffPLP_newFormat/NA",
    ]:
    for i in tqdm([str(x)+" PLP" for x in [1,2,3,4,5,6,8,10]] + ["NegCtrl 10 PLP"]):
        plp = i.split(" ")[0] + "PLP"
        segment = dir.split("/")[1]
        if i == "5 PLP" and segment in ["HA", "NA"]:
            continue
        rep = "rep0"
        outputdir = f"data/spacetx/plp_cumulative/{rep}/{segment}/{plp}"
        os.makedirs(outputdir, exist_ok=True)
        fulldir = f"data/raw/{dir}/{i}"

        primary_tile_fetcher = PrimaryTileFetcherOld(fulldir)
        nuclei_tile_fetcher = AuxTileFetcherOld(fulldir, 0)
        bf_tile_fetcher = AuxTileFetcherOld(fulldir, 2)
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

for dir in [
    "2023.10.19_increasingPLPs_PR8_6hpi_1MOI/PB1/r1",
    "2023.10.19_increasingPLPs_PR8_6hpi_1MOI/PB1/r2",
    "2023.10.19_increasingPLPs_PR8_6hpi_1MOI/HA/r1",
    "2023.10.19_increasingPLPs_PR8_6hpi_1MOI/HA/r2",
    "2023.11.10_increasingPLPs_PR8_6hpi_1MOI/NA/r1",
    "2023.11.10_increasingPLPs_PR8_6hpi_1MOI/NA/r2",
    "2023.11.16_increasingPLPs_PR8_6hpi_1MOI/HA/r3",
    "2023.11.16_increasingPLPs_PR8_6hpi_1MOI/HA/r4",
    ]:
    for i in tqdm([str(x)+" PLP" for x in [1,2,3,4,5,6,8,10]] + ["NegCtrl 10PLP"]):
        plp = i.split(" ")[0] + "PLP"
        _, segment, rep = dir.split("/")
        fovs = 10
        outputdir = f"data/spacetx/plp_cumulative/{rep.replace('r','rep')}/{segment}/{plp}"
        os.makedirs(outputdir, exist_ok=True)
        fulldir = f"data/raw/{dir}/{i}"

        primary_tile_fetcher = PrimaryTileFetcher(fulldir)
        nuclei_tile_fetcher = AuxTileFetcher(fulldir, 1)
        bf_tile_fetcher = AuxTileFetcher(fulldir, 0)
        write_experiment_json(
            path=outputdir,
            fov_count=fovs,
            tile_format=ImageFormat.TIFF,
            primary_image_dimensions=primary_image_dimensions,
            aux_name_to_dimensions=aux_images_dimensions,
            primary_tile_fetcher=primary_tile_fetcher,
            aux_tile_fetcher={"nuclei": nuclei_tile_fetcher, "brightfield": bf_tile_fetcher},
            dimension_order=(Axes.ROUND, Axes.CH, Axes.ZPLANE)
        )
        shutil.copy("code/data_formatting/codebooks/codebook_spot.json", outputdir+"/codebook.json")
