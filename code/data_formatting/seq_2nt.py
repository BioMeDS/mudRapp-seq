import os
import shutil
from tqdm import tqdm
from slicedimage import ImageFormat
from starfish.experiment.builder import write_experiment_json
from starfish.types import Axes
from helpers.seq_2nt import PrimaryTileFetcher, AuxTileFetcher, get_fovs, primary_image_dimensions, aux_images_dimensions

for rep in tqdm(range(3)):
    for moi in tqdm([0.3, 1], leave=False):
        for hpi in tqdm(range(9), leave=False):
            dir = "2022.05.09" if rep == 0 else "2023.10.19"
            dir += "_PR8_AllSegments_vRNA_mRNA_8hTimepoints_0.3_1_MOI"
            outputdir = f"data/spacetx/seq_2nt/rep{rep}/{moi:.1f}MOI/{hpi}hpi"
            os.makedirs(outputdir, exist_ok=True)
            fulldir = f"data/raw/{dir}/r{rep}_{moi}MOI"

            primary_fovs = get_fovs(fulldir, rep, moi, hpi, 1)
            primary_tile_fetcher = PrimaryTileFetcher(primary_fovs)
            nuclei_tile_fetcher = AuxTileFetcher(get_fovs(fulldir, rep, moi, hpi, 0))
            bf_tile_fetcher = AuxTileFetcher(get_fovs(fulldir, rep, moi, hpi, 5))
            fovs = len(primary_fovs)
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
            shutil.copy(f"code/data_formatting/codebooks/codebook_2nt.json", outputdir+"/codebook.json")
