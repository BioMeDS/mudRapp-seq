library(tidyverse)

spot_files <- dir("analysis/spot_detection/cDNA_vRNA/", "*_spots.csv", full.names = TRUE, recursive = TRUE)
names(spot_files) <- spot_files
spots <- spot_files %>% map_df(read_csv, col_types = cols(), .id = "file")

spots <- spots %>%
        mutate(
                file = str_remove(file, "analysis/spot_detection/cDNA_vRNA//"),
                file = str_remove(file, "_spots.csv"),
        ) %>%
        separate(file, c("rep", "molecule", "sample", "fov"), sep = "/")

outfile <- "analysis/spot_detection/cDNA_vRNA/all_spots.tsv.xz"
write_tsv(spots, outfile)
print(str_c("Successfully written ", nrow(spots), " spots to ", outfile))
