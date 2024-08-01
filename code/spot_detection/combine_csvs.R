library(tidyverse)

spots_from_path <- function(path, components) {
  spot_files <- dir(path, "*_spots.csv", full.names = TRUE, recursive = TRUE)
  names(spot_files) <- spot_files
  spots <- spot_files %>% map_df(read_csv, col_types = cols(), na = "", .id = "file")
  spots <- spots %>%
    mutate(
      file = str_remove(file, str_c(path, "/")),
      file = str_remove(file, "_spots.csv"),
    ) %>%
    separate(file, components, sep = "/") %>%
    group_by(across(all_of(components))) %>%
    mutate(n_cells = n_distinct(cell)) %>%
    ungroup()
  return(spots)
}

write_with_status <- function(df, path) {
  outfile <- str_c(path, "all_spots.tsv.xz")
  write_tsv(df, outfile)
  message(str_c("Wrote ", nrow(df), " spots to file: ", outfile))
}

# cDNA_vRNA
path <- "analysis/spot_detection/cDNA_vRNA/"
spots <- spots_from_path(path, c("rep", "molecule", "sample", "fov"))
write_with_status(spots, path)

# plp_individual
path <- "analysis/spot_detection/plp_individual/"
spots <- spots_from_path(path, c("rep", "segment", "sample", "fov"))
write_with_status(spots, path)

# plp_cumulative
path <- "analysis/spot_detection/plp_cumulative/"
spots <- spots_from_path(path, c("rep", "segment", "sample", "fov"))
write_with_status(spots, path)

# specificity
path <- "analysis/spot_detection/specificity/"
spots <- spots_from_path(path, c("rep", "sample", "fov"))
write_with_status(spots, path)

# seq_2nt
path <- "analysis/spot_detection/seq_2nt/"
spots <- spots_from_path(path, c("rep", "moi", "hpi", "fov")) %>%
  mutate(across(c(rep, moi, hpi, fov), parse_number))
write_with_status(spots, path)
