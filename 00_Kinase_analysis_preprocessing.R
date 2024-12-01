library(readr)
library(readxl)
library(tidyverse)
library(Biostrings)

source("config.R")

# ---------------------------------------------------------------------------- #

dat_orig <- readxl::read_excel(path_to_file)

dat_preproc <- dat_orig %>%
  dplyr::filter(Change == "Decreased") %>%
  dplyr::mutate(Accession = sub("\\|.*", "", Accession))

# ---------------------------------------------------------------------------- #

dat_proteome <- Biostrings::readAAStringSet(path_to_proteome)

dat_proteome <- data.frame(
  Accession = sub("^[^|]*\\|([^|]*)\\|.*", "\\1", names(dat_proteome)),
  `Fasta name` = sub(" .*", "", names(dat_proteome)),
  check.names = FALSE,
  Sequence = as.character(dat_proteome)
)

dat_proc <- dat_preproc %>%
  dplyr::left_join(dat_proteome, by = "Accession") %>%
  dplyr::mutate(Position = as.numeric(sub("^[A-Z]", "", Phosphosite)))
# ---------------------------------------------------------------------------- #

window_size <- 33
half_window <- (window_size - 1) / 2

dat_proc <- dat_proc %>%
  dplyr::mutate(
    `Residue window` = purrr::map2_chr(
      Sequence, Position,
      ~ substr(.x, max(.y - half_window, 1), min(.y + half_window, nchar(.x)))
    )
  ) %>%
  # IMPORTANT the sites too close 
  dplyr::filter(nchar(`Residue window`) == window_size) %>%
  dplyr::mutate(
    `Window start` = Position - half_window,
    `Window end` = Position + half_window,
    `Fasta name` = paste0(">", `Fasta name`, "%", `Window start`, "%", `Window end`),
    `Fasta line` = paste0(`Fasta name`, "\n", `Residue window`)
  )

# Save the full preprocessed dataset
readr::write_csv(dat_proc, paste0("dat_preproc/", condition, "_full_table.csv"))

# ---------------------------------------------------------------------------- #

# length(dat_proc$Phosphopeptide)
# length(unique(dat_proc$Phosphopeptide))
# length(unique(dat_proc$Gene_phosphosite))

print(paste0(
  length(unique(dat_proc$Phosphopeptide)), 
  " unique phosphopeptides (peptide sequences with phosphorylation) are present in the dataset. These phosphopeptides lie further than ", half_window, " residues from the protein start or end. These include a total of ",
  length(dat_proc$Phosphopeptide),
  " phosphosites (some phosphopeptides have multiple sites and different phosphopeptides might have the same sites). These map to ",
  length(unique(dat_proc$Gene_phosphosite)), 
  " unique phosphosites."
))

# Only unique sites are left to prevent duplicate predictions
dat_proc <- dat_proc %>%
  distinct(Gene_phosphosite, .keep_all = TRUE)

# ---------------------------------------------------------------------------- #

fasta_lines_y <- dat_proc %>%
  dplyr::filter(stringr::str_starts(Phosphosite, "Y")) %>%
  dplyr::pull(`Fasta line`)

fasta_lines_st <- dat_proc %>%
  dplyr::filter(stringr::str_starts(Phosphosite, "[ST]")) %>%
  dplyr::pull(`Fasta line`)

writeLines(fasta_lines_y, paste0("dataset/", condition, "_output_Y.fasta"))
writeLines(fasta_lines_st, paste0("dataset/", condition, "_output_ST.fasta"))
