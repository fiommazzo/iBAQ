## relative iBAQ: normalize iBAQ values across differente experiments

# We determined the relative iBAQ (riBAQ), which is the iBAQ for 
# a protein or protein group (calculated by MaxQuant) divided by all non-contaminant, 
# non-reversed iBAQ values for a replicate. 
# riBAQ is equivalent to normalized molar intensity30 (i), and is an example of a 
# transformation we call “fractional normalization.” 

library(tidyverse)


read_pg_ibaq <- function(folder){
  pg <- read_delim(folder,
                      comment = "#",
                      show_col_types = F) |>
    filter(is.na(`Potential contaminant`),
           is.na(Reverse),
           is.na(`Only identified by site`)) |>
    dplyr::select(matches("iBAQ"), `Gene names`, `Majority protein IDs`, Peptides) |>
    mutate(Protein = vapply(strsplit(`Majority protein IDs`,";"), `[`, 1, FUN.VALUE=character(1)),
           Gene = vapply(strsplit(`Gene names`,";"), `[`, 1, FUN.VALUE=character(1)),
           Identifier = paste0(Gene, "_", Protein)) |>
    dplyr::select(-Gene, -`Gene names`, -`Majority protein IDs`, -`Gene names`, -Protein, Peptides) |>
    column_to_rownames("Identifier") |>
    mutate(across(everything(), ~replace(.x, which(is.na(.x) ), 0)))
  
  return(pg)
}

pg_p1 <- read_pg_ibaq("proteinGroups_01.txt")
pg_p2 <- read_pg_ibaq("proteinGroups_02.txt") |>
  dplyr::select(Peptides, `iBAQ peptides`, `iBAQ`) |>
  rownames_to_column("ID") |>
  filter(!grepl("REV", ID),
         `iBAQ peptides` > 1) |>
  column_to_rownames("ID")


#calculate normalization factor
pg_p1$riBAQ <- (pg_p1$iBAQ / sum(pg_p1$iBAQ) *100)
pg_p2$riBAQ <- (pg_p2$iBAQ / sum(pg_p2$iBAQ) *100)

#bait normalization
pg_p1$bait_norm <- pg_p1$riBAQ / pg_p1$riBAQ[grepl("Pcgf1", rownames(pg_p1))]
pg_p2$bait_norm <- pg_p2$riBAQ / pg_p2$riBAQ[grepl("Pcgf2", rownames(pg_p2))]


#merge tables
tbl_CBXs <- pg_p1 |>
  rownames_to_column("ID") |>
  mutate(Gene = vapply(strsplit(ID,"_"), `[`, 1, FUN.VALUE=character(1))) |>
  filter(Gene %in% CBXs) |>
  dplyr::select(Gene, Peptides) |> full_join(
    pg_p2 |>
      rownames_to_column("ID") |>
      mutate(Gene = vapply(strsplit(ID,"_"), `[`, 1, FUN.VALUE=character(1))) |>
      filter(Gene %in% CBXs) |>
      dplyr::select(Gene, Peptides), by = "Gene")
