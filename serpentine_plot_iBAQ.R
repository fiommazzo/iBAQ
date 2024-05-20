library(tidyverse)

#first upload a protein group file

# known complex to label in the plot
PRC11 <- c("Pcgf1", "RNF2", "RING1", "RYBP", "YAF2", "KDM2B", "BCOR", "BCOL1", "USP7", "SKP1")
Other <- c("HIST4", "HISTH2A")

serpentine_plot <- function(proteinGroup, Complex, additional_category){
  proteinGroups |>
    filter(Potential.contaminant != "+",
           Reverse != "+",
           Only.identified.by.site != "+",
           Unique.peptides > 0,
           Peptides > 1) |>
    dplyr::select(Majority.protein.IDs, Gene.names, iBAQ) |>
    mutate(Majority.protein.IDs = vapply(strsplit(Majority.protein.IDs, ";"), "[", 1, FUN.VALUE = character(1)),
           Gene.names = vapply(strsplit(Gene.names, ";"), "[", 1, FUN.VALUE = character(1)),
           log10_iBAQ = log10(iBAQ)) |>
    arrange(-log10_iBAQ) |>
    mutate(complex = case_when(Gene.names %in% Complex ~ "Known Complex",
                               Gene.names %in% additional_category ~ "Other",
                               TRUE ~ "-"),
           Label = if_else(complex != "-", Gene.names, NA),
           rank = 1:n()) |> 
    filter(log10_iBAQ != "-Inf") |>
    ggplot(aes(x = rank, y = log10_iBAQ, col = complex)) +
    geom_point(alpha = 0.5) +
    geom_text(aes(label = Label), hjust = -0.7) +
    scale_color_manual(values = c(`-`="black",`Other`="orange", `Known Complex` = "purple")) +
    theme_classic()
}

serpentine_plot(proteinGroups, PRC11, Other)


