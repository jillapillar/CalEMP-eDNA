# CalEMP eDNA data - for 2027 report
# all samples from F24 - S26
# author: Jill Tupitza
# date: 10/5/2026

# open libraries
library(here)
library(ggplot2)
library(vegan)
library(dplyr)
library(stringr)
library(tidyr)
library(tidyverse)
library(cowplot)
library(forcats)

# read in data
#asvs.091726 <- read.csv(here("CalEMP-091726-data", "exported", "ASV_CalEMP_091726.csv"))
#asvs.testplate <- read.csv(here("CalEMP-Fall2024-testplate", "exported", "table_mifish_EMPA_test.csv")) # these samples were combined into 070726 dataset
#asvs.070726 <- read.csv(here("CalEMP-070726-data", "exported", "table_mifish_CalEMP_070726.csv"))
#asvs.093026 <- read.csv(here("CalEMP-093026-data", "exported_12S", "ASV_CalEMP_093026.csv")) # taxonomy assigned using normal 12S db
#asvs.093026_newref <- read.csv(here("CalEMP-093026-data", "exported_newref", "table_mifish_CalEMP_093026_newref.csv")) #taxonomy assigned using eDNA explorer 12S db

# read in re-run combined data 
#asvs.combined.mitohelper <- read.csv(here("CalEMP-combined-data", "exported_mitohelper", "table_mitohelper.csv")) # taxonomy assigned with MitoHelper db
#asvs.combined.rcrux <- read.csv(here("CalEMP-combined-data", "exported_rcrux", "table_rcrux.csv")) #taxonomy assigned with 12S RCrux (Patin et al 2022)
asvs.consensus <- read.csv(here("CalEMP-combined-data", "asv_table_all.csv")) #taxonomy run through BLAST, compared, and consensus reached
  
# new inventory - excel tracking program

# sample inventory

sample.list <- names(asvs.consensus)
sample.list 



################################################
#  percent total reads assigned to some level for entire dataset ----

# collapse identical ASVs - reuse this function - for the converged Consensus Lineage
reads.per.CL <- asvs.consensus %>%
  group_by(ConsensusLineage_converge) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

## reads across whole dataset
total.reads <- sum(reads.per.CL %>% select(where(is.numeric)), na.rm = TRUE)

# assigned reads 
assigned.reads <-reads.per.CL %>%
  filter(ConsensusLineage_converge != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)

# percent unassigned
perc.assigned <- 100*(assigned.reads)/total.reads

# print stats - 96.86%
print(perc.assigned)

# perc assigned to species
assigned.to.species.conv <- comb %>% filter(!is.na(Species_conv)) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>% sum()

perc.assigned.to.sp.conv <- 100*(assigned.to.species.conv/assigned.reads)
perc.assigned.to.sp.conv # 65.49%


###################################################

#### Compare three 12S reference databases:  RCrux, MitoHelper and eDNA explorer----

#  what percentage of ASVs were assigned to species by each database

# join two datasets and collapse Unassigned ASVs into one row
comb <- asvs.consensus %>% 
  rename(ConsensusLineage.rcrux = ConsensusLineage_rcrux) %>%
  rename(ConsensusLineage.EE = ConsensusLineage_explore) %>%
  rename(ConsensusLineage.mito = ConsensusLineage_mitohelper) %>%
  rename(ConsensusLineage.conv = ConsensusLineage_converge) %>%
  separate(
    ConsensusLineage.rcrux,
    into = c("Domain_rcrux", "Phylum_rcrux", "Class_rcrux", "Order_rcrux",
             "Family_rcrux", "Genus_rcrux", "Species_rcrux"),
    sep = ";\\s*",
    remove = FALSE
  ) %>%
  mutate(
    across(
      c(Domain_rcrux, Phylum_rcrux, Class_rcrux, Order_rcrux, Family_rcrux, Genus_rcrux, Species_rcrux),
      ~ sub("^[a-z]__", "", .x)
    )
  ) %>%
  separate(
    ConsensusLineage.EE,
    into = c("Domain_EE", "Phylum_EE", "Class_EE", "Order_EE",
             "Family_EE", "Genus_EE", "Species_EE"),
    sep = ";\\s*",
    remove = FALSE
  ) %>%
  mutate(
    across(
      c(Domain_EE, Phylum_EE, Class_EE, Order_EE, Family_EE, Genus_EE, Species_EE),
      ~ sub("^[a-z]__", "", .x)
    )
  ) %>%
  separate(
    ConsensusLineage.mito,
    into = c("Domain_mito", "Phylum_mito", "Class_mito", "Order_mito",
             "Family_mito", "Genus_mito", "Species_mito"),
    sep = ";\\s*",
    remove = FALSE
  ) %>%
  mutate(
    across(
      c(Domain_mito, Phylum_mito, Class_mito, Order_mito, Family_mito, Genus_mito, Species_mito),
      ~ sub("^[a-z]__", "", .x)
    )
  ) %>%
  separate(
    ConsensusLineage.conv,
    into = c("Domain_conv", "Phylum_conv", "Class_conv", "Order_conv",
             "Family_conv", "Genus_conv", "Species_conv"),
    sep = ";\\s*",
    remove = FALSE
  ) %>%
  mutate(
    across(
      c(Domain_conv, Phylum_conv, Class_conv, Order_conv, Family_conv, Genus_conv, Species_conv),
      ~ sub("^[a-z]__", "", .x)
    )
  )



#### Compare MitoHelper, eDNA Explorer and RCrux ref db ----
#  all samples - both assigned at a 97% identity threshold #

#  what percentage of ASVs were assigned to species by each database
# assigned to species by each database
assigned.to.species.rcrux <- comb %>% filter(!is.na(Species_rcrux)) %>% 
  summarise(across(where(is.numeric), sum), .groups = "drop") %>% sum()

assigned.to.species.EE <- comb %>% filter(!is.na(Species_EE)) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>% sum()

assigned.to.species.mito <- comb %>% filter(!is.na(Species_mito)) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>% sum()


# assigned reads from each dataset with reads collapsed per consensus lineage
reads.per.CL.rcrux <- comb %>%
  group_by(ConsensusLineage.rcrux) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

reads.per.CL.EE <- comb %>%
  group_by(ConsensusLineage.EE) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

reads.per.CL.mito <- comb %>%
  group_by(ConsensusLineage.mito) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")


#total assigned reads by each dataset
assigned.reads.rcrux <-reads.per.CL.rcrux %>%
  filter(ConsensusLineage.rcrux != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)

assigned.reads.EE <-reads.per.CL.EE %>%
  filter(ConsensusLineage.EE != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)

assigned.reads.mito <-reads.per.CL.mito %>%
  filter(ConsensusLineage.mito != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)


# percent assigned to species for each database
perc.assigned.to.sp.mito <- 100*(assigned.to.species.mito/assigned.reads.mito)
perc.assigned.to.sp.mito # 56.02%

perc.assigned.to.sp.EE <- 100*(assigned.to.species.EE/assigned.reads.EE)
perc.assigned.to.sp.EE # 42.79%

perc.assigned.to.sp.rcrux <- 100*(assigned.to.species.rcrux/assigned.reads.rcrux)
perc.assigned.to.sp.rcrux # 58.42%

comb_final <- comb %>%
  select(-c(ConsensusLineage.EE, ConsensusLineage.rcrux, ConsensusLineage.mito, Domain_EE, Phylum_EE, Class_EE, Order_EE, Family_EE, Genus_EE, Species_EE,
            Domain_mito, Phylum_mito, Class_mito, Order_mito, Family_mito, Genus_mito, Species_mito, Domain_rcrux, Phylum_rcrux, Class_rcrux,
            Order_rcrux, Family_rcrux, Genus_rcrux, Species_rcrux, total_reads)) %>%
  relocate(ConsensusLineage.conv, Domain_conv, Phylum_conv, Class_conv, Order_conv, Family_conv, Genus_conv, Species_conv) %>%
  filter(Class_conv == "Actinopteri" | is.na(Class_conv)) %>%
  group_by(ConsensusLineage.conv) %>%
  summarise(
    across(
      c(Domain_conv, Phylum_conv, Class_conv, Order_conv,
        Family_conv, Genus_conv, Species_conv),
      ~ first(.x)
    ),
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# write out final table with converged ConsensusLineage only
write.csv(asvs.consensus.sep, here("processed", "asv.table.FINAL.csv"))

# visualize using a Sankey diagram at each taxonomic level: ASVs 

# what percent of total ASVs were assigned to species?
# for the ASVs assigned to different levels by different databases:
  # (a) how many agreed?
  # (b) how many were assigned lower by (i) 12S usual and (ii) eDNA explorer db?
  # (c) list of taxa coverage done better by each db

# how many made it to family level

total_reads <- comb_final %>%
  select(where(is.numeric)) %>%
  unlist() %>%
  sum(na.rm = TRUE)

assigned_reads <- comb_final %>%
  filter(!is.na(Family_conv), Family_conv != "") %>%
  select(where(is.numeric)) %>%
  unlist() %>%
  sum(na.rm = TRUE)

perc_assigned_fam <- (assigned_reads/total_reads) *100 # 95.89% of assigned reads made it to family level

# what percent of reads made it to species level
total_reads <- comb_final %>%
  select(where(is.numeric)) %>%
  unlist() %>%
  sum(na.rm = TRUE)

assigned_reads_sp <- comb_final %>%
  filter(!is.na(Species_conv), Species_conv != "") %>%
  select(where(is.numeric)) %>%
  unlist() %>%
  sum(na.rm = TRUE)

perc_assigned_spec <- (assigned_reads_sp/total_reads) *100 # 63.09% of assigned reads made it to species

# pivot long 
abundance_long <- comb_final %>% 
  pivot_longer(
    cols = where(is.numeric),
    names_to = "sample_id",
    values_to = "reads")


  
# extract metadata from sample-id
sample_metadata <- abundance_long %>%
  distinct(sample_id) %>%
  mutate(sample_id = sample_id) %>%
  separate(sample_id, into = c("region_estuary", "station", "date", "filter", "replicate"),
           sep = "_", remove = FALSE) %>%
  separate(region_estuary, into = c("region", "estuary"),
           sep = "\\.") %>%
  mutate(date = as.Date(date, format = "%y%m%d")) %>%
  mutate(season = case_when(
    month(date) >= 8 ~ paste("Fall", year(date)),
    month(date) < 8  ~ paste("Spring", year(date)),
    TRUE ~ NA_character_
  )) %>%
  mutate(estuary = case_when(region == "NC" & estuary == "BOL" ~ "BOLINAS",
                             region == "SC" & estuary == "BOL" ~ "BOLSA",
                             TRUE ~ estuary))

# attach metadata to abundance table
abundance_table <- abundance_long %>%
  left_join(sample_metadata, by = "sample_id") %>%
  group_by(sample_id) %>%
  mutate(total_reads = sum(reads, na.rm = TRUE),
         relative_reads = reads/total_reads) %>%
  ungroup()


#### Taxonomic barplot ----
# function to loop through estuaries

##  color palette

Cxcx_minimal <- c("#33FF33", "#000033", "#FF0033","#330033","#663300", "#FFCCFF",
                  "#3333FF", "#FFCC66","#9999FF", "#99CC00", "darkmagenta")
            

Cxcx <- c("#FFCC99", "#333300", "#000000", "#663300", "#FFCCFF",
          "#CC0000", "#663333", "#996600", "#660000", "#99CC00", 
          "#003366", "#000066", "#330099", "#9999FF", "grey80",
          "grey20", "#FFFFFF", "#003300", "#33FF33", "#FF0033",
          "#330033", "#FF3300", "#3333FF", "#000033", "#FFCC66")

Cxcx_ext <- c("#FFCC99","green3", "chartreuse3", "#333300", "#000000",
              "#663300", "#FFCCFF","hotpink4", "hotpink", "#33FF33",
              "deeppink", "mediumorchid", "plum", "darkmagenta", "darkolivegreen",
              "lemonchiffon","darkgreen", "peachpuff", "lavenderblush3", "limegreen",
              "khaki","#CC0000", "#663333", "#996600", "#660000", "#99CC00", 
              "#003366", "#000066", "#330099", "#9999FF", "grey80",
              "grey20", "#FF33FF", "#003300", "#FF0033","dodgerblue",
              "#330033", "#FF3300", "#3333FF", "#000033", "#FFCC66")

# assign colors to species

all_species <- abundance_long %>%
  filter(
    Domain_conv == "Unassigned" |
      (!is.na(Species_conv) & Species_conv != "")
  ) %>%
  mutate(
    Species_conv = if_else(
      Domain_conv == "Unassigned",
      "Unassigned",
      Species_conv
    )
  ) %>%
  distinct(Species_conv) %>%
  arrange(Species_conv) %>%
  pull(Species_conv)
species_only <- all_species[all_species != "Unassigned"]

species_colors <- setNames(
  rep(Cxcx_ext, length.out = length(all_species)),
  all_species)

species_colors["Unassigned"] <- "grey50"
species_colors["Not assigned to species"] <- "grey25"

make_species_barplot <- function(estuary_name) {
  
  # Filter to one estuary
  estuary_data <- abundance_table %>%
    filter(estuary == estuary_name, reads > 0) %>%
    mutate(
      Species = case_when(Domain_conv == "Unassigned" ~ "Unassigned",
                          is.na(Species_conv) | Species_conv == "" ~ "Not assigned to species",
                          TRUE ~ Species_conv))
  
  # Order species by total reads within this estuary
  estuary_barplot <- estuary_data %>%
    group_by(Species) %>%
    mutate(
      total_reads = sum(reads, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      Species = forcats::fct_reorder(
        Species,
        total_reads,
        .desc = TRUE
      ),
      sample_label = paste0(
        "station: ", station,
        "\nreplicate: ", replicate
      )
    )
  
  # Make plot
  p <- ggplot(
    estuary_barplot,
    aes(
      x = sample_id,
      y = relative_reads,
      fill = Species
    )
  ) +
    geom_col() +
    scale_x_discrete(
      labels = setNames(
        estuary_barplot$sample_label,
        estuary_barplot$sample_id
      )
    ) +
    facet_wrap(
      ~ season,
      scales = "free_x"
    ) +
    scale_fill_manual(
      values = species_colors,
      drop = TRUE
    ) +
    labs(
      x = NULL,
      y = "Relative Reads",
      fill = "Species"
    ) +
    ggtitle(
      paste0(estuary_name, " Fish (12S) eDNA")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 0.9
      )
    )
  
  return(p)
}

# list of estuaries
estuaries <- abundance_table %>%
  distinct(estuary) %>%
  filter(!is.na(estuary)) %>%
  pull(estuary)

species_plots <- lapply(estuaries, make_species_barplot)
names(species_plots) <- estuaries

# save all estuary taxonomic species barplots automatically

dir.create("species_barplots", showWarnings = FALSE)

for (e in estuaries) {
  
  ggsave(
    filename = paste0(
      "species_barplots/",
      e,
      "_species_barplot.png"
    ),
    plot = species_plots[[e]],
    height = 6,
    width = 8,
    dpi = 500
  )
}



####### 
# family level barplots for each estuary
all_families <- abundance_long %>%
  filter(
    Domain_conv == "Unassigned" |
      (!is.na(Family_conv) & Family_conv != "")
  ) %>%
  mutate(
    Family = if_else(
      Domain_conv == "Unassigned",
      "Unassigned",
      Family_conv
    )
  ) %>%
  distinct(Family_conv) %>%
  arrange(Family_conv) %>%
  pull(Family_conv)

# connect each family to a color
family_only <- all_families[all_families != "Unassigned"]

family_colors <- setNames(
  rep(Cxcx_ext, length.out = length(all_families)),
  all_families)

family_colors["Unassigned"] <- "grey50"
family_colors["Not assigned to family"] <- "grey25"

make_family_barplot <- function(estuary_name) {
  
  # Filter to one estuary
  estuary_data_fam <- abundance_table %>%
    filter(estuary == estuary_name, reads > 0) %>%
    mutate(
      Family = case_when(Domain_conv == "Unassigned" ~ "Unassigned",
                          is.na(Family_conv) | Family_conv == "" ~ "Not assigned to family",
                          TRUE ~ Family_conv))
  
  # Order species by total reads within this estuary
  estuary_barplot_fam <- estuary_data_fam %>%
    group_by(Family) %>%
    mutate(
      total_reads = sum(reads, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      Family = forcats::fct_reorder(
        Family,
        total_reads,
        .desc = TRUE
      ),
      sample_label = paste0(
        "station: ", station,
        "\nreplicate: ", replicate
      )
    )
  
  # Make plot
  p <- ggplot(
    estuary_barplot_fam,
    aes(
      x = sample_id,
      y = relative_reads,
      fill = Family
    )
  ) +
    geom_col() +
    scale_x_discrete(
      labels = setNames(
        estuary_barplot_fam$sample_label,
        estuary_barplot_fam$sample_id
      )
    ) +
    facet_wrap(
      ~ season,
      scales = "free"
    ) +
    scale_fill_manual(
      values = family_colors,
      drop = TRUE
    ) +
    labs(
      x = NULL,
      y = " Relative Reads",
      fill = "Family"
    ) +
    ggtitle(
      paste0(estuary_name, " Fish (12S) eDNA")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 0.9
      )
    )
  
  return(p)
}

# list of estuaries
estuaries <- abundance_table %>%
  distinct(estuary) %>%
  filter(!is.na(estuary)) %>%
  pull(estuary)

family_plots <- lapply(estuaries, make_family_barplot)
names(family_plots) <- estuaries

# save all estuary taxonomic species barplots automatically

dir.create("family_barplots", showWarnings = FALSE)

for (e in estuaries) {
  ggsave(filename = paste0("family_barplots/", e, "_family_barplot.png"), plot = family_plots[[e]],
         height = 6, width = 8, dpi = 400)}



############################################
# nMDS for each region for each region
# color by estuary, shape by station

community_pa <- abundance_table %>%
  filter(!region %in% c("EB", "EB9")) %>%
  filter(!is.na(Species),
         Species != "",
         Species != "Unassigned",
         Domain != "Unassigned") %>%
  mutate(presence = if_else(reads > 0, 1, 0))

# matrix for nmds
community_matrix <- community_pa %>%
  select(sample_id, Species, presence) %>%
  group_by(sample_id, Species) %>%
  summarise(
    presence = max(presence),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Species,
    values_from = presence,
    values_fill = 0
  )
# snip off metadata
sample_metadata <- abundance_table %>%
  filter(!region %in% c("EB", "EB9")) %>%
  select(sample_id, region, estuary, station, replicate, filter, date, season) %>% distinct()

  
# make sure this snip still lines up with the community matrix
sample_metadata <- sample_metadata %>%
  filter(sample_id %in% community_matrix$sample_id) %>%
  arrange(match(sample_id, community_matrix$sample_id))

# construct nMDS with Jaccard distance - best for presence-absence data
make_nmds_plot <- function(region_name, n_species = 8) {
  
  # --------------------------------
  # 1. Samples belonging to region
  # --------------------------------
  
  region_metadata <- sample_metadata %>%
    filter(region == region_name)
  
  region_samples <- region_metadata$sample_id
  
  
  # --------------------------------
  # 2. Community matrix for region
  # --------------------------------
  
  region_matrix <- community_matrix %>%
    filter(sample_id %in% region_samples) %>%
    select(-sample_id)
  
  rownames(region_matrix) <- region_samples
  
  
  # --------------------------------
  # 3. Run nMDS
  # --------------------------------
  
  set.seed(123)
  
  nmds <- metaMDS(
    region_matrix,
    distance = "jaccard",
    binary = TRUE,
    k = 2,
    trymax = 100,
    trace = FALSE
  )
  
  
  # --------------------------------
  # 4. Extract sample coordinates
  # --------------------------------
  
  nmds_points <- as.data.frame(
    scores(nmds, display = "sites")
  )
  
  nmds_points$sample_id <- rownames(nmds_points)
  
  
  # Add metadata
  nmds_points <- nmds_points %>%
    left_join(
      region_metadata,
      by = "sample_id"
    )
  
  
  # --------------------------------
  # 5. Fit species to nMDS
  # --------------------------------
  
  species_fit <- envfit(
    nmds,
    region_matrix,
    permutations = 999
  )
  
  
  # Extract species vectors
  species_vectors <- as.data.frame(
    scores(
      species_fit,
      display = "vectors"
    )
  )
  
  species_vectors$Species <- rownames(species_vectors)
  
  # Add significance and r2
  species_vectors$p_value <- species_fit$vectors$pvals
  species_vectors$r2 <- species_fit$vectors$r
  
  
  # --------------------------------
  # 6. Select top species
  # --------------------------------
  
  top_species <- species_vectors %>%
    filter(p_value < 0.05) %>%
    arrange(desc(r2)) %>%
    slice_head(n = n_species)
  
  
  # --------------------------------
  # 7. Plot
  # --------------------------------
  
  p <- ggplot(
    nmds_points,
    aes(
      x = NMDS1,
      y = NMDS2,
      fill = estuary,
      shape = season
    )
  ) +
    
    # Sample points
    geom_point(
      size = 3,
      alpha = 0.8,
      color = "black"
    ) +
    
    # Season ellipses
   # stat_ellipse(
    #  aes(group = season),
     # linewidth = 0.8,
    #  alpha = 0.6,
    #  color = "grey20"
    #) +
    
    # Species vectors
    geom_segment(
      data = top_species,
      aes(
        x = 0,
        y = 0,
        xend = NMDS1,
        yend = NMDS2
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.5,
      arrow = arrow(
        length = unit(0.15, "cm")
      )
    ) +
    
    # Species labels
    geom_text_repel(
      data = top_species,
      aes(
        x = NMDS1,
        y = NMDS2,
        label = Species
      ),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    ) +
    
    scale_fill_manual(
      values = Cxcx_minimal,
      guide = guide_legend(
        override.aes = list(
          shape = 21,
          color = "black",
          alpha = 1
        )
      )
    ) +
    
    scale_shape_manual(
      values = c(21, 24)
    ) +
    
    theme_minimal() +
    
    labs(
      x = "NMDS1",
      y = "NMDS2",
      fill = "Estuary",
      shape = "Season",
      title = paste0(
        "nMDS - Region ",
        region_name
      )
    )
  
  
  return(p)
}


regions <- sort(unique(sample_metadata$region))
nmds_sp_plots <- lapply(regions, make_nmds_plot, n_species = 8)
names(nmds_sp_plots) <- regions

# peek
nmds_sp_plots[["CC"]]

# create new directory and save plot
dir.create("nmds_region_species", showWarnings = FALSE)

for (r in regions) {
  ggsave(filename = paste0("nmds_region_species/", r, "_nmds_sp.png"), plot = nmds_sp_plots[[r]],
         height = 6, width = 8, dpi = 400)}


############################################
# SEINE COMPARISON

# pull data from CalEMP portal for F24 + F25
SOP_8_fish_data <- read.csv("https://nexus.sccwrp.org/empachecker/export?tablename=tbl_fish_abundance_data")


# wrangle to match eDNA data
fish_seine <- SOP_8_fish_data %>%
  filter(lubridate::year(samplecollectiondate) %in% c(2024, 2025),
         projectid %in% c("EMPA", "Bight-Estuary-23")) %>%
  mutate(szsn = case_when(lubridate::month(samplecollectiondate) %in% c(1,2,3,4,5) ~ "Spring",
                            lubridate::month(samplecollectiondate) %in% c(8,9,10,11,12) ~ "Fall",
                            TRUE ~ NA_character_)) %>%
  mutate(year = case_when(lubridate::year(samplecollectiondate) %in% 2022 ~ "2022",
                          lubridate::year(samplecollectiondate) %in% 2023 ~ "2023",
                          lubridate::year(samplecollectiondate) %in% 2024 ~ "2024",
                          lubridate::year(samplecollectiondate) %in% 2025 ~ "2025",
                          lubridate::year(samplecollectiondate) %in% 2026 ~ "2026",
                          TRUE ~ NA_character_)) %>%
  mutate(season = paste(szsn, year, sep = " ")) %>%
  filter(season %in% c("Fall 2024", "Fall 2025"))
  
write.csv(fish_seine, here("processed", "fish_seine_calemp.csv"))

fish_seine_1 <- fish_seine %>%
  tidyr::separate(siteid, into = c("region", "estuary"),
                  sep = "-", remove = FALSE) %>%
  rename(Species = scientificname,
         station = stationno) %>%
  mutate(estuary = case_when(region == "NC" & estuary == "BOL" ~ "BOLINAS",
                             region == "SC" & estuary == "BOL" ~ "BOLSA",
                             TRUE ~ estuary)) %>%
  filter(!is.na(Species), Species != "Not recorded")


# Species lookup lists with standardized names - estuary, region, seine and eDNA taxonomic crosswalk
# use taxize package and base crosswalk on most updated WoRMS taxonomic ids
library(taxize)
library(worrms)
library(dplyr)
library(purrr)
library(tidyr)

# make species lookup table by region/ estuary - eDNA
species.lookup.est <- abundance_table %>%
  filter(!is.na(Species)) %>%
  distinct(region, estuary, Species) %>%
  arrange(region, estuary, Species)

edna_species <- abundance_table %>%
  filter(!is.na(Species_conv), Species_conv != "", Species_conv != "NA") %>%
  distinct(Species_conv) %>%
  arrange(Species_conv) %>%
  pull(Species_conv) 

# species lookup table - 
seine_species <- fish_seine_1 %>%
  distinct(Species) %>%
  pull(Species) %>% sort()

# get correct names, AphiaID, and acceptance status
worms_results <- map(
  edna_species,
  ~ tryCatch(
    get_wormsid_(
      sci_com = .x,
      marine_only = FALSE,
      fuzzy = FALSE,
      accepted = FALSE
    ),
    error = function(e) NULL
  )
)

library(tibble)

worms_lookup <- map2_dfr(
  worms_results,
  edna_species,
  function(result, original_name) {
    
    # Handle completely missing results
    if (is.null(result) || length(result) == 0) {
      return(tibble(
        original_name = original_name,
        AphiaID = NA_integer_,
        scientificname = NA_character_,
        authority = NA_character_,
        status = NA_character_
      ))
    }
    
    # Extract nested tibble(s), removing NULL elements
    result_tables <- result %>%
      purrr::compact()
    
    # Handle results containing only NULL
    if (length(result_tables) == 0) {
      return(tibble(
        original_name = original_name,
        AphiaID = NA_integer_,
        scientificname = NA_character_,
        authority = NA_character_,
        status = NA_character_
      ))
    }
    
    # Combine all returned matches
    bind_rows(result_tables) %>%
      mutate(original_name = original_name, .before = 1)
  }
)
# Inspect name matching - eDNA
# How many species were queried?
length(edna_species)

# How many original names are represented?
n_distinct(worms_lookup$original_name)

# How many matched?
sum(!is.na(worms_lookup$AphiaID))

# How many didn't match?
sum(is.na(worms_lookup$AphiaID))

# Status breakdown
worms_lookup %>%
  count(status) %>%
  print() 

# problem children
problem_taxa_edna <- worms_lookup %>%
  filter(status %in% c("unaccepted", "superseded combination", "junior homonym"))

problem_taxa_edna

# full record for problem children
problem_record_edna <- map(
  problem_taxa_edna$AphiaID,
  ~ wm_record(.x)
)

names(problem_record_edna) <- problem_taxa_edna$original_name

# build a corrections table 
problem_corrections <- map2_dfr(
  problem_record_edna,
  names(problem_record_edna),
  function(record, original_name) {
    
    if (is.null(record)) {
      return(tibble(
        original_name = original_name,
        valid_AphiaID = NA_integer_,
        valid_name = NA_character_,
        record_status = NA_character_
      ))
    }
    
    tibble(
      original_name = original_name,
      valid_AphiaID = record$valid_AphiaID,
      valid_name = record$valid_name,
      record_status = record$status
    )
  }
)

problem_corrections

# now build out the same type of list for seine fish names
# get correct names, AphiaID, and acceptance status
# SEINE ----
worms_results_seine <- map(
  seine_species,
  ~ tryCatch(
    get_wormsid_(
      sci_com = .x,
      marine_only = FALSE,
      fuzzy = FALSE,
      accepted = FALSE
    ),
    error = function(e) NULL
  )
)

library(tibble)

worms_lookup_seine <- map2_dfr(
  worms_results_seine,
  seine_species,
  function(result, original_name) {
    
    # Handle completely missing results
    if (is.null(result) || length(result) == 0) {
      return(tibble(
        original_name = original_name,
        AphiaID = NA_integer_,
        scientificname = NA_character_,
        authority = NA_character_,
        status = NA_character_
      ))
    }
    
    # Extract nested tibble(s), removing NULL elements
    result_tables <- result %>%
      purrr::compact()
    
    # Handle results containing only NULL
    if (length(result_tables) == 0) {
      return(tibble(
        original_name = original_name,
        AphiaID = NA_integer_,
        scientificname = NA_character_,
        authority = NA_character_,
        status = NA_character_
      ))
    }
    
    # Combine all returned matches
    bind_rows(result_tables) %>%
      mutate(original_name = original_name, .before = 1)
  }
)
# Inspect name matching - seine
# How many species were queried?
length(seine_species)

# How many original names are represented?
n_distinct(worms_lookup_seine$original_name)

# How many matched?
sum(!is.na(worms_lookup_seine$AphiaID))

# How many didn't match?
sum(is.na(worms_lookup_seine$AphiaID))

# Status breakdown
worms_lookup_seine %>%
  count(status) %>%
  print() 

# problem children - seine
problem_taxa_seine <- worms_lookup_seine %>%
  filter(status %in% c("unaccepted", "superseded combination", "junior homonym", "junior subjective synonym"))

problem_taxa_seine

# full record for problem children
problem_record_seine <- map(
  problem_taxa_seine$AphiaID,
  ~ wm_record(.x)
)

names(problem_record_seine) <- problem_taxa_seine$original_name

# build seine problem corrections table
problem_corrections_seine <- map2_dfr(
  problem_record_seine,
  names(problem_record_seine),
  function(record, original_name) {
    
    if (is.null(record)) {
      return(tibble(
        original_name = original_name,
        valid_AphiaID = NA_integer_,
        valid_name = NA_character_,
        record_status = NA_character_
      ))
    }
    
    tibble(
      original_name = original_name,
      valid_AphiaID = record$valid_AphiaID,
      valid_name = record$valid_name,
      record_status = record$status
    )
  }
)

# correct taxonomic names in each respective dataset
# start with eDNA
edna_corrections <- problem_corrections %>%
  filter(!is.na(valid_name)) %>%
  distinct(original_name, .keep_all = TRUE) %>%
  select(
    Species_conv = original_name,
    valid_name
  )

abundance_table <- abundance_table %>%
  mutate(Species_original = Species_conv) %>%
  left_join(
    edna_corrections,
    by = "Species_conv"
  ) %>%
  mutate(
    Species_conv = coalesce(valid_name, Species_conv)
  ) %>%
  select(-valid_name)

# now correct names with accepted WoRMS status in seine data

seine_corrections <- problem_corrections_seine %>%
  filter(!is.na(valid_name)) %>%
  distinct(original_name, .keep_all = TRUE) %>%
  select(
    Species = original_name,
    valid_name
  )

fish_seine_1 <- fish_seine_1 %>%
  mutate(Species_original = Species) %>%
  left_join(
    seine_corrections,
    by = "Species"
  ) %>%
  mutate(
    Species = coalesce(valid_name, Species)
  ) %>%
  select(-valid_name)



# convert fish to presence-absence
fish_presence <- fish_seine_1 %>%
  group_by(region, estuary, Species) %>%
  summarise(fish_detected = 1, .groups = "drop") %>%
  mutate(across(c(region, estuary, Species), str_trim))

# fish presence at genus level
fish_presence_genus <- fish_presence %>%
  mutate(Genus = sub(" .*", "", Species)) %>%
  group_by(region, estuary, Genus) %>%
  summarise(fish_detected = max(fish_detected), .groups = "drop")
  

# convert eDNA data to presence-absence
eDNA_presence <- abundance_table %>%
  filter(reads > 0, !is.na(Species_conv), Species_conv !="", Species_conv != "Unassigned") %>%
  distinct(region, estuary, Species_conv) %>%
  mutate(eDNA_detected = 1) %>%
  mutate(across(c(region, estuary, Species_conv), str_trim)) %>% rename(Species = Species_conv)

# convert eDNA to presence-absence at genus level
eDNA_presence_genus <- abundance_table %>%
  filter(reads > 0, !is.na(Genus_conv), Genus_conv !="", Genus_conv != "Unassigned") %>%
  distinct(region, estuary, Genus_conv) %>%
  mutate(eDNA_detected = 1) %>%
  mutate(across(c(region, estuary, Genus_conv), str_trim)) %>% rename(Genus = Genus_conv)

# what fish were ONLY detected using eDNA (were not found in a seine anywhere)?
fish_lookup <- fish_seine_1 %>% distinct(Species) #72
eDNA_lookup <- abundance_table %>% distinct(Species_conv) %>% arrange(Species_conv)

# what fish were detected by methods by region and estuary
seine_lookup_est <- fish_seine_1 %>% group_by(region, estuary, Species) %>% distinct(Species) %>% arrange(estuary)
seine_lookup_est

eDNA_lookup_est <- abundance_table %>% group_by(region, estuary, Species_conv) %>% filter(!region %in% c("EB", "EB9")) %>%
  distinct(Species_conv) %>% arrange(estuary) %>% rename(Species = Species_conv)
eDNA_lookup_est

# write out csvs with both seine and eDNA lookup lists 
write.csv(seine_lookup_est, here("processed", "seine_lookup_byestuary.csv"))
write.csv(eDNA_lookup_est, here("processed", "eDNA_lookup_byestuary.csv"))

# fish collection timepoints
fish_seine_2 <- fish_seine_1 %>%
  group_by(region, estuary, Species) %>%
  summarise(
    n_detections = n(),
    n_dates = n_distinct(samplecollectiondate),
    .groups = "drop"
  ) %>%
  arrange(desc(n_detections)) %>%
  filter(Species == "Atherinops affinis")

fish_seine_2

# combine datasets for both methods
comparison <- full_join(eDNA_presence, fish_presence, 
                        by = c("region", "estuary", "Species")) %>%
  mutate(eDNA_detected = replace_na(eDNA_detected, 0),
         fish_detected = replace_na(fish_detected, 0)) %>%
  filter(!str_detect(Species, regex("unknown|juvenile|\\bsp\\b|Argopecten|Navanax|Canis|Hemigrapsus|Homo|Bos tarurus|Cerithidea|Laevicardium|Macoma|Venerupis|NA|Mytilus|Carcinus|Melampus", 
                                    ignore_case = TRUE)))

# 
comparison_1 <- comparison %>%
  mutate(method = case_when(
    eDNA_detected == 1 & fish_detected ==1 ~ "Both",
    eDNA_detected == 1 & fish_detected == 0 ~ "eDNA only",
    eDNA_detected == 0 & fish_detected == 1 ~ "Fish seine only",
    TRUE ~ NA_character_))
  
# write out .csv of this comparative df for team
write.csv(comparison_1, here("processed", "edna_vs_seine_comparison_total.csv"))


### genus- level comparison
comparison_genus <- full_join(eDNA_presence_genus, fish_presence_genus, 
                                        by = c("region", "estuary", "Genus")) %>%
  mutate(eDNA_detected = replace_na(eDNA_detected, 0),
         fish_detected = replace_na(fish_detected, 0)) %>%
  filter(!str_detect(Genus, regex("unknown|juvenile|\\bsp\\b|Argopecten|Navanax|Canis|Homo|Bos|Cerithidea|Laevicardium|Macoma|Venerupis|NA|Mytilus|Carcinus|Melampus", 
                                    ignore_case = TRUE)))

comparison_1_genus <- comparison_genus %>%
  mutate(method = case_when(
    eDNA_detected == 1 & fish_detected ==1 ~ "Both",
    eDNA_detected == 1 & fish_detected == 0 ~ "eDNA only",
    eDNA_detected == 0 & fish_detected == 1 ~ "Fish seine only",
    TRUE ~ NA_character_))

# heatmap of presence-absence in each season in each region
region_data <- comparison_1 %>% filter(estuary %in% c("GOL", "VEN", "MAL"))
region_data <- comparison_1 %>% filter(estuary %in% c("SEAL", "BOL", "NEW", "BAT"))
region_data <- comparison_1 %>% filter(region == "NC")
region_data <- comparison_1_genus %>% filter(region == "CC")

heatmap_plot <- ggplot(region_data, aes( x = estuary, y = Species, fill = method)) +
         geom_tile(color = "white") + 
  scale_fill_manual(
    values = c(
      "eDNA only" = "#3333FF",
      "Fish seine only" = "#33FF33",
      "Both" = "#000033"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  ) +
  labs(x = "Estuary",
    y = "Species",
    fill = "Method",
    title = paste0("Fish community method comparison - ", "SC"))
heatmap_plot

ggsave(here("figures", "heatmap_SC1.jpg"), heatmap_plot,  
       height = 8, width = 6, dpi = 1000)


 # make heatmap better

region_data <- comparison_1 %>% filter(estuary %in% c("GOL", "VEN", "MAL"))
region_data <- comparison_1 %>% filter(estuary %in% c("SEAL", "BOL", "NEW", "BAT"))
region_data <- comparison_1 %>% filter(region == "CC")
region_data <- comparison_1_genus %>% filter(region == "CC")
region_data <- comparison_1 %>% filter(region == "SC")

library(ComplexHeatmap)
library(circlize)
# toggle back between Genus and Species - don't overwrite the plot

heatmap_data <- region_data %>%
  select(Species, estuary, method) %>%
  pivot_wider(
    names_from = estuary,
    values_from = method,
    values_fill = "Not detected"
  )

mat <- as.matrix(heatmap_data[, -1])
rownames(mat) <- heatmap_data$Species


# Convert categorical methods into numeric values
mat_num <- matrix(
  match(mat, c("Not detected", "eDNA only",
               "Fish seine only", "Both")),
  nrow = nrow(mat),
  dimnames = dimnames(mat)
)

heatmap_betterplot <- Heatmap(
  mat_num,
  name = "Detection",
  col = c(
    "1" = "white",
    "2" = "#3333FF",
    "3" = "#33FF33",
    "4" = "#000033"
  ),
  heatmap_legend_param = list(
    at = 1:4,
    labels = c(
      "Not detected",
      "eDNA only",
      "Fish seine only",
      "Both methods"
    ),
    title = "Detection method"
  ),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_names_gp = grid::gpar(fontsize = 10),
  rect_gp = grid::gpar(col = "grey90", lwd = 0.5)
)

png(
  here("figures", "heatmaps", "heatmap_CC.png"),
  width = 6,
  height = 8,
  units = "in",
  res = 1000
)

draw(heatmap_betterplot)

dev.off()

#############################################
# TOTAL READS vs ABUNDANCE OF MOST COMMON SPP
# for species picked up widely with both methods
# Gasterosteus aculeatus, Menida beryllina, Quietula y-cauda, Atherinops affinis, 
# Mugil cephalus,Gillichthys mirabilis, Acanthogobius, Clevelandia ios, Cymatogaster aggregata,
# Eucyclogobius newberryi, Gambusia affinis, Luciana parva, Engraulis mordax, Paralabrax maculofasciatus

seine_abun <- fish_seine_1 %>%
  select(region, estuary, station, season, Species, commonname, status, abundance) %>%
  filter(abundance != -88) 
seine_abun

# only select fish picked up by both methods - combine net replicates for each station
# what we want is the total abundance of each species recorded at each station
seine_abun_trim <- seine_abun %>%
  filter(Species %in% c("Gasterosteus aculeatus", "Menida beryllina", "Quietula y-cauda",
                         "Atherinops affinis", "Mugil cephalus", "Gillichthys mirabilis",
                        "Acanthogobius flavimanus", "Clevelandia ios", "Cymatogaster aggregata",
                         "Eucyclogobius newberryi", "Gambusia affinis", "Luciana parva", "Engraulis mordax",
                        "Paralabrax maculofasciatus")) %>%
  group_by(region, estuary, station, season, Species) %>%
  summarise(Abundance = sum(abundance, na.rm = TRUE),
  across(c(commonname, status), first), .groups = "drop")

seine_abun_trim

seine_abun_GA <- seine_abun_trim %>% filter(Species == "Gastrosteus aculeatus")

# eDNA was not taken side-by-side with seines - combine replicates to get total reads of each species
# recorded at each station
eDNA_reads <- abundance_table %>%
  filter(reads > 0, !is.na(Species), Species !="", Species != "Unassigned") %>%
  distinct(region, estuary, season, station, replicate, Species, reads) %>%
  group_by(region, estuary, season, station, Species) %>%
  summarise(Reads = sum(reads, na.rm = TRUE), .groups = "drop")

eDNA_reads  

eDNA_reads_GA <- eDNA_reads %>% filter(Species == "Gastrosteus aculeatus")

eDNA_reads$station <- as.integer(eDNA_reads$station)

# combine datasets
read_abun_data <- seine_abun_trim %>% left_join(eDNA_reads, by = c("region", "estuary", "station", "season", "Species"))
read_abun_data

library(ggpubr)

# plot abundance vs reads for each taxa
abun_reads_plot <- read_abun_data %>%
  filter(Species %in% c("Gasterosteus aculeatus", "Eucyclogobius newberryi", "Clevelandia ios")) %>%
  ggplot(aes(x = Abundance, y = Reads)) +
  geom_point(alpha = 0.7, size = 4, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  facet_wrap(~commonname, scales = "free", ncol = 4) +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  stat_cor(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "bottom",
    r.digits = 3,
    p.accuracy = 0.0001,
    size = 4,
    aes(label = paste(after_stat(rr.label),
                      after_stat(p.label),
                      sep = "~`,`~"))
  ) +
  theme_bw() +
  theme(strip.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 14, hjust = 1, angle = 45),
        axis.text.y = element_text(size = 16)) +
  labs(
    x = "Fish abundance (log 1 + x)",
    y = "eDNA reads (log 1 + x)"
  )
abun_reads_plot

ggsave(here("figures", "reads_vs_abundance_three.jpg"), abun_reads_plot, 
       height = 8, width = 14, dpi = 500)

#ggsave(here("figures", "reads_vs_abundance.jpg"), abun_reads_plot, 
 #      height = 8, width = 8, dpi = 500)

#############################################
# INVASIVES - map 




##############################################
# Venn Diagrams
library(eulerr)
library(ggplotify)

comparison_2 <- comparison_1 %>% filter(region %in% c("SC", "CC", "NC")) 


edna_set <- comparison_2 %>%
  filter(region =="CC", eDNA_detected == 1) %>%
  pull(Species) %>%
  unique()

seine_set <- comparison_2 %>%
  filter(region =="CC", fish_detected == 1) %>%
  pull(Species) %>%
  unique()

# Check the actual numbers
length(edna_set)
length(seine_set)
length(intersect(edna_set, seine_set))

fit <- euler(list(
  "eDNA" = edna_set,
  "Fish seine" = seine_set
))

venn_diagram <- plot(
  fit,
  fills = list(
    fill = c("#2878B5", "#D98B35"),
    alpha = 0.65
  ),
  edges = list(col = "grey30"),
  labels = list(font = 2, fontsize = 20),
  quantities = list(
    type = "counts",
    fontsize = 20,
    font = 2
  ),
  main = "12S eDNA metabarcoding vs. fish seine - Central Coast"
)

venn_diagram

png(here("figures", "venn_diagram_CC.png"), width = 2000, height = 1000, res = 200)
plot(venn_diagram)
dev.off()


# heatmap of ADLC only



heatmap_ADLC <- comparison_2 %>%
  filter(estuary == "ADLC") %>%
  ggplot(aes(x = method, y = reorder(Species, method), fill = method)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c(
    "eDNA only" = "#3B9AB2",
    "Both" = "#F21A00",
    "Fish seine only" = "#78B7C5"
  )) +
  labs(
    title = "Species Detection by Method",
    subtitle = "ADLC Estuary",
    x = "Detection Method",
    y = "Species",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "italic", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
heatmap_ADLC

##############################################
# Estuary type (lagoonal vs open)
estuary_info <- fish_seine_1 %>%
  select(region, estuary, estuaryname, estuaryclass, mpastatus, estuarytype) %>%
  mutate(estuary = case_when(
    region == "SC" & estuary == "BOL" ~ "BOLSA",
    region == "NC" & estuary == "BOL" ~ "BOLINAS",
    TRUE ~ estuary
  )) %>%
  distinct()

estuary_info

# fish seine nMDS and PERMANOVA
fish_abun_table <- fish_seine_1 %>%
  mutate(Genus = sub(" .*", "", Species)) %>%
  mutate(sample_id = cur_group_id(),
         .by = c(estuary, season, station)) 

# give each sample id a unique number so we get rid of duplicates
fish_abun_table <- fish_abun_table %>%
  group_by(region, estuary, season, station) %>%
  mutate(sample_id = cur_group_id()) %>%
  ungroup()


community_fish <- fish_abun_table %>%
  filter(!is.na(Genus)) %>%
  filter(!region %in% c("EB", "EB9")) %>%
  group_by(sample_id, Genus) %>%
  summarise(
    abundance = sum(abundance[abundance >= 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(presence = if_else(abundance > 0, 1, 0))

fish_metadata <- fish_abun_table %>%
  filter(!region %in% c("EB", "EB9")) %>%
  select(
    sample_id,
    region,
    estuary,
    estuaryclass,
    mpastatus,
    estuarytype,
    season
  ) %>%
  distinct()

# fish community matrix and nMDS
comm_matrix_fish <- community_fish %>%
  group_by(sample_id, Genus) %>%
  summarise(
    abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Genus,
    values_from = abundance,
    values_fill = 0
  ) %>%
  column_to_rownames("sample_id") %>%
  as.data.frame()

# remove zero rows (shouldn't be any)
comm_matrix_fish <- comm_matrix_fish[
  rowSums(comm_matrix_fish) > 0,
  ,
  drop = FALSE
]

#make sure all those sample ids are of the same class
metadata_nmds <- fish_metadata %>%
  mutate(sample_id = as.character(sample_id)) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  filter(sample_id %in% rownames(comm_matrix_fish)) %>%
  slice(match(rownames(comm_matrix_fish), sample_id))

# run nmds - stress = 0.15 (acceptable)
set.seed(123)

nmds_fish <- metaMDS(
  comm_matrix_fish,
  distance = "bray",
  k = 2,
  trymax = 100,
  autotransform = FALSE,
  wascores = TRUE
)

nmds_scores <- as.data.frame(
  scores(nmds_fish, display = "sites")
) %>%
  rownames_to_column("sample_id") %>%
  left_join(metadata_nmds, by = "sample_id")

cols <- c(
  "NC" = "darkgreen",
  "CC" = "cornflowerblue",
  "SC" = "deeppink"
)

nmds_plot <- ggplot(
  nmds_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    colour = region,
    shape = estuarytype
  )
) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = cols) +
  theme_classic() +
  labs(
    colour = "Region",
    shape = "Estuary class"
  )

print(nmds_plot)

ggsave(here("figures", "nmds_estuarytype_seine.jpg"), 
       nmds_plot, height = 6, width = 8, dpi = 400)

# PERMANOVA
bray_dist <- vegdist(comm_matrix_fish, method = "bray")


permanova <- adonis2(
  comm_matrix_fish ~ estuaryclass + estuarytype + region,
  data = fish_metadata,
  method = "bray",
  permutations = 999,
  by = "margin"
)

permanova

# eDNA nMDS and PERMANOVA - pool replicates at each station
abun_table_est <- abundance_table %>% left_join(estuary_info, by = "estuary")

#############################
pooled <- abun_table_est %>%
  filter(!is.na(Genus_conv)) %>%
  filter(!region.x %in% c("EB", "EB9")) %>%
  group_by(region.x, estuary, station, season) %>%
  mutate(sample_id = cur_group_id()) %>%
  ungroup()


community <- pooled %>%
  group_by(sample_id, region.x, estuary, station, season, Genus_conv) %>%
  summarise(reads = sum(reads), .groups = "drop") %>%
  mutate(presence = as.integer(reads > 0)) %>%
  select(-reads) %>%
  pivot_wider(
    names_from = Genus_conv,
    values_from = presence,
    values_fill = 0
  )


metadata <- pooled %>%
  select(sample_id, region.x, estuary, station, season,
         estuaryclass, mpastatus, estuarytype) %>%
  distinct(sample_id, .keep_all = TRUE)


comm_matrix <- community %>%
  select(-sample_id, -region.x, -estuary, -station, -season) %>%
  as.data.frame()

rownames(comm_matrix) <- community$sample_id


nmds <- metaMDS(
  comm_matrix,
  distance = "jaccard",
  binary = TRUE,
  k = 2,
  trymax = 100,
  autotransform = FALSE
)


nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("sample_id") %>%
  mutate(sample_id = as.integer(sample_id)) %>%
  left_join(metadata, by = "sample_id")

cols <- c(
  "NC" = "darkgreen",
  "CC" = "cornflowerblue",
  "SC" = "deeppink"
)

nmds_plot <- ggplot(
  nmds_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    colour = region.x,
    shape = estuarytype
  )
) +
  geom_point(size = 4, alpha = 0.8) +
  theme_classic() +
  scale_color_manual(values = cols) +
  labs(
    title = "eDNA community composition",
    colour = "Region",
    shape = "Estuary Open/Closed status"
  )

print(nmds_plot)



ggsave(here("figures", "nmds_estuarytype_pooled.jpg"), 
       nmds_plot, height = 6, width = 8, dpi = 400)

# PERMANOVA
bray_dist <- vegdist(comm_matrix, method = "jaccard")

metadata_clean <- metadata %>%
  filter(
    !is.na(estuaryclass),
    !is.na(mpastatus),
    !is.na(region.x),
    !is.na(estuarytype)
  )

comm_matrix_clean <- comm_matrix[
  metadata_clean$sample_id,
]

permanova <- adonis2(
  comm_matrix_clean ~ estuaryclass + estuarytype + region.x,
  data = metadata_clean,
  method = "bray",
  permutations = 999,
  by = "margin"
)

permanova

##############################################
# Species of special concern
fish_info <- read.csv(here("processed", "Fish-index-info.csv")) # info from expert panel on fish life and guild info

trait_comparison <- comparison_1 %>%
  left_join(fish_info, by = c("Species" = "Scientific.name"))

species_summary <- trait_comparison %>%
  group_by(Species) %>%
  summarise(
    eDNA_detected = max(eDNA_detected, na.rm = TRUE),
    fish_detected = max(fish_detected, na.rm = TRUE),
    Ecological.category = first(Ecological.categories),
    .groups = "drop"
  )

ecological_summary <- species_summary %>%
  group_by(Ecological.category) %>%
  summarise(
    Total_species = n(),
    eDNA_detections = sum(eDNA_detected),
    Seine_detections = sum(fish_detected),
    Both_methods = sum(eDNA_detected == 1 & fish_detected == 1),
    eDNA_only = sum(eDNA_detected == 1 & fish_detected == 0),
    Seine_only = sum(eDNA_detected == 0 & fish_detected == 1),
    .groups = "drop"
  ) %>%
  mutate(
    eDNA_percent = round(eDNA_detections / Total_species * 100, 1),
    Seine_percent = round(Seine_detections / Total_species * 100, 1)
  ) %>%
  arrange(desc(Total_species))


write.csv(ecological_summary, here("processed", "ecological.summary.csv"))
# picks up more marine (70% vs 45%) and freshwater taxa

# create table of species with microhabitat associations
species_summary_microhabitat <- trait_comparison %>%
  group_by(Species) %>%
  summarise(
    eDNA_detected = max(eDNA_detected, na.rm = TRUE),
    fish_detected = max(fish_detected, na.rm = TRUE),
    Micro.habitat.association = first(Micro.habitat.association),
    .groups = "drop"
  )

ecological_summary_microhabitat <- species_summary_microhabitat %>%
  group_by(Micro.habitat.association) %>%
  summarise(
    Total_species = n(),
    eDNA_detections = sum(eDNA_detected),
    Seine_detections = sum(fish_detected),
    Both_methods = sum(eDNA_detected == 1 & fish_detected == 1),
    eDNA_only = sum(eDNA_detected == 1 & fish_detected == 0),
    Seine_only = sum(eDNA_detected == 0 & fish_detected == 1),
    .groups = "drop"
  ) %>%
  mutate(
    eDNA_percent = round(eDNA_detections / Total_species * 100, 1),
    Seine_percent = round(Seine_detections / Total_species * 100, 1)
  ) %>%
  arrange(desc(Total_species))

# picks up more hard structure assoicated taxa

# let's look at trophic level
species_summary_trophic <- trait_comparison %>%
  group_by(Species) %>%
  summarise(
    eDNA_detected = max(eDNA_detected, na.rm = TRUE),
    fish_detected = max(fish_detected, na.rm = TRUE),
    Trophic.level = first(Trophic.level),
    .groups = "drop"
  )

ecological_summary_trophic <- species_summary_trophic %>%
  group_by(Trophic.level) %>%
  summarise(
    Total_species = n(),
    eDNA_detections = sum(eDNA_detected),
    Seine_detections = sum(fish_detected),
    Both_methods = sum(eDNA_detected == 1 & fish_detected == 1),
    eDNA_only = sum(eDNA_detected == 1 & fish_detected == 0),
    Seine_only = sum(eDNA_detected == 0 & fish_detected == 1),
    .groups = "drop"
  ) %>%
  mutate(
    eDNA_percent = round(eDNA_detections / Total_species * 100, 1),
    Seine_percent = round(Seine_detections / Total_species * 100, 1)
  ) %>%
  arrange(desc(Total_species))

View(ecological_summary_trophic)

##############################################

# overview heatmap of taxa by ecological category
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

heatmap_data <- trait_comparison %>%
  mutate(
    Detection = case_when(
      eDNA_detected == 1 & fish_detected == 1 ~ "Both",
      eDNA_detected == 1 & fish_detected == 0 ~ "eDNA only",
      eDNA_detected == 0 & fish_detected == 1 ~ "Seine only",
      TRUE ~ "Neither"
    )
  )


heatmap_data <- heatmap_data %>%
  mutate(
    detection_value = case_when(
      Detection == "Neither" ~ 0,
      Detection == "eDNA only" ~ 1,
      Detection == "Seine only" ~ 2,
      Detection == "Both" ~ 3
    )
  )

heatmap_matrix <- heatmap_data %>%
  select(Species, estuary, detection_value) %>%
  distinct() %>%
  pivot_wider(
    names_from = estuary,
    values_from = detection_value,
    values_fill = 0
  ) %>%
  column_to_rownames("Species") %>%
  as.matrix()

species_metadata <- heatmap_data %>%
  select(
    Species,
    Detection,
    Ecological.categories,
    Trophic.level,
    Micro.habitat.association
  ) %>%
  distinct() %>%
  group_by(Species) %>%
  slice(1) %>%
  ungroup()

species_metadata <- species_metadata %>%
  slice(match(rownames(heatmap_matrix), Species))

rownames(species_metadata) <- species_metadata$Species

estuary_metadata <- heatmap_data %>%
  select(estuary, region) %>%
  distinct() %>%
  arrange(region, estuary) %>%
  as.data.frame()

rownames(estuary_metadata) <- estuary_metadata$estuary

estuary_metadata <- estuary_metadata[colnames(heatmap_matrix), ]

detection_colours <- c(
  "0" = "white",
  "1" = "#3333FF",
  "2" = "#33FF33",
  "3" = "#000033"
)


region_colours <- c(
  "NC" = "darkgreen",
  "CC" = "cornflowerblue",
  "SC" = "deeppink"
)



eco_levels <- unique(na.omit(species_metadata$Ecological.categories))

eco_colours <- setNames(
  colorRampPalette(
    c("grey20", "darkolivegreen", "#E34A33", "#B30000")
  )(length(eco_levels)),
  eco_levels
)

top_ha <- HeatmapAnnotation(
  Region = estuary_metadata$region,
  col = list(
    Region = region_colours
  )
)

right_ha <- rowAnnotation(
  Ecological.category = species_metadata$Ecological.categories,
  col = list(
    Ecological.category = eco_colours
  )
)



ht <- Heatmap(
  heatmap_matrix ,
  
  name = "Detection",
  
  col = detection_colours,
  
  #top_annotation = top_ha,
  right_annotation = right_ha,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  row_split = species_metadata$Ecological.categories,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 9),
  
  column_names_rot = 45,
  
  column_title = "Fish detections across California estuaries",
  row_title = "Fish species",
  
  border = TRUE,
  
  heatmap_legend_param = list(
    at = 0:3,
    labels = c(
      "Neither",
      "eDNA only",
      "Seine only",
      "Both"
    )
  )
)

draw(ht)

png(
  here("figures", "heatmaps", "heatmap_EcoCat.png"),
  width = 10,
  height = 12,
  units = "in",
  res = 1000
)

draw(ht)

dev.off()


###############################################
# stats
head(abundance_table)

sample_summary <- abundance_table %>%
  distinct(sample_id, estuary, station, season) %>%
  group_by(estuary, station, season) %>%
  summarise(n_distinct(sample_id),
            .groups = "drop") %>%
  arrange(season, estuary, station)
sample_summary
View(sample_summary)
