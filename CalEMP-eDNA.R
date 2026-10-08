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
asvs.091726 <- read.csv(here("CalEMP-091726-data", "exported", "ASV_CalEMP_091726.csv"))
asvs.testplate <- read.csv(here("CalEMP-Fall2024-testplate", "exported", "table_mifish_EMPA_test.csv")) # these samples were combined into 070726 dataset
asvs.070726 <- read.csv(here("CalEMP-070726-data", "exported", "table_mifish_CalEMP_070726.csv"))
asvs.093026 <- read.csv(here("CalEMP-093026-data", "exported_12S", "ASV_CalEMP_093026.csv")) # taxonomy assigned using normal 12S db
asvs.093026_newref <- read.csv(here("CalEMP-093026-data", "exported_newref", "table_mifish_CalEMP_093026_newref.csv")) #taxonomy assigned using eDNA explorer 12S db

# sample inventory
sample.inventory <- read.csv(here("samples_received_9.23.26.csv"))
samples.received <- n_distinct(sample.inventory$samples.received)
total.samples <- n_distinct(sample.inventory$masterInventory)
samples.remaining <- total.samples - samples.received

print(samples.received)
print(total.samples)

################################################
# stats: percent total reads assigned per dataset ----
## testplate dataset

# collapse identical ASVs - reuse this function
reads.per.CL.testplate <- asvs.testplate %>%
  group_by(ConsensusLineage) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

## reads across whole dataset
total.reads.testplate <- sum(reads.per.CL.testplate %>% select(where(is.numeric)), na.rm = TRUE)

# assigned reads 
assigned.reads.testplate <-reads.per.CL.testplate %>%
  filter(ConsensusLineage != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)

# percent unassigned
perc.assigned.testplate <- 100*(assigned.reads.testplate)/total.reads.testplate

# print stats - testplate - 87.91%
print(perc.assigned.testplate)

## 070726 dataset

# collapse identical ASVs - reuse this function
reads.per.CL.070726 <- asvs.070726 %>%
  group_by(ConsensusLineage) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

## reads across whole dataset
total.reads.070726 <- sum(reads.per.CL.070726 %>% select(where(is.numeric)), na.rm = TRUE)

# unassigned reads 
assigned.reads.070726 <-reads.per.CL.070726 %>%
  filter(ConsensusLineage != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)

# percent unassigned
perc.assigned.070726 <- 100*(assigned.reads.070726)/total.reads.070726

# print stats - 070726 - 96.03%
print(perc.assigned.070726)

##################
## 091726 dataset

# collapse identical ASVs - reuse this function
reads.per.CL.091726 <- asvs.091726 %>%
  group_by(ConsensusLineage) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

## reads across whole dataset
total.reads.091726 <- sum(reads.per.CL.091726 %>% select(where(is.numeric)), na.rm = TRUE)

# unassigned reads 
assigned.reads.091726 <-reads.per.CL.091726 %>%
  filter(ConsensusLineage != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)

# percent unassigned
perc.assigned.091726 <- 100*(assigned.reads.091726)/total.reads.091726

# print stats - 091726 - 92.42%
print(perc.assigned.091726)

################################
## 093026 dataset

# collapse identical ASVs - reuse this function
reads.per.CL.093026 <- asvs.093026 %>%
  group_by(ConsensusLineage) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

## reads across whole dataset
total.reads.093026 <- sum(reads.per.CL.093026 %>% select(where(is.numeric)), na.rm = TRUE)

# unassigned reads 
assigned.reads.093026 <-reads.per.CL.093026 %>%
  filter(ConsensusLineage != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)

# percent unassigned
perc.assigned.093026 <- 100*(assigned.reads.093026)/total.reads.093026

# print stats - 093026 - 96.16 %
print(perc.assigned.093026)

################################
## 093026 - NEWREF dataset - 

# collapse identical ASVs - reuse this function
reads.per.CL.093026_newref <- asvs.093026_newref %>%
  group_by(ConsensusLineage) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

## reads across whole dataset
total.reads.093026_newref <- sum(reads.per.CL.093026_newref %>% select(where(is.numeric)), na.rm = TRUE)

# unassigned reads 
assigned.reads.093026_newref <-reads.per.CL.093026_newref %>%
  filter(ConsensusLineage != "Unassigned") %>%
  select(where(is.numeric)) %>%
  sum(na.rm = TRUE)

# percent unassigned
perc.assigned.093026_newref <- 100*(assigned.reads.093026_newref)/total.reads.093026_newref

# print stats - 093026 NEWREF - 73.54 %
print(perc.assigned.093026_newref)

# almost identical percent reads assigned to some level for the two different 12S databases


###################################################

#### Compare standard 12S ref database to eDNA explorer ref db ----
# use 093017 samples as a test for this - both assigned at a 95% identity threshold #

#  what percentage of ASVs were assigned to species by each database

# join two datasets and collapse Unassigned ASVs into one row
comb_093026 <- asvs.093026 %>% full_join(asvs.093026_newref %>% select(OTU.ID, ConsensusLineage), by = "OTU.ID") %>%
  rename(ConsensusLineage.12S = ConsensusLineage.x) %>%
  rename(ConsensusLineage.EE = ConsensusLineage.y) %>%
  separate(
    ConsensusLineage.12S,
    into = c("Domain_12S", "Phylum_12S", "Class_12S", "Order_12S",
             "Family_12S", "Genus_12S", "Species_12S"),
    sep = ";\\s*",
    remove = FALSE
  ) %>%
  mutate(
    across(
      c(Domain_12S, Phylum_12S, Class_12S, Order_12S, Family_12S, Genus_12S, Species_12S),
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
  )

# assigned to species by each database
assigned.to.species.12S <- comb_093026 %>% filter(!is.na(Species_12S)) %>% 
  summarise(across(where(is.numeric), sum), .groups = "drop") %>% sum()

assigned.to.species.EE <- comb_093026 %>% filter(!is.na(Species_EE)) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>% sum()

#total assigned reads (from above)
assigned.reads.093026
assigned.reads.093026_newref

# percent 12S assigned to species - 44.56%
perc.assigned.to.sp.12S <- 100*(assigned.to.species.12S/assigned.reads.093026)
perc.assigned.to.sp.12S

# percent EE assigned to species - 69.54% (out of an initial 73% of total receiving any assignment)
per.assigned.to.sp.EE <- 100*(assigned.to.species.EE/ assigned.reads.093026_newref)
per.assigned.to.sp.EE

########## we'll stop analysis here since there was a clear winner in which database assigned more reads overall
  
# what percent of total ASVs were assigned to species?
# for the ASVs assigned to different levels by different databases:
  # (a) how many agreed?
  # (b) how many were assigned lower by (i) 12S usual and (ii) eDNA explorer db?
  # (c) list of taxa coverage done better by each db

####
# wrangle datasets together ----
# combine ASV tables
combined <- asvs.070726 %>%
  full_join(asvs.091726 %>% select(-ConsensusLineage),
            by = "OTU.ID") %>%
  full_join(asvs.093026 %>% select(-ConsensusLineage),
            by = "OTU.ID")

sample.bind <- bind_rows(asvs.093026, asvs.091726, asvs.070726) %>%
  mutate(across(where(is.numeric), ~ replace(.x, is.na(.x), 0)))


# collapse ASVs by ConsensusLineage
collapsed <- sample.bind %>%
  group_by(ConsensusLineage) %>%
  summarise(across(where(is.numeric), sum),
            .groups = "drop")

# clean taxonomic hierarchy
asv.calemp <- collapsed %>%
  separate(
  ConsensusLineage,
  into = c("Domain", "Phylum", "Class", "Order",
           "Family", "Genus", "Species"),
  sep = ";\\s*",
  remove = FALSE
) %>%
  mutate(
    across(
      c(Domain, Phylum, Class, Order, Family, Genus, Species),
      ~ sub("^[a-z]__", "", .x)
    )
  )

# pivot long 
abundance_long <- asv.calemp %>% 
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
  ))

# attach metadata to abundance table
abundance_table <- abundance_long %>%
  left_join(sample_metadata, by = "sample_id")

# taxonomic barplots
# let's start with one estuary
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
              "grey20", "#FFFFFF", "#003300", "#FF0033","dodgerblue",
              "#330033", "#FF3300", "#3333FF", "#000033", "#FFCC66")

# assign colors to species

all_species <- abundance_long %>%
  filter(
    Domain == "Unassigned" |
      (!is.na(Species) & Species != "")
  ) %>%
  mutate(
    Species = if_else(
      Domain == "Unassigned",
      "Unassigned",
      Species
    )
  ) %>%
  distinct(Species) %>%
  arrange(Species) %>%
  pull(Species)
species_only <- all_species[all_species != "Unassigned"]

species_colors <- setNames(
  rep(Cxcx_ext, length.out = length(all_species)),
  all_species)

species_colors["Unassigned"] <- "grey50"

make_species_barplot <- function(estuary_name) {
  
  # Filter to one estuary
  estuary_data <- abundance_table %>%
    filter(estuary == estuary_name) %>%
    filter(reads > 0) %>%
    filter(
      Domain == "Unassigned" |
        (!is.na(Species) & Species != "")
    ) %>%
    mutate(
      Species = if_else(
        Domain == "Unassigned",
        "Unassigned",
        Species
      )
    )
  
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
      y = reads,
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
      y = "Reads",
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
    dpi = 400
  )
}



####### 
# family level barplots for each estuary
all_families <- abundance_long %>%
  filter(
    Domain == "Unassigned" |
      (!is.na(Family) & Family != "")
  ) %>%
  mutate(
    Species = if_else(
      Domain == "Unassigned",
      "Unassigned",
      Family
    )
  ) %>%
  distinct(Family) %>%
  arrange(Family) %>%
  pull(Family)

# connect each family to a color
family_only <- all_families[all_families != "Unassigned"]

family_colors <- setNames(
  rep(Cxcx_ext, length.out = length(all_families)),
  all_families)

family_colors["Unassigned"] <- "grey50"

make_family_barplot <- function(estuary_name) {
  
  # Filter to one estuary
  estuary_data_fam <- abundance_table %>%
    filter(estuary == estuary_name) %>%
    filter(reads > 0) %>%
    filter(
      Domain == "Unassigned" |
        (!is.na(Family) & Family != "")
    ) %>%
    mutate(
      Species = if_else(
        Domain == "Unassigned",
        "Unassigned",
        Family
      )
    )
  
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
      y = reads,
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
      scales = "free_x"
    ) +
    scale_fill_manual(
      values = family_colors,
      drop = TRUE
    ) +
    labs(
      x = NULL,
      y = "Reads",
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

# write out species list 
write.csv(all_species, here("processed", "species_lookup.csv"))

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

# correct errors and remove EB from nmds or it will get fucked
sample_metadata <- sample_metadata %>%
  mutate(region = if_else(region == "C", "CC", region))
  

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
         projectid %in% c("EMPA", "Bight-Estuary-23"),
         method == "count") %>%
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
  

fish_seine_1 <- fish_seine %>%
  tidyr::separate(siteid, into = c("region", "estuary"),
                  sep = "-", remove = FALSE) %>%
  rename(Species = scientificname,
         station = stationno) %>%
  filter(!is.na(Species), Species != "Not recorded")

# convert fish to presence-absence
fish_presence <- fish_seine_1 %>%
  group_by(region, estuary, station, season, Species) %>%
  summarise(fish_detected = 1, .groups = "drop")

# convert eDNA data to presence-absence
eDNA_presence <- abundance_table %>%
  filter(reads > 0, !is.na(Species), Species !="", Species != "Unassigned") %>%
  distinct(region, estuary, season, station, Species, reads) %>%
  mutate(edna_detected = 1)

eDNA_presence$station <- as.integer(eDNA_presence$station)

# combine datasets for both methods
comparison <- full_join(eDNA_presence, fish_presence, by = c("region", "estuary", "station", "season", "Species")) %>%
  mutate(eDNA_detected = replace_na(edna_detected, 0),
         fish_detected = replace_na(fish_detected, 0))


comparison_1 <- comparison %>%
  mutate(method = case_when(
    eDNA_detected == 1 & fish_detected ==1 ~ "Both",
    eDNA_detected == 1 & fish_detected == 0 ~ "eDNA only",
    eDNA_detected == 0 & fish_detected == 1 ~ "Fish seine only"))
  
# write out .csv of this comparative df for team
write.csv(comparison_1, here("processed", "edna_vs_seine_comparison.csv"))

# heatmap of presence-absence in each season in each region
region_data <- comparison_1 %>% filter(estuary %in% c("GOL", "VEN", "MAL"))
region_data <- comparison_1 %>% filter(estuary %in% c("SEAL", "BOL", "NEW", "BAT"))
region_data <- comparison_1 %>% filter(region == "NC")
region_data <- comparison_1 %>% filter(region == "CC")

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

ggsave(here("figures", "heatmap_CC.jpg"), heatmap_plot,  
       height = 10, width = 8, dpi = 500)

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

# eDNA was not taken side-by-side with seines - combine replicates to get total reads of each species
# recorded at each station
eDNA_reads <- abundance_table %>%
  filter(reads > 0, !is.na(Species), Species !="", Species != "Unassigned") %>%
  distinct(region, estuary, season, station, replicate, Species, reads) %>%
  group_by(region, estuary, season, station, Species) %>%
  summarise(Reads = sum(reads, na.rm = TRUE), .groups = "drop")

eDNA_reads  

eDNA_reads$station <- as.integer(eDNA_reads$station)

# combine datasets
read_abun_data <- seine_abun_trim %>% left_join(eDNA_reads, by = c("region", "estuary", "station", "season", "Species"))
read_abun_data



# plot abundance vs reads for each taxa
abun_reads_plot <- ggplot(read_abun_data, aes(x = Abundance, y = Reads)) +
  geom_point(alpha = 0.7, na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, na.rm = TRUE) +
  facet_wrap(~commonname, scales = "free", ncol = 4) +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  theme_bw() +
  labs(
    x = "Fish abundance (log 1 + x)",
    y = "eDNA reads (log 1 + x)"
  )
abun_reads_plot
ggsave(here("figures", "reads_vs_abundance.jpg"), abun_reads_plot, 
       height = 8, width = 8, dpi = 500)

#############################################
# INVASIVES - map 



