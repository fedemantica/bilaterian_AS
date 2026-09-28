########################################################################
########### Script to filter TS exons by minimum coverage ##############
########################################################################

##### Upload libraries
library(argparser)
library(tidyverse)

##### Define arguments
#Create a parser
p = arg_parser("Require that TS exons have at least one noVLOW or higher-quality sample in the tissue of interest")

#Add command line arguments
p = add_argument(p, arg="--ts_exons", help="File with the modified output of get_tisAS.pl, which contains the TS exons to filter")
p = add_argument(p, arg="--inclusion_table", help="Inclusion table for the species of interest")
p = add_argument(p, arg="--sample_tissue", help="File including the sample - tissue correspondence for all species of interest")
p = add_argument(p, arg="--taus", help="File with event and corresponding Tau")
p = add_argument(p, arg="--exons_gene_dict", help="File with correspondence between vastID and geneID")
p = add_argument(p, arg="--gene_orthogroups", help="Gene orthogroup file containing bilaterian conserved genes, with col1=OG_ID, col2=Species, col3=GeneID")
p = add_argument(p, arg="--output", help="Output file with panAS exons")
p = add_argument(p, arg="--min_dPSI_ts", help="minimum delta PSI to select one tissue, two tissues or no tissues", type="numeric")
p = add_argument(p, arg="--specialization_PSI_cutoff", help="maximum PSI in other tissues for an exon to be considered TS", type="numeric")

#Parse the command line arguments
argv = parse_args(p)

#Define arguments
ts_exons_file = argv$ts_exons
inclusion_table_file = argv$inclusion_table
sample_tissue_file = argv$sample_tissue
taus_file = argv$taus
exon_to_gene_file = argv$exons_gene_dict
gene_orthogroups_file = argv$gene_orthogroups
output_file = argv$output
min_dPSI_ts = argv$min_dPSI_ts
specialization_PSI_cutoff = argv$specialization_PSI_cutoff

#### Header of the ts_exons_file in input
# GENE
# EventID
# COORD
# LENGTH
# FullCO
# COMPLEX
# Tissue
# Direction
# N_Rep
# Total_N_others
# Av_Tis
# Av_Others
# dPSI
# SD_Tis
# SD_Others
# Min_Others
# Max_Others
# Valid_Tis
# dPSI_First_Second 
# dPSI_Second_Third 
# dPSI_First_Third
# Second_Tis

#########################################
########## MAIN #########################
#########################################

#### Read in files
vast_tools_df = read.delim(inclusion_table_file, header=TRUE)
sample_tissue_df = read.delim(sample_tissue_file, header=TRUE, col.names=c("Sample", "Tissue"))
ts_exons_df = read.delim(ts_exons_file, header=TRUE)
taus_df =  read.delim(taus_file, header=TRUE)
exon_to_gene_df = read.delim(exon_to_gene_file, header=TRUE)
gene_orthogroups_df = read.delim(gene_orthogroups_file, header=FALSE, col.names=c("OG_ID", "Species", "GeneID"))

#### Set variables
#Samples
species_samples = colnames(vast_tools_df)[grep("\\.Q", colnames(vast_tools_df), invert=TRUE)]
species_samples = species_samples[!(species_samples %in% c("GENE", "EVENT", "COORD", "LENGTH", "FullCO", "COMPLEX"))]

#Quality of samples
species_samples_quality = colnames(vast_tools_df)[grep("\\.Q", colnames(vast_tools_df), invert=FALSE)]

#########################################
########## GENERATE COVERAGE INFO #######
#########################################

#### For each event, assign to each tisuse the NO_COV or YES_COV labels depending of all samples have quality N/VLOW or not.
vast_tools_coverage_long_df = vast_tools_df %>%
  dplyr::select("EVENT", all_of(species_samples_quality)) %>%
  filter(grepl("EX", EVENT)) %>% #select only cassette exon events
  pivot_longer(cols=species_samples_quality, names_to = "Sample", values_to = "Quality_scores") %>%
  mutate(Sample = sub("\\.Q", "", Sample),
         First_quality_score = sub(",.*", "", Quality_scores)) %>%
  ### Add tissue
  left_join(sample_tissue_df, by="Sample") %>%
  ### Add YES_COV or NO_COV labels
  mutate(Coverage = ifelse(First_quality_score %in% c("N", "VLOW"), "NO_COV", "YES_COV")) %>%
  group_by(EVENT, Tissue) %>%
  summarize(Tissue_coverage = ifelse(all(Coverage=="NO_COV"), "NO_COV", "YES_COV")) %>%
  ungroup()

###############################################
######## ADD ALL EXTRA INFO ###################
###############################################

ts_exons_extra_info_df = ts_exons_df %>%
  ### Associate Tissue with tissue specificity based on delta PSI between top tissues
  mutate(Associated_Tis = case_when(dPSI_First_Second >= min_dPSI_ts ~ Tissue,
                                        dPSI_First_Second < min_dPSI_ts & dPSI_Second_Third >= min_dPSI_ts ~ paste0(Tissue, ";", Second_Tis),
                                        TRUE ~ "NO_tissue")) %>%
  ### Add coverage information
  left_join(vast_tools_coverage_long_df, by=c("EventID" = "EVENT", "Tissue" = "Tissue")) %>%
  rename(Tissue_coverage_First = Tissue_coverage) %>%
    left_join(vast_tools_coverage_long_df, by=c("EventID" = "EVENT", "Second_Tis" = "Tissue")) %>%
  rename(Tissue_coverage_Second = Tissue_coverage) %>%
  #### Add label based on the PSI in Av_Others
  mutate(Specialized_pattern = ifelse(Av_Others <= 15, "YES_SPEC", "NO_SPEC")) %>%
  #### Modify associated tissue based on covarage information
  mutate(Associated_Tis_coverage = case_when(Associated_Tis == Tissue & Tissue_coverage_First == "YES_COV" ~ Associated_Tis,
                                             Associated_Tis == Tissue & Tissue_coverage_First == "NO_COV" ~ "NO_tissue_cov",
                                             grepl(";", Associated_Tis) & Tissue_coverage_First == "YES_COV" & Tissue_coverage_Second == "YES_COV" ~ Associated_Tis,
                                             grepl(";", Associated_Tis) & Tissue_coverage_First == "NO_COV" & Tissue_coverage_Second == "NO_COV" ~ "NO_tissue_cov",
                                             grepl(";", Associated_Tis) & Tissue_coverage_First == "NO_COV" & Tissue_coverage_Second == "YES_COV" ~ Second_Tis,
                                             grepl(";", Associated_Tis) & Tissue_coverage_First == "YES_COV" & Tissue_coverage_Second == "NO_COV" ~ Tissue,
                                             TRUE ~ Associated_Tis)) %>% #These are for the NO_tissue cases (this does not respect the dPSI).
  #### Add Taus
  left_join(taus_df, by=c("EventID"="EVENT")) %>%
  #### Add gene orthogroups
  left_join(exon_to_gene_df %>% distinct(), by="EventID") %>%
  left_join(gene_orthogroups_df, by="GeneID") %>%
  dplyr::select(-Species)

#########################################
########## SAVE TO OUTPUT ###############
#########################################
write.table(ts_exons_extra_info_df, output_file, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
