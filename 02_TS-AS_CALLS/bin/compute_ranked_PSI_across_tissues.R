###################################################################################
##### Script to compute the Tau of exon events based on the PSI by tissue #########
###################################################################################

##### Upload libraries
library(argparser)
library(tidyverse)

##### Define arguments
#Create parser
p = arg_parser("Compute PSIs and PSI ranks across tissues")

#Add command line arguments
p = add_argument(p, arg="inclusion_table", help="Inclusion table for the species of interest")
p = add_argument(p, arg="sample_tissue", help="File including the sample - tissue correspondence for all species of interest")
p = add_argument(p, arg="output", help="Output file")

#Parse the command line arguments
argv = parse_args(p)

#Define arguments
inclusion_table_file = argv$inclusion_table
sample_tissue_file = argv$sample_tissue
output_file = argv$output

#########################################
########## FUNCTIONS ####################
#########################################

##### Define a function to compute tau
compute_psi_tau = function(psi_vector) {
  tau = sum(1-psi_vector/max(psi_vector, na.rm=TRUE), na.rm=TRUE)/(length(psi_vector)-1)
  return(tau)
}

#########################################
########## MAIN #########################
#########################################

#### Read in files
vast_tools_df = read.delim(inclusion_table_file, header=TRUE)
sample_tissue_correspondence_df = read.delim(sample_tissue_file, header=TRUE)
  
#### Get sample tissue correspondence data.frame
sample_tissue_dict_df = sample_tissue_correspondence_df %>%
  dplyr::select(Group, Tissue) %>%
  distinct() %>%
  rename(Sample = Group)
  
#### Set variables
#Samples
species_samples = colnames(vast_tools_df)[grep("\\.Q", colnames(vast_tools_df), invert=TRUE)]
species_samples = species_samples[!(species_samples %in% c("GENE", "EVENT", "COORD", "LENGTH", "FullCO", "COMPLEX"))]
tot_species_samples = length(species_samples)
#Quality of samples
species_samples_quality = colnames(vast_tools_df)[grep("\\.Q", colnames(vast_tools_df), invert=FALSE)]
  
vast_tools_quality_long_df = vast_tools_df %>%
  dplyr::select("EVENT", species_samples_quality) %>%
  ### Select only cassette exon events
  filter(grepl("EX", EVENT)) %>% 
  pivot_longer(cols=species_samples_quality, names_to = "Sample", values_to = "Quality_scores") %>%
  mutate(Sample = sub("\\.Q", "", Sample),
         First_quality_score = sub(",.*", "", Quality_scores)) %>%
  filter(First_quality_score != "N") %>%
  unite(EVENT_Sample, EVENT, Sample)

#### Select only the event and the samples that have coverage (no N)
EVENT_Samples_with_coverage_vector = as.vector(vast_tools_quality_long_df$EVENT_Sample) 
  
exons_tau_df = vast_tools_df %>%
   dplyr::select("EVENT", species_samples) %>%
   filter(grepl("EX", EVENT)) %>% #select only cassette exon events
   pivot_longer(cols=species_samples, names_to = "Sample", values_to = "PSI") %>%
   filter(!(is.na(PSI))) %>%
   unite(EVENT_Sample, EVENT, Sample, remove=FALSE) %>%
   #### Select only the event and the samples that have coverage (no N)
   filter(EVENT_Sample %in% EVENT_Samples_with_coverage_vector) %>%
   left_join(sample_tissue_dict_df, by="Sample") %>%
   group_by(EVENT, Tissue) %>%
   ### Compute average PSI among samples from the same tissue
   summarise(PSI=mean(PSI)) %>% 
   ungroup() %>%
   group_by(EVENT) %>% mutate(Rank = rank(-PSI, ties.method = "first"))
 
###### Save to file
write.table(exons_tau_df, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
