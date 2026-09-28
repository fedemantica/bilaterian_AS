########################################################################
########### Script to filter TS exons by minimum coverage ##############
########################################################################

##### Upload libraries
.libPaths("/users/mirimia/fmantica/R/x86_64-pc-linux-gnu-library/4.2")

library(argparser, lib.loc="/users/mirimia/fmantica/R/x86_64-pc-linux-gnu-library/4.2")
library(tidyverse, lib.loc="/users/mirimia/fmantica/R/x86_64-pc-linux-gnu-library/4.2")

##### Define arguments
#Create a parser
p = arg_parser("Script to define TS and PATAS exon events based on coverage, specialization profile and Tau")

#Add command line arguments
p = add_argument(p, arg="--input", help="File with the modified output of get_tisAS.pl, which contains the TS exons to filter associated with the extra info")
p = add_argument(p, arg="--output_all", help="Output file with all classes of TS exons")
p = add_argument(p, arg="--output_ts", help="Output file with selected TS exons")
p = add_argument(p, arg="--output_patas", help="Output file with PATAS exons")
p = add_argument(p, arg="--output_discarded", help="Output file with Discarded exons")
p = add_argument(p, arg="--tau_cutoff", help="Minimum Tau required for tissue-specificity", type="numeric")

#Parse the command line arguments
argv = parse_args(p)

#Define arguments
input_ts_exons_file = argv$input
output_all_file = argv$output_all
output_ts_file = argv$output_ts
output_patas_file = argv$output_patas
output_discarded_file = argv$output_discarded
tau_cutoff = argv$tau_cutoff

#########################################
########## MAIN #########################
#########################################

#### Read in files
input_ts_exons_df = read.delim(input_ts_exons_file, header=TRUE)

all_final_exons_df = input_ts_exons_df %>%
  mutate(Exon_type = case_when(Associated_Tis_coverage != "NO_tissue" & Associated_Tis_coverage != "NO_tissue_cov" & Specialized_pattern=="YES_SPEC" & Tau >= tau_cutoff ~ "TS_AS",
                               Associated_Tis_coverage != "NO_tissue" & Associated_Tis_coverage != "NO_tissue_cov" ~ "PATAS",
                               TRUE ~ "Discarded"))

###############################
#### Select TS exons ##########
###############################
#### TS exons are defined are those exons that have:
## Associated_Tis_coverage != "NO_tissue" and "NO_tissue_cov"
## Specialized_pattern == "YES_SPEC"
## Tau ≥ Tau_cutoff

ts_exons_df = all_final_exons_df %>%
  filter(Exon_type == "TS_AS")

###############################
#### Select PATAS exons #######
###############################
#### PATAS exons are defined as those that have:
## Associated_Tis_coverage != "NO_tissue" and "NO_tissue_cov"

patas_exons_df = all_final_exons_df %>%
  filter(Exon_type == "PATAS")

###############################
#### Select DISCARDED exons ###
###############################
#### Discarded exons are all the ones that are neither TS nor PATAS. This is just to keep track of what we are losing
discarded_exons_df = all_final_exons_df %>%
  filter(Exon_type == "Discarded")

#########################################
########## SAVE TO OUTPUT ###############
#########################################
write.table(all_final_exons_df, output_all_file, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(ts_exons_df, output_ts_file, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(patas_exons_df, output_patas_file, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(discarded_exons_df, output_discarded_file, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
