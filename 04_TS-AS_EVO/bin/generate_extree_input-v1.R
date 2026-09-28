#####################################################
########### UPDATE LIBRARY PATH #####################
#####################################################

#.libPaths(c("/software/mi/Rlib3.6/", "/software/as/el7.2/EasyBuild/CRG/software/R/3.6.0-foss-2019a/lib64/R/library", .libPaths()))

.libPaths(c("/users/mirimia/fmantica/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

#Upload libraries
library(tidyverse)
library(dplyr)

packageVersion("dplyr")

#####################################################
########### READ EXTERNAL ARGUMENTS #################
#####################################################
Args = commandArgs(trailingOnly=TRUE);
exon_orthogroups_file = Args[1]
gene_orthogroups_file = Args[2]
output_file = Args[3]
ordered_tissues = unlist(strsplit(Args[4], ","))
ordered_species = unlist(strsplit(Args[5], ","))
ts_call_dir = Args[6]
dPSI_prefix = Args[7]
all_reference_annot_dir = Args[8]


#####################################################
#################### MAIN ###########################
#####################################################

### Upload exon orthogroups
#Header: ExCluster_ID	GeneID	Coordinate	Species	Membership_score	VastID	Overlap_groupID	OG_ID
exon_orthogroups_df = read.delim(exon_orthogroups_file, header=TRUE, col.names=c("ExOG_ID", "GeneID",	"Coordinate",	"Species",	"Membership_score",	"VastID",	"Overlap_groupID",	"OG_ID")) %>%
  ### Add column with info
  mutate(info = "info")

### Upload gene orthogroups with gene names
gene_orthogroups_gene_names_df = read.delim(gene_orthogroups_file, header=FALSE, col.names=c("OG_ID", "Species", "GeneID", "GeneName"))

#### Get PSI tables
all_species_PSI_df = data.frame()
for (species in ordered_species) {
  species_PSI_df = data.frame() 
  for (tissue in ordered_tissues) {
    ts_as_file = paste0(ts_call_dir, "/", species, "/", dPSI_prefix, tissue, ".tab")
    ts_df = read.delim(ts_as_file, header=TRUE, sep="\t")
    ### In some cases these files are empty
    if (nrow(ts_df) == 0) {
      species_PSI_df = species_PSI_df %>%
        mutate(!!paste0("Av_", tissue) := NA)
    } else {
      ### Select column with PSI info
      filtered_ts_df = ts_df %>%
        dplyr::select(EventID, LENGTH, paste0("Av_", tissue))
      ### Join to species dataframe
      if (nrow(species_PSI_df) == 0) {species_PSI_df = filtered_ts_df} else {
        species_PSI_df = species_PSI_df %>%
          left_join(filtered_ts_df, by=c("EventID", "LENGTH"))
        }
      }
    }
  species_PSI_df = species_PSI_df %>%
    mutate(Species = species)
  all_species_PSI_df = all_species_PSI_df %>%
    dplyr::bind_rows(species_PSI_df)
}

### Add correct length based REFERENCE tables
all_species_annot_df = data.frame()
for (species in ordered_species) {
  ### Upload file with reference annotations
  annot_file = paste0(all_reference_annot_dir, "/REFERENCE-ALL_ANNOT-", species, ".tab")
  annot_df = read.delim(annot_file, header=TRUE) %>%
    ### Filter for the actual length
    dplyr::select(EVENT, LE_n) %>% 
    mutate(Species = species)
  ### Join to final file
  all_species_annot_df = all_species_annot_df %>%
    dplyr::bind_rows(annot_df)
}

all_species_annot_df = all_species_annot_df %>%
  dplyr::rename("EventID"="EVENT")

### Replace info in the final table
final_all_species_PSI_df = all_species_PSI_df %>%
  left_join(all_species_annot_df, by=c("Species", "EventID")) %>%
  dplyr::select(-LENGTH) %>%
  dplyr::rename("LENGTH"="LE_n") %>%
  dplyr::select(colnames(all_species_PSI_df))


### Modify colnames
colnames(final_all_species_PSI_df)[grep("Av_", colnames(final_all_species_PSI_df))] = paste0(colnames(final_all_species_PSI_df)[grep("Av_", colnames(final_all_species_PSI_df))], "_PSI")
colnames(final_all_species_PSI_df) = gsub("Av_", "", colnames(final_all_species_PSI_df))
colnames(final_all_species_PSI_df)[colnames(final_all_species_PSI_df) == "LENGTH"] = "exon_length"


################################
## INTEGRATE INFO WITH EX OGS ##
################################
final_exon_orthogroups_df = exon_orthogroups_df %>%
  ### Add PSI
  left_join(final_all_species_PSI_df, by=c("VastID"="EventID", "Species"="Species")) %>%
  ### Add gene names
  left_join(gene_orthogroups_gene_names_df, by=c("OG_ID", "Species", "GeneID"))


### Save to file
write.table(final_exon_orthogroups_df, file=output_file, row.names = FALSE, col.names=TRUE, sep="\t", quote=FALSE)
