#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import os

parser = argparse.ArgumentParser(description="Script to add extra information to the output of ExOrthist")
parser.add_argument("--ex_orthogroups", "-ex", required=True, metavar="ex_orthogroups", help="EX_clusters.tab file returned by an ExOrthist run")
parser.add_argument("--vastID_coords", "-vc", required=True, metavar="vastID_coords", help="File containing the correspondence bewteen vastID and exon coords (first and third column, respectively)", nargs="+")
parser.add_argument("--overlapping_groups", "-ov", required=True, metavar="output", help="overlapping_EXs_by_species.tab from the same ExOrthist run as the EX_clusters.tab")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

############################
##### DEFINE ARGUMENTS #####
############################

args = parser.parse_args()
ex_orthogroups_file = args.ex_orthogroups
vastID_coords_files = args.vastID_coords
overlapping_groups_file = args.overlapping_groups
output_file = args.output

############################
##### READ IN FILES ########
############################

#Header: ExCluster_ID, GeneID, Coordinate, Species, Membership_score
ex_orthogroups_df = pd.read_table(ex_orthogroups_file, sep="\t", index_col=None, header=0)

#Header: ExonID, GeneID, Exon_coords_A, C1_Ref, C2_Ref
all_species_vastID_coords_df = pd.DataFrame() 
for vastID_coords_file in vastID_coords_files:
  vastID_coords_df = pd.read_table(vastID_coords_file, sep="\t", index_col=None, header=0)
  vastID_coords_df["Species"] = re.sub("_.*", "", os.path.basename(vastID_coords_file))
  all_species_vastID_coords_df = pd.concat([all_species_vastID_coords_df, vastID_coords_df])

#Header: OV_EX_Aae_1	AAEL000001	chr3:102845303-102845712:+	Aae	10 (no real header)
overlapping_groups_df = pd.read_table(overlapping_groups_file, sep="\t", index_col=None, header=None, names=["Overlap_groupID", "GeneID", "Coordinate", "Species", "Number_hits"])

############################
##### FORMAT INPUTS ########
############################

#Add Species_Coordinate entries
ex_orthogroups_df["Coordinate"] = [re.sub(":\-", "", re.sub(":\+", "", element)) for element in list(ex_orthogroups_df["Coordinate"])]
ex_orthogroups_df["Species_Coordinate"] = [element[0]+"_"+element[1] for element in zip(list(ex_orthogroups_df["Species"]), list(ex_orthogroups_df["Coordinate"]))]

all_species_vastID_coords_df["Species_Coordinate"] = [element[0]+"_"+element[1] for element in zip(list(all_species_vastID_coords_df["Species"]), list(all_species_vastID_coords_df["Exon_coords_A"]))]

overlapping_groups_df["Coordinate"] = [re.sub(":\-", "", re.sub(":\+", "", element)) for element in list(overlapping_groups_df["Coordinate"])] #Remove strand from coordinate
overlapping_groups_df["Species_Coordinate"] = [element[0]+"_"+element[1] for element in zip(list(overlapping_groups_df["Species"]), list(overlapping_groups_df["Coordinate"]))]

############################
##### MAIN #################
############################

#### Generate a dictionary with key=Coordinate, value=vastID
SpeciesCoord_vastID_dict = pd.Series(all_species_vastID_coords_df.ExonID.values, index=all_species_vastID_coords_df.Species_Coordinate).to_dict()
SpeciesCoord_overlapping_group_dict = pd.Series(overlapping_groups_df.Overlap_groupID.values, index=overlapping_groups_df.Species_Coordinate).to_dict()
overlapping_groups_df["ExonID"] = overlapping_groups_df["Species_Coordinate"].map(SpeciesCoord_vastID_dict)

#### Generate a dictionary with key=Overlap_GroupID, value = list of associated vastIDs
overlapping_groups_df = overlapping_groups_df.dropna(subset=["ExonID"])
collapsedIDs_overlapping_groups_df = overlapping_groups_df.groupby(["Species", "Overlap_groupID"])["ExonID"].apply(list).reset_index(name="ExonID")
overlapping_group_vastID_dict = pd.Series(collapsedIDs_overlapping_groups_df.ExonID.values, index=collapsedIDs_overlapping_groups_df.Overlap_groupID)

#### Add overlapping groups and vastID to exon orthogroups df.
#If the coords have vastID, translate directly;
ex_orthogroups_df["VastID"] = ex_orthogroups_df["Species_Coordinate"].map(SpeciesCoord_vastID_dict)

#If not, and there is only one vastID associated with that overlap group, select that vastID.
ex_orthogroups_df["Overlap_groupID"] = ex_orthogroups_df["Species_Coordinate"].map(SpeciesCoord_overlapping_group_dict) 
ex_orthogroups_df["Overlap_group_VastID"] = ex_orthogroups_df["Overlap_groupID"].map(overlapping_group_vastID_dict)
ex_orthogroups_df["Category"] = ["NO_VASTID" if isinstance(element[0], float) and isinstance(element[1], list) else "VASTID" for element in zip(list(ex_orthogroups_df["VastID"]), list(ex_orthogroups_df["Overlap_group_VastID"]))]

#Isolate for later use
clear_vastID_df = ex_orthogroups_df.loc[ex_orthogroups_df["Category"]=="VASTID"]

unclear_vastID_df = ex_orthogroups_df.loc[ex_orthogroups_df["Category"]=="NO_VASTID"]
unclear_vastID_df["VastID_num"] = [len(element) for element in list(unclear_vastID_df["Overlap_group_VastID"])]
single_unclear_vastID_df = unclear_vastID_df.loc[unclear_vastID_df["VastID_num"]==1]
single_unclear_vastID_df["VastID"] = [element[0] for element in list(single_unclear_vastID_df["Overlap_group_VastID"])]

#If there are more vastIDs associated with the same overlap group, select the one with the highest overlap with the original exon.
multiple_unclear_vastID_df = unclear_vastID_df.loc[unclear_vastID_df["VastID_num"]>1] 

#For each exon in each overlapping group, get the exons with VastID which has the highest overlap
######## NB: This code is not efficient at all if there are many entries

SpeciesCoords_selected_exon_dict = {}
for entry in list(multiple_unclear_vastID_df["Species_Coordinate"]):
  #### Select the Species and the Exon coordinates
  exon_coords = list(multiple_unclear_vastID_df.loc[multiple_unclear_vastID_df["Species_Coordinate"]==entry]["Coordinate"])[0]
  exon_coords_start = re.sub("-.*", "", re.sub(".*:", "", exon_coords))
  exon_coords_stop = re.sub(".*-", "", exon_coords)
  exon_range = range(int(exon_coords_start), int(exon_coords_stop))
  #### Select the coordinates of the exons to be compared
  other_exons = list(multiple_unclear_vastID_df.loc[multiple_unclear_vastID_df["Species_Coordinate"]==entry]["Overlap_group_VastID"])[0]
  all_other_exons_overlaps = []
  for other_exon in other_exons:
    other_exon_coords = list(all_species_vastID_coords_df.loc[all_species_vastID_coords_df["ExonID"]==other_exon]["Exon_coords_A"])[0]
    other_exon_start = re.sub("-.*", "", re.sub(".*:", "", other_exon_coords))
    other_exon_stop = re.sub(".*-", "", other_exon_coords)
    other_exon_range = range(int(other_exon_start), int(other_exon_stop))
    other_exon_overlap = len(range(max(exon_range[0], other_exon_range[0]), min(exon_range[-1], other_exon_range[-1])+1))
    all_other_exons_overlaps = all_other_exons_overlaps + [other_exon_overlap]
  ### Select the exon with the highest overlap
  max_overlap_exon = other_exons[all_other_exons_overlaps.index(max(all_other_exons_overlaps))]
  SpeciesCoords_selected_exon_dict[entry] = max_overlap_exon  

### Add the selected exon to the dataframe
multiple_unclear_vastID_df["VastID"] = multiple_unclear_vastID_df["Species_Coordinate"].map(SpeciesCoords_selected_exon_dict)

### Combine all dataframes and save to file
#Remove unused columns
single_unclear_vastID_df = single_unclear_vastID_df.drop(columns=["VastID_num"])
multiple_unclear_vastID_df = multiple_unclear_vastID_df.drop(columns=["VastID_num"])

final_df = pd.concat([clear_vastID_df, single_unclear_vastID_df, multiple_unclear_vastID_df])
final_df = final_df.drop(columns=["Category", "Species_Coordinate", "Overlap_group_VastID"])
final_df = final_df.sort_values(by=["ExCluster_ID", "Species", "Membership_score"])

#### Save to output file
final_df.to_csv(output_file, sep="\t", header=True, index=False, na_rep="NA")
