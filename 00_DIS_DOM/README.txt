README
27/01/2025 - Irene Zapata Bódalo

Explanation of the Snakemake file to obtain the disorder rate and overlapping functional domains (PFAM) of EX events for (20) species.
* All needed scripts are in the folder: bin/
* Paths to all files and variables in: config.yaml
* Resources needed to launch the jobs for each rule in the cluster (Slurm) in: new_cluster.json
* To run the script: activate snakemake using a conda environment: snakemake-V5  
	/users/mirimia/izapata/.conda/envs/snakemake-V5

* Then, run the bash script: run_snakemake_newC.sh
	You need to specify (go) to the snakemake folder (cd path...). 
	This script contains the order to run the snakemake file, using the information in "new_cluster.json"
	to assign the required cluster resources for each rule. The language is a mix of old and new cluster
	managers systems (Univa Grid Engine, then Slurm). There is quite a lot of "unused" code in there.

! This snakemake script is ment to be run obtaining the final LOG file, generated with the last "CLEANUP" step. 
It also works if the target files are one of the 2 "CLEANUP" LOG files, one for the analysis of functional domains,
one for the disorder rates (see explanation of step (E)). Since these steps eliminate and compress all intermediate 
and final files (.gz), running intermediate steps after these compressions could be problematic: the script is not 
prepared to unzip (or even consider) compressed versions of the files. Thus, if you want to run the script to obtain 
intermediate files, it is better to select intermediate targets, without reaching the "CLEANUP" steps that 
compress and eliminate files.  If you want the Snakemake to detect some of the compressed files in the 
intermediate rules, decompress them. 

* Steps of the file:
For the explanations, "$sp" refers to the "species".
The data for each species is stored in a separate folder, with the name of the species. 

A) General data preparation:
1- Incorporate "extra exons" ("Fake Transcripts") into the annotaion file (GTF) of each species.
	- Rules:
		* Incorporate_FakeTranscripts1
		* Incorporate_FakeTranscripts2
	- Scripts:
		* A1_generate_annotations_copy.pl

2- Obtain genetic (DNA, CDS) and protein sequences. 
First, Information of each protein/CDS in the header relative to its
structure (exons, introns, coordinades, strand...) followed by the sequence in fasta format.
Second, parsed file with information per exon (geneID, proteinID, sequence, position...).

	- Rules:
		* GetCDS_annot_from_GTF
		* CDS_prot_parse
	- Scripts:
		* 1-GetCDS_annot_from_GTF.pl
		* translate.pl
		* 2-Parse_CDS_UTR.pl
		* 35-Parse_CDS_prot_annot.pl
		
Input files:
	> REFERENCE-ALL_ANNOT-$sp.tab
	> $sp.Event-Gene.IDs.txt
	> Annotations GTF
	> Genome sequence (gDNA) FASTA
		
Main outputs:		
	< Annotations GTF with extra annotations ($sp_annot_fake.gtf.gz)
	< $sp_CDS_annot.fasta
	< $sp_prot_annot.fasta
	< $sp_CDS_annot_parsed.tab
	< $sp_prot_annot_parsed.tab
	
B) Preparation of data for analysis.
Protein sequences are split in individual files and stored in folders to be analyzed by batches.
Intermediate step to optimize the analysis and parallelized it.
	- Rules:
		* FileStruct_Parts
		* PFAMinput
		
	- Scritps:
		* StructureDataInput_Dis.pl

C) PFAM domains: analysis of protein sequences (multiexon).  
1- Obtention of PFAM domains: name/id/clan of the domain detected in each exon 
(with information of the gene, protein, position/coordinades of the exon, kind of overlap...)

	- Rules: 
		* Run_PFAM_By_Batch
		* Combine_Parsed_PFAM

	- Scripts/tools:
		* HMMER module: module load HMMER/3.3.2-gompi-2022a
		* pfam_scan.py
		* pfam libraries (in bin/pfam_scan/pfam_dir)
		* 3-Parse_parsed_pfamout_copy.pl

2- Finally, "REFERENCE" file with PFAM domains, as in VastDB.
	- Rules: 
		* RefDomTab_VastDB
		
	- Scripts: 
		* 7-PutDomAndEv_C1AC2-v2_IZB.pl
	
Input files:
	> fasta files of protein sequences (generated and organized in step B)

Main output files:
	< $sp_pfam_parsed_sorted.tab: detected domains, information relative to each CD
	< $sp_prot_annot-PFAM.tab: all exons, genomic/protein coordinades, information on the overlap
	< REFERENCE-ALL_ANNOT-$sp-PFAM-v2.tab



D) Disorder rate of exons: analysis of protein sequences (multiexon).
1- Obtention of the estimate probability of being disordered per residue (iupred) and parsing
to obtain average disorder rate of EX exons (A) and their up/downstream sequences (C1, C2),
computed with and without the binarization of the raw probability estimate (cutoff).
	- Rules:
		* Run_Iupred_By_Batch
		* Combine_Parsed_Iupred
		
	- Scripts:
		* iupred2a.py
		* Parse_IUPRED2_rec_NoCutOff.pl: 

2- Finally, "REFERENCE" file with average disorder rates of EX exons, as in VastDB. 
Computed with (previous versions of VastDB) and without cutoff (new alternative). 
	- Rules:
		* RefDisTab_VastDB
		
	- Scripts:
		* 6-CalculateAndPutDis_C1AC2-v2_IZB.pl
		
Input files:
	> fasta files of protein sequences (generated and organized in step B)

Main output files:
	< $sp-iupred_result.txt: estimate per residue (with and without binarization with cutoff)
	< REFERENCE-ALL_ANNOT-$sp-DisI-v2.tab
	
E) Clean-up: elimination of unnecessary intermediate files, compression of the rest of files.
Also, all "LOG" files are combined into a single LOG file, summarizing the analysis and the generated docs.
	- Rules:
		* Clean_Domain_Analysis
		* Clean_Disorder_Analysis
		* Clean_Final
		
Input files:
	> Uncompressed Reference files
	> Intermediate uncompressed files (to compress or eliminate)
	> Intermediate LOG files

Main output files:
	< LOG-$sp_DISDOM.txt
	< Compressed versions of the files









	
