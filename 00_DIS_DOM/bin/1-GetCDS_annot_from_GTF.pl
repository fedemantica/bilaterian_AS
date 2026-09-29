#!/usr/bin/perl

# Script edited by Irene (Sep 2024)
# provided a GTF and gDNA.fasta creates a BioMart-like file with annot CDS (then translate.pl)
# 26/11/16: The CDS column 8 (t[7]) is the offset in the translation
#           if 1 => it starts in phase 2; if 2 => starts in phase 1

die "\nUsage: GetCDS_annot_from_GTF.pl GTF gDNA.fasta [or only Sp]\n\n" if !$ARGV[0];

#if ($ARGV[2] =~ /(Spu|Sp2|Bla)/) {
#    $stop_included = 1;
#} else {
#    $stop_included = 0;
#}

#$stop_included = 0; # 1 = Spu, Bla; 0 = Ensembl
#print "\nMake sure that the 'shift' in the seq is right, dependending on whether or not the STOP is encoded in the CDS coord.";
#print "\n(now shift is: $stop_included because the specie is $ARGV[2]).\n";
#print "I think this script is NOT prepared to handle a 'shift' different than 0, revise it if necessary!!.\n\n";

#sleep(2);  # wait a couple of seconds (to read the print?)


# Output directory
$outdir = $ARGV[3];
$sp = $ARGV[2];
$LOGfile=$ARGV[4];
#$LOGfile="$outdir"."/"."$sp"."_LOG.txt";

#$date = localtime();  # Get current date and time
open (LOG, ">>$LOGfile");
#print LOG "Input files: \n\t$ARGV[0]\n\t$ARGV[1]\n";
#print LOG "The stop shift is $stop_included because the specie is $ARGV[2].\n";
#print LOG "I think this script is NOT prepared to handle a 'shift' different than 0, revise it if necessary!!.\n\n";


# Open GTF and Fasta files of the specie of interest

$gtf_file= $ARGV[0];

    if ($gtf_file =~ /.gz$/) {
	open(GTF, "gunzip -c $gtf_file |") || die "It cannot open pipe to $gtf_file";
    }else {
	open(GTF, $gtf_file) || die "It cannot open $gtf_file";
    }
#open (GTF, $ARGV[0]);
open (DNA, $ARGV[1]);
($sp)=$ARGV[1]=~/.+\/(.+?)_gDNA/;
die "Can't identify Sp in $ARGV[1]\n" if !$sp;

# OLD CODE, not used in SNAKEMAKE (Irene)
#$dir="/users/mirimia/mirimia/XPIPE";

#$pre_version=""; # if only provides a species, changes the name
#if ($ARGV[1]){
#    open (GTF, $ARGV[0]);
#    open (DNA, $ARGV[1]);
#    ($sp)=$ARGV[1]=~/.+\/(.+?)_gDNA/;
#    die "Can't identify Sp in $ARGV[1]\n" if !$sp;
#}
#elsif (!$ARGV[1] && length($ARGV[0])==3){
#    $sp=$ARGV[0];
#    open (GTF, "$dir/GTF/$sp"."_annot.gtf") || die "Can't open GTF";
#    open (DNA, "$dir/GENOMES/$sp"."_gDNA.fasta") || die "Can't open gDNA";
#    $pre_version=1;
#}

# first, load sequences
while (<DNA>){
    chomp;
    if (/\>/){
	($chr)=/\>(.+)/;
    }
    else {
#	print "Parsing $chr\n";
	$_=~tr/a-z/A-Z/;
	$seq{$chr}.=$_ if $chr;
    }
    # Dictionary of sequence per chromosome (key=$chr)
}
close DNA;

print "\nParsing GTF\n";
while (<GTF>){
    chomp;
    @t=split(/\t/);
    # splits the columns in an array
    ($g)=/gene_id \"(.+?)\"/;
    # looks for "gene_id" and stores content
    ($tr)=/transcript_id \"(.+?)\"/;
    # content of transcript_id
    ($ex)=/exon_number \"(.+?)\"/;
    # content of exon number
    $tr_g{$tr}=$g;
    $chr{$tr}=$t[0];
    $str{$tr}="1" if $t[6] eq "+";
    $str{$tr}="-1" if $t[6] eq "-";
    # These are dictionaries. Each of them store different info in keys that are the "tr", transcripts
    # values of the dictionaries are the gene id (tr_g), chromosome, and strand

	# If the "feature" (3rd col) is "exon"
    if ($t[2] eq "exon"){
	print LOG "*** Parsing GTF: No transcript detected\n" if !$tr;  # there is no transcript_id
	print LOG  "*** Parsing GTF: No exon number detected for $tr\n" if !$ex;  # there is no exon number
	#Dictionaries for positions of the exons:
	$exon_ini{$tr}{$ex}=$t[3];	# key: transcript_id - value: key:exon_number - value:  "start"
	$exon_end{$tr}{$ex}=$t[4];	# key: transcript_id - value: key:exon_number - value:  "end"
	$tot_lines_ex++;
    }
    # If the "feature" is "CDS"
    elsif ($t[2] eq "CDS"){
	($prot)=/protein_id \"(.+?)\"/;
	print LOG "*** Parsging GTF: No protein_ID detected\n" if !$prot;
	$prot=$tr if !$prot;
	print LOG "*** Parsing GTF: No exon number detected for $tr,$prot\n" if !$ex;
	$CDS_ini{$tr}{$ex}=$t[3];  # position "start" of CDS (key: transcript- value=key:exon-value="start")
	$CDS_end{$tr}{$ex}=$t[4]; # position "end" of CDS (key: transcript- value=key:exon-value="end")
	$offset{$tr}{$ex}=$t[7];	# value of "frame" (key: transcript- value=key:exon-value="frame")
	$tr_prot{$tr}=$prot;	# key: transcript_id -  value: protein_id
	$tot_lines_CDS++;
    }
    
    # If the "feature" is "CDS"
    elsif ($t[2] eq "stop_codon"){
	$stop{$tr}=$t[3];	# position "start" of the STOP codon in the transcript
	$stop_ex{$tr}=$ex;	# exon that corresponds to the transcript of the annotated STOP codon
	$tot_lines_STOP++;
    }
}
close GTF;

while ($tot_lines_ex=~s/(\d+)(\d{3})/\1,\2/){}
while ($tot_lines_CDS=~s/(\d+)(\d{3})/\1,\2/){}
while ($tot_lines_STOP=~s/(\d+)(\d{3})/\1,\2/){}

print "Parsed: Exons ($tot_lines_ex); CDS ($tot_lines_CDS); STOP ($tot_lines_STOP)\n\n";
print LOG "\nParsed: Exons ($tot_lines_ex); CDS ($tot_lines_CDS); STOP ($tot_lines_STOP)\n";
print "Creating sequence files\n";

#### OLD CODE, now there is no "pre_version" (Irene)
#if ($pre_version){
#    $output_file="$sp"."_CDS_annot-Vpre.fasta";
#    $output_file_off="$sp"."_CDS_annot-Vpre-off.fasta"; # 26/11/16: to correct for offset
#    $output_file2="$sp"."_prot_annot-Vpre.fasta";
#    $pro="$sp"."_CDS_annot-Vpre-off.pro";
#}
#else {
#    $output_file="$outdir"."/"."$sp"."_CDS_annot-vB.fasta";
#    $output_file_off="$outdir"."/"."$sp"."_CDS_annot-vB-off.fasta"; # 26/11/16: to correct for offset
#    #$output_file2="$sp"."_prot_annot-vB.fasta";
#    #$pro="$sp"."_CDS_annot-vB-off.pro";
#}

#######

# Creating output files
$output_file="$outdir"."/"."$sp"."_CDS_annot.fasta";
$output_file_off="$outdir"."/"."$sp"."_CDS_annot-off.fasta"; 

open (O, ">$output_file");
open (O_OFF, ">$output_file_off");

######


# only looping through coding transcripts
# Uses the protein_ID-transcript_ID dictionary
# sorts all keys ($tr) and iterates them 
foreach $tr (sort keys %tr_prot){
    $exon_starts="";
    $exon_stops="";
    $CDS_starts="";
    $CDS_stops="";
    $CDS_end=0;
    $CDS_ini=0;
    $CDS_seq="";
    $offset="";

	# goes to the $CDS_ini dictionary and gets the keys that correspond to each $tr
	# These "keys" are the exons (exon number) --> now in $ex 
	# The value is the start position 
	# Dictionary: (key: transcript- value:(key:exon-value="start"))
	# The sort thing is to order the keys ($ex) numericalle, from less to more
    foreach $ex (sort {$a<=>$b} (keys %{$CDS_ini{$tr}})){
	
	# For first exon, it starts in the CDS of the first exon.
	# Gets start position of the CDS:
	$base_CDS=$CDS_ini{$tr}{$ex} if (!$exon_starts); # sets co = 1 for CDS
	# Gets the offset:
	$offset=$offset{$tr}{$ex} if (!$exon_starts); # the first exon of CDS
	
	# Goes to the exon start/stop positions, and stores the values (as a string $) adding a ";",
	# it appends the positions of each exon ($ex) when iterating them (foreach $ex)
	$exon_starts.="$exon_ini{$tr}{$ex};";
	$exon_stops.="$exon_end{$tr}{$ex};";
	
	$le_CDS_ex=$CDS_end{$tr}{$ex}-$CDS_ini{$tr}{$ex}+1;  # length of the CDS of the exon: CDS_end - CDS_in
	$CDS_ini=$CDS_end+1;  # Now the start of the CDS is 1 (or just 1 after the end of the last CDS analyzed??)
	$CDS_end=$CDS_ini+$le_CDS_ex-1;  # End of the CDS: CDS_ini + length of exon (-1)
	# If there is a $stop for this transcript 
	# and the exon of the iteration ($ex) is the same where there is a stop ($stop_ex{$tr}=$ex)
	# And there is no "$stop_included" (or $stop_included=0)
	$CDS_end+=3 if $stop{$tr} && $stop_ex{$tr} eq $ex && !$stop_included;
	# CDS_ends advances 3, so incorportates a STOP codon
	# But the way the sequence is obtained is NOT with CDS_end, so this is just about getting the position
	# of the end of the exon right (if contains stop, advance 3 so "CDS_init" starts after the exon)

	$CDS_starts.="$CDS_ini;"; # Gets CDS initial position (it's +1 the end)and appends it 
	$CDS_stops.="$CDS_end;"; # Gets CDS end position and appends it 
	
	# If the strand of the transcript is +
	if ($str{$tr} eq "1"){
		# gets the sequence of the chromosome where the transcript is:
		# $seq{$chr{$tr}}: sequence of all chromosome
		# From the start position (dictionary): $CDS_ini{$tr}{$ex}-1 (-1 because now 0 is 1)
		# and advances the whole length: $le_CDS_ex
	    $seq_bit=substr($seq{$chr{$tr}},$CDS_ini{$tr}{$ex}-1,$le_CDS_ex);
	    $CDS_seq.=$seq_bit;
	}
	# strand is -:
	elsif ($str{$tr} eq "-1"){
		# gets sequence
	    $seq_bit=substr($seq{$chr{$tr}},$CDS_ini{$tr}{$ex}-1,$le_CDS_ex);
	    # reverses sequence (order of letters and then A->T etc)
	    # so obtains the complementary chain, inverted. So the proper sequence
	    $seq_bit=join("", reverse split (//, $seq_bit));
	    $seq_bit=~tr/ACGTacgt/TGCAtgca/;
	    $CDS_seq.=$seq_bit;	    
	}
    }  # We stop iterating the exons!!!!
    
    # If there is a stop codon for the transcript
    # and there is no "$stop_included"
    if ($stop{$tr} && !$stop_included){
	$stop_codon=substr($seq{$chr{$tr}},$stop{$tr}-1,3);
	# stop_codon sequence: seq of chromosome of the $tr, go to the position, 3 letters
	# But if the strand is "-", reverse it
	if ($str{$tr} eq "-1"){
	    $stop_codon=join("", reverse split (//, $stop_codon));
            $stop_codon=~tr/ACGTacgt/TGCAtgca/;
	}
	# Add the STOP codon to the sequence
	$CDS_seq.=$stop_codon if !$stop_included;	
    }
	
	# chop removes last character (the last ";" added)
    chop($exon_starts);
    chop($exon_stops);
    chop($CDS_starts);
    chop($CDS_stops);
    
    print O ">$tr_g{$tr}|$tr_prot{$tr}|$exon_starts|$exon_stops|$CDS_starts|$CDS_stops|$str{$tr}\n$CDS_seq\n";

    ### added to correct for wrong offsets (i.e. not start codon)
    $CDS_seq_off="NN$CDS_seq" if $offset==1;
    $CDS_seq_off="N$CDS_seq" if $offset==2;
    $CDS_seq_off=$CDS_seq if $offset==0;
    
    print O_OFF ">$tr_g{$tr}|$tr_prot{$tr}|$exon_starts|$exon_stops|$CDS_starts|$CDS_stops|$str{$tr}\n$CDS_seq_off\n";
}
close O;
close O_OFF;

# 30/08/2024-Irene: I skip these orders to run them separately in snakemake
# print "\nBy now, I should have the files: $output_file and $output_file_off\n and the plan is to get $pro and $output_file2.";
#print "\nNow the CDS sequences will be translated\n";  # Irene's commenting

#print "Translating CDS sequences\n";
#system "translate.pl $output_file_off";

#system "mv $pro $output_file2";
#system "rm $output_file_off";
#print "Done\n\n";
close LOG;
