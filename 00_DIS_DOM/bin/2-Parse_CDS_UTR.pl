#!/usr/bin/perl

# Crucial file to get which exon are actually annotated as CDS.


# Irene's modification (Sep2024)
die "\nUsage: 2-Parse_CDS_UTR.pl GTF Sp version (without \"v\")\n\n" if !$ARGV[1];

$sp=$ARGV[1]; # new from 26/11/16
#$v=$ARGV[2]; # e.g. 84 or 80preB
#$v="V" if !$v; # so then can be renamed

my $outDir = $ARGV[2]."/".$sp."/";
die "Needs version code for output\n" if !$sp;

my $date = localtime();  # Get current date and time
$LOGfile=$ARGV[3];

open (LOG, ">>$LOGfile");
#print LOG "\n\t2-Parse_CDS_UTR.pl on $date, specie: $sp. \n";

#open (I, "Hsa_Gene-Tr-Prot.tab") || die "IDs\n";
#while (<I>){
#    chomp;
#    @t=split(/\t/);
#    $t_p{$t[1]}=$t[2];
#}

open (O, ">"."$outDir"."$sp"."_OKex.tab"); # only really done for ANNOT, so "pre" (14/11/16); not anymore 26/11/16 (all novel exons as GTF now)

# Adapted to open also compressed GTF
#open (IN, $ARGV[0]); #GTF from ensembl or novel
$gtf_file= $ARGV[0];
    if ($gtf_file =~ /.gz$/) {
	open(IN, "gunzip -c $gtf_file |") || die "It cannot open pipe to $gtf_file";
    }else {
	open(IN, $gtf_file) || die "It cannot open $gtf_file";
    }

while (<IN>){
    chomp;
    @t=split(/\t/);
    
    # Checks the ATTRIBUTE column
    # simplified just in case:
    ($g)=$t[8]=~/gene_id \"(.+?)\"/;
    ($tr)=$t[8]=~/transcript_id \"(.+?)\"/;
    ($n)=$t[8]=~/exon_number \"(\d+?)\"/;
    ($pr)=$t[8]=~/protein_id \"(.+?)\"/;
    
    # If "feature" is "exon", $ex: length exon
    if ($t[2] eq "exon"){
	$ex="$t[3]-$t[4]";  # writes it as a string, like: "1567-1584" 
    }
    # If "feature" is "CDS", and there is no "protein_id", $pr: transcript ID
    elsif ($t[2] eq "CDS"){
	$pr=$tr if !$pr;
	print O "$g\t$pr\t$ex\n"; # this prints the exon sequence when it's followed by a CDS entry (from which picks names, etc).
	$a++;
    }
}

print "Total number of CDS lines: $a\n";
print LOG "Total number of CDS lines: $a\n";
#print LOG "Output: $outDir $sp _OKex.tab\n";
