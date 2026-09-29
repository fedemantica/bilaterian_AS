#!/usr/bin/perl


# Arguments
$sp = $ARGV[0];
$batch = $ARGV[1];
$inDir = $ARGV[2]; # path to input files: IUPRED_PARTS 
$outDir = $ARGV[3];
$cutoff=0.5; # to call it DISORDER or NOT


# Load all iupred files of the batch
@files=glob($inDir."/".$sp."_".$batch."/*.iupred");

# Creates a file Hs2_batch-IUPRED.txt in "IUPRED_batches" folder
system "mkdir $outDir" unless (-e "$outDir");
open (O, ">".$outDir."/".$sp."_".$batch."_parsed.txt");
foreach $f (@files){
    ($prot)=$f=~/.+\/(.+)\./;
    open (I, $f);
    $processed++;

    while (<I>){
	chomp;
	@t=split(/\t/);
	if (/\#/){}
	else {
	    $pos=$t[0];
	    $aa=$t[1];
	    print "ISSUE: $prot\t$pos\t$aa\t$dis\n" if $t[2]!~/\d/;
	    $rawDis=$t[2];
	    $dis=1 if $t[2] >= $cutoff;
	    $dis=0 if $t[2] < $cutoff;
	    
	    if ($aa=~/[A-Z]/){
		print O "$prot\t$pos\t$aa\t$dis\t$rawDis\n";
	    }
	}
    }
    close I;
}

print "Procesed files: $processed\n";


