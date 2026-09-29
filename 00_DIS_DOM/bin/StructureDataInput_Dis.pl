use English;
# POST and PREMATCH

# expects: ARGV[0]: input, ARGV[1]: sp ,ARGV[2]: outdir

die "\nNeeds fasta of protein sequences (parsed)\n\n" if !$ARGV[0];

$sp = $ARGV[1];
$dir=$ARGV[2]."/".$sp;
$outDir= $dir."/DISDOM";
$PartsDir = $outDir."/PARTS";
$batches = $ARGV[3];

my $date = localtime();  # Get current date and time
$LOGfile=$ARGV[4];



# Open Reference File
#print LOG ">".$outDir.$sp."-locationFiles.tab";

system "mkdir $outDir" unless (-e "$outDir");
system "mkdir $PartsDir" unless (-e "$PartsDir");
#system "mkdir $outDir/LOG_IUPRED" unless (-e "$outDir/LOG_IUPRED");
#open (ProtPlace, ">".$outDir."/LOG_IUPRED/".$sp."-locationFiles.tab");
$LocationFile=$outDir."/".$sp."-locationFiles.tab";
open (ProtPlace, ">$LocationFile");

# Open file with sequences
open (I, $ARGV[0]) || die "Cannot open the input\n"; 

$/="\>";
<I>;

while (<I>){
    /\n/;
    $name=$PREMATCH;
    $seq=$POSTMATCH;
    ($pr)=$name=~/\|([^\|]+)\|/; 
    
    next if $name!~/\;/;
    
    
    $seq=~s/[\n\>]//g;
    
    
    #my $folder_part = int($a++ / $fXfolder) + 1;
    my $folder_part = ($a++ % $batches) + 1;
    my $folderName = $sp."_".$folder_part;
    my $folderPath = $PartsDir."/".$folderName;
    
    system "mkdir $folderPath" unless (-e "$folderPath");
    
    my $myFile = $PartsDir."/".$folderName."/".$pr.".fasta";
    #print LOG "\t$myFile\n";
    
    
    if (-e "$myFile"){}
    else {
		open (OUTS, ">".$myFile) or die "Cannot open file:\n $myFile\n";
		print OUTS ">$name\n$seq";
		close OUTS;
	    $z++;
	    }
    
    $RefFile{$pr}=$folderName;
}

for my $p (keys %RefFile) {
	$fold = $RefFile{$p};
	print ProtPlace "$p\t$fold\n";
	$b++;
}


close ProtPlace;

open (LOG, ">>$LOGfile");
#print LOG "\nProtein sequences have been split into individual fasta files and stored in $batches folders in:\n\t $PartsDir\n";
#print LOG "\nProtein sequences have been split into individual fasta files and stored in $batches folders in, a fasta file per protein.";
#print LOG "\nThe location of individual protein fastas in:\n\t $LocationFile \nInput file with protein fasta sequences:\n\t $ARGV[0]\n";
#print LOG "\nInput file (protein fasta sequences):\n\t $ARGV[0]";
print LOG "\n Control log (multiexon proteins): $a\t$z\t$b\n"; # This was an internal control (I don't remember exactly of what, I want to see it (Irene Oct2024))
close LOG;

print "$a\t$z\t$b\n";

##

