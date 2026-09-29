#!/usr/bin/perl
use English;

# Modified Irene Sep2024: combination of 2 of Manu's scripts:
# 3-Parse_CDS_annot.pl
# 5-Parse_prot_annot.pl

#($sp,$v)=$ARGV[0]=~/(.{3})\_CDS_annot-(.+?)\.fasta/;
# I think, for "Hs2_CDS_annot-v136.fasta":
# $sp=Hs2 and $v=v136

#$arg0 = $ARGV[0];
#$arg1 = $ARGV[1];
#$arg2 = $ARGV[2];
#$arg3 = $ARGV[3];
#$arg4 = $ARGV[4];

$sp=$ARGV[3];

my $outDir = $ARGV[4]."/".$sp."/";

my $date = localtime();  # Get current date and time
$LOGfile=$ARGV[5];

open (LOG, ">>$LOGfile");
#print LOG "\n\t3-Parse_CDS_annot.pl on $date, specie: $sp \n";
#print LOG "Arguments are:\n\t$arg0 \n\t$arg1 \n\t$arg2  \n\t$arg3  \n\t$arg4 \n";
#print LOG "Outdir: $outDir \n";


#die "\nUsage: 3-Parse_CDS_annot.pl Sp_CDS_annot-vVpre.fasta\n\n" if !$ARGV[0];

open (OK, $ARGV[0]) || die "Needs OK exons\n";
while (<OK>){
    chomp;
    @t=split(/\t/);
    # creates the OK dictionary key: prot_id - value:key: "ex_start-ex_end" - value: 1
    $OK{$t[1]}{$t[2]}=1;
    $pre++;
}
close OK;


# Open the CDS annotated fasta and gets name and sequence

open (O_CDS, ">"."$outDir"."$sp"."_CDS_annot_parsed.tab");
open (I_CDS, $ARGV[1]); #$sp_CDS_annot-$version.fasta
$/="\>";
<I_CDS>;
while (<I_CDS>){
    /\n/;
    $name=$PREMATCH;
    $seq=$POSTMATCH;
    $seq=~s/\n//g;
    # eliminate ">" and start sequence wiht "x"
    $seq=~s/\>//g;
    $seq="x$seq";
    # Separate fields of the name
    @t=split(/\|/,$name);
    $g=$t[0];  # gene ID
    $prot=$t[1];	# protein ID
    $str=$t[6];		# strand

	# If seqence is NOT "unavailable"
	# split it in an array
    if ($seq!~/unavailable/){
	@S=split("",$seq);
	
	# It takes the fields of positions, splits them (;) and stores in arrays
	# Then, it reoders "absolute exon positions" ("start" and "end") if the strand is "-"
	# !! I don't understand the reordering of start and end! I think its because now, if str is (-)
	# the first pair "start-end" is "more advanced" than the second (of the exons)
	# it's like Okazaki fragments!! (replication)
	if ($str==1){ 
	    @coI=sort{$a<=>$b}(split(/\;/,$t[2]));  # exon start 
	    @coF=sort{$a<=>$b}(split(/\;/,$t[3]));  # exon end
	    @exI=sort{$a<=>$b}(split(/\;/,$t[4]));  # exon start relative to CDS
	    @exF=sort{$a<=>$b}(split(/\;/,$t[5]));  # exon end relative to CDS
	}
	elsif ($str==-1){ 
	    @coI=sort{$b<=>$a}(split(/\;/,$t[2]));  # exon start --> REVERSE order if strand is ("-")
	    @coF=sort{$b<=>$a}(split(/\;/,$t[3]));  # exon end --> REVERSE order if strand is ("-")
	    @exI=sort{$a<=>$b}(split(/\;/,$t[4]));  # exon start relative to CDS
	    @exF=sort{$a<=>$b}(split(/\;/,$t[5]));  # exon end relative to CDS
	}
	
	$b=0;
	for $i (0..$#coI){ # the 5'UTR exon may also include coding (I is always UTR)
		# $#coI: last index of the array (length -1)
		# It iterates the $coI array from 0 to the last
	    $ex="$coI[$i]-$coF[$i]";  # string: "exon_start-exon_end"
	    
	    # Dictionary "OK" --> key: prot_iD - value:key:  "exon_start-exon_end" - value: 1
	    # If this exon (these coordinades) exist in the protein (CDS)
	    
	    if ($OK{$prot}{$ex}){
		$coIn[$b]=$coI[$i];  # now coIn[0] is that "exon start"
		$coFn[$b++]=$coF[$i];  # now coFn[0] is that "exon end". After this, $b increases!!!
	    } # but i thik all exons should be in the protein, no???
	    else {
			print LOG "An exon not in the OK dictionary! Is that expected?\n";
		} 
	} # so, for each exon that is annotated in the CDS, its coordinades are added to coIn/coFn variables
	
	for $i (0..$#exI){
		# $exI: values of the exon start relative to the CDS (original, from input file)
		
#	    $pI=sprintf("%.0f",$exI[$i]/3)+1 if ($exI[$i]+2)%3==0;
#	    $pI=sprintf("%.0f",$exI[$i]/3) if ($exI[$i])%3==0 || ($exI[$i]+1)%3==0;
#	    $pF=sprintf("%.0f",$exF[$i]/3)+1 if ($exF[$i]+2)%3==0;
#	    $pF=sprintf("%.0f",$exF[$i]/3) if $exF[$i]%3==0 || ($exF[$i]+1)%3==0;
	    
	    # @S is an array of the sequence of the protein
	    $exS=join("",@S[$exI[$i]..$exF[$i]]);  # sequence of the exon ($i exon)
	    
	    # I assume ALL exons where in the OK dictionary (previous if)
	    print O_CDS "$g\t$str\t$prot\t$i\t$exI[$i]-$exF[$i]\t$coIn[$i]-$coFn[$i]\t$exS\n";
	    $a++;
	}
    }
}

print "Parsing CDS. Total exons: $pre/$a\n";
print LOG "Parsing CDS. Total exons: $pre/$a\n";
#print LOG "Output: $outDir $sp _CDS_annot_parsed.tab\n";


##### Now process prot_annot! Old "5-Parse_prot_annot.pl" #####
my $date = localtime();  # Get current date and time
#print LOG "\n\t5-Parse_prot_annot on $date, specie: $sp \n";
open (O_prot, ">"."$outDir"."$sp"."_prot_annot_parsed.tab");
open (I_prot, $ARGV[2]); #$sp_prot_annot-$version.fasta
$/="\>";
<I_prot>;
while (<I_prot>){
    /\n/;
    $name=$PREMATCH;
    $seq=$POSTMATCH;
    $seq=~s/\n//g;
    $seq=~s/\>//g;

    $seq="x$seq";
    @t=split(/\|/,$name);
    $g=$t[0];
    $prot=$t[1];
    $str=$t[6];

    if ($seq!~/unavailable/){
	@S=split("",$seq);
	
	if ($str eq "1"){ 
	    @coI=sort{$a<=>$b}(split(/\;/,$t[2]));
	    @coF=sort{$a<=>$b}(split(/\;/,$t[3]));
	    @exI=sort{$a<=>$b}(split(/\;/,$t[4]));
	    @exF=sort{$a<=>$b}(split(/\;/,$t[5]));
	}
	elsif ($str eq "-1"){ 
	    @coI=sort{$b<=>$a}(split(/\;/,$t[2]));
	    @coF=sort{$b<=>$a}(split(/\;/,$t[3]));
	    @exI=sort{$a<=>$b}(split(/\;/,$t[4]));
	    @exF=sort{$a<=>$b}(split(/\;/,$t[5]));
	}

	#for $i (0..$#exI){
	 #   $exI[$i]=$exI[$i]+$Ns{$prot};
	#}
	#for $i (0..$#exF){
	#    $exF[$i]=$exF[$i]+$Ns{$prot};
	#}

#	$total=abs($exF[$#exF]-$exI[0])+1;
#	$test3n="3n" if $total%3==0;
#	$test3n="NO" if $total%3!=0;
#	print "$prot\t$total\t$test3n\n";

	$b=0;
	for $i (0..$#coI){ #This is done to remove all fully UTR exons
	    $ex="$coI[$i]-$coF[$i]";
	    if ($OK{$prot}{$ex}){
		$coIn[$b]=$coI[$i];
		$coFn[$b++]=$coF[$i];
	    }
	}

	for $i (0..$#exI){
		# GET the aminoacid number of the exon!!!! pI: initial, pF: final

		# Takes the exon position relative to the CDS (either the initial or the end)
		# Divides by 3 the position of the initiation codon. 
		# If its a "start" of a codon (like 1) --> 1+2/3 = 3 --> divisible by 3 (%3)
		# In this case, it divides by 3 and rounds: 0.33-->0 and adds 1
		# pI=1+ start position (4/3=1.3 -> 1+1 --> pI:2) (19/3=6.3 -> 6+1 --> pI: 7)
		
		# If it's the 2nd or 3rd of a codon (like 2 or 3): --> 3/3, 2+1/3 --> divisible by 3 
		# Then just divides by 3 and rounds
		# (5/3=1.6 --> pI=2)(6/3=2 --> pI=2) (2/3=0.6 --> pI=1)
	    $pI=sprintf("%.0f",$exI[$i]/3)+1 if ($exI[$i]+2)%3==0;
	    $pI=sprintf("%.0f",$exI[$i]/3) if ($exI[$i])%3==0 || ($exI[$i]+1)%3==0;
	    $pF=sprintf("%.0f",$exF[$i]/3)+1 if ($exF[$i]+2)%3==0;
	    $pF=sprintf("%.0f",$exF[$i]/3) if $exF[$i]%3==0 || ($exF[$i]+1)%3==0;
	    
	    $exS=join("",@S[$pI..$pF]);
	    print O_prot "$g\t$str\t$prot\t$i\t$pI-$pF\t$coIn[$i]-$coFn[$i]\t$exS\n";
	    $z++;
	}
    }
}

print "Parsing protein. Total exons: $pre/$z\n";
print LOG "Parsing protein. Total exons: $pre/$z\n";



