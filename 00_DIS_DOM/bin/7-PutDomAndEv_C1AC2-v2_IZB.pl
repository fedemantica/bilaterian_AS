#!/usr/bin/perl

# IreneJuly 2025: Important note about the "foreach" loop:
# Each event has different domains depending on the isoform in input PFAM file,
# Only the "longest" stays. So if there are 2 domains in one isoforma and 1 
# domain in another, the domains listed for A are those 2.
# If only 1 domain in both, the longest domain name. 
# If same name but overlap 1 is 15=81 and overlap 2 is 9=81, then 15=81 stays
# just because it has more characters, longest string.
# Totally possible that the domains of C1, A, C2 are not part of same isoform.
# (also, C2 maybe it's not really after A in a real isoform)
# Look at the plot vs the annotation: https://vastdb.crg.eu/event/HsaEX6064755@hg38
# Domain in C2 (plot, speciphic isoform) but not the reference
# Also, look at the tracks to understand it.


# Script adapted by Irene, Nov 2024
die "\nUsage: PutDomAndEv_C1AC2-v2.pl Sp_prot_annot-v-DomType.tab REFERENCE-ALL_ANNOT.tab\n\n" if !$ARGV[1];

$PfamRes=$ARGV[0];
$RefFile=$ARGV[1];
$EventGene=$ARGV[2];
$outFile=$ARGV[3];
$logWarnings=$ARGV[4];

open (LOGW, ">>".$logWarnings);
open (DOM,$PfamRes) || die "Needs Domain parsed file (e.g. Hsa_prot_annot-v71-PFAM.tab)\n";
<DOM>; # added 26/08/19 (old ones do not have header)
while (<DOM>){
    chomp;
    @t=split(/\t/);
    $ex="$t[0]=$t[5]";
    ($i,$f)=$t[5]=~/(.+?)\-(.+)/;
    $coA="$t[0]=$i";
    $coB="$t[0]=$f";
    
    $D{$ex}.="$t[6];";
    $D{$coA}.="$t[6];";
    $D{$coB}.="$t[6];";
}
close DOM;


open (O, ">".$outFile);

open (IDs, $EventGene);
<IDs>;
while (<IDs>){
    chomp;
    @t=split(/\t/);
    $id_g{$t[0]}=$t[1];
}
close IDs;

#open (LOGW, ">>".$logWarnings);
open (I, $RefFile) || die "Needs the REFERENCE\n";
$head=<I>;
chomp($head);
@head=split(/\t/,$head);
$headF=join("\t",@head[0..7]);
print O "$headF\tC1_DOM\tA_DOM\tC2_DOM\n";
while (<I>){
    chomp;
    @t=split(/\t/);
    $ev=$t[1];
    $g="";
    $g=$id_g{$t[1]} if !$g;

    print LOGW "\tWarning: No gene ID for $ev\n" if !$g;
    next if !$g;

#### Mathing C1    
    $ex=$coA=$coB=$ex2="";
    if ($t[5] !~ /Alt5/){ # different for Alt5
	($start,$end)=$t[8]=~/\:(\d+?)\-(\d+)/;
	$ex="$g=$&";
	$coA="$g=$start";
	$coB="$g=$end";
	
	$allC1="";
	$allC1=$D{$ex};
	$allC1=$D{$coA} if !$allC1;
	$allC1=$D{$coB} if !$allC1;
    }
    elsif ($t[5] =~ /Alt5/){ # different for Alt5
	($temp_A,$temp_B)=$t[4]=~/\:(.*?)\-(.*?)\,/;
	($str)=$t[6]=~/.+\:([\+\-])/;
	$allC1="";

	if ($str eq "+" && $temp_B =~/\+/){
	    @temp_B=split(/\+/,$temp_B);
	    foreach $t_B (@temp_B){
		$coB="$g=$t_B";
		$allC1=$D{$coB} if !$allC1 && $D{$coB};
	    }
	}
	elsif ($str eq "-" && $temp_A =~/\+/){
	    @temp_A=split(/\+/,$temp_A);
	    foreach $t_A (@temp_A){
		$coA="$g=$t_A";
		$allC1=$D{$coA} if !$allC1 && $D{$coA};
	    }
	}
	else {
	    print LOGW "\nUnexpected pattern of + and strand for $t[1]\n";
	}
    }
    
    $allC1="NO" if ($allC1=~/\;/ && $allC1!~/\=/);
    $allC1="na" if !$allC1;
    
    %done=();
    @An=();
    $zz=0;
    if ($allC1=~/\;.+?\;/){
	@A=split(/\;/,$allC1);
	$keep="";
	foreach $z (@A){
	    $keep=$z if length($z)>length($keep);
	}
	$allC1=$keep;
    }
    elsif ($allC1=~/\;/) {
	chop($allC1);
#	$allC1=~s/\;$//g;
    }
##AS
    $ex=$coA=$coB=$ex2="";
    unless ($t[9] eq "NA" || $t[1]=~/ALT/){ # temporarily ignoring the A seq ALTA/D
	($start,$end)=$t[9]=~/\:(\d+?)\-(\d+)/;
	$ex="$g=$&";
	$coA="$g=$start";
	$coB="$g=$end";
	unless ($t[1]=~/ALT/){
	    ($c2)=$t[2]=~/\:(.+)/; # second attempt to get it
	    $ex2="$g=$c2";
	}
	
	$allA="";
	$allA=$D{$ex};
	$allA=$D{$ex2} if !$allA && $t[1]!~/ALT/;
	$allA=$D{$coA} if !$allA;
	$allA=$D{$coB} if !$allA;
	
	$allA="NO" if ($allA=~/\;/ && $allA!~/\=/);
	$allA="na" if !$allA;

	%done=();
	@An=();
	$zz=0;

#	# Irene Debugging
	# print("What is ex: $ex\n");
	# if ($ex eq "ENSG00000111641=:6563081-6563170"){
	# 	print("\nallA before the splitting: $allA\n");
	# 	}
	####
	if ($allA=~/\;.+?\;/){
#		# Irene Debugging		
		# if ($ex eq "ENSG00000111641=:6563081-6563170"){
		# 	print("\nallA entering the keep loop: $allA\n");
		# }
		#
		@A=split(/\;/,$allA);
		$keep="";
		foreach $z (@A){
#			#	# Irene Debugging
			#	if ($ex eq "ENSG00000111641=:6563081-6563170"){
			#	print("\nInside the foreach, z is $z \nand keep is $keep\n");
			#	print("The length of z is ". length($z) ." and the length of keep is ".length($keep) ."\n\n");
			#	}
		####
			$keep=$z if length($z)>length($keep);

		}
		$allA=$keep;
#		#		# Irene Debugging
		#if ($ex eq "ENSG00000111641=:6563081-6563170"){
		#	print("\nallA after the keep: $allA\n");
		#}
		####

    	}
	elsif ($allA=~/\;/) {
		chop($allA);

#		# Irene Debugging
		#		if ($ex eq "ENSG00000111641=:6563081-6563170"){
		#			print("\nallA chopped: $allA\n");
		#		}
		####		
#	    $allA=~s/\;$//g;
	}
    }
    else {
	$allA="na";
    }
##C2
    $ex=$coA=$coB=$ex2="";
    if ($t[5] !~ /Alt3/){
	($start,$end)=$t[10]=~/\:(\d+?)\-(\d+)/;
	$ex="$g=$&";
	$coA="$g=$start";
	$coB="$g=$end";
	
	$allC2="";
	$allC2=$D{$ex};
	$allC2=$D{$coA} if !$allC2;
	$allC2=$D{$coB} if !$allC2;
    }
    elsif ($t[5] =~ /Alt3/){ # different for Alt3
	($temp_A,$temp_B)=$t[4]=~/\,(.*?)\-(.*)/;
	($str)=$t[6]=~/.+\:([\+\-])/;
	$allC2="";

	if ($str eq "+" && $temp_A =~ /\+/){
	    @temp_A=split(/\+/,$temp_A);
	    foreach $t_A (@temp_A){
		$coA="$g=$t_A";
		$allC2=$D{$coA} if !$allC2 && $D{$coA};
	    }
	}
	elsif ($str eq "-" && $temp_B =~ /\+/){
	    @temp_B=split(/\+/,$temp_B);
	    foreach $t_B (@temp_B){
		$coB="$g=$t_B";
		$allC2=$D{$coB} if !$allC2 && $D{$coB};
	    }
	}
	else {
	    print LOGW "\nUnexpected pattern of + and strand for $t[1] ($str and $temp_B/$temp_A)\n";
	}
    }

    $allC2="NO" if ($allC2=~/\;/ && $allC2!~/\=/);
    $allC2="na" if !$allC2;

    %done=();
    @An=();
    $zz=0;
    if ($allC2=~/\;.+?\;/){
	@A=split(/\;/,$allC2);
	$keep="";
	foreach $z (@A){
	    $keep=$z if length($z)>length($keep);
	}
	$allC2=$keep;
    }
    elsif ($allC2=~/\;/) {
	chop($allC2);
#	$allC2=~s/\;$//g;
    }
    
    ##
    $precol=join("\t",@t[0..7]);
    while ($allC1=~s/\;$//){}
    while ($allA=~s/\;$//){}
    while ($allC2=~s/\;$//){}	
    print O "$precol\t$allC1\t$allA\t$allC2\n";
}

close LOGW;
