#!/usr/bin/perl

# ARGV[0]= disorder, position, by position
# ARGV[1]= parsed annotated protein file (tab)
# ARGV[2]= REFERENCE with C1 A C2

# Adaptation of the script: Irene Oct 2024

# Arguments
$DisFile=$ARGV[0]; # DIS
$ParsProtFile=$ARGV[1]; # I
$RefFile=$ARGV[2];  # REF
$EvGenIDs=$ARGV[3]; # IDs
$sp=$ARGV[4];  
$OutDis=$ARGV[5];
$logWarnings=$ARGV[6];

#$OutDis=$datadir."/".$sp."/REFERENCE-ALL_ANNOT-".$sp."-DisI-v2.tab";  # O

# updated on 03/02/17: comments

die "\nUsage: CalculateAndPutDis_C1AC2.pl Sp-DisorderN-v.txt Sp_prot_annot_parsed-v.tab REFERENCE.tab\n\n" if !$ARGV[0];


#### parse the disorder file
open (DIS, $DisFile) || die "Needs Disorder\n"; # IUPRED result
while (<DIS>){
    chomp;
    @t=split(/\t/);
    $Dr{$t[0]}{$t[1]}=$t[4];
    $Dc{$t[0]}{$t[1]}=$t[3];
    $seq{$t[0]}{$t[1]}=$t[2];
}
close DIS;

#### parse annotated protein file
open (I,$ParsProtFile ) || die "Parsed annotated protein file\n";  # prot_annot_parsed tab
while (<I>){
    chomp;
    @t=split(/\t/);
    $g=$t[0];
    $pr=$t[2];
    #    $aa=$t[4];

    ($start,$end)=$t[4]=~/(\d+?)\-(\d+)/;
    $ex="$g=$t[5]";
    ($ei,$ef)=$t[5]=~/(\d+?)\-(\d+)/;
    $coA="$g=A$ei";
    $coB="$g=B$ef";

    for $i ($start..$end){
	if ($seq{$pr}{$i}){
	    if ($total{$ex} && $pr=~/f\d.+?[ABC]/){}
	    else {
		#$ProtID{$ex}={$pr};
		#$AAnum{$ex}={$aa};
		#$GeneID{$ex}={$g};
		$total{$ex}++;
		$tally{$ex}+=$Dc{$pr}{$i};
		$total{$coA}++;
		$tally{$coA}+=$Dc{$pr}{$i};
		$total{$coB}++;
		$tally{$coB}+=$Dc{$pr}{$i};
		
		$tallyR{$ex}+=$Dr{$pr}{$i};
		$tallyR{$coA}+=$Dr{$pr}{$i};
		$tallyR{$coB}+=$Dr{$pr}{$i};
	    }
	}
    }
}   
   
close I;


#### gets GeneIDs for each event
open (IDs, $EvGenIDs) || die "Can't find the Event-Gene.IDs key";
#<IDs>;
while (<IDs>){
    chomp;
    @t=split(/\t/);
    $id_g{$t[0]}=$t[1];
}
close IDs;


### Open log for the warnings
open (LOGW, ">>".$logWarnings);

#### Starts analyzing the events
open (O, ">".$OutDis); # output

open (REF, $RefFile) || die "Reference file with C1 and C2\n";  # REFERENCE-ALL tab
$head=<REF>;
chomp($head);
@head=split(/\t/,$head);
$headF=join("\t",@head[0..7]);

print O "$headF\tC1_Dis\tA_Dis\tC2_Dis\tC1_rDis\tA_rDis\tC2_rDis\n";
while (<REF>){
    chomp;
    @t=split(/\t/);

    $ev=$t[1];
    $g="";
    $g=$id_g{$t[1]} if !$g;
    
    print LOGW "\tWarning: No gene ID for $ev\n" if !$g;
    next if !$g;

    $pre=join("\t",@t[0..7]);
    #$PIDa=$GIDa=$PaaA="na";
    $disC1=$disA=$disC2=$RdisC1=$RdisA=$RdisC2="na";
## C1
    $ex=$coA=$coB="";
    if ($t[5]!~/Alt5/){
	($start,$end)=$t[8]=~/(\d+?)\-(\d+)/;
	$ex="$g=$&";
	$coA="$g=A$start";
	$coB="$g=B$end";
	$disC1=sprintf("%.3f",($tally{$ex}+$tally{$coA}+$tally{$coB})/($total{$ex}+$total{$coA}+$total{$coB})) if ($total{$ex}+$total{$coA}+$total{$coB})>0;
	$RdisC1=sprintf("%.3f",($tallyR{$ex}+$tallyR{$coA}+$tallyR{$coB})/($total{$ex}+$total{$coA}+$total{$coB})) if ($total{$ex}+$total{$coA}+$total{$coB})>0;
    }
    elsif ($t[5]=~/Alt5/){
	($temp_A,$temp_B)=$t[4]=~/\:(.*?)\-(.*?)\,/;
	($str)=$t[6]=~/.+\:([\+\-])/;
	$temp_sum=$temp_sumR=$temp_tot=0;

	if ($str eq "+" && $temp_B =~/\+/){
	    @temp_B=split(/\+/,$temp_B);
	    foreach $t_B (@temp_B){
		$coB="$g=B$t_B";
		$temp_sum+=$tally{$coB};
		$temp_sumR+=$tallyR{$coB};
		$temp_tot+=$total{$coB};
	    }
	}
	elsif ($str eq "-" && $temp_A =~/\+/){
	    @temp_A=split(/\+/,$temp_A);
	    foreach $t_A (@temp_A){
		$coA="$g=A$t_A";
		$temp_sum+=$tally{$coA};
		$temp_sumR+=$tallyR{$coA};
		$temp_tot+=$total{$coA};
	    }
	}
	else {
	    print "Unexpected pattern of + and strand for $t[1]\n";
	}
	$disC1=sprintf("%.3f",$temp_sum/$temp_tot) if $temp_tot>0;
	$RdisC1=sprintf("%.3f",$temp_sumR/$temp_tot) if $temp_tot>0;
    }
## AS
    $ex=$coA=$coB="";
    unless ($t[9] eq "NA" || $t[1]=~/ALT/){
	($start,$end)=$t[9]=~/(\d+?)\-(\d+)/;
	$ex="$g=$&";
	$coA="$g=A$start";
	$coB="$g=B$end";
	$disA=sprintf("%.3f",($tally{$ex}+$tally{$coA}+$tally{$coB})/($total{$ex}+$total{$coA}+$total{$coB})) if ($total{$ex}+$total{$coA}+$total{$coB})>0;
	$RdisA=sprintf("%.3f",($tallyR{$ex}+$tallyR{$coA}+$tallyR{$coB})/($total{$ex}+$total{$coA}+$total{$coB})) if ($total{$ex}+$total{$coA}+$total{$coB})>0;
	#$GIDa=$GeneID{$ex};
	#$PIDa=$ProtID{$ex};
	#$PaaA=$AAnum{$ex};
    }
## C2
    $ex=$coA=$coB="";
    if ($t[5]!~/Alt3/){
	($start,$end)=$t[10]=~/(\d+?)\-(\d+)/;
	$ex="$g=$&";
	$coA="$g=A$start";
	$coB="$g=B$end";
	$disC2=sprintf("%.3f",($tally{$ex}+$tally{$coA}+$tally{$coB})/($total{$ex}+$total{$coA}+$total{$coB})) if ($total{$ex}+$total{$coA}+$total{$coB})>0;
	$RdisC2=sprintf("%.3f",($tallyR{$ex}+$tallyR{$coA}+$tallyR{$coB})/($total{$ex}+$total{$coA}+$total{$coB})) if ($total{$ex}+$total{$coA}+$total{$coB})>0;
    }
    elsif ($t[5]=~/Alt3/){
	($temp_A,$temp_B)=$t[4]=~/\,(.*?)\-(.*)/;
	($str)=$t[6]=~/.+\:([\+\-])/;
	$temp_sum=$temp_sumR=$temp_tot=0;

	if ($str eq "-" && $temp_B =~/\+/){
	    @temp_B=split(/\+/,$temp_B);
	    foreach $t_B (@temp_B){
		$coB="$g=B$t_B";
		$temp_sum+=$tally{$coB};
		$temp_sumR+=$tallyR{$coB};
		$temp_tot+=$total{$coB};
	    }
	}
	elsif ($str eq "+" && $temp_A =~/\+/){
	    @temp_A=split(/\+/,$temp_A);
	    foreach $t_A (@temp_A){
		$coA="$g=A$t_A";
		$temp_sum+=$tally{$coA};
		$temp_sumR+=$tallyR{$coA};
		$temp_tot+=$total{$coA};
	    }
	}
	else {
	    print LOGW "Unexpected pattern of + and strand for $t[1]\n";
	}
	$disC2=sprintf("%.3f",$temp_sum/$temp_tot) if $temp_tot>0;
	$RdisC2=sprintf("%.3f",$temp_sumR/$temp_tot) if $temp_tot>0;
    }
##    
  
    print O "$pre\t$disC1\t$disA\t$disC2\t$RdisC1\t$RdisA\t$RdisC2\n";
}
