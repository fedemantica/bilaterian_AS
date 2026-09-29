#!/usr/bin/perl

# it makes the DOMAINS file ONLY for actual Domains. The "All" makes a gigantic file (many TB).
# Adapted: Irene Oct2024

$pfamR=$ARGV[0];
$aprotp=$ARGV[1];
$outfile=$ARGV[2];



open (PS, $pfamR) || die "Can't find Parsed ps_out\n";
open (PROT, $aprotp) || die "Can't find the parsed prots\n";
open (O, ">".$outfile);

print O "GeneID\tStrand\tProteinID\tExon_Rank\tProt_Coord\tGenome_Coord\tDomain_hit\n";

<PS>;
while (<PS>){
    chomp;
    @t=split(/\t/);
    
    $a=0 if $t[1] ne $prB;

    $t[4]=~s/\.+//;
    $dom{$t[1]}[$a]="$t[2]-$t[3]";
    $name{$t[1]}[$a++]="$t[4]=$t[5]";
    
    $prB=$t[1];
}
close PS;

while (<PROT>){
    chomp;
    @t=split(/\t/);
    $pr=$t[2];
    $pre=join("\t",@t[0..5]);
    ($ei,$ef)=$t[4]=~/(\d+?)\-(\d+)/;
    $type="";
    
    for $z (0..$#{$dom{$pr}}){
	($di,$df)=$dom{$pr}[$z]=~/(.+?)\-(\d+)/;
	$ledom=abs($df-$di)+1;
	$leex=abs($ef-$ei)+1;

	if ($ei>$df || $ef<$di){} #nothing

	elsif ($ei<=$di && $ef>=$di){
	    if ($ef<$df){ #partial overlap of both on the left side
		$dif=$ef-$di;
		$perc=sprintf("%.1f",100*$dif/$ledom);
		$perc="0.1" if $perc eq "0.0";
		$percE=sprintf("%.1f",100*$dif/$leex);
		$type.="$name{$pr}[$z]=PU($perc=$percE),";
	    }
	    elsif ($ef>=$df){ #exon fully contains the domain
		$percE=sprintf("%.1f",100*$ledom/$leex);
		$type.="$name{$pr}[$z]=WD(100=$percE),";
	    }
	}
	elsif ($ei>$di && $ei<=$df){ #partial overlap of both on the right side
	    if ($ef>$df){
		$dif=$df-$ei;
		$perc=sprintf("%.1f",100*$dif/$ledom);
		$perc="0.1" if $perc eq "0.0";
		$percE=sprintf("%.1f",100*$dif/$leex);
		$type.="$name{$pr}[$z]=PD($perc=$percE),";
	    }
	    elsif ($ef<=$df){ #the whole exon is contained in the domain
		$dif=$ef-$ei;
		$perc=sprintf("%.1f",100*$dif/$ledom);
		$perc="0.1" if $perc eq "0.0";
		$type.="$name{$pr}[$z]=FE($perc=100),";
	    }
	}
    }
    chop($type);

    print O "$pre\t$type\n";
}
