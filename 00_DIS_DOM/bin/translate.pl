#!/usr/bin/perl
$t{"TTT"} = "F";
$t{"TTC"} = "F";
$t{"TTA"} = "L";
$t{"TTG"} = "L";
$t{"TCT"} = "S";
$t{"TCC"} = "S";
$t{"TCA"} = "S";
$t{"TCG"} = "S";
$t{"TAT"} = "Y";
$t{"TAC"} = "Y";
$t{"TAA"} = "*";
$t{"TAG"} = "*";
$t{"TGT"} = "C";
$t{"TGC"} = "C";
$t{"TGA"} = "*";
$t{"TGG"} = "W";
$t{"CTT"} = "L";
$t{"CTC"} = "L";
$t{"CTA"} = "L";
$t{"CTG"} = "L";
$t{"CCT"} = "P";
$t{"CCC"} = "P";
$t{"CCA"} = "P";
$t{"CCG"} = "P"; 
$t{"CAT"} = "H";
$t{"CAC"} = "H";
$t{"CAA"} = "Q";
$t{"CAG"} = "Q";
$t{"CGT"} = "R";
$t{"CGC"} = "R";
$t{"CGA"} = "R";
$t{"CGG"} = "R";
$t{"ATT"} = "I";
$t{"ATC"} = "I";
$t{"ATA"} = "I";
$t{"ATG"} = "M";
$t{"ACT"} = "T";
$t{"ACC"} = "T";
$t{"ACA"} = "T";
$t{"ACG"} = "T";
$t{"AAT"} = "N";
$t{"AAC"} = "N";
$t{"AAA"} = "K";
$t{"AAG"} = "K";
$t{"AGT"} = "S";
$t{"AGC"} = "S";
$t{"AGA"} = "R";
$t{"AGG"} = "R";
$t{"GTT"} = "V";
$t{"GTC"} = "V";
$t{"GTA"} = "V";
$t{"GTG"} = "V";
$t{"GCT"} = "A";
$t{"GCC"} = "A";
$t{"GCA"} = "A";
$t{"GCG"} = "A";
$t{"GAT"} = "D";
$t{"GAC"} = "D";
$t{"GAA"} = "E";
$t{"GAG"} = "E";
$t{"GGT"} = "G";
$t{"GGC"} = "G";
$t{"GGA"} = "G";
$t{"GGG"} = "G";

@c = (0,2,1);

($root) = $ARGV[0] =~ /([^\.]+)/;
# Take name of the file: {specie}_CDS_annot-off.fasta
# keep only until dot: {specie}_CDS_annot-off

print "Translating CDS sequences\n";  # Irene commenting

# OLD: Output file: {specie}_CDS_annot-off.pro	
#open (OUT, ">$root.pro");

# Output file now: 
$root =~ s/CDS/prot/;  # replace "CDS" with "prot"
$root =~ s/-off//;  # remove "-off"
# Output file: {specie}_prot_annot.fasta	
open (OUT, ">$root.fasta");

$/ = ">";
$_ = <>;
while (<>) {
	/\n/;
	$head = $`;
	$seq = $';
	$seq =~ tr/a-z/A-Z/;
	#($phase) = /\((\d+)/;
	#$head =~ s#\((\d+)\)#sprintf "(%d)", ($1-$c[$phase])/3#e;
	$seq =~ s/[^A-Za-z]+//g;

	#$seq =~ s/.{$c[$phase]}//;
$seq =~ tr/a-z/A-Z/;
	$seq =~ s/[a-z]+//g;
	print OUT ">$head\n";
	while ($seq =~ s/([A-Z\-]{3})//){
		if ($t{$1}) {print OUT $t{$1}}
		else {print OUT "X"}}
	print OUT "\n";
}
print "Translation done\n\n"; # Irene commenting	








