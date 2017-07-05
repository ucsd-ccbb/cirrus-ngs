#!/usr/bin/perl

$in2 = "/shared/workspace/RNASeqPipeline/kallisto_deseq_workflow/scripts/gencode.v23.metadata.EntrezGene"; 
open (IN2, $in2) or die "cannot open ";
$in5 = "/shared/workspace/RNASeqPipeline/kallisto_deseq_workflow/scripts/Hsa_gene_symbol_description.txt";
open (IN5, $in5) or die "cannot open ";

while($line = <IN2>) {
  chomp $line;
  ($t, $g) = (split /\t/, $line)[0,1];
  $genehash{$t} = $g;
  $mash{$g} = 0;
}
while($line = <IN5>) {
  chomp $line;
  ($gene,$sym,$descr) = (split /\t/, $line)[0,1,2];
  $gene2descr{$gene} = "$sym\t$descr";
}

$sample = $ARGV[0];

open(SF, "$sample\/abundance\.txt");
$line = <SF>; #skip header
while(<SF>)
{
	chomp($_);
	($transcr, $count) = (split /\t/,$_)[0,3];
	$transcr = (split /\|/,$transcr)[0];
	$mash{$genehash{$transcr}}+=$count;
}

close SF;

$outfile = ">".$sample."_counts.txt";
open (OUTG, $outfile);
print OUTG "gene\tsymbol\tdescription\t$sample\n";

foreach $key (sort (keys(%mash))) {
	print OUTG "$key\t$gene2descr{$key}\t$mash{$key}\n";
}
close(OUTG);

