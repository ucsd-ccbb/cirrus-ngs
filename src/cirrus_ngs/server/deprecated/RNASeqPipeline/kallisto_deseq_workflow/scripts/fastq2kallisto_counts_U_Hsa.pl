#does read mapping assuming PAIRED reads

#!/usr/bin/perl

my $inputfile = $ARGV[0];

# input file should contain these lines: 
# path to fastq (or fastq.gz) files
# list of sample prefixes to process into counts and sam, for example (without the "#"):

#/media/rs_EBS/RNA-seq/Olefsky
#LXS26ER
#OXS01ER

if (!open(IN, "$inputfile")) {
	print  ("Could not open $inputfile\n");
	die;
}

$dir = <IN>;
chomp($dir);
$dir=~s/^\s+//;
$dir=~s/\s+$//;
if ($dir =~ / /) { print "Path may not contain spaces, sorry.";
   die
}
if ($dir !~ /\/$/) { $dir.="/"} #append / if needed


opendir(DIR, $dir) or die $!;
my @files = grep{ /\.fastq/ } readdir(DIR); #only the .fastq files
closedir(DIR);

while($sample = <IN>) {
   chomp($sample);
   print "processing sample $sample\n";

   system("gunzip $dir$sample\*.gz"); #unzip fastq files that were gzipped
   system("cat $dir$sample\*_R1_*fastq > $dir$sample\_R1.fastq");
#   system("cat $dir$sample\*_R2_*fastq > $dir$sample\_R2.fastq");

#  system("/media/rs_EBS/kallisto/kallisto_linux-v0.42.1/kallisto quant -i /media/rs_EBS/ENSEMBL.homo_sapiens.release-75/gencode23/kallisto_index -o $sample $dir$sample\_R1.fastq $dir$sample\_R2.fastq");
   system("/shared/workspace/software/kallisto_linux-v0.42.1/kallisto quant -i /shared/workspace/software/kallisto_linux-v0.42.1/kallisto_index/kallisto_index -o $sample --single -l 50 $dir$sample\_R1.fastq");
#system("awk -F\"\t\" '\$1~/^@/ || \$7~/=/ {print \$0}' $dir$sample\_Aligned.out.sam > $dir$sample.sam");

   system("wc -l $dir$sample\_R1.fastq >> raw_counts_times_4.txt");
   system("rm $dir$sample*fastq");

}
