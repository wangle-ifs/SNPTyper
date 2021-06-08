#!/usr/bin/perl -w
use strict;
=head1 Program Description

It works for ilmn sequencing platform.

SNP calling & Genotyping

=head1 Contact & Version

  Author: YanDong, yandong.cao@analyses.cn
  Version: 0.1,  Date: 2018-10-20

=head1 Command-line Option
	--snpcsv	require, SNP.csv file, format of csv: chrom,position,rs
	--fastq		require, fastq file, PE fq1.gz,fq2.gz; SE fq.gz
	--target	require, target region file in bed format
	--index		optional, hg19 or hg38, default = hg19
	--num		optional, max number of reads after downsampling, default = 0 and skip this step of downsampling
	--prefix	optional, prefix of all output files, default= WL
	--outdir	optional, directory of output files, default= WL

=head1 Usage Exmples

  perl analyses-extractSnp.pl  --snpcsv SNP.csv  --fastq fq1.gz,fq2.gz  --target target.bed  --index hg19  --prefix WL  --outdir WL

=cut

use strict;
use FindBin qw($Bin);
use Getopt::Long;
use File::Basename;
use List::Util;

my ($help, $snpcsv, $fastq, $index, $prefix, $outdir, $target, $num) = ("", ".", ".", "hg19", "WL", "WL", ".", 0);
GetOptions(
    "snpcsv=s"	=>	\$snpcsv,
    "fastq=s"	=>	\$fastq,
    "target=s"	=>	\$target,
	"index:s"	=>	\$index,
	"outdir:s"	=>	\$outdir,
	"prefix:s"	=>	\$prefix,
	"num:i"		=>	\$num,
	);

die `pod2text $0` if $fastq eq "." || $target eq "." || $index eq ".";

system("mkdir -p $outdir");
my $srtbam = "$outdir/$prefix.srt.bam";
my $yourvcf = "$outdir/$prefix.vcf";
my $result = "$outdir/$prefix.xls";
my @fastq = split /,/, $fastq;
my $bwa = "/home/wangle/bwa-0.7.12/bwa";
my $samtools = "/home/wangle/samtools-1.5/samtools";
my $varscan2 = "/home/wangle/varscan-master/VarScan.v2.4.2.jar";
my $hg19 = "/home/wangle/hg19/hg19.fasta";
my $hg38 = "/home/wangle/hg38/hg38.fasta";
my $seqtk = "/home/wangle/seqtk/seqtk";
my $reference = $hg19;
$reference = $hg38 if $index eq "hg38";

if ($num > 0) {
	for (my $i = 0; $i < @fastq; $i++) {
		system("$seqtk sample $fastq[$i] $num > $fastq[$i].$num.fq");
		$fastq[$i] = "$fastq[$i].$num.fq";
	}
}
#print "$bwa mem -t 12 -M -R \"\@RG\tID:Default\tLB:Library\tPL:ILLUMINA\tSM:$prefix\" $reference @fastq | $samtools view -Shb - | $samtools sort - > $srtbam\n";
#print "$samtools mpileup -l $target -f $reference $srtbam | java -jar $varscan2 mpileup2cns --output-vcf 1 > $yourvcf\n";
system("$bwa mem -t 12 -M -R \"\@RG\tID:Default\tLB:Library\tPL:ILLUMINA\tSM:$prefix\" $reference @fastq | $samtools view -Shb - | $samtools sort - > $srtbam");
system("$samtools mpileup -d 999999 -l $target -f $reference $srtbam | java -jar $varscan2 mpileup2cns --output-vcf 1 > $yourvcf");

my %SNP;
open IN, $snpcsv or die "$! $snpcsv\n";
while (<IN>) {
	chomp;
	my @array = split /,/, $_;
	next unless exists $array[2] && $array[2] =~/^rs/;
	$SNP{"$array[0]:$array[1]"} = $array[2];
}
close IN;

my %variants;
open IN, $yourvcf or die "$! $yourvcf\n";
while (<IN>) {
	chomp;
	my @array = split /\t/, $_;
	next if $array[0] =~/\#/;
	#chr9	3584543	.	C	A	.	PASS	ADP=1872;WT=0;HET=1;HOM=0;NC=0	GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR	0/1:255:1872:1872:863:1004:53.63%:0E0:37:36:400:463:468:536
	if (exists$SNP{"$array[0]:$array[1]"}) {
		my @value = split /:/, $array[-1];
		my ($genotype,$freq) = @value[0,6];
		next if $genotype eq "./.";
		my @alt = split /,/, $array[4];
		@alt=($array[3],@alt);
		my @depth = split /,/, $value[5];
		@depth = ($value[4], @depth);
		my $refDep = $depth[0];

		my @genotype = split /\//, $genotype;
		$genotype = "$alt[$genotype[0]]/$alt[$genotype[1]]";
		
		my $depth = 0;
		if ($genotype[0] eq $genotype[1]) {
			$depth = $depth[$genotype[0]];
		}
		else {
			$depth = $depth[$genotype[0]] + $depth[$genotype[1]];
		}
		$refDep = 0 if $genotype[0] * $genotype[1] > 0;
		my $altDep = $depth - $refDep;
		
		my ($A, $T, $C, $G, $other) = (0, 0, 0, 0, "no ohter");
		$A = $depth[$genotype[0]] if $alt[$genotype[0]] eq "A";
		$A = $depth[$genotype[1]] if $alt[$genotype[1]] eq "A";
		$T = $depth[$genotype[0]] if $alt[$genotype[0]] eq "T";
		$T = $depth[$genotype[1]] if $alt[$genotype[1]] eq "T";
		$C = $depth[$genotype[0]] if $alt[$genotype[0]] eq "C";
		$C = $depth[$genotype[1]] if $alt[$genotype[1]] eq "C";
		$G = $depth[$genotype[0]] if $alt[$genotype[0]] eq "G";
		$G = $depth[$genotype[1]] if $alt[$genotype[1]] eq "G";
		$other = "$alt[$genotype[0]]:$depth[$genotype[0]]\;" if (length$alt[$genotype[0]]) > 1;
		$other .= "$alt[$genotype[1]]:$depth[$genotype[1]]\;" if (length$alt[$genotype[1]]) > 1;
		$variants{"$array[0]:$array[1]"} = {depth=>$depth, freq=>$freq, genotype=>$genotype, ref=>$array[3], refDep=>$refDep, altDep=>$altDep, A=>$A, T=>$T, C=>$C, G=>$G, other=>$other};
	}
}
open OUT, ">$result" or die "$! $result\n";
print OUT "#location\tSNP_ID\tgenotype\treference allele\tdepth of all\tdepth of reference\tdepth of variants\tdepth of allele A\tdepth of allele T\tdepth of allele C\tdepth of allele G\tother\n";
foreach my $loci (sort keys%SNP) {
	if (exists$variants{$loci}) {
		print OUT "$loci\t$SNP{$loci}\t$variants{$loci}{'genotype'}\t$variants{$loci}{'ref'}\t$variants{$loci}{'depth'}\t$variants{$loci}{'refDep'}\t$variants{$loci}{'altDep'}\t$variants{$loci}{'A'}\t$variants{$loci}{'T'}\t$variants{$loci}{'C'}\t$variants{$loci}{'G'}\t$variants{$loci}{'other'}\n";
	}
	else {
		print OUT "$loci\t$SNP{$loci}\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n";
	}
}
close OUT;

__END__
#bwa mem index fastq1 fastq2 | samtools sort - | samtools mpileup | java -jar VarScan.jar mpileup2cns
