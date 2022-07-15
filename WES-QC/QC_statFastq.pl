#!/home/fanyucai/software/perl/perl-v5.24.1/bin/perl -w
use strict;
use warnings;
use JSON;
use File::Spec;
use File::Basename qw(basename dirname);
use Getopt::Long;

my $version="v2.0";

my $path_curf = File::Spec->rel2abs($0);
my $work_dir||=dirname($path_curf);

my $Rscript="module purge && module load mro/3.5.1 && Rscript";
my $R_cont="$work_dir/QC.content.r";
my $R_qual="$work_dir/QC.qual.r";
my $R_percent="$work_dir/QC.CleanReadsPercent.r";


##############################################   Usage   ############################################
my ($input,$sample,$outdir,$help);
GetOptions(
        "i:s"=>\$input,
        "s:s"=>\$sample,
        "o:s"=>\$outdir,
        "h|help!"  => \$help,
);


my $Usage=<<USAGE;

Program : Capture $version 
Tuthor : Yincong.gu 2018-6-19
Discriptions : This script will carry out statistics on QC.

options:
	-i              fastp json file        [ force ]
	-s              sample name            [ force ]
	-o              outdir                 [ force ]
	-help|h         print help information
	
Example:
	perl $0 -w [work dir] -o [outdir]
USAGE
die $Usage if( ! $input || !$sample || ! $outdir ||$help );



my $dir = "$outdir/$sample";
system("mkdir -p $dir") if ( ! -e "$dir" );

my $QC="$dir/$sample.QC.stat.xls";
my $quality="$dir/$sample.quality";
my $acgtn="$dir/$sample.acgtn";
my $filter="$dir/$sample.filter";
open(ACTG,">$acgtn");
open(QUAL,">$quality");
open(FILT,">$filter");
print ACTG "pos\tA\tT\tC\tG\tN\tGC\n";
print QUAL "pos\tmean\n";
print FILT "mode\tnumber\tperc\n";

open(OUT,">$QC");
print OUT "Samples\tRaw_Reads\tClean_Reads\tClean_Reads_Percent\tRaw_Base\tClean_Base\tClean_Base_Percent\tGC_Content\t>Q20\t>Q30\n";


open( my $JSON, "<", "$input") or die "Cannot open the json output of sample $sample from fastp$!\n";
my $jstring = "";
while ( defined(my $line = <$JSON> ) ){
	chomp $line;
	$jstring.= $line;
}

my $js_object = decode_json($jstring);
my $RawR   = $js_object->{summary}{before_filtering}{total_reads};
my $CleanR = $js_object->{summary}{after_filtering}{total_reads};
my $percR = sprintf("%.2f",($CleanR*100/$RawR) );

my $RawB   = $js_object->{summary}{before_filtering}{total_bases};
my $CleanB = $js_object->{summary}{after_filtering}{total_bases};
my $percB = sprintf("%.2f",($CleanB*100/$RawB) );

my $shortRB=&unit($RawB);
my $shortCB=&unit($CleanB);

my $Q20    = sprintf("%.2f",$js_object->{summary}{after_filtering}{q20_rate}*100);
my $Q30    = sprintf("%.2f",$js_object->{summary}{after_filtering}{q30_rate}*100);
my $GC     = sprintf("%.2f",$js_object->{summary}{after_filtering}{gc_content}*100);
my $LQ     = $js_object->{filtering_result}{low_quality_reads};
my $N      = $js_object->{filtering_result}{too_many_N_reads};
my $short  = $js_object->{filtering_result}{too_short_reads};
my @A=(@{$js_object->{read1_after_filtering}{content_curves}{A}},@{$js_object->{read2_after_filtering}{content_curves}{A}});
my @T=(@{$js_object->{read1_after_filtering}{content_curves}{T}},@{$js_object->{read2_after_filtering}{content_curves}{T}});
my @C=(@{$js_object->{read1_after_filtering}{content_curves}{C}},@{$js_object->{read2_after_filtering}{content_curves}{C}});
my @G=(@{$js_object->{read1_after_filtering}{content_curves}{G}},@{$js_object->{read2_after_filtering}{content_curves}{G}});
my @N=(@{$js_object->{read1_after_filtering}{content_curves}{N}},@{$js_object->{read2_after_filtering}{content_curves}{N}});
my @GCs=(@{$js_object->{read1_after_filtering}{content_curves}{GC}},@{$js_object->{read2_after_filtering}{content_curves}{GC}});
my @mean = (@{$js_object->{read1_after_filtering}{quality_curves}{mean}},@{$js_object->{read2_after_filtering}{quality_curves}{mean}});
my $length1 = $js_object->{read1_after_filtering}{total_cycles};
my $length2 = $js_object->{read2_after_filtering}{total_cycles};



### print  ###
printf FILT "high_quality\t$CleanR\t%s\n",($CleanR*100/$RawR);
printf FILT "low_quality\t$LQ\t%s\n",($LQ*100/$RawR);
printf FILT "too_many_N\t$N\t%s\n",($N*100/$RawR);
printf FILT "too_short_reads\t$short\t%s\n",($short*100/$RawR);
foreach my $i (0..($length1-1)) {
	printf ACTG "%s\t$A[$i]\t$T[$i]\t$C[$i]\t$G[$i]\t$N[$i]\t$GCs[$i]\n",($i+1);
	printf QUAL "%s\t$mean[$i]\n",($i+1);
}
foreach my $i (0..($length2-1)) {
	my $j=$i+$length1;
	printf ACTG "%s\t$A[$j]\t$T[$j]\t$C[$j]\t$G[$j]\t$N[$j]\t$GCs[$j]\n",($i+1);
	printf QUAL "%s\t$mean[$j]\n",($i+1);
}

close QUAL;close ACTG;close FILT;

print OUT "$sample\t$RawR\t$CleanR\t$percR%\t$RawB\($shortRB\)\t$CleanB\($shortCB\)\t$percB%\t$GC%\t$Q20%\t$Q30%\n";
system("$Rscript $R_qual infile=$outdir/$sample/$sample.quality outfile=$outdir/$sample/$sample.quality.png");
system("$Rscript $R_cont infile=$outdir/$sample/$sample.acgtn outfile=$outdir/$sample/$sample.acgtn.png");
system("$Rscript $R_percent $outdir/$sample/$sample.filter $outdir/$sample/$sample.QCsummary.png");



########### function #############
sub unit {
	my $input=shift;
	my $value;
	if($input > 100000000 ) {
		$value=sprintf("%.2f",($input/1000000000) )."G";
	} elsif ($input > 100000 ) { 
		$value=sprintf("%.2f",($input/1000000) )."M";
	} elsif ($input > 10 ) { 
		$value=sprintf("%.2f",($input/1000) )."K";
	}
	return($value);
}
