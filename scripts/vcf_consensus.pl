#!/usr/bin/perl 
# 
# Quick hack to convert affy files to semi VCF files
# 
# 
# Kim Brugger (12 Mar 2014), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib "/software/packages/vcftools_0.1.11/lib/perl5/site_perl";
use Vcf;

my $MIN_QUAL = 500;

my @singles;

my $vcf_file = shift;

my $vcf = Vcf->new(file=>$vcf_file);
$vcf->parse_header();


my %counts;

my $exit_counter = 50000;

while (my $entry = $vcf->next_data_hash()) {
#  last if ( ! $exit_counter--);

  my %GTs;

  foreach my $gtype ( keys %{$$entry{ gtypes }} ) {
    next if $gtype eq "AFFY";
    next if ( $$entry{ gtypes }{ $gtype}{GT } eq "." || $$entry{ gtypes }{ $gtype}{GT } eq "./.");
    $$entry{ gtypes }{ $gtype}{GT } = "0/1" if ($$entry{ gtypes }{ $gtype}{GT } eq "1/0");

    $GTs{ $$entry{ gtypes }{ $gtype}{ GT } }++;

  }


  if ( int(keys %GTs) > 1 ) {
    my $v;
#    map { $v .= "\t$_:$GTs{$_}"} keys %GTs;
#    print "$$entry{CHROM}:$$entry{POS}$v\n";
#    print "$$entry{CHROM}:$$entry{POS} " . join("\t", keys %GTs) ."\n";
  }


  foreach my $gt ( keys %GTs ) {
    $counts{ $GTs{ $gt } }++;

    if (  $GTs{ $gt } >= 4) {
#      print "$$entry{CHROM}:$$entry{POS} $gt -- $GTs{ $gt }\n";
#      print Dumper (\%GTs);
    }
  }


#  print Dumper( \%GTs );

#  die Dumper( $entry );

  


}

foreach my $count ( sort { $a <=> $b } keys %counts ) {
  print "$count:$counts{$count}\t";
}

print "\n";

#print Dumper( \%counts );
