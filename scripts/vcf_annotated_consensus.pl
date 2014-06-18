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

#print Dumper( $vcf->_read_column_names());

my %counts;

my $exit_counter = 5000;

while (my $entry = $vcf->next_data_hash()) {

  my $CSQ_line = $$entry{INFO}{'CSQ'};
  my @CSQs;
    
  map {push @CSQs, [split(/\|/, $_)] } split(",", $CSQ_line);
  
    
  next if ( $$entry{QUAL} eq "." || $$entry{QUAL} < $MIN_QUAL );
    
  my $interesting_variant = 0;

  foreach my $CSQ ( @CSQs) {
    my ($Allele,$ENS_gene, $HGNC,$RefSeq,$feature,$effects,$CDS_position,$Protein_position,$Amino_acid,$Existing_variation,$SIFT,$PolyPhen,$HGVSc,$Distance) = @$CSQ;
      
      
    if ( $effects =~ /stream/ || 
	 $effects =~ /intron/ || 
	 $effects =~ /UTR/  ) {
      next;
    }
    
    next if ( ! $HGNC || ! $RefSeq );
    
#    print "$effects $HGNC $RefSeq $HGVSc\n";
    
      
    $interesting_variant++;
  }
  
  next if ( ! $interesting_variant );


  my %GTs;

  foreach my $gtype ( keys %{$$entry{ gtypes }} ) {
    next if ( $$entry{ gtypes }{ $gtype}{GT } eq ".");

    $GTs{ $$entry{ gtypes }{ $gtype}{GT } }++;

  }



  foreach my $gt ( keys %GTs ) {
    $counts{ $GTs{ $gt } }++;

    if (  $GTs{ $gt } == 1) {
      print "$$entry{CHROM}:$$entry{POS} $gt\n";
#      print Dumper (\%GTs);
    }
  }


#  print Dumper( \%GTs );

#  die Dumper( $entry );

  last if ( ! $exit_counter--);
  


}

foreach my $count ( sort { $a <=> $b } keys %counts ) {
  print "$count:$counts{$count}\t";
}

print "\n";

#print Dumper( \%counts );
