#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (18 Jun 2014), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib "/software/packages/vcftools_0.1.12a/lib/perl5/site_perl";
use Vcf;



my $vcf_file = shift;
my $vcf = Vcf->new(file=>$vcf_file);
$vcf->parse_header();

#print $vcf->format_header();

my $vcf_out = Vcf->new( version=>"4.0" );
add_vcf_out_header( $vcf_out );

my @columns = $vcf->get_samples();

$vcf_out->add_columns(  @columns );


#$vcf_out->add_columns('test');
print $vcf_out->format_header();

while (my $entry = $vcf->next_data_hash()) {

  next if ( $$entry{CHROM} =~ /M/);

  $$entry{CHROM} =~ s/chr//;


  my $filter = join(",", @{$$entry{FILTER}} );
  next if ( $filter !~ /PASS/);

  my $alt = join(",", @{$$entry{ALT}});
  next if ( $alt !~/^[ACGTN]+$/ );

#  print Dumper( $entry );
  my @line = ($$entry{CHROM}, $$entry{POS}, $$entry{ID}, $$entry{REF}, $alt, 1000, $filter, "","GT:GQ", );
  my @sample_GTs;
  foreach my $column ( @columns ) {
#    print Dumper( $column);
    $$entry{ gtypes }{ $column }{GT } = "0/1" if ($$entry{ gtypes }{ $column}{GT } eq "1/0");

    push @sample_GTs, "$$entry{gtypes}{$column}{GT}:".int($$entry{QUAL});
  }
  push @line, join("\t", @sample_GTs);

#  print Dumper( \@line );

  print $vcf_out->format_line(\@line);

#  exit;


}


# 
# 
# 
# Kim Brugger (18 Jun 2014)
sub add_vcf_out_header {
  my ( $vcf ) = @_;

  $vcf_out->add_header_line({key=>'FORMAT', ID=>'GT',Number=>-1,Type=>'String',Description=>'Genotypes'});
  $vcf_out->add_header_line({key=>'FORMAT', ID=>'GQ',Number=>-1,Type=>'Integer',Description=>'Genotype quality'});
  
}
