#!/usr/bin/perl 
# 
# 
# 
# 
# Kim Brugger (13 Jun 2014), contact: kim.brugger@easih.ac.uk

use strict;
use warnings;
use Data::Dumper;

use lib "/software/packages/vcftools_0.1.11/lib/perl5/site_perl";
use Vcf;

my $vcf_file1 = shift;
my $res = readin_vcf( $vcf_file1);

my ( $id, $non_id ) = (0,0);

my $vcf_file2 = shift;
my $vcf = Vcf->new(file=>$vcf_file2);
$vcf->parse_header();
while (my $entry = $vcf->next_data_hash()) {
  if ( $$res{"$$entry{CHROM}:$$entry{POS}"} ) {
    my $ref_alt = "$$entry{REF}:". join(",",@{$$entry{ALT}});
    if ( $$res{"$$entry{CHROM}:$$entry{POS}"} eq $ref_alt) {
#      print "$$entry{CHROM}:$$entry{POS} :: IDENTICAL\n";
      $id++;
    }
    else {
      print "$$entry{CHROM}:$$entry{POS} :: NON-IDENTICAL ". $$res{"$$entry{CHROM}:$$entry{POS}"} . " ne $ref_alt\n";
      $non_id++;
    }
  }
}

print "ID: $id, NON-ID: $non_id\n";



# 
# 
# 
# Kim Brugger (13 Jun 2014)
sub readin_vcf {
  my ($infile) = @_;
  
  my %res;

  my $vcf = Vcf->new(file=>$infile);
  $vcf->parse_header();
  while (my $entry = $vcf->next_data_hash()) {
    $res{"$$entry{CHROM}:$$entry{POS}"} = "$$entry{REF}:". join(",",@{$$entry{ALT}});
#    $res{"$$entry{CHROM}:$$entry{POS}"} = "$$entry{REF}:". join(",",@{$$entry{ALT}});
  }

  return \%res;
}

