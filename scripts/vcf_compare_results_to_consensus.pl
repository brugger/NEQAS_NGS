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

my $vcf_file1     = shift;
my $vcf_reference = shift;


my $vcf = Vcf->new(file=>$vcf_file1);
$vcf->parse_header();

my @columns = $vcf->get_samples();

my ($agree, $disagree, $unknown) = (0,0,0);

while (my $entry = $vcf->next_data_hash()) {

  foreach my $column ( @columns ) {
    $$entry{ gtypes }{ $column }{GT } = "0/1" if ($$entry{ gtypes }{ $column }{ GT } eq "1/0");

    $$entry{ CHROM } =~ s/chr//;
    
    ref_data( $$entry{ CHROM }, $$entry{ POS}, $$entry{gtypes}{$column}{GT});
    #print "$column: $$entry{gtypes}{$column}{GT}\t";
  }



#  exit;
}


print "score Agree: $agree Disagree: $disagree Unknown $unknown\n";




# 
# 
# 
# Kim Brugger (18 Jun 2014)
sub ref_data {
  my ( $chr, $pos, $GT) = @_;

#  open(my $i, "tabix /software/packages/easih-toolbox/local_dbsnp/00-All.vcf.gz $chr:$pos-$pos | ");
#  print "tabix $vcf_reference $chr:$pos-$pos \n";
  open(my $i, "tabix $vcf_reference $chr:$pos-$pos | ");
  my $l = <$i>;
  if ( ! $l ) {
#    print "$chr:$pos\t$GT\tUnknown\n";
    $unknown++;
    return;
  }

  chomp $l;
  my @F = split("\t", $l);
  my ($similar, $diff) = (0,0);
  my $ref_GTs = "";
  for(my $i=9;$i<@F;$i++) {
    my ( $rGT,$rQC) = split(":", $F[$i]);
    $rGT = "0/1" if ( $rGT eq "1/0");
    next if ( $rGT =~ /\./);
    if ( $GT eq $rGT ) {
      $similar++;
    }
    else { 
      $ref_GTs .= "$rGT ";
      $diff++;
    }
#    print "$F[$i]\n";

  }

  

  if ( $similar > $diff ) {
    $agree++;
  }
  elsif ($diff == 1 || $similar == $diff ) {
    $unknown++;
  }
  else { 
    print "$chr:$pos\t$GT\t+$similar\t-$diff $ref_GTs\n";
    $disagree++;
  }


}
