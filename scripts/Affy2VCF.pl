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

my $vcf_out = Vcf->new( version=>"4.1");

$vcf_out->add_header_line({key=>'INFO', ID=>'AC',Number=>-1,Type=>'Integer',Description=>'Allele count in genotypes'});
$vcf_out->add_header_line({key=>'FORMAT', ID=>'GT',Number=>-1,Type=>'String',Description=>'Genotypes'});
$vcf_out->add_header_line({key=>'INFO', ID=>'FQ',Number=>-1,Type=>'Float',Description=>'Allele frequency'});
$vcf_out->add_columns('test');

print $vcf_out->format_header();

my $exit_counter = 200;

my $infile = shift;
open(my $in, $infile) || die "Either could not open infile '$infile': $!\n";
my $data = 0;
while(<$in>) {

  last if ( !$exit_counter--);
  
  if ( ! $data ) {
    $data = 1 if ( /\[data\]/i);
    next;
  }
  
  next if (/^SNP/);
  next if (!/rs\d+/);
  chomp;
  $_ =~ s/\r//;
  my @F = split("\t", $_);
  my ( $chr, $pos, $allele1, $allele2, $score, $freq ) = ( $F[5], $F[6], $F[2], $F[3], $F[7], $F[8]);

  next if ( !$chr );

  $chr = "X" if ($chr eq "XY");

  next if ( $allele1 !~ /[ACGT]/);

#  my $ref = ref_base("$chr:$pos-$pos");

  my ( $ref, $ALT, $minus_strand ) = dbSNP_info( "$chr:$pos-$pos" );
  if ( ! $ref ) {
    print STDERR "Uknown region $chr:$pos-$pos\n";
    next;
  }


  if ( $minus_strand ) {
    $allele1 = revDNA( $allele1);
    $allele2 = revDNA( $allele2);
  }
  print "$ref $ALT == $allele1 $allele2 $minus_strand $chr:$pos\n";
  next;

  my $alt = "";
  my $AC = 1;
  my $GT = "0/0";
  if ( $allele1 eq $allele2 && $allele1 eq $ref ) {
    $alt = $ref;
    next;
  }
  elsif ( $allele1 eq $allele2 ) {
    $alt = $allele1;
    $AC = 2;
    $GT = "1/1";
  }
  else {
    $AC = 1;
    $GT = "0/1";
    if ( $allele1 ne $ref ) {
      $alt = $allele1;
    }
    else {
      $alt = $allele2;
    }
  }

  my @line = [$chr, $pos, ".", $ref, $alt, $score, "PASS", "AC=$AC;FQ=$freq", "GT", $GT];

  print $vcf_out->format_line(@line);

#  die Dumper(\@F);

  
}




# 
# 
# 
# Kim Brugger (17 Mar 2014)
sub revDNA {
   my ($dna) = @_;

  $dna =~ tr/[ACGT]/[TGCA]/;
  $dna = reverse($dna);
  return $dna;

}



# 
# 
# 
# Kim Brugger (17 Mar 2014)
sub dbSNP_info {
  my ($position) = @_;

  open(my $i, "tabix /software/packages/easih-toolbox/local_dbsnp/00-All.vcf.gz $position | ");
  my $l = <$i>;
  close( $i );
  return if ( ! $l );
  chomp $l;
  my @F = split("\t", $l);

  my $reverse_SNP = 0;
  $reverse_SNP = 1 if ( $F[7] =~ /;RV;/);

  return ($F[3], $F[4], $reverse_SNP);
}




# 
# 
# 
# Kim Brugger (13 Mar 2014)
sub ref_base {
  my ( $position ) = @_;

  open (my $i, "/software/bin/samtools faidx /refs/human_1kg/human_g1k_v37.fasta $position | ");
  <$i>;
  my $base = <$i>;
  close( $i );

  return undef if ( ! $base );

  chomp $base;

  return $base;
  
}
