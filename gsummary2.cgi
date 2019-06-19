#!/usr/local/bin/perl

 use strict;
 use warnings;

 #######################################################################
 # USER INPUT

 # Take file and threshold
   my $tablefile = $ARGV[0];
   if (!$ARGV[0]){
    die "Useage: perl gsummary.cgi tablefile\n";
   }
   
#######################################################################

# Name the output file
  my ($m, $c) = split(/\./, $tablefile);
  my $out = join ('', "$m", '.out');
  my $log = join ('', "$m", '.log');
  print "Progress will be logged to $log\n";
  open (LOG, '>'."$log");
  print LOG "GENE SYNTENY SUMMARY RUN\n\n";
  print "Program will summarize attributes of the tablular data to be input to the synteny analysis from $tablefile \n";

# Open the input and output files
  open (M, $tablefile) || die "File $tablefile not found\n";
  print "tablefile file $tablefile opened \n";
  print "NOTE: Assumed that the input file is sorted on Gallus chromosome ID and then gene location \n";
  print LOG "tablefile file $tablefile opened \n";
  open (OUT, '>'."$out");
  print "Output file $out opened \n";
  print LOG "Output file $out opened \n";
  
# Initialise scalars and arrays

  my (@gal_gene_id, @gal_chr_id, @gal_gene_order, @gal_strand, @gal_start, @gal_end, @pog_gene_id, @pog_chr_id, @pog_gene_order, @pog_gene_reverse, @pog_strand, @pog_start, @pog_end, @gal_align_rate, @pog_align_rate, @percent_identity, @synteny_inf, @gal_gene, @pog_gene );

  my $i;
  my $linecount;

 # Read in records from the table and decompose to arrays

   $linecount=0;
   while (<M>){
    $linecount = $linecount + 1;
    $i = $linecount - 2;  # Perl array counter starts at 0
    if ($linecount > 1) {    # Skip header
      # Split the line into variables
        ($gal_gene_id[$i], $gal_chr_id[$i], $gal_gene_order[$i], $gal_strand[$i], $gal_start[$i], $gal_end[$i], $pog_gene_id[$i], $pog_chr_id[$i], $pog_gene_order[$i], $pog_gene_reverse[$i], $pog_strand[$i], $pog_start[$i], $pog_end[$i], $gal_align_rate[$i], $pog_align_rate[$i], $percent_identity[$i], $synteny_inf[$i], $gal_gene[$i], $pog_gene[$i]) = split("\t", $_);
    }
   }
   print "$linecount records from $tablefile processed\n";
   print "Header line discarded\n";
   my $records = $linecount - 1;

# INSERT HERE A CHECK OF THE ORDER OF THE INPUT FILE

# Calculate the number of unique scaffolds (use hash feature, unique key);

  my %gal_scaffold_count;
  my %pog_scaffold_count;
  my $a;
  my $b;
  my $gal_size; 
  my $pog_size;

  $i = 0;
  while ($i < $records-1) {
    ($a,$b) = split(/_/,$gal_chr_id[$i]);
    if (!$b || $b ne 'random'){
      $gal_scaffold_count{$gal_chr_id[$i]} = $i;
    }
    if ($pog_chr_id[$i] ne 'NA') {
      $pog_scaffold_count{$pog_chr_id[$i]} = $i;
    }
    $i = $i +1;
  }
  $gal_size = scalar(keys %gal_scaffold_count);
  print "$gal_size well defined scaffolds in table for reference genome\n";
  $pog_size = scalar(keys %pog_scaffold_count);
  print "$pog_size well defined scaffolds in table for Pogona genome\n";

# Calculate the number of runs in each of reference and Pogona

  $i = 0;
  my $pog_hold = 'xxx';
  my $pog_count = 0;
  my $gal_hold = 'xxx';
  my $gal_count = 0;
  while ($i < $records-1) {
    if ($pog_chr_id[$i] ne 'NA' && $pog_chr_id[$i] ne $pog_hold) {
      $pog_count= $pog_count + 1;
      $pog_hold = $pog_chr_id[$i];
    }
    ($a,$b) = split(/_/,$gal_chr_id[$i]);
    if (!$b || $b ne 'random'){
      if ($gal_chr_id[$i] ne 'NA' && $gal_chr_id[$i] ne $gal_hold) {
        $gal_count= $gal_count + 1;
        $gal_hold = $gal_chr_id[$i];
      }
    } 
    $i = $i +1;
  }

  print "$gal_count gene runs in table for reference genome\n";
  print "$pog_count gene runs in table for Pogona genome\n";
  if ($gal_count != $gal_size) {
    print "Warning: File may not be ordered on Reference scaffolds and gene order";
  }

# Remove all records with NA
  my @pog_id;
  my @gal_id;
  my @pog_order;
  my @pog_reverse;
  $i = 0;
  while ($i < $records-1) {
    if ($pog_chr_id[$i] eq 'NA') {
      $gal_chr_id[$i] = 'NA';
      $pog_gene_order[$i] = 'NA';
      $pog_gene_reverse[$i] = 'NA';
    }
  $i = $i + 1;
  }
  @gal_id = grep { $_ ne 'NA' } @gal_chr_id;
  @pog_id = grep { $_ ne 'NA' } @pog_chr_id;
  @pog_order = grep { $_ ne 'NA' } @pog_gene_order;
  @pog_reverse = grep { $_ ne 'NA' } @pog_gene_reverse;

  $records = @gal_id;

# Parse through the file to identify where scaffolds abut (at the gene level of resolution) with homology to a contiguous region in the reference
  my $scafcount;

  $i = 1; # start with the second line
  while ($i < $records-1) {
      if ($pog_id[$i] ne $pog_id[$i-1]) {
      # Reset the counter
        $scafcount = 1;
      # Check to see if the scaffolds abut
        if (($pog_order[$i] == 1 && $pog_order[$i-1] == 1) ||
            ($pog_order[$i] == 1 && $pog_reverse[$i-1] == 1) ||
            ($pog_reverse[$i] == 1 && $pog_reverse[$i-1] == 1) ||
            ($pog_reverse[$i] == 1 && $pog_order[$i-1] == 1) &&
            ($gal_id[$i] eq $gal_id[$i-1])) {
          print "$pog_id[$i-1] \+ $pog_id[$i] \~ $gal_id[$i] \n";
          print OUT "$pog_id[$i-1] $pog_id[$i] $gal_id[$i] \n";
        }
      }
    $i = $i +1;
  }


#  print "\n", 'Line 2', $gal_gene_id[1], $gal_chr_id[1], $gal_gene_order[1], $gal_strand[1], $gal_start[1], $gal_end[1], $pog_gene_id[1], $pog_chr_id[1], $pog_gene_order[1], $pog_strand[1], $pog_start[1], $pog_end[1], $gal_align_rate[1], $pog_align_rate[1], $percent_identity[1], $synteny_inf[1], $gal_gene[1], $pog_gene[1];
#  print "\n", 'Line 3', $gal_gene_id[2];
#  print "\n", 'Line 1', $gal_gene_id[0];

 close (OUT);
 print "No records written to $out , closed \n";
 close (M);
 print "tablefile file $tablefile closed \n";

 exit;
