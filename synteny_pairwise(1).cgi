#!/usr/local/bin/perl

 use strict;
 use warnings;

# This program will take tabular data comparing a reference taxon (reference)with a target taxon (target).
#
# Useage: perl synteny_pairwise.cgi inputfile.txt
#
# The input files are
#
#   (a) A file cross referencing the target taxon scaffolds to those of the reference taxon.
#       Refer to a copy of the file for the format. Files are of the form synteny_Ref-Tgt.txt. They are tab delimited.
#   (b) A file containing the scaffold lengths. This file has a header line, then pairs of values
#       -- scaffold name and length, separated by a blank. Fixed file and filename of scaffold_lengths.txt.
#
# The program looks for target scaffolds that abut end to end (forward and reverse), face to face, back to back
#   based on synteny with the reference. These abutting scaffolds are thus putatively joined into an interim superscaffold,
#   and form part of the output.
#
# The output is a table of these interim 2-value superscaffolds, together with the name of the reference taxon scaffold that spans the
#   junction between the two, and the length of the new superscaffold. The output file is called <inputfile>.out .
#
# Manual inspection may be required to add junctions that are obvious by eye, but missed by the strict algoritm applied here.
#
# Author: Arthur Georges, georges@aerg.canberra.edu.au
#
#######################################################################
# CONFIGURATION

  my $script = 'synteny_pairwise.cgi';
  my $version = '21-May-2019';
  my $slfile = 'scaffold_lengths.txt';
  my $out_ext = '.out';
  my $log_ext = '.log';

  print "Running $script\n";
  print "Version $version\n";
  
 #######################################################################
 # USER INPUT

 # Take file
   my $inputfile = $ARGV[0];
   if (!$ARGV[0]){
    die "Useage: perl $script inputfile\n";
   }
   
#######################################################################
# I/O

# Name the output file
  my ($m, $c) = split(/\./, $inputfile);
  my $out = join ('', "$m", "$out_ext");
  my $log = join ('', "$m", "$log_ext");
  print "Progress will be logged to $log\n";
  open (LOG, '>'."$log");
  print LOG "GENE SYNTENY SUMMARY RUN\n\n";
  print LOG "$script $version\n";

# Open the input and output files

  open (IN, $inputfile) || die "File $inputfile not found\n";
  print "Input file $inputfile opened \n";
  print "NOTE: Assumed that the input file is sorted on reference taxon chromosome ID and then gene location, tab delimited \n";
  print LOG "Input file $inputfile opened \n";

  open (SL, $slfile) || die "File $slfile not found\n";
  print "Input file $slfile opened, contains scaffold lengths \n";
  print LOG "Input file $slfile opened, contains scaffold lengths \n";

  open (OUT, '>'."$out");
  print "Output file $out opened \n\n";
  print LOG "Output file $out opened \n\n";
  
#######################################################################
# INITIALIZE SCALERS, ARRAYS AND HASHES

# Input variables
  my (
      @ref_gene_id, @ref_chr_id, @ref_gene_order, @ref_strand, @ref_start, @ref_end, @target_gene_id, @target_chr_id,
      @target_gene_order, @target_gene_reverse, @target_strand, @target_start, @target_end
     );
# Loop index variables
  my $i;
# Counters
  my $linecount;
  my $c1;
  my $c2;
  my $c3;
  my $c4;
  my $c5;
  my $c6;
  my $c7;
  my $c8;
  my $c9;
  my $orphaned_count;
# Flags
  my ($flag1, $flag2, $flag3, $flag4);
  my $orphaned;
# Hashes
  my %max_hash;
  my %min_hash;
  my %gene_count;
  my %length_hash;
# Arrays
  my @target_id;
  my @ref_id;
  my @target_order;
  my @target_reverse;
# Scalars
  my $scaffold;
  my $length;
  my $order;
  my $genes;

#######################################################################
# DATA INPUT

# First the records provided by Hrdip

 # Read in records from the table and decompose to arrays

   $linecount=0;
   while (<IN>){
    $linecount++;
    $i = $linecount - 2;  # Perl array counter starts at 0
    if ($linecount > 1) {    # Skip header
      # Split the line into variables
        ($ref_gene_id[$i], $ref_chr_id[$i], $ref_gene_order[$i], $ref_strand[$i], $ref_start[$i], $ref_end[$i],
         $target_gene_id[$i], $target_chr_id[$i], $target_gene_order[$i], $target_gene_reverse[$i], $target_strand[$i],
         $target_start[$i], $target_end[$i]) = split("\t", $_);
    }
   }
   print "$linecount records from $inputfile processed\n";
   print "  Header line discarded\n";
   print LOG "$linecount records from $inputfile processed\n";
   print LOG "  Header line discarded\n";
   
  $linecount = $linecount - 1;

# Then the scaffold lengths into a hash

  while (<SL>) {
    ($scaffold, $length) = split(' ', $_);
    $length_hash{$scaffold} = $length;
  }

#######################################################################
# PRELIMINARY DATA MANIPULATION
  
# Remove all records with NA or unknown scaffolds
  $i = 0;
  while ($i < $linecount-1) {
#    if ($target_chr_id[$i] eq 'NA' || $ref_chr_id[$i] eq 'Un_random' || $ref_chr_id[$i] eq 'W_random' || $ref_chr_id[$i] eq 'Z_random') {
    if ($target_chr_id[$i] eq 'NA') {
      $ref_chr_id[$i] = 'NA';
      $target_gene_order[$i] = 'NA';
      $target_gene_reverse[$i] = 'NA';
    }
  $i++;
  }
  @ref_id = grep { $_ ne 'NA' } @ref_chr_id;
  @target_id = grep { $_ ne 'NA' } @target_chr_id;
  @target_order = grep { $_ ne 'NA' } @target_gene_order;
  @target_reverse = grep { $_ ne 'NA' } @target_gene_reverse;

  $linecount = @target_id;

# Calculate minimum and maximum gene order number for each Koala scaffold

  # Initialize the hashes
  $i = 0;
  while ($i < $linecount) {
    $max_hash{$target_id[$i]} = 0;
    $min_hash{$target_id[$i]} = 9999;
    $gene_count{$target_id[$i]} = 0;
    $i++;
  }
  # Calculate the minumum and maxiumum gene order score for each scaffold
  $i = 0;
  while ($i < $linecount) {
    $gene_count{$target_id[$i]}++;
    if (  $max_hash{$target_id[$i]} < $target_order[$i] ) {
      $max_hash{$target_id[$i]} = $target_order[$i];
    }
    if ( $min_hash{$target_id[$i]} > $target_order[$i] ) {
      $min_hash{$target_id[$i]} = $target_order[$i];
    }
    $i++;
  }

#######################################################################
# FIND ABUTTING SCAFFOLDS

# Parse through the file to identify where scaffolds abut (at the gene level of resolution)
# whn ordered against their homology to a contiguous region in the reference

  $c1=0;
  $c2=0;
  $c3=0;
  $c4=0;
  $c5=0;
  $c6=0;
  $c7=0;
  $c8=0;
  $c9=0;
  $orphaned_count=0;

  $i = 1; # start with the second line
  while ($i < $linecount-1) {
    if ($target_id[$i] ne $target_id[$i-1]) {
      # Reset the flags
        $flag1=0;
        $flag2=0;
        $flag3=0;
        $flag4=0;

      # Check for orphaned terminal genes, and skip if found
        $orphaned = 0;
        if ( defined($target_id[$i+1]) ) {

        # Looking for orphans. If a scaffold contains more than one gene, and if the focus gene is terminal, and if the flanking genes are different

          if ( ($max_hash{$target_id[$i]} > $min_hash{$target_id[$i]}) && (( $target_order[$i] == $max_hash{$target_id[$i]} ) || ( $target_order[$i] == $min_hash{$target_id[$i]} )) && ( ($target_id[$i] ne $target_id[$i+1]) ) ) {

#               if ($target_id[$i] eq 'scf000165') {print "\n\n/[ ID: $pog_id[$i-1] -- $target_id[$i] -- $target_id[$i+1]\n  Order: $min_hash{$target_id[$i-1]} - $max_hash{$target_id[$i-1]} : $min_hash{$target_id[$i]} - $max_hash{$target_id[$i]} : $min_hash{$target_id[$i+1]} - $max_hash{$target_id[$i+1]} /] \n\n";}

               $orphaned = 1;
               $orphaned_count++;
          }
        }
        if ( $orphaned == 1 ) {
#         print "Skipping orphaned gene: $target_id[$i] $target_order[$i] [ $min_hash{$target_id[$i]} - $max_hash{$target_id[$i]} ] \n";
         print LOG "Skipping orphaned gene: $target_id[$i] $target_order[$i] [ $min_hash{$target_id[$i]} - $max_hash{$target_id[$i]} ] \n";
         $i++;
         next;
        }
        
      # Check to see if the scaffolds abut, and flag if they do
      # Back to back
      # If the gene order number at the junction is the minimum for both scaffold 1 and 2 ...
        if (($target_order[$i] == $min_hash{$target_id[$i]} && $target_order[$i-1] == $min_hash{$target_id[$i-1]}) && ($ref_id[$i] eq $ref_id[$i-1])) {
          $flag1 =  1;
        }
      # End to end
        if (($target_order[$i] == $min_hash{$target_id[$i]} && $target_order[$i-1] == $max_hash{$target_id[$i-1]}) && ($ref_id[$i] eq $ref_id[$i-1])) {
          $flag2 = 1;
        }
      # Reverse end to end
        if (($target_order[$i] == $max_hash{$target_id[$i]} && $target_order[$i-1] == $min_hash{$target_id[$i-1]}) && ($ref_id[$i] eq $ref_id[$i-1])) {
          $flag3 = 1;
        }
      # Face to face
        if (($target_order[$i] == $max_hash{$target_id[$i]} && $target_order[$i-1] == $max_hash{$target_id[$i-1]}) && ($ref_id[$i] eq $ref_id[$i-1])) {
          $flag4 = 1;
        }
      # When 2 or more flags are raised, it means that one of the scaffolds has only one gene,
      # and so its beginning and end are both listed at gene order 1. The combinations are
      #
      # One flag raised
      #   flag 1, back to back, -- [N ... 1][1 ... N]
      #   flag 2, end to end, +-   [1 ... N][1 ... N]
      #   flag 3, end to end, reverse order -+ [N ... 1][N ... 1]
      #   flag 4, face to face, ++ [1 ... N][N ... 1]
      # Two flags raised
      #   flags 1 and 2, first scaffold with one gene only, 0-
      #     [1][1 ... N] ~ --
      #     [N][1 ... N] ~ +-
      #   flags 2 and 4, second scaffold with one gene only, +0
      #     [1 ... N][1] ~ +-
      #     [1 ... N][N] ~ ++
      #   flags 3 and 4, first scaffold with one gene only, 0+
      #     [1][N ... 1] ~ -+
      #     [N][N ... 1] ~ ++
      #   flags 1 and 3, second scaffold with one gene only, -0
      #     [N ... 1][1] ~ --
      #     [N ... 1][N] ~ -+
      # Three flags raised
      #   Error
      # All four flags raised
      #   Both first and second scaffold with one gene only, 00
      #

        $order = '  ';
      # Back to back
        if ($flag1 == 1 && $flag2 == 0 && $flag3 == 0 && $flag4 == 0) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = '--';
            $c1++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }
      # End to end
        if ($flag1 == 0 && $flag2 == 1 && $flag3 == 0 && $flag4 == 0) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = '+-';
            $c2++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }
      # Reverse end to end
        if ($flag1 == 0 && $flag2 == 0 && $flag3 == 1 && $flag4 == 0) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = '-+';
            $c3++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }
      # Face to face
        if ($flag1 == 0 && $flag2 == 0 && $flag3 == 0 && $flag4 == 1) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = '++';
            $c4++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }

        if ($flag1 == 1 && $flag2 == 1 && $flag3 == 0 && $flag4 == 0) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = 'O-';
            $c5++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }
        if ($flag1 == 0 && $flag2 == 1 && $flag3 == 0 && $flag4 == 1) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = '+O';
            $c6++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }
        if ($flag1 == 0 && $flag2 == 0 && $flag3 == 1 && $flag4 == 1) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = 'O+';
            $c7++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }
        if ($flag1 == 1 && $flag2 == 0 && $flag3 == 1 && $flag4 == 0) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = '-O';
            $c8++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }

        if ($flag1 == 1 && $flag2 == 1 && $flag3 == 1 && $flag4 == 1) {
        # Take out those scaffolds that are linked by an unassigned scaffold in the reference genome
          if ($ref_id[$i] ne 'Un_random' && $ref_id[$i] ne 'W_random' && $ref_id[$i] ne 'Z_random') {
            $order = 'OO';
            $c9++;
            $length = $length_hash{$target_id[$i-1]} + $length_hash{$target_id[$i]};
            $genes = $max_hash{$target_id[$i-1]} + $max_hash{$target_id[$i]};
#            print "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print LOG "$target_id[$i-1] \+ $target_id[$i] \~ $ref_id[$i] $order min = $min_hash{$target_id[$i-1]}, $min_hash{$target_id[$i]} max = $max_hash{$target_id[$i-1]}, $max_hash{$target_id[$i]} genecount = $gene_count{$target_id[$i-1]}, $gene_count{$target_id[$i]}\n";
            print OUT "$target_id[$i-1] $target_id[$i] $ref_id[$i] $order $length $genes\n";
          }
        }

        if ($flag1 + $flag2 + $flag3 + $flag4 == 3) { print "ERROR: Junction between $target_id[$i-1] and $target_id[$i] confused\n" };
      }
    $i++;
  }
  print "\nPutative junctions based on $orphaned_count orphaned terminal genes excluded\n\n";
  if ($c3 > 0 ) { print "$c3 scaffold pairs joined end to end\n"; }
  if ($c2 > 0 ) { print "$c2 scaffold pairs joined, end to end, but in reverse order\n"; }
  if ($c1 > 0 ) { print "$c1 scaffold pairs joined back to back\n"; }
  if ($c4 > 0 ) { print "$c4 scaffold pairs joined face to face\n"; }
  if ($c5 > 0 ) { print "$c5 scaffold 1 with one gene only, joined to start of scaffold 2\n"; }
  if ($c6 > 0 ) { print "$c6 scaffold 2 with one gene only, joined to end of scaffold 1\n"; }
  if ($c7 > 0 ) { print "$c7 scaffold 1 with one gene only, joined to end of scaffold 2\n"; }
  if ($c8 > 0 ) { print "$c8 scaffold 2 with one gene only, joined to start of scaffold 1\n"; }
  if ($c9 > 0 ) { print "$c9 scaffold 1 and 2 with one gene only, joined\n"; }

  print LOG "\nPutative junctions based on $orphaned_count orphaned terminal genes excluded\n\n";
  print LOG "$c3 scaffold pairs joined end to end\n";
  print LOG "$c2 scaffold pairs joined, end to end, but in reverse order\n";
  print LOG "$c1 scaffold pairs joined back to back\n";
  print LOG "$c4 scaffold pairs joined face to face\n";
  print LOG "$c5 scaffold 1 with one gene only, joined to start of scaffold 2\n";
  print LOG "$c6 scaffold 2 with one gene only, joined to end of scaffold 1\n";
  print LOG "$c7 scaffold 1 with one gene only, joined to end of scaffold 2\n";
  print LOG "$c8 scaffold 2 with one gene only, joined to start of scaffold 1\n";
  print LOG "$c9 scaffold 1 and 2 with one gene only, joined\n";

  my $count = $c1 + $c2 + $c3 + $c4 + $c5 + $c6 + $c7 + $c8 + $c9;

 close (OUT);
 print "\n$count joined pairs of scaffolds written to $out , closed \n";
 print LOG "\n$count records written to $out , closed \n";
 close (IN);
 print "Input file $inputfile closed \n";
 close (IN);
 print "Scaffold lengths file $slfile closed \n";
 print LOG "Scaffold lengths file $slfile closed \n";
 close (LOG);
 print "Log file closed \n";


 exit;
