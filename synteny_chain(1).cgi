#!/usr/local/bin/perl

 use strict;
 use warnings;

 # This program will take the output from synteny_combine.cgi and attempt to create superscaffolds from the pairwise
 # data of the output file from synteny_combine.cgi (synteny.out).
 #
 # The output is a table of putative target super-scaffolds together with their support drawn from consensus across the reference taxa.
 #
 # Author: Arthur Georges, georges@aerg.canberra.edu.au
 #

#######################################################################
# CONFIGURATION

  my $script = 'synteny_chain.cgi';
  my $version = '21-May-2019';
  my $out_ext = '.out';
  my $log_ext = '.log';
# Input files
  my $input = 'synteny.out';
  my $slfile = 'scaffold_lengths.txt';
# Taxon list
  my $ref1_flag = 1;
  my $ref2_flag = 1;
  my $condition;
  my $trailer;
# Create expression for selecting reference taxa
      if ($ref1_flag == 1) {
        $condition = $condition.' && $ref1 ne \'NA\'';
        $trailer = $trailer.'_ref1';
      }
      if ($ref2_flag == 1) {
        $condition = $condition.' && $ref2 ne \'NA\'';
        $trailer = $trailer.'_ref2';
      }
      $condition =~ s/ && //;
      $trailer =~ s/_//;

# Output files
  my $output = 'synteny_chain_'."$trailer"."$out_ext";
  my $logfile = 'synteny_chain_'."$trailer"."$log_ext";

  print "Running $script\n";
  print "Version $version\n";

  print "Considering gene synteny across the following taxa:\n";
  print "    $condition\n\n";

 
#######################################################################

# INTIALIZE SCALARS AND ARRAYS

  my $linecount;
  my $used_count;
  my ($scaffold, $scaffold1, $scaffold2);
  my ($key, $value, $ref1, $ref2);
  my @stack;
  my $stacksize;
  my $a;
  my %results;
  my $direction;
  my $i;
  my $j;
  my $flag1;
  my $flag2;
  my $SS;
  my $SS_length;
  my %length_hash;
  my $length;

 #######################################################################
 # OPEN THE OUTPUT FILE
  open (OUT, '>'."$output");
  print "Output file $output opened \n";

 # OPEN THE LOG FILE
  open (LOG, '>'."$logfile");
  print "Log file $logfile opened \n";
  print LOG "SUMMARY OF ANALYSIS TO CHAIN INTERIM 2-SCAFFOLDS IN TARGET INTO SUPERSCAFFOLDS \n\n";

 # OPEN THE INPUT FILES
  open (IN, $input) || die "File $input not found\n";
  print "Input file $input opened \n";
  print LOG "Input file $input opened \n";

  open (SL, $slfile) || die "File $slfile not found\n";
  print "Input file $slfile opened, contains scaffold lengths \n";
  print LOG "Input file $slfile opened, contains scaffold lengths \n";

 # Read the scaffold lengths into a hash

  while (<SL>) {
    ($scaffold, $length) = split(' ', $_);
    #$scaffold =~ s/NW_//;
	#$scaffold =~ s/.1//;
    $length_hash{$scaffold} = $length;
  }

# READ THE JOINED SCAFFOLDS INTO A HASH TREE
# ASSUME THAT EACH SCAFFOLD HAS A RIGHT AND A LEFT SEQUENCE (WHICH MAY BE MISSING)

  my %left;
  my %right;

# PARSE THROUGH THE INPUT FILE, ASSIGNING JOINED SCAFFOLDS TO LEFT AND RIGHT
# WHERE THEY EXIST. IF LEFT IS OCCUPIED, POPULATE RIGHT. IF LEFT AND RIGHT
# ARE OCCUPIED, FLAG THE CONFLICT.

# USE ONLY DATA SUPPORTED BY ALL TAXA IN THE FIRST INSTANCE

   $linecount=0;
   $used_count=0;
   while (<IN>){
    $linecount = $linecount + 1;
    # Split the line into variables
      chomp($_);
      ($scaffold, $ref1, $ref2) = split(", ", $_);
      ($a, $scaffold, $direction) = split("__", $scaffold);
      ($scaffold1, $scaffold2) = split("-", $scaffold, 2);
    # Assign the scaffolds to left and right
    # Select records with desired support
       if ( eval $condition ) {
#      if ($gal ne 'NA'
#           && $ref1 ne 'NA'
#           && $ref2 ne 'NA'
#          ) {
        $used_count++;
        if (!$left{$scaffold1}) {
          $left{$scaffold1} = $scaffold2;
        } elsif (!$right{$scaffold1}) {
          $right{$scaffold1} = $scaffold2;
        } else {
          print "CONFLICT: $scaffold2 not placed -- Left and right of $scaffold1 already occupied by scaffolds\n";
        }
        if (!$left{$scaffold2}) {
          $left{$scaffold2} = $scaffold1;
        } elsif (!$right{$scaffold2}) {
          $right{$scaffold2} = $scaffold1;
        } else {
          print "CONFLICT: $scaffold1 not placed -- Left and right of $scaffold2 already occupied by scaffolds\n";
        }
      }
    }
    print "$linecount records read in from $input\n";
    print "$used_count records showing synteny across selected taxa\n";
    close (IN);
    open (IN, $input) || die "File $input not found\n";

# PRINT OUT THE RESULTS
   while (<IN>){
    # Split the line into variables
      chomp($_);
      ($scaffold, $ref1, $ref2) = split(", ", $_);
      ($a, $scaffold, $direction) = split("__", $scaffold);
      ($scaffold1, $scaffold2) = split("-", $scaffold, 2);
    # Select records with desired support
       if ( eval $condition ) {
#     if ($gal ne 'NA'
#           && $ref1 ne 'NA'
#           && $ref2 ne 'NA'
#          ) {
    # Set up an initial stack
      @stack = ();
      $stack[0] = $scaffold1;
      $stack[1] = $scaffold2;
    # Add to the stack, first extending the bottom, then the top
      $flag1 = 0;
      $flag2 = 0;
      do {
        $stacksize = @stack;
        if ( $right{$stack[$stacksize-1]} && $right{$stack[$stacksize-1]} ne $stack[$stacksize-2]){
          print LOG "Pushed scaffold $right{$stack[$stacksize-1]} onto the bottom of the stack\n";
          push ( @stack, $right{$stack[$stacksize-1]} );
        } elsif ( $left{$stack[$stacksize-1]} && $left{$stack[$stacksize-1]} ne $stack[$stacksize-2] ) {
          print LOG "Pushed scaffold $left{$stack[$stacksize-1]} onto the bottom of the stack\n";
          push ( @stack, $left{$stack[$stacksize-1]} );
        } else {
          $flag1 = 1;
        }
        if ( $right{$stack[0]} && $right{$stack[0]} ne $stack[1]){
          print LOG "Pushed scaffold $right{$stack[0]} onto the top of the stack\n";
          unshift ( @stack, $right{$stack[0]} );
        } elsif ( $left{$stack[0]} && $left{$stack[0]} ne $stack[1] ) {
          print LOG "Pushed scaffold $left{$stack[0]} onto the top of the stack\n";
          unshift ( @stack, $left{$stack[0]} );
        } else {
          $flag2 = 1;
        }
      } while ($flag1*$flag2 == 0);
      # print "\n";

    # Push the results into a hash
    # First, order them to remove mirror images, chain them, record length
      $stacksize = @stack;
      $SS_length = 0;
      $SS = " ";
      if ($stack[$stacksize-1] gt $stack[0]) {
        $j = 0;
        while ($j < $stacksize) {
          $SS = "$SS"."-"."$stack[$j]";
          $SS_length = $SS_length + $length_hash{$stack[$j]};
          $j = $j + 1;
        }
      } else {
        $j = $stacksize-1;
        while ($j > -1) {
          $SS = "$SS"."-"."$stack[$j]";
          $SS_length = $SS_length + $length_hash{$stack[$j]};
          $j = $j - 1;
        }
      }

    # Then remove niusense characters
      $SS =~ s/ -0/ -/;
      $SS =~ s/ -0/ -/;
      $SS =~ s/ -0/ -/;
      $SS =~ s/ -0/ -/;
      $SS =~ s/ -0/ -/;
      $SS =~ s/ -//;
      $SS =~ s/-0/-/g;
      $SS =~ s/-0/-/g;
      $SS =~ s/-0/-/g;
      $SS =~ s/-0/-/g;
      $SS =~ s/-0/-/g;
      $SS =~ s/ //g;

    # Then shove it in the hash, unique records only retained
      $results{$SS} = $SS_length;
    }
    }

 # FOR EACH SUPER-SCAFFOLD IN THE COMBINED HASH, PRINT

   while ( ($key,$value) = each %results ) {
     print OUT "$value $key \n";
   }

 close (SL);
 print "\nFile of scaffold lengths closed \n";
 close (IN);
 print "Input file synteny.out closed \n";
 close (OUT);
 print "Output file $output closed \n";
 close (LOG);
 print "Log file $logfile closed \n\n";

 exit;
