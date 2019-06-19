#!/usr/local/bin/perl

 use strict;
 use warnings;

 # This program will take the output from synteny_pairwise.cgi for reference taxa (flanking taxa)
 # against the target taxon.
 #
 # The output is a table of putative junctions between target taxon scaffolds together 
 # with the support for these junctions drawn from consensus across the reference taxa.
 #
 # Briefly, two target taxon scaffolds are joined to form a 2-scaffold when they align 
 # consecutively to the reference genome AND when they abut, end to end.
 #
 # Manual inspection may be required to add junctions that are obvious by eye, but missed by the strict algoritm applied
 # here.
 #
 # Author: Arthur Georges, georges@aerg.canberra.edu.au
 #
#######################################################################
# CONFIGURATION

  my $script = 'synteny_combine.cgi';
  my $version = '21-May-2019';
  my $out_ext = '.out';
  my $log_ext = '.log';
# Input files
  my $ref1_tablefile = 'MONDO.out';
  my $ref2_tablefile = 'PHACI.out';
# Output file
  my $output = 'synteny.out';
  my $logfile = 'synteny.log';

  print "Running $script\n";
  print "Version $version\n";

#######################################################################

# INTIALIZE SCALARS AND ARRAYS

  my $linecount;
  my ($scaffold, $scaffold1, $scaffold2);
  my $reference;
  my $order;
  my %master_hash;
  my %ref1_hash;
  my %ref2_hash;
  my $size;
  my ($key, $ref1, $ref2);

 #######################################################################
 # OPEN THE OUTPUT FILE
  open (OUT, '>'."$output");
  print "Output file $output opened \n";

 # OPEN THE LOG FILE
  open (LOG, '>'."$logfile");
  print LOG "SUMMARY OF ANALYSIS TO LINK SCAFFOLDS IN target BY PAIRWISE CONSENSUS \n\n";

 # OPEN THE INPUT FILES
  open (WAL, $ref1_tablefile) || die "File $ref1_tablefile not found\n";
  print "tablefile file $ref1_tablefile opened \n";
  print LOG "tablefile file $ref1_tablefile opened \n";
  print "NOTE: Assumed that the input file is sorted on reference chromosome ID and then gene location \n";
  open (OPO, $ref2_tablefile) || die "File $ref2_tablefile not found\n";
  print "tablefile file $ref2_tablefile opened \n";
  print LOG "tablefile file $ref2_tablefile opened \n\n";

# READ THE JOINED SCAFFOLDS INTO A HASH FOR EACH REFERENCE SPECIES

   $linecount=0;
   while (<WAL>){
    $linecount = $linecount + 1;
    # Split the line into variables
      ($scaffold1, $scaffold2, $reference, $order) = split(" ", $_);
    # Add to the hash
      $scaffold = 'SS__'.$scaffold1.'-'.$scaffold2.'__'.$order;
      $scaffold =~ s/scf//g;
      if ($scaffold1 gt $scaffold2) {
        $order = reverse($order);
        $scaffold = 'SS__'.$scaffold2.'-'.$scaffold1.'__'.$order;
        $scaffold =~ s/scf//g;
      }
      $ref1_hash{$scaffold} = $reference;
      $master_hash{$scaffold} = 1;
   }
   print "$linecount records from $ref1_tablefile processed\n";
   print LOG "$linecount records from $ref1_tablefile processed\n";
   $size = keys %master_hash;
   print "  $size records in the combined hash tables\n";

   $linecount=0;
   while (<OPO>){
    $linecount = $linecount + 1;
    # Split the line into variables
      ($scaffold1, $scaffold2, $reference, $order) = split(" ", $_);
    # Add to the hash
      $scaffold = 'SS__'.$scaffold1.'-'.$scaffold2.'__'.$order;
      $scaffold =~ s/scf//g;
      if ($scaffold1 gt $scaffold2) {
        $order = reverse($order);
        $scaffold = 'SS__'.$scaffold2.'-'.$scaffold1.'__'.$order;
        $scaffold =~ s/scf//g;
      }
      $ref2_hash{$scaffold} = $reference;
      $master_hash{$scaffold} = 1;
   }
   print "$linecount records from $ref2_tablefile processed\n";
   print LOG "$linecount records from $ref2_tablefile processed\n";
   $size = keys %master_hash;
   print "  $size records in the combined hash tables\n";

 # FOR EACH SCAFFOLD IN THE COMBINED HASH, CHECK FOR VALUES IN THE REFERENCE HASHES

   foreach my $key ( keys %master_hash ){
     $ref1='NA'; $ref2='NA'; 
     if ($ref1_hash{$key}) { $ref1 = $ref1_hash{$key} };
     if ($ref2_hash{$key}) { $ref2 = $ref2_hash{$key} };
     print OUT "$key, $ref1, $ref2\n";
   }

 close (OUT);
 print "Putative super-scaffold name, and spanning scaffolds for";
 print "  ref1, ref2 written to output file $output\n";
 print "  NA for no overlap\n\n";
 print "outputfile file $output closed \n";
 print LOG "Putative super-scaffold name, and spanning scaffolds for";
 print LOG "  ref1, ref2 written to output file\n";
 print LOG "  NA for no overlap\n";
 print LOG "Output file $output closed \n\n";


 close (WAL);
 print "tablefile file $ref1_tablefile closed \n";
 close (OPO);
 print "tablefile file $ref2_tablefile closed \n";

 close (LOG);
 print "Log file $logfile closed \n\n";

 exit;
