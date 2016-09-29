###############################################################################
# readMapping.pl
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This is a wrapper funtion to map metagenomics reads to reference genomes
# using BWA. It was developed to map reads from OMD-TOIL to our SAGs and MAGs.
################################################################################
#!/usr/bin/perl
use strict;
use warnings;

# Import packages
use File::Path qw(make_path);

################################################################################
### User-defined files and folder structure
################################################################################

my $mtFolder = '../../data/metatranscriptomes';
my $genomeFolder = '../../data/refGenomes/concat';
my $mapFolder = '../../data/mapping/bamFiles';
my $numThreads = 1;
my $numChildren = 10;

# Check if the output directory exists and create it if necessary
if (-d $mapFolder) {
  print "Output directory exists \n";
}
# If not, create it
else {
  print "Creating output directory \n";
  make_path($mapFolder);
}

my @genomeList = glob($genomeFolder.'/*.fna');
my @sampleList = glob($mtFolder.'/*.fastq');

################################################################################
### Step 2: Simultaneously Index and Map Metatranscriptomes to Reference Genomes
################################################################################

my @commandArray;

foreach my $samplePath (@sampleList) {
  foreach my $genomePath (@genomeList) {
    if ($genomePath =~ /.+\/(.+).fna/) {
      my $genome = $1;
      if ($samplePath =~ /.+\/(.+).fastq/) {
	my $sample = $1;
	system("bbmap.sh t=".$numThreads." in=".$mtFolder."/".$sample.".fastq outm=".$mapFolder."/".$sample."-".$genome.".sam ref=".$genomeFolder."/".$genome.".fna ambig=all minid=0.85 nodisk sam=1.3");
      }
    }
  }
}
