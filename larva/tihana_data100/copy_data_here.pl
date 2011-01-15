#!/usr/bin/perl

use strict;

my $prefix = "/groups/zlatic/zlaticlab/Data/results";
my $listname = "FilesToCopy.txt";
my $dstdir = "train";
my $nFilesPerDir = 1;
my $DEBUG = 0;

# read file list
if(! -e $listname){
    print "File list $listname does not exist.\n";
    exit(1);
}
open(IN,"$listname") or die("Error: Could not open file $listname for reading\n");
open(OUT,">FilesCopied.txt");

# loop over lines in the input file
my $srcdir;
my $line;
while( $line = <IN> ){

    $line =~ s/^\s+//;
    $line =~ s/\s+$//;

    $srcdir = $prefix . "/" . $line;

    # check to make sure this is a directory
    if(! -e $srcdir){
        print "\n***WARNING: Input directory >>$srcdir<< does not exist. Skipping.***\n\n";
	next;
    }

    # get all outlines
    my @outlines = glob("$srcdir/*.outline");
    my $noutlines = @outlines;

    # choose $nFilesPerDir outlines from this directory
    my @samples = ();
    my $nsamples = 0;
    if($noutlines < $nFilesPerDir){
	print "\n\n***WARNING: Input directory $srcdir only contains $noutlines < $nFilesPerDir outlines. Only choosing $noutlines files.***\n\n";
	for(my $i = 0; $i < $noutlines; $i++){
	    push(@samples,$i);
	    $nsamples++;
	}
    }
    else{

	while(1){
	    my $randidx = int(rand($noutlines));
	    my $isselected = 0;
	    foreach my $sample(@samples){
		if($randidx == $sample){
		    $isselected = 1;
		    last;
		}
	    }
	    if(!$isselected){
		push(@samples,$randidx);
		$nsamples++;
	    }
	    if($nsamples >= $nFilesPerDir){
		last;
	    }
	}
    }

    # copy these files -- assuming that the name will be unique
    foreach my $sample(@samples){
	my $srcfile = $outlines[$sample];
	my $cmd = "cp \"$srcfile\" \"$dstdir\"/.";
	print OUT "$srcfile\n";
	if($DEBUG){
	    print "$cmd\n";
	}
	`$cmd`;
    }
}

close(IN);
close(OUT);
