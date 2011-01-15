#!/usr/bin/perl

use strict;
use File::Spec;

my $numArgs = $#ARGV + 1;
if($numArgs < 1){
    print "Usage: check_train.pl <datadir>\n";
    print "numArgs = $numArgs\n";
    exit 1;
}

my $datadir = $ARGV[0];
my $trainfile = File::Spec->catfile($datadir,"train.txt");
if(! -e $trainfile){
    print "Train file $trainfile does not exist.\n";
    exit 1;
}
print "Checking train file $trainfile...\n";

open(FILE,$trainfile);

while(my $line = <FILE>){
    chomp($line);
    if(length($line) == 0){
	next;
    }
    if(-e $line){
	print "GOOD: $line\n";
    }
    else{
	print " BAD: $line\n";
    }
}

close(FILE);
