#!/usr/bin/perl

use strict;

my $filename = $ARGV[0];

open(FILE,"<$filename");
my @childs = ();

my $count = 1;
while(my $line = <FILE>){
    chomp $line;
    if(length($line) == 0){
	next;
    }
    if($line =~ /^#/){
	next;
    }
    my $expdir = $line;
    my $cmd1 = "./register_cmd.sh $expdir";
    my $cmd2 = "./classifysex_cmd.sh $expdir";
    my $cmd3 = "./computeperframefeatures_cmd.sh $expdir";
    my $pid = fork();
    if($pid){
	push(@childs, $pid);
	$count++;
    }
    else{
	sleep($count);
	print "$cmd1\n";
	`$cmd1`;
	print "$cmd2\n";
	`$cmd2`;
	print "$cmd3\n";
	`$cmd3`;
	exit(0);
    }
	
}

close(FILE);

foreach (@childs) {
    waitpid($_, 0);
    print "child $_ finished\n";
}
