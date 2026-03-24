#!/usr/bin/env perl

use warnings;
use strict;

my ($adn, $l, $header);

while ( <> ) { 
    chomp;

    ## First line is known, a header, so print it and process next one.
    if ( $. == 1 ) { 
        printf qq|%s_%s\n|, $_, q|200-300|;
        next;
    }   

    ## Concat adn while not found a header.
    if ( '>' ne substr $_, 0, 1 ) { 
        if ( ! $l ) { $l = length }
        $adn .= $_; 
        if ( ! eof ) { next }
    }   
    else {
        $header = sprintf qq|%s_%s\n|, $_, q|200-300|;
    }   

    ## Extract range 200-300 and insert newlines to set same length of 
    ## line than before.
    my $s = substr $adn, 10937, 100;
    $s =~ s/(.{$l})/$1\n/g;
    printf qq|%s\n|, $s; 
    undef $adn;

    ## If not end of file, print the header of the following adn.
    if ( ! eof ) { print $header }
}