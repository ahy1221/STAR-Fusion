#!/usr/bin/env perl

use strict;
use warnings;
use v5.20;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use Process_cmd;
use DelimParser;
use Data::Dumper;
use SAM_entry;
use JSON::XS;
use Mojo::DOM;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);  
my $HTML_TEMPLATE = "./util/fuse.html";
my $json = "./test.json";

open(FILE, $HTML_TEMPLATE) or die "Can not open $HTML_TEMPLATE";
my $content = do {local $/=undef; <FILE> };
close(FILE);

open(JSON, $json) or die "Can not open $json file";
my $json_string = do {local $/=undef; <JSON> };
close(JSON);
my $dom = Mojo::DOM->new($content);
my $json_hash = decode_json($json_string);
#print Dumper \$json_hash;
foreach my $vts (keys %$json_hash) {
    say $vts;
    my $left_side_seq = $json_hash->{$vts}->{left_side_ref}->{seq};
    my $right_side_seq = $json_hash->{$vts}->{right_side_ref}->{seq};
    my $leftpart_rightmost_coord = $json_hash->{$vts}->{left_side_ref}->{rightmost_coord};
    my $rightpart_leftmost_coord = $json_hash->{$vts}->{right_side_ref}->{leftmost_coord};
    my $left = $json_hash->{$vts}->{left};
    my $right = $json_hash->{$vts}->{right};
    my $refseq_row = <<REF;
<tr>
<td class="alignright" colspan= "3"> $left_side_seq </td>
<td class="alignleft"> $right_side_seq </td>
</tr>
REF
    $dom = $dom->find('div.refseq ')->last->append("$refseq_row")->root;
    foreach my $read (keys %{$json_hash->{$vts}->{alignment}}) {
        #print Dumper \$json_hash->{$vts}->{alignment};
        #--- Handle left seq
        my $lseq = '';
        my $rseq = '';
        if ($json_hash->{$vts}->{alignment}->{$read}->{$left}->{genome_max_coord} != $leftpart_rightmost_coord)  {
            $lseq = reverse($json_hash->{$vts}->{alignment}->{$read}->{$left}->{alignment}); 
        } else {
            $lseq = $json_hash->{$vts}->{alignment}->{$read}->{$left}->{alignment}; 
        } 
            
        
        if ($json_hash->{$vts}->{alignment}->{$read}->{$right}->{genome_min_coord} != $rightpart_leftmost_coord)  {
         $rseq = reverse($json_hash->{$vts}->{alignment}->{$read}->{$right}->{alignment}); 
     } else {
         $rseq = $json_hash->{$vts}->{alignment}->{$read}->{$right}->{alignment}; 
         
     } 
         
     
        
        my $table_row = <<READS;
<tr>
<td> 0001 </td>
<td> 0|0 </td>
<td class="alignright"> $lseq <\td>
<td class="alignleft">  $rseq <\td>
</tr>
READS
        $dom = $dom->find('div.reads ')->last->append("$table_row")->root;
    }
    print $dom->to_string;
    last;
}



__DATA__
