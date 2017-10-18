#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use List::Util qw(min max);
use Process_cmd;
use DelimParser;
use Data::Dumper;
use SAM_entry;
use JSON::XS;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);  


my $usage = <<__EOUSAGE__;

######################################################################################
#
#  --vts <string>              vts predictions 'preliminary' output file (not abridged)
#  --junc_file <string>        chimeric junction file 
#  --chimeric_sam <string>     chimeric_sam file #
#  --output_prefix <string>     output prefix
#
######################################################################################

__EOUSAGE__

    ;


my $help_flag;

my $vts_results_file;
my $junc_file;
my $chimeric_sam;
my $output_prefix;

my $genome_fa = "/media/heyao/YaoDataS/reference/GenomeLib/GRCh37_gencode_v19_CTAT_lib_July192017_plus_hepatitas/ref_genome.fa";


&GetOptions( 'help|h' => \$help_flag,

             'vts=s' => \$vts_results_file,
             
             'junc_file=s' => \$junc_file,
             'chimeric_sam=s' => \$chimeric_sam,
             'output_prefix=s' => \$output_prefix,
    );

if ($help_flag) {
    die $usage;
}

unless ($vts_results_file && $chimeric_sam && $output_prefix) {
    die $usage;
}


 main: {


     #### Index VTS to reads
     my %vts_to_junction_reads; #{vts_complex_name}->(j1, j2, j3, ... jn)

     open (my $fh, $vts_results_file) or die "Error, cannot open file $vts_results_file";
     my $tab_reader = new DelimParser::Reader($fh, "\t");

     while (my $row = $tab_reader->get_row()) {
         my $vts_name = $row->{'#VirusIntegrationName'} or die "Error, cannot get fusion name from " . Dumper($row);
         $vts_name = $vts_name . "|" . $row->{'LeftBreakpoint'} . "|" . $row->{'RightBreakpoint'};
         #print $vts_name, "\n";
         my @junction_reads_list = split(/,/, $row->{JunctionReads});


         foreach my $read_name (@junction_reads_list) {
             push @{$vts_to_junction_reads{$vts_name}}, $read_name;
         }
             
     }
     close $fh;

     #### Get Chimeric Alignment from junction file

     open(my $junc_fh, $junc_file) or die "Error, cannot open file $junc_file";

     my %read_to_alignment; # {read_name}->{left|right}->(first base, CIGAR);
     while (<$junc_fh>) {
         chomp;
         my @F = split /\t/;
         my $read_name = $F[9];
         
         my ($leftFirstBase, $leftAlignment) = ($F[10], $F[11]);
         my ($rightFirstBase, $rightAlignment) = ($F[12], $F[13]);

         push @{$read_to_alignment{$read_name}{left}}, ($leftFirstBase, $leftAlignment);
         push @{$read_to_alignment{$read_name}{right}}, ($rightFirstBase, $rightAlignment);

     }         
     close $junc_fh;

     ####  Get sam entry from the sam file
     my %read_to_sam; #Assume a core read should have 3 sam entry (Chimeric)
     open(my $sam_fh, $chimeric_sam) or die "Error, cannot open file $chimeric_sam";

     while (<$sam_fh>) {
         my $aln = new SAM_entry($_);
         my $read_name = &SAM_entry::get_core_read_name($aln);
         push @{$read_to_sam{$read_name}}, $aln;
         
     }

     close($sam_fh);
     ## Selecting junction reads we need 
     my %vts_alignment_json;
     foreach my $vts (keys %vts_to_junction_reads) { # For each integration sites
         #--- Get reference
         my ($vts_name, $l_bkpt, $r_bkpt) = split("\\|", $vts);
         #print "$vts_name\t$l_bkpt\t$r_bkpt\n";
         my ($leftside_name, $rightside_name) = split("--", $vts_name); 
         if ($leftside_name =~ /NC_00/) {
             $vts_alignment_json{$vts}{left} = 'virus';
             $vts_alignment_json{$vts}{right} = 'host';
         } else {
             $vts_alignment_json{$vts}{right} = 'virus';
             $vts_alignment_json{$vts}{left} = 'host';
         }
         my @read_names = @{$vts_to_junction_reads{$vts}};

         my %target_aln;
         my %target_sam;

         my %alignment_vis_json; #{read_name}->{'host','virus'}->{alignmentinfo} 

         foreach my $read_name (@read_names) {
            $target_aln{$read_name} =  $read_to_alignment{$read_name};
            $target_sam{$read_name} = $read_to_sam{$read_name};
            foreach my $aln (@{$target_sam{$read_name}}) { # A pair of reads should have at least three part
                
                my %aln_temp;
                my $ref_pos = &SAM_entry::get_scaffold_name($aln);
                my ($genome_coord, $align_coord) = &SAM_entry::get_alignment_coords($aln);
                my ($genome_min_coord, $genome_max_coord) = &SAM_entry::get_genome_span($aln);
                my $sequence = &SAM_entry::get_sequence($aln);
                my $cigar = &SAM_entry::get_cigar_alignment($aln);
                # Skip non-chimeric reads (if an alignment has neither S cigar or S length less than 12 )
                if ($cigar !~ /(1\d\dS)|(1[2-9]S)|([2-9][0-9]S)/ ) {
                    next;
                }

                my $rl = &reconstruct_from_cigar($sequence, $cigar);
                # Assign chimeric part 
                if ($ref_pos =~ /NC_00/) {
                    $aln_temp{ref} = $ref_pos;
                    $aln_temp{alignment} = $rl;
                    $aln_temp{genome_min_coord} = $genome_min_coord;
                    $aln_temp{genome_max_coord} = $genome_max_coord;
                    $aln_temp{cigar} = $cigar;
                    $alignment_vis_json{$read_name}{virus} = \%aln_temp;
                } else {
                    $aln_temp{ref} = $ref_pos;
                    $aln_temp{alignment} = $rl;
                    $aln_temp{genome_min_coord} = $genome_min_coord;
                    $aln_temp{genome_max_coord} = $genome_max_coord;
                    $aln_temp{cigar} = $cigar;
                    $alignment_vis_json{$read_name}{host} = \%aln_temp;
                } 
               
            }

         }
         $vts_alignment_json{$vts}{'alignment'} = \%alignment_vis_json;
         my %left_side_ref = &get_reference_seq($l_bkpt, 'left', $vts_alignment_json{$vts}{'left'}, \%alignment_vis_json);
         my %right_side_ref = &get_reference_seq($r_bkpt, 'right', $vts_alignment_json{$vts}{'right'}, \%alignment_vis_json);

         $vts_alignment_json{$vts}{'left_side_ref'} = \%left_side_ref;
         $vts_alignment_json{$vts}{'right_side_ref'} = \%right_side_ref;

     }
     
     #print Dumper \%vts_alignment_json;
     # Output json
     my $j = JSON::XS->new->utf8->pretty(1);
     my $output = $j->encode(\%vts_alignment_json);
     print $output; 
     
}

sub get_reference_seq {
    my ($bkpt, $which_part, $virus_or_host, $alignment_href) = @_;
    #print "$bkpt,$which_part,$virus_or_host\n";
    my ($ref, $bkpos, $strand) = split(":", $bkpt);


    my $leftmost ; 
    my $rightmost;
    if (($which_part eq 'left') & ($virus_or_host eq 'host')) {
        # check strand for identifying 5'end
        if ( $strand eq '+') {
            $leftmost = $bkpos - 150;
        } elsif ($strand eq '-') {
            $leftmost = $bkpos + 150;
        } else {
            die "This $which_part is virus but NOT $virus_or_host";
        }

        $rightmost = $bkpos; 

    } elsif (($which_part eq 'left') & ($virus_or_host eq 'virus')) {
        
        $rightmost = $bkpos; 

        if ($strand ne '*') {die "Not Virus for lefthand side"};
        foreach my $k (keys %$alignment_href) {
           my $coord = $alignment_href->{$k}->{virus}->{genome_max_coord};
           if ($coord == $bkpos) { # ------------// cis breakpoint
                $leftmost = $bkpos - 150;
           } else {
                $leftmost = $bkpos + 150;
           }
           last;
        }



    } elsif (($which_part eq 'right') & ($virus_or_host eq 'host')) {
        
        if ( $strand eq '+') {
            $rightmost = $bkpos + 150;
        } elsif ($strand eq '-') {
            $rightmost = $bkpos - 150;
        } else {
            die "This $which_part is virus but NOT $virus_or_host";
        }
        $leftmost = $bkpos; 


    } elsif (($which_part eq 'right') & ($virus_or_host eq 'virus')) {
        $leftmost = $bkpos;
        if ($strand ne '*') {die "Not Virus for righthand side"};
 #       print Dumper \$alignment_href;
        foreach my $k (keys %$alignment_href) {
            my $coord = $alignment_href->{$k}->{virus}->{genome_min_coord};
            if ($coord == $bkpos) { # //------------ cis breakpoint
                $rightmost = $bkpos + 150;
            } else {
                $rightmost = $bkpos - 150;
            }
            last;
        }

    }


    #-- use samtools to get sequence
    my %final_seq;
    $final_seq{leftmost_coord} = $leftmost;
    $final_seq{rightmost_coord} = $rightmost;
   
    if ($leftmost < $rightmost) {
 #       print "samtools faidx $genome_fa $ref:$leftmost-$rightmost\n";
        open(LFA, "samtools faidx $genome_fa $ref:$leftmost-$rightmost|");
        my @aux = undef;
        while ((my $name, my $seq) = readfq(\*LFA, \@aux)) {
            $final_seq{seq} = $seq;
#            print "[$l_ref:$l_leftmost-$l_rightmost]:", $seq, "\n";
   #         $seq =~ tr/ACGTacgt/TGCAtgca/;
   #         $seq = reverse($seq);

        }
        close(LFA);
    } else {
        
    } {
 #       print "samtools faidx $genome_fa $ref:$rightmost-$leftmost\n";
        open(RFA, "samtools faidx $genome_fa $ref:$rightmost-$leftmost|");
        my @aux = undef;
        while ((my $name, my $seq) = readfq(\*RFA, \@aux)) {
            $final_seq{seq} = reverse($seq);

        }
        close(RFA);
        
    } 
        
    
    return %final_seq;

}

sub reconstruct_from_cigar {
    # $seq is the read sequence
    # $coord is the alignment coord on read
    # $cigar is the read cigar string
    my ($seq, $cigar) = @_;

    my $read_alignment = '';
    my $r_lend = 0;
    my $r_rend = 0;
    while ($cigar =~ /(\d+)([A-Z])/g) {
        my $len =  $1;
        my $code = $2;


        if ($code eq 'M' || $code eq 'I') { #aligned bases match or mismatch
            $r_rend = $r_lend + $len;
            $read_alignment = $read_alignment . substr($seq, $r_lend, $r_rend);
            $r_lend = $r_rend;
        } elsif ($code eq 'D' || $code eq 'N') { #insertion in the genome or gap in read
            $read_alignment = $read_alignment . '-' x $len;
        } elsif ($code eq 'S' || $code eq 'H') {
            $r_rend = $r_lend + $len;
            #$read_alignment = $read_alignment . lc(substr($seq, $r_lend, $r_rend));
            $r_lend = $r_rend;
        }
    }

    return($read_alignment);
}
     

         


sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (! @$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my $name = /^.(\S+)/? $1 : '';
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}         
