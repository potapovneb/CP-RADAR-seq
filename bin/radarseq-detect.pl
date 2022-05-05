#!/usr/bin/perl -w

##################################################################################
# Utilizing RADAR-seq to detect and quantitate DNA damage on a genome-wide scale #
# Copyright (C) 2022 New England Biolabs, Inc.                                   #
#                                                                                #
# This program is free software: you can redistribute it and/or modify           #
# it under the terms of the GNU General Public License as published by           #
# the Free Software Foundation, either version 3 of the License, or              #
# (at your option) any later version.                                            #
#                                                                                #
# This program is distributed in the hope that it will be useful,                #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                 #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  #
# GNU General Public License for more details.                                   #
#                                                                                #
# You should have received a copy of the GNU General Public License              #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.         #
##################################################################################

use strict;
use Getopt::Long;

### command-line options
my $o_edist  = 30;
my $o_ndist  = 30;
my $o_nsites =  4;
my $o_mtotal =  5;
my @o_bases  = ("A;C");
my %o_bases  = ();
my $o_nsfile = "";
my $o_msfile = "";

GetOptions(
    "bases=s"  => \@o_bases,
    "edist=f"  => \$o_edist,
    "ndist=f"  => \$o_ndist,
    "nsites=f" => \$o_nsites,
    "mtotal=f" => \$o_mtotal,
    "nsfile=s" => \$o_nsfile,
    "msfile=s" => \$o_msfile,
    );

for my $x (split(/;/, join(";", @o_bases))) {$o_bases{$x} = 1;}

### command-line usage
if (@ARGV == 0)
{
    print "usage: $0 [options] reference.fasta modifications.csv.gz\n\n";
    print "options:\n";
    print "  --bases\tA or C or A,C\n";
    print "  --edist\textension distance\n";
    print "  --ndist\tnicking site distance\n";
    print "  --nsites\tmodifiable sites tolerated in extension\n";
    print "  --mtotal\tminimum number of modified bases\n";
    print "  --nsfile\tnicking site file\n";
    print "  --msfile\tmethylation site file\n";
    exit;
}

my $reference_file = shift @ARGV;
my $modlist = shift @ARGV;

print STDERR "Loading reference file '$reference_file'...\n";
my $refs = load_references($reference_file);

print STDERR "Loading methylation sites '$o_msfile'...\n";
my $methylation_sites = ( $o_msfile ) ? load_nicking_sites($o_msfile) : {};

print STDERR "Loading nicking sites '$o_nsfile'...\n";
my $nicking_sites = ( $o_nsfile ) ? load_nicking_sites($o_nsfile) : {};

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print STDERR "Calculating cumulative counts in the reference...\n";

### this is used to determine number of "modifiable" bases in an interval
my %data = ();
my @bases = qw/A C G T/;

for my $rname (keys %$refs)
{
    ### init cumulative sum
    my %csum = map { $_ => 0 } @bases;

    ### store cumulative sum
    $data{$rname}{0} = { %csum };

    for (my $i = 0; $i < length($$refs{$rname}); $i++)
    {
        ### extract the base
        my $b = uc(substr($$refs{$rname}, $i, 1));

        ### accumulate counts
        $csum{$b}++;

        ### store cumulative sum
        $data{$rname}{$i + 1} = { %csum };
    }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my %cov = ();   # coverage
my @glist = (); # global patch list
my $countt = 0;

open(LIST, $modlist) || die "Can't open '$modlist'";

while (my $modfile = <LIST>)
{
    chomp($modfile);
    
    print STDERR "  Reading '$modfile'...\n";

    my %llist = (); ## local patch list

    open(IPD, "gzip -cd $modfile |") || die "Can't open '$modfile'";

    while (my $line = <IPD>)
    {
        ### skip CSV header line
        if (substr($line, 0, 7) eq "refName")
        {
            next;
        }

        chomp($line);

        my ($rname,
            $pos,
            $strand,
            $base,
            $score,
            $tMean,
            $tErr,
            $modelPrediction,
            $ipdRatio,
            $coverage,
            $frac,
            $fracLow,
            $fracUpp) = split(",", $line);
        
        # next if ( $strand != 0);
        
        $rname =~ s/^\"|\"$//g;
        
        $cov{$rname}{$pos}{$strand}++;
        
        ### focus on specified bases only (eg. m6A, m4C)
        next if (%o_bases && ! exists $o_bases{$base});
        
        ### skip known methylation sites
        next if (exists $$methylation_sites{$rname}{$strand}{$pos});

        ### start and extend patches
        if ($ipdRatio > 2.0)
        {
            ### determine extension start
            my $start = ( exists $llist{$rname}{$strand} ) ? $llist{$rname}{$strand}[-1]{"end"} : $pos;
            
            ### determine number of sites where modification is possible
            my $nsites = 0;

            for my $b (keys %o_bases)
            {
                if ($strand == 0)
                {
                    $nsites += $data{$rname}{$pos}{$b} - $data{$rname}{$start}{$b};
                }
                else
                {
                    my $bc = complement($b);
                    
                    $nsites += $data{$rname}{$pos}{$bc} - $data{$rname}{$start}{$bc};
                }
            }
            
            ### patch extension
            if (! exists $llist{$rname}{$strand} || ($pos - $start) > $o_edist || $nsites > $o_nsites)
            {
                ### start a new patch if:
                ###  1. no previous patches defined
                ###  2. modified base too far from the last known patch
                ###  3. modified base is within 'radius of attraction' but there
                ###     are too many non-modified bases between patch start and
                ###     current modified base

                my @path = split(/\//, $modfile);
                my @name = split(/\./, $path[-1]);

                # my $p = index($modfile, "/1-strata/");
                # my $stratum = substr($modfile, $p + 10, 5);

                push @{$llist{$rname}{$strand}},
                {
                    "modfile"=> $name[0],
                    "rname"  => $rname,
                    "strand" => $strand,
                    "start"  => $pos,
                    "end"    => $pos,
                    "mA"     => 0,
                    "mC"     => 0,
                    "mG"     => 0,
                    "mT"     => 0,
                    "mTotal" => 0,
                    "iA"     => 0,
                    "iC"     => 0,
                    "iG"     => 0,
                    "iT"     => 0,
                    "iTotal" => 0,
                    "ipds"   => {}
                };
            }
            else
            {
                ###  extend previous patch
                $llist{$rname}{$strand}[-1]{"end"} = $pos;
            }

            ### count modified bases
            $llist{$rname}{$strand}[-1]{"m$base"}++;
            $llist{$rname}{$strand}[-1]{"mTotal"}++;

            ### sum up IPD values
            $llist{$rname}{$strand}[-1]{"i$base"} += $ipdRatio;
            $llist{$rname}{$strand}[-1]{"iTotal"} += $ipdRatio;

            ### keep track of all IPDs for future reference
            $llist{$rname}{$strand}[-1]{"ipds"}{$pos} = $ipdRatio;
        }
    }

    ### push found patches to the global list
    for my $rname (keys %llist)
    {
        for my $strand (keys %{$llist{$rname}})
        {
            push @glist, @{$llist{$rname}{$strand}};
        }
    }

    $countt++;
}

close(LIST);

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print STDERR "Writing patches file...\n";

open(PATCHES, ">", "patches.csv") || die "Can't write 'patches.csv'";

my @extra_headers = ();

if ($o_nsfile)
{
    @extra_headers = qw/Site Distance Location Status/;
}

### print patch info
print PATCHES join(",",
    qw/Modfile Reference Strand Start End Length mA mC mG mT mTotal A C G T Total iA iC iG iT iTotal Modifiable Fraction/,
    @extra_headers,
    qw/Sequence Ipds/), "\n";

for my $patch (sort { $$a{"modfile"}   cmp $$b{"modfile"}
               || $$a{"rname"}  cmp $$b{"rname"} 
               || $$a{"strand"} <=> $$b{"strand"} 
               || $$a{"start"}  <=> $$b{"start"} } @glist
               )
{
    my $rname  = $$patch{"rname"};
    my $strand = $$patch{"strand"};
    my $start  = $$patch{"start"};
    my $end    = $$patch{"end"};
    my $length = $$patch{"end"} - $$patch{"start"} + 1;

    ### number of bases
    my $nA = 0;
    my $nC = 0;
    my $nG = 0;
    my $nT = 0;

    if ($strand == 0)
    {
        $nA = $data{$rname}{$end}{"A"} - $data{$rname}{$start-1}{"A"};
        $nC = $data{$rname}{$end}{"C"} - $data{$rname}{$start-1}{"C"};
        $nG = $data{$rname}{$end}{"G"} - $data{$rname}{$start-1}{"G"};
        $nT = $data{$rname}{$end}{"T"} - $data{$rname}{$start-1}{"T"};
    }
    else
    {
        $nA = $data{$rname}{$end}{"T"} - $data{$rname}{$start-1}{"T"};
        $nC = $data{$rname}{$end}{"G"} - $data{$rname}{$start-1}{"G"};
        $nG = $data{$rname}{$end}{"C"} - $data{$rname}{$start-1}{"C"};
        $nT = $data{$rname}{$end}{"A"} - $data{$rname}{$start-1}{"A"};
    }

    ### number of modifiable bases
    my $modifiable = 0;
    $modifiable += $nA if ( exists $o_bases{"A"} );
    $modifiable += $nC if ( exists $o_bases{"C"} );
    $modifiable += $nG if ( exists $o_bases{"G"} );
    $modifiable += $nT if ( exists $o_bases{"T"} );

    ### check if patch is close to a known nicking site
    my @extra_values = ();

    if (%$nicking_sites)
    {
        my $min = 1e9; ## distance to the closest nicking site
        my $hit = "NA"; ## closest nicking site ID
        my $loc = "OUTSIDE"; ## is nicking side OUTSIDE or INSIDE patch

        for my $pos (keys %{$$nicking_sites{$rname}{$strand}})
        {
            my $dist = ($strand == 0) ? abs($start - $pos) : abs($end - $pos);

            if ($dist < $min)
            {
                $min = $dist;
                $hit = $pos;
            }

            if ($pos >= $start && $pos <= $end)
            {
                $loc = "INSIDE";
            }
        }

        my $status = ($min < $o_ndist || $loc eq "INSIDE") ? "TRUE" : "FALSE";

        @extra_values = ($hit, $min, $loc, $status);
    }

    ### combine all IPDs to a parseable string
    my $ipds = join(";" , map { sprintf("%i=%.3f", $_ - $start, $$patch{"ipds"}{$_}) } (sort { $a <=> $b } keys %{$$patch{"ipds"}}) );

    ### extract patch sequence
    my $seq = substr($$refs{$rname}, $start - 1, $length);

    if ($strand == 1)
    {
        $seq = complement($seq);
    }

    ### present modified bases as CAPITAL letters
    my $nseq = "";

    for (my $i = 0; $i < $length; $i++)
    {
        if (exists $$patch{"ipds"}{$i + $start})
        {
            $nseq .= uc(substr($seq, $i, 1));
        }
        else
        {
            $nseq .= lc(substr($seq, $i, 1));
        }
    }

    if ($$patch{"mTotal"} < $o_mtotal)
    {
        next;
    }

    ### print all info
    print PATCHES join(",",
        $$patch{"modfile"},
        $$patch{"rname"},
        $$patch{"strand"},
        $$patch{"start"},
        $$patch{"end"},
        $length,
        $$patch{"mA"},
        $$patch{"mC"},
        $$patch{"mG"},
        $$patch{"mT"},
        $$patch{"mTotal"},
        $nA,
        $nC,
        $nG,
        $nT,
        $nA + $nC + $nG + $nT,
        ( $$patch{"mA"} > 0 ) ? $$patch{"iA"} / $$patch{"mA"} : 0,
        ( $$patch{"mC"} > 0 ) ? $$patch{"iC"} / $$patch{"mC"} : 0,
        ( $$patch{"mG"} > 0 ) ? $$patch{"iG"} / $$patch{"mG"} : 0,
        ( $$patch{"mT"} > 0 ) ? $$patch{"iT"} / $$patch{"mT"} : 0,
        ( $$patch{"mTotal"} > 0 ) ? $$patch{"iTotal"} / $$patch{"mTotal"} : 0,
        $modifiable,
        sprintf("%.4f", $$patch{"mTotal"} / $modifiable),
        @extra_values,
        $nseq,
        $ipds,
    ), "\n";
}

close(PATCHES);
print STDERR "done.\n";

print STDERR "Writing coverage file...\n";
open(COV, ">", "coverage.csv") || die "Can't open 'coverage.csv'";
print COV join(",", qw/Reference Position N0 N1/), "\n";

for my $rname (keys %cov)
{
    for my $pos (sort { $a <=> $b } keys %{$cov{$rname}})
    {
        my @values = ($rname,$pos);

        for my $strand (0, 1)
        {
            if (exists $cov{$rname}{$pos}{$strand})
            {
                push @values, $cov{$rname}{$pos}{$strand};
            }
            else
            {
                push @values, 0;
            }
        }

        print COV join(",", @values), "\n";
    }
}

close(COV) ;
print STDERR "done.\n";

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub load_references {
    my ($file) = @_;
    
    my %refs = ();
    my $name = "";
    
    open(FA, $file) || die "Can't open '$file'";
    
    while (my $line = <FA>)
    {
        chomp($line);
    
        if (substr($line, 0, 1) eq ">")
        {
            $name = substr($line, 1);
        }
        else
        {
            $refs{$name} .= $line;
        }
    }
    
    close(FA);
    
    return \%refs;
}

sub load_nicking_sites {
    my ($file) = @_;

    my %data = ();
    my @head = ();

    open(IN, $file) || die "Can't open $file'";

    while (my $line = <IN>)
    {
        chomp($line);
        
        my @tokens = split(/,/, $line);
        
        if (! @head)
        {
            @head = @tokens;
        }
        else
        {
            my %entry = ();
            
            for(my $i = 0; $i < @tokens; $i++)
            {
                $entry{$head[$i]} = $tokens[$i];
            }

            $data{$entry{"Reference"}}{$entry{"Strand"}}{$entry{"Position"}} = 1;
        }
    }
    
    close(IN);

    return \%data;
}

sub complement {
    my ($seq) = @_;

    my %bases = (
    "A" => "T",
    "C" => "G",
    "G" => "C",
    "T" => "A",
    "a" => "t",
    "c" => "g",
    "g" => "c",
    "t" => "a",
    "-" => "-",
    );

    my $complement = join("", map { $bases{$_} } split(//, $seq));

    return $complement;
}
