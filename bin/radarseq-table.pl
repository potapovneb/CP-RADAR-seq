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

if (@ARGV == 0)
{
    print "usage: $0 patches_scored.csv\n";
    print "\n\n";
    exit;
}

my $file = shift @ARGV;

my @head = ();

print "Reference,Strand,Start,End,Sequence\n";

open(CSV, $file) || die "Can't open '$file'";

while (my $line = <CSV>)
{
    chomp($line);

    my @tokens = split(/,/, $line);

    if (@head == 0)
    {
        @head = @tokens;
    }
    else
    {
        my %entry = ();

        for (my $i = 0; $i < @tokens; $i++)
        {
            $entry{$head[$i]} = $tokens[$i];
        }

        if ($entry{"SVM_Prediction"} == 1)
        {
            print join(",", $entry{"Reference"}, $entry{"Strand"}, $entry{"Start"}, $entry{"End"}, $entry{"Sequence"}), "\n";
        }
    }
}

close(CSV);
