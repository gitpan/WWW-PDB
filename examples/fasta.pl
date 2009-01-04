#!/usr/bin/perl

use warnings;
use strict;

use WWW::PDB;

my $pdb = new WWW::PDB;

my $seq = q(
SHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKG
GPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVD
);

printf("%s\n", $_) for $pdb->fasta($seq, 100);

#print "$_\n" for $pdb->fasta('2ili', 'A', 100);

