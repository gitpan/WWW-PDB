#!/usr/bin/perl

use warnings;
use strict;

use WWW::PDB;

my $pdb = new WWW::PDB;

my $seq = <<'END_SEQ';
SHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKG
GPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVD
END_SEQ

print $pdb->blast($seq, 10.0, 'BLOSUM62', 'HTML');

print $pdb->blast('2ili', 'A', 10.0, 'BLOSUM62', 'HTML');

print $pdb->blast($seq, 10.0);

print $pdb->blast('2ili', 'A', 10.0);

