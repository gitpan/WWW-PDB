use Test::More tests => 8;
BEGIN { use_ok('WWW::PDB') };

my $pdb;
ok($pdb = new WWW::PDB, "new() test");

ok($pdb->get_structure('2ili'), "get_structure() test");
ok($pdb->get_structure_factors('2ili'), "get_structure_factors() test");

my $seq = <<'END_SEQ';
SHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKG
GPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVD
END_SEQ

ok($pdb->blast($seq, 10.0, 'BLOSUM62', 'HTML'), "4-arg blast() test");
ok($pdb->blast('2ili', 'A', 10.0, 'BLOSUM62', 'HTML'), "5-arg blast() test");
ok($pdb->blast($seq, 10.0), "2-arg blast() test");
ok($pdb->blast('2ili', 'A', 10.0), "3-arg blast() test");

