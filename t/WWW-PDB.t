use Test::More tests => 7;
BEGIN { use_ok('WWW::PDB') };

SKIP: {
    my $conn_ok;
    eval 'use Net::Config qw(%NetConfig); $conn_ok = $NetConfig{test_hosts}';
    skip 'Net::Config needed for network-related tests', 6 if $@;
    skip 'No network connection', 6 unless $conn_ok;

    WWW::PDB->import(':file');

    ok(get_structure('2ili'), "get_structure() test");
    ok(get_structure_factors('2ili'), "get_structure_factors() test");

my $seq = <<'END_SEQ';
SHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKG
GPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVD
END_SEQ

ok(WWW::PDB->blast($seq, 10.0, 'BLOSUM62', 'HTML'), "4-arg blast() test");
ok(WWW::PDB->blast('2ili', 'A', 10.0, 'BLOSUM62', 'HTML'), "5-arg blast() test");
ok(WWW::PDB->blast($seq, 10.0), "2-arg blast() test");
ok(WWW::PDB->blast('2ili', 'A', 10.0), "3-arg blast() test");

}
