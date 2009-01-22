package WWW::PDB;

=head1 NAME

WWW::PDB - Perl interface to the Protein Data Bank

=head1 SYNOPSIS

  use WWW::PDB qw(:all);

  # set directory for caching downloads
  WWW::PDB->cache('/foo/bar');
  
  my $fh = get_structure('2ili');
  print while <$fh>;
  
  my @pdbids = WWW::PDB->keyword_query('carbonic anhydrase');
  for(@pdbids) {
      my $citation = WWW::PDB->get_primary_citation_title($_),
      my @chains   = WWW::PDB->get_chains($_);
      printf("%s\t%s\t[%s]\n", $_, $citation, join(', ', @chains));
  }

  my $seq = q(
      VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK
      TYFPHFDLSHGSAQVKGHGKKVADALTAVAHVDDMPNAL
  );
  print WWW::PDB->blast($seq, 10.0, 'BLOSUM62', 'HTML');

=head1 DESCRIPTION

The Protein Data Bank (PDB) was established in 1971 as a repository of the
atomic coordinates of protein structures (Bernstein I<et al.>, 1997).  It
has since outgrown that role, proving invaluable not only to the research
community but also to students and educators (Berman I<et al.>, 2002).

L<WWW::PDB> is a Perl interface to the Protein Data Bank.  It provides
functions for retrieving files, optionally caching them locally.
Additionally, it wraps the functionality of the PDB's SOAP web services.

=cut

use 5.006;
use strict;
use warnings;

use Carp;
use Exporter;
use Fcntl;
use File::Path;
use File::Spec;
use IO::File;
use IO::Uncompress::Gunzip;
use Net::FTP;
use SOAP::Lite;

use constant {
    BOOLEAN => 0,
    DOUBLE  => 1,
    INT     => 2,
    SELF    => 3,
    STRING  => 4,
};

our @ISA = qw(Exporter);
our $VERSION = '0.00_02';
$VERSION = eval $VERSION;

our %EXPORT_TAGS = (
    file   => [qw(get_structure get_structure_factors)],
    status => [qw(get_status is_current is_obsolete is_unreleased
                  is_model is_unknown)],
);
our @EXPORT_OK = map {@$_} values %EXPORT_TAGS;
$EXPORT_TAGS{all} = \@EXPORT_OK;

my($uri, $proxy, $ftp, $cache, $soap);

=head1 FUNCTIONS

=head2 CUSTOMIZATION

Let's start with some functions that let you customize how the module does
its job.  You probably won't play with any of these very often (if at all)
except for C<cache>, which is recommended for anyone that expects to do
extensive work with a set of files: that way you don't waste resources
downloading them each time.

=over 4

=item WWW::PDB->ftp( [ $HOST ] )

Returns the host name for the PDB FTP archive, first setting it to $FTP if
it's specified.  Default value is F<ftp.wwpdb.org>.

=cut

sub ftp {
    return $ftp = $_[1] ? $_[1] : $ftp || 'ftp.wwpdb.org';
}

=item WWW::PDB->cache( [ $DIR ] )

Returns the local cache directory, first setting it to $DIR if it's
specified.  If C<defined>, the module will look for files here first and
also use the directory to store any downloads.

=cut

sub cache {
    $cache = $_[1] if $_[1];
    return $cache;
}

=item WWW::PDB->ns( [ $URI ] )

Returns the namespace URI for the PDB web services, first setting it to $URI
if it's specified.  Default value is http://www.pdb.org/pdb/services/pdbws.

=cut

sub ns {
    my $tmp = $uri;
    $uri = $_[1] ? $_[1] : $uri || 'http://www.pdb.org/pdb/services/pdbws';
    $_[0]->soap->ns($uri) unless $tmp && $tmp eq $uri;
    return $uri;
}

=item WWW::PDB->proxy( [ $URI ] )

Returns the proxy for the PDB web services, first setting it to $URI if it's
specified.  Default value is http://www.pdb.org/pdb/services/pdbws.

=cut

sub proxy {
    my $tmp = $proxy;
    $proxy = $_[1] ? $_[1] : $proxy || 'http://www.pdb.org/pdb/services/pdbws';
    $_[0]->soap->proxy($proxy) unless $tmp && $tmp eq $proxy;
    return $proxy;
}

=item WWW::PDB->soap( [ $CLIENT ] )

Returns the client L<SOAP::Lite> object used by this module to talk to
the PDB's SOAP interface, first setting it to $CLIENT if it's specified.
It's best not to access it directly, but if you must, this is how.

=cut

sub soap {
    return $soap = $_[1] ? $_[1] : $soap ||
        SOAP::Lite->ns($_[0]->ns)->proxy($_[0]->proxy);
}

=back

=head2 FILE RETRIEVAL

Each of the following functions takes a PDB ID as input and returns a file
handle (or C<undef> on failure).  You can import these into your namespace
with the C<file> tag, as in C<use WWW::PDB qw(:file)>.

=over 4

=item get_structure( $PDBID )

Retrieves the structure in PDB format.

=cut

sub get_structure {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $class = shift;
    my $pdbid = lc(shift);
    return $pdbid =~ /^.(..).$/
        ? $class->_get_file(qw(pub pdb data structures divided pdb), $1,
        "pdb${pdbid}.ent.gz") : undef;
}

=item get_structure_factors( $PDBID )

Retrieves the structure factors file.

=cut

sub get_structure_factors {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $class = shift;
    my $pdbid = lc(shift);
    return $pdbid =~ /^.(..).$/
        ? $class->_get_file(qw(pub pdb data structures divided
        structure_factors), $1, "r${pdbid}sf.ent.gz") : undef;
}

=back

=head2 PDB ID STATUS

The following functions deal with the status of PDB IDs.  You can import
them into your namespace with the C<status> tag:
C<use WWW::PDB qw(:status)>.

=over 4

=item get_status( $PDBID )

Finds the status of the structure with the given $PDBID.  Return is in
C<qw(CURRENT OBSOLETE UNRELEASED MODEL UNKNOWN)>.

=cut

sub get_status {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    return 'UNKNOWN' if length($_[1]) != 4;
    my $class = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getIdStatus', $pdbid
    );
    return $ret;
}

=item is_current( $PDBID )

Checks whether or not the specified $PDBID corresponds to a current
structure.  Implemented for orthogonality, all this does is check
if C<get_status> returns C<CURRENT>.

=cut

sub is_current {
    my $class = UNIVERSAL::isa($_[0], __PACKAGE__) ? shift : __PACKAGE__;
    return $class->get_status(@_) eq 'CURRENT';
}

=item is_obsolete( $PDBID )

Checks whether or not the specified $PDBID corresponds to an obsolete
structure.  This is actually defined by the PDB web services interface.

=cut

sub is_obsolete {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'isStructureIdObsolete', $pdbid
    );
    return $ret;
}

=item is_unreleased( $PDBID )

Checks whether or not the specified $PDBID corresponds to an unreleased
structure.  Implemented for orthogonality, all this does is check
if C<get_status> returns C<UNRELEASED>.

=cut

sub is_unreleased {
    my $class = UNIVERSAL::isa($_[0], __PACKAGE__) ? shift : __PACKAGE__;
    return $class->get_status(@_) eq 'UNRELEASED';
}

=item is_model( $PDBID )

Checks whether or not the specified $PDBID corresponds to a model
structure.  Implemented for orthogonality, all this does is check
if C<get_status> returns C<MODEL>.

=cut

sub is_model {
    my $class = UNIVERSAL::isa($_[0], __PACKAGE__) ? shift : __PACKAGE__;
    return $class->get_status(@_) eq 'MODEL';
}

=item is_unknown( $PDBID )

Checks whether or not the specified $PDBID is unknown.  Implemented
for orthogonality, all this does is check if C<get_status> returns
C<UNKNOWN>.

=cut

sub is_unknown {
    my $class = UNIVERSAL::isa($_[0], __PACKAGE__) ? shift : __PACKAGE__;
    return $class->get_status(@_) eq 'UNKNOWN';
}

=back

=head2 PDB WEB SERVICES

The following methods are the interface to the PDB web services.

=over 4

=item blast( $SEQUENCE , $CUTOFF , $MATRIX , $OUTPUT_FORMAT )

=item blast( $PDBID , $CHAINID, $CUTOFF , $MATRIX , $OUTPUT_FORMAT )

=item blast( $SEQUENCE , $CUTOFF )

=item blast( $PDBID , $CHAINID , $CUTOFF )

Performs a BLAST against sequences in the PDB and returns the output of
the BLAST program. XML is used if the output format is unspecified.

=cut

sub _blast_pdb {
    my($self, $sequence, $cutoff, $matrix, $output_format) = 
        _wrap(\@_ => [SELF, STRING, DOUBLE, STRING, STRING]);

    my $ret = $self->_call(
        'blastPDB', $sequence, $cutoff, $matrix, $output_format
    );
    return $ret;
}

sub _blast_structure_id_pdb {
    my $self = shift;

    # I keep getting "ERROR: No Results Found" using the PDB's 5 argument
    # form of blastPDB. Here's a workaround:    
    my $seq = $self->get_sequence(shift, shift);
    return $self->blast($seq, @_);

#   my($self, $pdbid, $chainid, $cutoff, $matrix, $output_format) =
#       _wrap(\@_ => [SELF, STRING, STRING, DOUBLE, STRING, STRING]);
#   my $ret = $self->_call(
#       'blastPDB', $pdbid, $chainid, $cutoff, $matrix, $output_format
#   );
#   return $ret;
}

sub _blast_query_xml {
    my($self, $sequence, $cutoff) = _wrap(\@_ => [SELF, STRING, DOUBLE]);
    my $ret      = $self->_call(
        'blastQueryXml', $sequence, $cutoff
    );
    return $ret;
}

sub _blast_structure_id_query_xml {
    my($self, $pdbid, $chainid, $cutoff)
        = _wrap(\@_ => [SELF, STRING, STRING, DOUBLE]);
    my $ret     = $self->_call(
        'blastStructureIdQueryXml', $pdbid, $chainid, $cutoff
    );
    return $ret;
}

sub blast {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self = shift;
    my $ret;
    my $c = scalar(@_);
    if   ($c == 4) { $ret = $self->_blast_pdb(@_) }
    elsif($c == 5) { $ret = $self->_blast_structure_id_pdb(@_) }
    elsif($c == 2) { $ret = $self->_blast_query_xml(@_) }
    elsif($c == 3) { $ret = $self->_blast_structure_id_query_xml(@_) }
    else { confess 'Called blast with unexpected number of arguments' }
    return $ret;
}

=item fasta( $SEQUENCE , $CUTOFF )

=item fasta( $PDBID , $CHAINID , $CUTOFF )

Takes a sequence or PDB ID and chain identifier and runs FASTA using the
specified cut-off. The results are overloaded to give PDB IDs when used
as strings, but they can also be explicitly probed for a C<pdbid> or
FASTA C<cutoff>:

  printf("%s %s %s\n", $_, $_->pdbid, $_->cutoff)
      for $pdb->fasta("2ili", "A");

=cut

sub _fasta_query {
    my $self     = shift;
    my $sequence = _to_string(shift);
    my $cutoff   = _to_double(shift);
    my $ret      = $self->_call(
        'fastaQuery', $sequence, $cutoff
    );
    return $ret;
}

sub _fasta_structure_id_query {
    my($self, $pdbid, $chainid, $cutoff)
        = _wrap(\@_ => [SELF, STRING, STRING, DOUBLE]);
    my $ret = $self->_call(
        'fastaStructureIdQuery', $pdbid, $chainid, $cutoff
    );
    return $ret;
}

sub fasta {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self = shift;
    my $c    = scalar(@_);
    my $ret;
    if   ($c == 2) { $ret = $self->_fasta_query(@_) }
    elsif($c == 3) { $ret = $self->_fasta_structure_id_query(@_) }
    else { confess 'Called fasta with unexpected number of arguments' }
    $_ = bless(\"$_", 'WWW::PDB::_FastaResult') for @$ret;
    return wantarray ? @$ret : $ret;
}

=item get_chain_length( $PDBID , $CHAINID )

Returns the length of the specified chain.

=cut

sub get_chain_length {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self    = shift;
    my $pdbid   = _to_string(shift);
    my $chainid = _to_string(shift);
    my $ret     = $self->_call(
        'getChainLength', $pdbid, $chainid
    );
    return $ret;
}

=item get_chains( $PDBID )

Returns a list of all the chain identifiers for a given structure, or a
reference to such a list in scalar context.

=cut

sub get_chains {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getChains', $pdbid
    );
    return wantarray ? @$ret : $ret;
}

=item get_cif_chain( $PDBID , $CHAINID )

Converts the specified author-assigned chain identifier to its mmCIF
equivalent.

=cut

sub get_cif_chain {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self    = shift;
    my $pdbid   = _to_string(shift);
    my $chainid = _to_string(shift);
    my $ret     = $self->_call(
        'getCifChain', $pdbid, $chainid
    );
    return $ret;
}

=item get_cif_chain_length( $PDBID , $CHAINID )

Returns the length of the specified chain, just like C<get_chain_length>,
except it expects the chain identifier to be the mmCIF version.

=cut

sub get_cif_chain_length {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self    = shift;
    my $pdbid   = _to_string(shift);
    my $chainid = _to_string(shift);
    my $ret     = $self->_call(
        'getCifChainLength', $pdbid, $chainid
    );
    return $ret;
}

=item get_cif_chains( $PDBID )

Returns a list of all the mmCIF chain identifiers for a given structure, or
a reference to such a list in scalar context.

=cut

sub get_cif_chains {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getCifChains', $pdbid
    );
    return wantarray ? @$ret : $ret;
}

=item get_cif_residue( $PDBID , $CHAINID , $RESIDUEID )

Converts the specified author-assigned residue identifier to its mmCIF
equivalent.

=cut

sub get_cif_residue {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self      = shift;
    my $pdbid     = _to_string(shift);
    my $chainid   = _to_string(shift);
    my $residueid = _to_string(shift);
    my $ret       = $self->_call(
        'getCifResidue', $pdbid, $chainid, $residueid
    );
    return $ret;
}

=item get_current_pdbids( )

Returns a list of the identifiers (PDB IDs) corresponding to "current"
structures (i.e. not obsolete, models, etc.), or a reference to such a
list in scalar context.

=cut

sub get_current_pdbids {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self = shift;
    my $ret  = $self->_call(
        'getCurrentPdbIds'
    );
    return wantarray ? @$ret : $ret;
}

=item get_ec_nums( @PDBIDS )

=item get_ec_nums( )

Retrieves the Enzyme Classification (EC) numbers associated with the
specified PDB IDs or with all PDB structures if called with no arguments. 

=cut

sub get_ec_nums {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self = shift;
    my $ret;
    if(@_) {
        my @pdbids = map(_to_string($_), map { ref($_) ? @$_ : $_ } @_);
        $ret = $self->_call(
            'getEcNumsForStructures', \@pdbids
        );
    }
    else {
    	$ret = $self->_call(
    	    'getEcNums'
    	);
    }
    $_ = bless(\"$_", 'WWW::PDB::_EcNumsResult') for @$ret;
    return wantarray ? @$ret : $ret;
}

=item get_entities( $PDBID )

Returns a list of the entity IDs for a given structure, or a reference
to such a list in scalar context.

=cut

sub get_entities {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getEntities', $pdbid
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item get_genome_details( )

Retrieves genome details for all PDB structures.

=cut

sub get_genome_details {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self = shift;
    my $ret  = $self->_call(
        'getGenomeDetails'
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item get_kabsch_sander( $PDBID , $CHAINID )

Finds secondary structure for the given chain.

=cut

sub get_kabsch_sander {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self    = shift;
    my $pdbid   = _to_string(shift);
    my $chainid = _to_string(shift);
    my $ret     = $self->_call(
        'getKabschSander', $pdbid, $chainid
    );
    return $ret;
}

=item get_obsolete_pdbids( )

Returns a list of the identifiers (PDB IDs) corresponding to obsolete
structures, or a reference to such a list in scalar context.

=cut

sub get_obsolete_pdbids {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self = shift;
    my $ret  = $self->_call(
        'getObsoletePdbIds'
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item get_primary_citation_title( $PDBID )

Finds the title of the specified structure's primary citation (if it has
one).

=cut

sub get_primary_citation_title {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getPrimaryCitationTitle', $pdbid
    );
    return $ret;
}

=item get_pubmed_ids( )

Retrieves the PubMed IDs associated with all PDB structures.

=cut

sub get_pubmed_ids {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self = shift;
    my $ret  = $self->_call(
        'getPubmedIdForAllStructures'
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item get_pubmed_id( $PDBID )

Retrieves the PubMed ID associated with the specified structure.

=cut

sub get_pubmed_id {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getPubmedIdForStructure', $pdbid
    );
    return $ret;
}

=item get_release_dates( @PDBIDS )

Maps the given PDB IDs to their release dates.

=cut

sub get_release_dates {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self   = shift;
    my @pdbids = map(_to_string($_), map { ref($_) ? @$_ : $_ } @_);
    my $ret    = $self->_call(
        'getReleaseDates', \@pdbids
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item get_sequence( $PDBID , $CHAINID )

Retrieves the sequence of the specified chain.

=cut

sub get_sequence {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self    = shift;
    my $pdbid   = _to_string(shift);
    my $chainid = _to_string(shift);
    my $ret = $self->_call(
        'getSequenceForStructureAndChain', $pdbid, $chainid
    );
    return $ret;
}

=item get_space_group( $PDBID )

Returns the space group of the specified structure (the
C<symmetry.space_group_name_H_M> field according to the mmCIF dictionary).

=cut

sub get_space_group {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getSpaceGroupForStructure', $pdbid
    );
    return $ret;
}

=item homology_reduction_query( @PDBIDS , $CUTOFF )

Reduces the set of PDB IDs given as input based on sequence homology.

=cut

sub homology_reduction_query {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self   = shift;
    my $cutoff = _to_int(int(pop));
    my @pdbids = map(_to_string($_), map { ref($_) ? @$_ : $_ } @_);
    my $ret    = $self->_call(
        'homologyReductionQuery', \@pdbids, $cutoff
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item keyword_query( $KEYWORD_EXPR [, $EXACT_MATCH [, $AUTHORS_ONLY ] ] )

Runs a keyword query with the specified expression. Search can be made
stricter by requiring an exact match or restricting the search to
authors. Both boolean arguments are optional and default to false. Returns
a list of PDB IDs or a reference to such a list in scalar context.

=cut

sub keyword_query {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self         = shift;
    my $keyword      = _to_string(shift);
    my $exact_match  = _to_boolean(shift || 0);
    my $authors_only = _to_boolean(shift || 0);
    my $ret          = $self->_call(
        'keywordQuery', $keyword, $exact_match, $authors_only
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item pubmed_abstract_query( $KEYWORD_EXPR )

Runs a keyword query on PubMed Abstracts. Returns a list of PDB IDs or
a reference to such a list in scalar context.

=cut

sub pubmed_abstract_query {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self    = shift;
    my $keyword = _to_string(shift);
    my $ret     = $self->_call(
        'pubmedAbstractQuery', $keyword
    );
    return $ret && wantarray ? @$ret : $ret;
}

=back

=head3 UNTESTED

The following methods are defined by the PDB web services interface, so
they are wrapped here, but they have not been tested.

=over 4

=item get_annotations( $STATE_FILE )

Given a string in the format of a ViewState object from Protein
Workshop, returns another ViewState object.

=cut

sub get_annotations {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self       = shift;
    my $state_file = _to_string(shift);
    my $ret        = $self->_call(
        'getAnnotations', $state_file
    );
    return $ret;
}

=item get_atom_site( $PDBID )

Returns the first atom site object for a structure.

=cut

sub get_atom_site {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getAtomSite', $pdbid
    );
    return $ret;
}

=item get_atom_sites( $PDBID )

Returns the atom site objects for a structure.

=cut

sub get_atom_sites {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getAtomSites', $pdbid
    );
    return $ret;
}

=item get_domain_fragments( $PDBID , $CHAINID , $METHOD )

Finds all structural protein domain fragments for a given structure.

=cut

sub get_domain_fragments {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self    = shift;
    my $pdbid   = _to_string(shift);
    my $chainid = _to_string(shift);
    my $method  = _to_string(shift);
    my $ret = $self->_call(
        'getDomainFragments', $pdbid, $chainid, $method
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item get_first_struct_conf( $PDBID )

Finds the first struct_conf for the given structure.

=cut

sub get_first_struct_conf {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getFirstStructConf', $pdbid
    );
    return $ret;
}

=item get_first_struct_sheet_range( $PDBID )

Finds the first struct_sheet_range for the given structure.

=cut

sub get_first_struct_sheet_range {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getFirstStructSheetRange', $pdbid
    );
    return $ret;
}

=item get_struct_confs( $PDBID )

Finds the struct_confs for the given structure.

=cut

sub get_struct_confs {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self  = shift;
    my $pdbid = _to_string(shift);
    my $ret   = $self->_call(
        'getStructConfs', $pdbid
    );
    return $ret;
}

=item get_struct_sheet_ranges( $PDBID )

Finds the struct_sheet_ranges for the given structure.

=cut

sub get_struct_sheet_ranges {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my($self, $pdbid) = _wrap(\@_ => [SELF, STRING]);
    my $ret   = $self->_call(
        'getStructSheetRanges', $pdbid
    );
    return $ret;
}

=item get_structural_genomics_pdbids( )

Finds info for structural genomics structures.

=cut

sub get_structural_genomics_pdbids {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my $self = _wrap(\@_ => [SELF]);
    my $ret  = $self->_call(
        'getStructureGenomicsPdbIds'
    );
    return $ret && wantarray ? @$ret : $ret;
}

=item xml_query( $XML )

Runs any query that can be constructed, pretty much.

=cut

sub xml_query {
    unshift @_, __PACKAGE__ # add the package name unless already there
        unless defined($_[0]) && UNIVERSAL::isa($_[0], __PACKAGE__);

    my($self, $xml) = _wrap(\@_ => [SELF, STRING]);
    my $ret  = $self->_call(
        'xmlQuery', $xml
    );
    return $ret && wantarray ? @$ret : $ret;
}

=back

=cut

################################################################################

sub _get_file {
    my($self, @dir) = @_;
    my $file        = pop @dir;
    my($dir, $local_path, $store, $fh);
    if($self->cache) {
        $dir        = File::Spec->catfile($self->cache, @dir);
        $local_path = File::Spec->catfile($dir, $file);
    }
    unless($self->cache && ($store = new IO::File($local_path))) {
        my $ftp;
        if(   ($ftp = new Net::FTP($self->ftp, Debug => 0)) # connect
            && $ftp->login(qw(anonymous -anonymous@))       # login
            && $ftp->cwd(join('', map("/$_", @dir)))        # chdir
        ) {
            # store in temporary file unless there's a cache
            $store = IO::File->new_tmpfile unless $self->cache # cache exists
                && File::Path::mkpath($dir)                    # mkdir
                && ($store = new IO::File($local_path, '+>')); # create file
            
            # seek to start if successful get otherwise delete file
            if($ftp->get($file => $store)) {
                seek($store, 0, SEEK_SET);
            }
            else {
                undef $store;
                $self->cache and unlink $local_path;
            }
            
            # clean up
            $ftp->quit;
        }
    }
    
    # if file stored, decompress it
    if($store) {
        $fh = IO::File->new_tmpfile;
        IO::Uncompress::Gunzip::gunzip($store => $fh);
        seek($fh, 0, SEEK_SET);
        close $store;
    }
    
    return $fh;
}

sub _call {
    my $self   = shift;
    my $result = $self->soap->call(@_);
    confess $result->faultstring if $result->fault;
    return $result->result;
}

sub _wrap {
    my @data = @{shift()};
    my @type = @{shift()};  
    return map {
        my $type = shift @type;
        if($type == BOOLEAN) {
            $_ = SOAP::Data->type(boolean => ($_ ? 1 : 0));
        }
        elsif($type == DOUBLE) {
            $_ = SOAP::Data->type(double => $_);
        }
        elsif($type == INT) {
            $_ = SOAP::Data->type('int' => $_);
        }
        elsif($type == STRING) {
            $_ = SOAP::Data->type(string => $_);
        }
    $_} @data;
}

sub _to_int {
    my $var = shift;
    return SOAP::Data->type('int' => $var);
}

sub _to_string {
    my $var = shift;
    return SOAP::Data->type(string => $var);
}

sub _to_boolean {
    my $var = shift;
    return SOAP::Data->type(boolean => ($var ? 1 : 0));
}

sub _to_double {
    my $var = shift;
    return SOAP::Data->type(double => $var);
}

################################################################################

package WWW::PDB::_FastaResult;

use overload '""' => sub { shift->pdbid };

sub pdbid  {
    return substr(${$_[0]}, 0, 4);
}

sub cutoff {
    return substr(${$_[0]}, 5);
}

################################################################################

package WWW::PDB::_EcNumsResult;

use overload '""' => sub { scalar shift->ec };

sub pdbid {
    return substr(${$_[0]}, 0, 4);
}

sub chainid {
    return substr(${$_[0]}, 5, 1);
}

sub ec {
    local $_ = substr(${$_[0]}, 7);
    return wantarray ? split(', ', $_) : $_;
}

################################################################################

1;

__END__

=head1 REFERENCES

=over 4

=item 1.

Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N.,
Weissig, H., Shindyalov, I. N. & Bourne, P. E. (2000).
I<Nucleic Acids Res.> B<28>(1), 235-242.

=item 2.

Bernstein, F. C., Koetzle, T. F., Williams, G. J. B., Meyer, Jr., E. F.,
Brice, M. D., Rodgers, J. R., Kennard, O., Shimanouchi, T. & Tasumi,
M. (1977). I<Eur. J. Biochem.> B<80>(2), 319-324.

=back

=head1 SEE ALSO

The PDB can be accessed via the web at L<http://www.pdb.org/>. The
Java API documentation for the PDB's web services is located at
L<http://www.rcsb.org/robohelp_f/webservices/pdbwebservice.html>.

=head1 BUGS

Please report them:
L<http://rt.cpan.org/Public/Dist/Display.html?Name=WWW-PDB>

=head1 AUTHOR

Miorel-Lucian Palii, E<lt>mlpalii@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008-2009 by Miorel-Lucian Palii

This library is free software; you can redistribute it and/or modify it
under the same terms as Perl itself, either Perl version 5.8.8 or, at
your option, any later version of Perl 5 you may have available.

=cut
