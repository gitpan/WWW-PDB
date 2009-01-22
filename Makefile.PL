use 5.006;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME              => 'WWW::PDB',
    VERSION_FROM      => 'lib/WWW/PDB.pm',
    LICENSE           => 'perl',
    PREREQ_PM         => {
                             Carp                   => 0,
                             Exporter               => 0,
                             Fcntl                  => 0,
                             File::Path             => 0,
                             File::Spec             => 0,
                             IO::File               => 0,
                             IO::Uncompress::Gunzip => 0,
                             Net::FTP               => 0,
                             SOAP::Lite             => 0,
                         },
    ($] >= 5.005 ?
      (ABSTRACT_FROM  => 'lib/WWW/PDB.pm',
       AUTHOR         => 'Miorel-Lucian Palii <mlpalii@gmail.com>') : ()),
);

