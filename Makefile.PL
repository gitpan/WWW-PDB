use 5.008;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME              => 'WWW::PDB',
    VERSION_FROM      => 'lib/WWW/PDB.pm',
    LICENSE           => 'perl',
    PREREQ_PM         => {
                             Carp                   => 0,
                             Fcntl                  => 0,
                             File::Path             => '1.08',
                             File::Spec             => '3.12',
                             IO::File               => '1.13',
                             IO::Uncompress::Gunzip => '2.015',
                             Net::FTP               => '2.77',
                             SOAP::Lite             => '0.710.08',
                         },
    ($] >= 5.005 ?
      (ABSTRACT_FROM  => 'lib/WWW/PDB.pm',
       AUTHOR         => 'Miorel-Lucian Palii <mlpalii@gmail.com>') : ()),
);

