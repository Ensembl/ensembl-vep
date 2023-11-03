requires 'DBI';
requires 'Set::IntervalTree';
requires 'JSON';
requires 'Text::CSV';
recommends 'DBD::mysql', '<= 4.050'; # newer versions do not support MySQL 5
recommends 'PerlIO::gzip';
recommends 'IO::Uncompress::Gunzip';
recommends 'Bio::DB::BigFile';
recommends 'Sereal';
recommends 'HTML::Lint';
recommends 'Capture::Tiny';
