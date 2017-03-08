=head1 LICENSE

Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig - BigWig annotation source

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::BigWig;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::IO::Parser::BigWig;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::File);

sub parser {
  my $self = shift;

  unless(exists($self->{parser})) {

    # This is something of a hack.
    # The BigWig parser uses the Kent source tree which calls exit() if the file doesn't exist.
    # We can't use -e $self->file as the file might be a remote one.
    # So, we fork, attempt to open the file, and check the return code of the forked process.
    my $pid = fork();

    # parent process, capture the return code and die nicely if it failed
    if($pid) {
      if(waitpid($pid, 0) > 0) {
        my $rc = $? >> 8;
        throw("ERROR: Failed to open ".$self->file."\n") if $rc > 0;
      }
    }

    # child process, attempt to open the file
    else {
      my $test = Bio::EnsEMBL::IO::Parser::BigWig->open($self->file);
      exit(0);
    }

    $self->{parser} = Bio::EnsEMBL::IO::Parser::BigWig->open($self->file);
  }
  return $self->{parser};
}

sub valid_chromosomes {
  my $self = shift;
  return $self->{valid_chromosomes} ||= [keys %{$self->parser->{cache}->{chromosomes}}];
}

sub _get_record_name {
  my $self = shift;
  my $parser = $self->parser;

  return $self->report_coords ?
    sprintf(
      '%s:%i-%i',
      $parser->get_seqname,
      $parser->get_start,
      $parser->get_end
    ) :
    $parser->get_score;
}

sub _record_overlaps_VF {
  my $self = shift;
  my $vf = shift;

  my $parser = $self->parser();
  my $type = $self->type();

  if($type eq 'overlap') {
    return overlap($parser->get_start + 1, $parser->get_end, $vf->{start}, $vf->{end});
  }
  elsif($type eq 'exact') {
    return $parser->get_start + 1 == $vf->{start} && $parser->get_end == $vf->{end};
  }
}

1;