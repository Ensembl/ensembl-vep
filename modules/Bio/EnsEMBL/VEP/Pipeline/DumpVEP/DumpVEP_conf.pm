=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::VEP::Pipeline::DumpVEP::DumpVEP_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
 # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly
use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

sub default_options {
  my ($self) = @_;

  # The hash returned from this function is used to configure the
  # pipeline, you can supply any of these options on the command
  # line to override these default values.
    
  # You shouldn't need to edit anything in this file other than
  # these values, if you find you do need to then we should probably
  # make it an option here, contact the variation team to discuss
  # this - patches are welcome!

  return {

    %{ $self->SUPER::default_options()
      },    # inherit other stuff from the base class
    hive_force_init => 1,
    hive_use_param_stack => 0,
    hive_use_triggers => 0,
    hive_auto_rebalance_semaphores => 0, 
    hive_no_init => 0,
    # a name for your pipeline (will also be used in the name of the hive database)    
    pipeline_name           => 'dump_vep',

    # a directory to keep hive output files and your registry file, you should
    # create this if it doesn't exist
    pipeline_dir            => '/hps/nobackup/production/ensembl/'.$ENV{'USER'}.'/'.$self->o('pipeline_name').'/'.$self->o('ensembl_release'),

    # contains frequency data
    data_dir                => '/nfs/production/panda/ensembl/variation/data/',
    dump_vep_data_dir       => $self->o('data_dir') . '/dump_vep',
        
    # dump databases of this version number
    ensembl_release => undef,
    eg_version => undef,
    
    # add refseq, merged dumps?
    refseq => 1,
    merged => 1,
    
    # tabix-convert species with var DBs?
    # this creates an extra tar.gz file for each of these species
    # the web interface and REST API use these in preference to the non-converted ones
    convert => 1,
   
    # Update since release 95 run the vep pipeline with the `-division vertebrates` flag. If you only want to run the pipeline for human, you will need these two flags `-species homo_sapiens` and `-division vertebrates`
 
    # include or exclude the following species from the dumps, run for a division or all the species on the server
    species => [],
    antispecies => [],
    division    => [],
    run_all     => 0,

    # include LRGs in dumps
    lrg => 1,

    # don't change this unless you know what you're doing!!!
    region_size => 1e6,
    
    # special flags apply to certain species
    species_flags => {
      
      # human has SIFT, PolyPhen, regulatory data and a file containing SNP frequencies
      homo_sapiens => {
        
        # assembly-specific stuff
        assembly_specific => {
          GRCh37 => {
            bam => $self->o('dump_vep_data_dir').'/GCF_000001405.25_GRCh37.p13_knownrefseq_alns.bam',
            freq_vcf => [
              {
                file => $self->o('dump_vep_data_dir').'/1KG.phase3.GRCh37.vcf.gz',
                pops => [qw(AFR AMR EAS EUR SAS)],
                name => '1000genomes',
                version => 'phase3'
              },
              {
                file => $self->o('dump_vep_data_dir').'/ESP6500SI-V2-SSA137.vcf.gz',
                pops => [qw(AA EA)],
                name => 'ESP',
              },
              # {
              #   file => $self->o('data_dir').'/ExAC.0.3.GRCh37.vcf.gz',
              #   pops => ['', qw(AFR AMR Adj EAS FIN NFE OTH SAS)],
              #   name => 'ExAC',
              #   prefix => 'ExAC',
              #   version => 0.3,
              # },
              {
                file => $self->o('data_dir').'/gnomAD/v2.1/grch37/exomes/gnomad.exomes.r2.1.sites.chr+++CHR+++_noVEP.vcf.gz',
                pops => ['', qw(afr amr asj eas fin nfe oth sas)],
                name => 'gnomAD',
                prefix => 'gnomAD',
                version => 'r2.1',
              },
            ],
          },
          GRCh38 => {
            bam => $self->o('dump_vep_data_dir').'/GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bam',
            freq_vcf => [
              {
                file => $self->o('dump_vep_data_dir').'/1KG.phase3.GRCh38_2018_02_26.vcf.gz',
                pops => [qw(AFR AMR EAS EUR SAS)],
                name => '1000genomes',
                version => 'phase3'
              },
              {
                file => $self->o('dump_vep_data_dir').'/ESP6500SI-V2-SSA137_GRCh38.vcf.gz',
                pops => [qw(AA EA)],
                name => 'ESP',
                version => 'V2-SSA137',
              },
              # {
              #   file => $self->o('data_dir').'/ExAC.0.3.GRCh38.vcf.gz',
              #   pops => ['', qw(AFR AMR Adj EAS FIN NFE OTH SAS)],
              #   name => 'ExAC',
              #   prefix => 'ExAC',
              #   version => 0.3,
              # },
              {
                file => $self->o('data_dir').'/gnomAD/v2.1.1/grch38/exomes/gnomad.exomes.r2.1.1.sites.+++CHR+++.liftover_grch38_no_VEP.vcf.gz',
                pops => ['', qw(afr amr asj eas fin nfe oth sas)],
                name => 'gnomAD',
                prefix => 'gnomAD',
                version => 'r2.1.1',
                use_chr_prefix => 1,
              },
            ],
          },
        }
      },
    },

    # configuration for the various resource options used in the pipeline
    # EBI farm users should either change these here, or override them on the
    # command line to suit the EBI farm. The names of each option hopefully
    # reflect their usage, but you may want to change the details (memory
    # requirements, queue parameters etc.) to suit your own data
        
    default_lsf_options => '-q production-rh74 -R"select[mem>4000] rusage[mem=4000]" -M4000',
    urgent_lsf_options  => '-q production-rh74 -R"select[mem>2000] rusage[mem=2000]" -M2000',
    highmem_lsf_options => '-q production-rh74 -R"select[mem>15000] rusage[mem=15000]" -M15000', # this is Sanger LSF speak for "give me 15GB of memory"
    long_lsf_options    => '-q production-rh74 -R"select[mem>2000] rusage[mem=2000]" -M2000',
    
    debug => 0,
    qc => 1,
  };
}


sub resource_classes {
  my ($self) = @_;
  return {
    'default' => { 'LSF' => $self->o('default_lsf_options') },
    'urgent'  => { 'LSF' => $self->o('urgent_lsf_options')  },
    'highmem' => { 'LSF' => $self->o('highmem_lsf_options') },
    'long'    => { 'LSF' => $self->o('long_lsf_options')    },
  };
}

# Override the default method, to force an automatic loading of the registry in all workers
sub beekeeper_extra_cmdline_options {
  my ($self) = @_;
  return
      ' -reg_conf ' . $self->o('registry'),
  ;
}

sub pipeline_wide_parameters {  # these parameter values are visible to all analyses, can be overridden by parameters{} and input_id{}
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class

        'refseq'     => $self->o('refseq'),
        'ensembl_release' => $self->o('ensembl_release'),
        'eg_version' => $self->o('eg_version'),
        'pipeline_dir' => $self->o('pipeline_dir'),
        'debug' => $self->o('debug'),
        'region_size' => $self->o('region_size'),
        'division' => $self->o('division'),
    };
}

sub pipeline_analyses {
  my ($self) = @_;
   
  my @analyses = (
    {
      -logic_name    => 'species_factory',
      -module         => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
      -parameters    => {
        species     => $self->o('species'),
        antispecies => $self->o('antispecies'),
        division    => $self->o('division'),
        run_all     => $self->o('run_all'),
      },
      -input_ids     => [{}],
      -rc_name       => 'default',
      -hive_capacity => 1,
      -max_retry_count => 0,
      -flow_into     => {
        '1' =>  $self->o('debug') ? [] : ['distribute_dumps' => {}],
        '2' => ['init_dump_vep_core'],
        '7' => WHEN(
            '#refseq#' => 'init_dump_vep_otherfeatures'
          )
      },
    },

    {
     -logic_name    => 'init_dump_vep_core',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::InitDump',
      -parameters    => {
        group           => 'core',
      },
      -rc_name       => 'default',
      -hive_capacity => 1,
      -max_retry_count => 0,
      -flow_into     => {
        '2' => ['create_dump_jobs'],
      },
    },
    
    {
     -logic_name    => 'init_dump_vep_otherfeatures',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::InitDump',
      -parameters    => {
        group           => 'otherfeatures',
      },
      -rc_name       => 'default',
      -hive_capacity => 1,
      -max_retry_count => 0,
      -can_be_empty   => 1,
      -flow_into     => {
        '2' => ['create_dump_jobs'],
      },
    },

    {
      -logic_name    => 'create_dump_jobs',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::CreateDumpJobs',
      -parameters    => {
        merged          => $self->o('merged'),
        lrg             => $self->o('lrg'),
      },
      -rc_name       => 'default',
      -analysis_capacity => 20,
      -max_retry_count => 0,
      -flow_into     => {
        '3' => ['dump_vep_core'],
        '4' => ['dump_vep_otherfeatures'],
        '5' => ['dump_vep_variation'],
        '6' => ['dump_vep_regulation'],

        '7' => ['merge_vep'],
        '8' => ['join_vep'],
        '9' => ['qc_vep'],
        '10' => $self->o('debug') ? [] : ['finish_dump'],
      },
      -wait_for      => ['init_dump_vep_core','init_dump_vep_otherfeatures'],
    },

    ## these analyses do the feature dumping
    {
      -logic_name    => 'dump_vep_core',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Core',
      -parameters    => {
        species_flags  => $self->o('species_flags'),
      },
      -rc_name       => 'default',
      -analysis_capacity => 10,
      -max_retry_count => 0,
    },
    {
      -logic_name    => 'dump_vep_otherfeatures',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Otherfeatures',
      -parameters    => {
        species_flags  => $self->o('species_flags'),
      },
      -rc_name       => 'default',
      -analysis_capacity => 10,
      -can_be_empty   => 1,
      -max_retry_count => 0,
      -wait_for      => ['create_dump_jobs'],
    },
    {
      -logic_name    => 'dump_vep_variation',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Variation',
      -parameters    => {
        species_flags  => $self->o('species_flags'),
        convert        => $self->o('convert'),
      },
      -rc_name       => 'default',
      -analysis_capacity => 10,
      -can_be_empty   => 1,
      -max_retry_count => 0,
    },
    {
      -logic_name    => 'dump_vep_regulation',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::Dumper::Regulation',
      -parameters    => {
        species_flags  => $self->o('species_flags'),
      },
      -rc_name       => 'default',
      -analysis_capacity => 10,
      -can_be_empty   => 1,
      -max_retry_count => 0,
    },

    {
      -logic_name    => 'merge_vep',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::MergeVEP',
      -parameters    => {
      },
      -wait_for      => ['dump_vep_core', 'dump_vep_otherfeatures'],
      -analysis_capacity => 20,
      -rc_name       => 'default',
      -failed_job_tolerance => 0,
      -can_be_empty   => 1,
    },

    {
      -logic_name    => 'join_vep',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::JoinVEP',
      -parameters    => {
        convert => $self->o('convert'),
      },
      -wait_for      => ['dump_vep_core', 'dump_vep_otherfeatures', 'dump_vep_variation', 'dump_vep_regulation', 'merge_vep'],
      -analysis_capacity => 20,
      -rc_name       => 'default',
      -failed_job_tolerance => 0,
      -max_retry_count => 0,
    },

    {
      -logic_name    => 'qc_vep',
      -module        => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::QCDump',
      -parameters    => {
        convert => $self->o('convert'),
        lrg     => $self->o('lrg'),
      },
      -wait_for      => ['join_vep'],
      -analysis_capacity => 50,
      -rc_name       => 'default',
      -failed_job_tolerance => 0,
      -max_retry_count => 0,
    },

    {
      -logic_name => 'finish_dump',
      -module     => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::FinishDump',
      -parameters => { },
      -wait_for   => ['join_vep'],
      -max_retry_count => 0,
    },

    {
      -logic_name => 'distribute_dumps',
      -module     => 'Bio::EnsEMBL::VEP::Pipeline::DumpVEP::DistributeDumps',
      -parameters => { },
      -wait_for   => ['join_vep'],
      -max_retry_count => 0,
    },
  );

  return \@analyses;
}

1;

