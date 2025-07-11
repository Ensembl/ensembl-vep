=head1 LICENSE

Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

# EnsEMBL module for Bio::EnsEMBL::VEP::Constants
#
#

=head1 NAME

Bio::EnsEMBL::VEP::Constants - Constants used by VEP classes

=head1 SYNOPSIS

Not invoked directly

=head1 DESCRIPTION

Contains version data, column headers and various other constants

=cut

package Bio::EnsEMBL::VEP::Constants;

use strict;
use warnings;

use base qw(Exporter);

our $VEP_VERSION     = 114;
our $VEP_SUB_VERSION = 2;

our @EXPORT_OK = qw(
  @FLAG_FIELDS
  %FIELD_DESCRIPTIONS
  @DEFAULT_OUTPUT_COLS
  $MAX_NOT_ORDERED_VARIANTS
);

# contains an ordered map between command line flags and output columns
our @FLAG_FIELDS = (

  # general
  { flag => 'individual',      fields => ['IND','ZYG'] },
  { flag => 'individual_zyg',  fields => ['ZYG'] },
  { flag => 'allele_number',   fields => ['ALLELE_NUM'] },
  { flag => 'show_ref_allele', fields => ['REF_ALLELE'] },
  { flag => 'uploaded_allele', fields => ['UPLOADED_ALLELE'] },
  { flag => 'user',            fields => ['IMPACT','DISTANCE','STRAND','FLAGS'] },
  { flag => 'flag_pick',       fields => ['PICK'] },
  { flag => 'flag_pick_allele',fields => ['PICK'] },
  { flag => 'flag_pick_allele_gene', fields => ['PICK'] },
  { flag => 'variant_class',   fields => ['VARIANT_CLASS']},
  { flag => 'minimal',         fields => ['MINIMISED']},
  { flag => 'spdi',            fields => ['SPDI']},
  { flag => 'ga4gh_vrs',       fields => ['GA4GH_VRS']},

  # gene-related
  { flag => 'symbol',          fields => ['SYMBOL','SYMBOL_SOURCE','HGNC_ID'] },
  { flag => 'biotype',         fields => ['BIOTYPE'] },
  { flag => 'canonical',       fields => ['CANONICAL'] },
  { flag => 'mane_select',     fields => ['MANE', 'MANE_SELECT']},
  { flag => 'mane',            fields => ['MANE', 'MANE_SELECT', 'MANE_PLUS_CLINICAL']},
  { flag => 'tsl',             fields => ['TSL']},
  { flag => 'flag_gencode_primary', fields => ['GENCODE_PRIMARY']},
  { flag => 'appris',          fields => ['APPRIS']},
  { flag => 'ccds',            fields => ['CCDS'] },
  { flag => 'protein',         fields => ['ENSP'] },
  { flag => 'uniprot',         fields => ['SWISSPROT', 'TREMBL', 'UNIPARC', 'UNIPROT_ISOFORM'] },
  { flag => 'xref_refseq',     fields => ['RefSeq'] },
  { flag => 'refseq',          fields => ['REFSEQ_MATCH', 'REFSEQ_OFFSET'] },
  { flag => 'merged',          fields => ['REFSEQ_MATCH', 'SOURCE', 'REFSEQ_OFFSET'] },
  { flag => 'use_transcript_ref', fields => ['GIVEN_REF', 'USED_REF']},
  { flag => 'bam_edited',      fields => ['BAM_EDIT']},
  { flag => 'custom',          fields => ['SOURCE'] },
  { flag => 'gene_phenotype',  fields => ['GENE_PHENO'] },
  { flag => 'nearest',         fields => ['NEAREST'] },
  { flag => 'ambiguity',       fields => ['AMBIGUITY'] },

  # non-synonymous predictions
  { flag => 'sift',            fields => ['SIFT'] },
  { flag => 'polyphen',        fields => ['PolyPhen'] },

  # transcript/protein stuff
  { flag => 'numbers',         fields => ['EXON','INTRON'] },
  { flag => 'domains',         fields => ['DOMAINS'] },
  { flag => 'mirna',           fields => ['miRNA'] },
  { flag => 'hgvs',            fields => ['HGVSc','HGVSp','HGVS_OFFSET'] },
  { flag => 'hgvsg',           fields => ['HGVSg'] },

  # frequency stuff
  { flag => 'af',              fields => ['AF'] },
  { flag => 'af_1kg',          fields => ['AF','AFR_AF','AMR_AF','EAS_AF','EUR_AF','SAS_AF'] },
  { flag => 'af_gnomad',       fields => ['gnomADe_AF','gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_ASJ_AF','gnomADe_EAS_AF','gnomADe_FIN_AF', 'gnomADe_MID_AF', 'gnomADe_NFE_AF','gnomADe_REMAINING_AF','gnomADe_SAS_AF'] },
  { flag => 'af_gnomade',      fields => ['gnomADe_AF','gnomADe_AFR_AF','gnomADe_AMR_AF','gnomADe_ASJ_AF','gnomADe_EAS_AF','gnomADe_FIN_AF', 'gnomADe_MID_AF', 'gnomADe_NFE_AF','gnomADe_REMAINING_AF','gnomADe_SAS_AF'] },
  { flag => 'af_gnomadg',      fields => ['gnomADg_AF','gnomADg_AFR_AF', ,'gnomADg_AMI_AF', 'gnomADg_AMR_AF','gnomADg_ASJ_AF','gnomADg_EAS_AF','gnomADg_FIN_AF','gnomADg_MID_AF','gnomADg_NFE_AF','gnomADg_REMAINING_AF','gnomADg_SAS_AF'] },
  { flag => 'max_af',          fields => ['MAX_AF', 'MAX_AF_POPS'] },
  { flag => 'check_frequency', fields => ['FREQS'] },

  # misc variation stuff
  { flag => 'check_existing',  fields => ['CLIN_SIG','SOMATIC','PHENO'] },
  { flag => 'pubmed',          fields => ['PUBMED'] },
  { flag => 'check_svs',       fields => ['SV'] },
  { flag => 'dont_skip',       fields => ['CHECK_REF'] },
  { flag => 'overlaps',        fields => ['OverlapBP', 'OverlapPC']},
  { flag => 'shift_length',    fields => ['SHIFT_LENGTH']},
  { flag => 'var_synonyms',    fields => ['VAR_SYNONYMS']},

  # regulatory
  { flag => 'regulatory',      fields => ['BIOTYPE','MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE', 'TRANSCRIPTION_FACTORS'] },
  { flag => 'cell_type',       fields => ['CELL_TYPE'] },
);


# field descriptions for output headers
our %FIELD_DESCRIPTIONS = (
  'Uploaded_variation' => 'Identifier of uploaded variant',
  'ID'                 => 'Identifier of uploaded variant',
  'Location'           => 'Location of variant in standard coordinate format (chr:start or chr:start-end)',
  'Allele'             => 'The variant allele used to calculate the consequence',
  'Gene'               => 'Stable ID of affected gene',
  'Feature'            => 'Stable ID of feature',
  'Feature_type'       => 'Type of feature - Transcript, RegulatoryFeature or MotifFeature',
  'Consequence'        => 'Consequence type',
  'cDNA_position'      => 'Relative position of base pair in cDNA sequence',
  'CDS_position'       => 'Relative position of base pair in coding sequence',
  'Protein_position'   => 'Relative position of amino acid in protein',
  'Amino_acids'        => 'Reference and variant amino acids',
  'Codons'             => 'Reference and variant codon sequence',
  'Existing_variation' => 'Identifier(s) of co-located known variants',
  'IMPACT'             => 'Subjective impact classification of consequence type',
  'CANONICAL'          => 'Indicates if transcript is canonical for this gene',
  'MANE'               => 'MANE (Matched Annotation from NCBI and EMBL-EBI) set(s) the transcript belongs to',
  'MANE_SELECT'        => 'MANE Select (Matched Annotation from NCBI and EMBL-EBI) Transcript',
  'MANE_PLUS_CLINICAL' => 'MANE Plus Clinical (Matched Annotation from NCBI and EMBL-EBI) Transcript',
  'TSL'                => 'Transcript support level',
  'APPRIS'             => 'Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods',
  'CCDS'               => 'Indicates if transcript is a CCDS transcript',
  'SYMBOL'             => 'Gene symbol (e.g. HGNC)',
  'SYMBOL_SOURCE'      => 'Source of gene symbol',
  'SOURCE'             => 'Source of transcript',
  'HGNC_ID'            => 'Stable identifer of HGNC gene symbol',
  'ENSP'               => 'Protein identifer',
  'FLAGS'              => 'Transcript quality flags',
  'SWISSPROT'          => 'UniProtKB/Swiss-Prot accession',
  'TREMBL'             => 'UniProtKB/TrEMBL accession',
  'UNIPARC'            => 'UniParc accession',
  'UNIPROT_ISOFORM'    => 'Direct mappings to UniProtKB isoforms',
  'NEAREST'            => 'Identifier(s) of nearest transcription start site',
  'miRNA'              => 'SO terms of overlapped miRNA secondary structure feature(s)',
  'HGVSc'              => 'HGVS coding sequence name',
  'HGVSp'              => 'HGVS protein sequence name',
  'HGVSg'              => 'HGVS genomic sequence name',
  'SPDI'               => 'Genomic SPDI notation',
  'GA4GH_VRS'          => 'GA4GH Variation Representation Specification',
  'SIFT'               => 'SIFT prediction and/or score',
  'PolyPhen'           => 'PolyPhen prediction and/or score',
  'EXON'               => 'Exon number(s) / total',
  'INTRON'             => 'Intron number(s) / total',
  'DOMAINS'            => 'The source and identifer of any overlapping protein domains',
  'MOTIF_NAME'         => 'The stable identifier of a transcription factor binding profile (TFBP) aligned at this position',
  'MOTIF_POS'          => 'The relative position of the variation in the aligned TFBP',
  'HIGH_INF_POS'       => 'A flag indicating if the variant falls in a high information position of the TFBP',
  'MOTIF_SCORE_CHANGE' => 'The difference in motif score of the reference and variant sequences for the TFBP',
  'TRANSCRIPTION_FACTORS' => 'List of transcription factors which bind to the transcription factor binding profile',
  'CELL_TYPE'          => 'List of cell types and classifications for regulatory feature',
  'IND'                => 'Individual name',
  'ZYG'                => 'Zygosity of individual genotype at this locus',
  'SV'                 => 'IDs of overlapping structural variants',
  'FREQS'              => 'Frequencies of overlapping variants used in filtering',
  'AF'                 => 'Frequency of existing variant in 1000 Genomes combined population',
  'AFR_AF'             => 'Frequency of existing variant in 1000 Genomes combined African population',
  'AMR_AF'             => 'Frequency of existing variant in 1000 Genomes combined American population',
  'ASN_AF'             => 'Frequency of existing variant in 1000 Genomes combined Asian population',
  'EAS_AF'             => 'Frequency of existing variant in 1000 Genomes combined East Asian population',
  'EUR_AF'             => 'Frequency of existing variant in 1000 Genomes combined European population',
  'SAS_AF'             => 'Frequency of existing variant in 1000 Genomes combined South Asian population',
  'AA_AF'              => 'Frequency of existing variant in NHLBI-ESP African American population',
  'EA_AF'              => 'Frequency of existing variant in NHLBI-ESP European American population',
  'ExAC_AF',           => 'Frequency of existing variant in ExAC combined population',
  'ExAC_Adj_AF',       => 'Adjusted frequency of existing variant in ExAC combined population',
  'ExAC_AFR_AF',       => 'Frequency of existing variant in ExAC African/American population',
  'ExAC_AMR_AF',       => 'Frequency of existing variant in ExAC American population',
  'ExAC_EAS_AF',       => 'Frequency of existing variant in ExAC East Asian population',
  'ExAC_FIN_AF',       => 'Frequency of existing variant in ExAC Finnish population',
  'ExAC_NFE_AF',       => 'Frequency of existing variant in ExAC Non-Finnish European population',
  'ExAC_OTH_AF',       => 'Frequency of existing variant in ExAC other combined populations',
  'ExAC_SAS_AF',       => 'Frequency of existing variant in ExAC South Asian population',
  'gnomADe_AF',         => 'Frequency of existing variant in gnomAD exomes combined population',
  'gnomADe_AFR_AF',     => 'Frequency of existing variant in gnomAD exomes African/American population',
  'gnomADe_AMR_AF',     => 'Frequency of existing variant in gnomAD exomes American population',
  'gnomADe_ASJ_AF',     => 'Frequency of existing variant in gnomAD exomes Ashkenazi Jewish population',
  'gnomADe_EAS_AF',     => 'Frequency of existing variant in gnomAD exomes East Asian population',
  'gnomADe_FIN_AF',     => 'Frequency of existing variant in gnomAD exomes Finnish population',
  'gnomADe_MID_AF',     => 'Frequency of existing variant in gnomAD exomes Mid-eastern population',
  'gnomADe_NFE_AF',     => 'Frequency of existing variant in gnomAD exomes Non-Finnish European population',
  'gnomADe_REMAINING_AF',     => 'Frequency of existing variant in gnomAD exomes remaining combined populations',
  'gnomADe_SAS_AF',     => 'Frequency of existing variant in gnomAD exomes South Asian population',
  'gnomADg_AF',         => 'Frequency of existing variant in gnomAD genomes combined population',
  'gnomADg_AFR_AF',     => 'Frequency of existing variant in gnomAD genomes African/American population',
  'gnomADg_AMI_AF',     => 'Frequency of existing variant in gnomAD genomes Amish population',
  'gnomADg_AMR_AF',     => 'Frequency of existing variant in gnomAD genomes American population',
  'gnomADg_ASJ_AF',     => 'Frequency of existing variant in gnomAD genomes Ashkenazi Jewish population',
  'gnomADg_EAS_AF',     => 'Frequency of existing variant in gnomAD genomes East Asian population',
  'gnomADg_FIN_AF',     => 'Frequency of existing variant in gnomAD genomes Finnish population',
  'gnomADg_MID_AF',     => 'Frequency of existing variant in gnomAD genomes Mid-eastern population',
  'gnomADg_NFE_AF',     => 'Frequency of existing variant in gnomAD genomes Non-Finnish European population',
  'gnomADg_REMAINING_AF',     => 'Frequency of existing variant in gnomAD genomes remaining combined populations',
  'gnomADg_SAS_AF',     => 'Frequency of existing variant in gnomAD genomes South Asian population',
  'MAX_AF',            => 'Maximum observed allele frequency in 1000 Genomes, ESP and ExAC/gnomAD',
  'MAX_AF_POPS'        => 'Populations in which maximum allele frequency was observed',
  'DISTANCE'           => 'Shortest distance from variant to transcript',
  'CLIN_SIG'           => 'ClinVar clinical significance of the dbSNP variant',
  'BIOTYPE'            => 'Biotype of transcript or regulatory feature',
  'PUBMED'             => 'Pubmed ID(s) of publications that cite existing variant',
  'ALLELE_NUM'         => 'Allele number from input; 0 is reference, 1 is first alternate etc',
  'VAR_SYNONYMS'       => 'List of known variation synonyms and their sources',
  'REF_ALLELE'         => 'Reference allele',
  'STRAND'             => 'Strand of the feature (1/-1)',
  'PICK'               => 'Indicates if this consequence has been picked as the most severe',
  'SOMATIC'            => 'Somatic status of existing variant',
  'REFSEQ_MATCH'       => 'RefSeq transcript match status',
  'REFSEQ_OFFSET'      => 'HGVS adjustment length required due to mismatch between RefSeq transcript and the reference genome',
  'BAM_EDIT'           => 'Indicates success or failure of edit using BAM file',
  'GIVEN_REF'          => 'Reference allele from input',
  'USED_REF'           => 'Reference allele as used to get consequences',
  'VARIANT_CLASS'      => 'SO variant class',
  'PHENO'              => 'Indicates if existing variant(s) is associated with a phenotype, disease or trait; multiple values correspond to multiple variants',
  'GENE_PHENO'         => 'Indicates if gene is associated with a phenotype, disease or trait',
  'MINIMISED'          => 'Alleles in this variant have been converted to minimal representation before consequence calculation',
  'HGVS_OFFSET'        => 'Indicates by how many bases the HGVS notations for this variant have been shifted',
  'AMBIGUITY'          => 'Allele ambiguity code',
  'OverlapBP'          => 'Number of base pairs overlapping with the corresponding structural variation feature',
  'OverlapPC'          => 'Percentage of corresponding structural variation feature overlapped by the given input',
  'CHECK_REF'	       => 'Reports variants where the input reference does not match the expected reference',
  'UPLOADED_ALLELE'    => 'The variant allele uploaded',
  'SHIFT_LENGTH'       => 'Reports the number of bases the insertion or deletion has been shifted relative to the underlying transcript due to right alignment before consequence calculation',
  'GENCODE_PRIMARY'    => 'Reports if transcript is GENCODE primary'
);

our @DEFAULT_OUTPUT_COLS = qw(
  Uploaded_variation
  Location
  Allele
  Gene
  Feature
  Feature_type
  Consequence
  cDNA_position
  CDS_position
  Protein_position
  Amino_acids
  Codons
  Existing_variation
);

our %TS_TV = (
  'A/G' => 'Ts',
  'G/A' => 'Ts',
  'C/T' => 'Ts',
  'T/C' => 'Ts',
  'A/C' => 'Tv',
  'C/A' => 'Tv',
  'G/T' => 'Tv',
  'T/G' => 'Tv',
  'C/G' => 'Tv',
  'G/C' => 'Tv',
  'A/T' => 'Tv',
  'T/A' => 'Tv',
);

our %COLOUR_KEYS = (
  'polyphen' => {
    'unknown' => 'blue',
    'benign' => 'green',
    'possibly damaging' => 'orange',
    'probably damaging' => 'red',
  },
  'sift' => {
    'tolerated' => 'green',
    'deleterious' => 'red',
  },
  
  # copied from COLOUR.ini in web code via browser to check colours
  'consequences' => {
    'intergenic_variant'                => 'gray',
    'intron_variant'                    => '#02599c',
    'upstream_gene_variant'             => '#a2b5cd',
    'downstream_gene_variant'           => '#a2b5cd',
    '5_prime_utr_variant'               => '#7ac5cd',
    '3_prime_utr_variant'               => '#7ac5cd',
    'splice_region_variant'             => '#ff7f50',
    'splice_donor_variant'              => '#ff7f50',
    'splice_acceptor_variant'           => '#ff7f50',
    'frameshift_variant'                => '#ff69b4',
    'transcript_ablation'               => '#ff0000',
    'transcript_amplification'          => '#ff69b4',
    'inframe_insertion'                 => '#ff69b4',
    'inframe_deletion'                  => '#ff69b4',
    'synonymous_variant'                => '#76ee00',
    'stop_retained_variant'             => '#76ee00',
    'missense_variant'                  => '#ffd700',
    'initiator_codon_variant'           => '#ffd700',
    'stop_gained'                       => '#ff0000',
    'stop_lost'                         => '#ff0000',
    'mature_mirna_variant'              => '#458b00',
    'non_coding_exon_variant'           => '#32cd32',
    'nc_transcript_variant'             => '#32cd32',
    'incomplete_terminal_codon_variant' => '#ff00ff',
    'nmd_transcript_variant'            => '#ff4500',
    'coding_sequence_variant'           => '#458b00',
    'tfbs_ablation'                     => 'brown',
    'tfbs_amplification'                => 'brown',
    'tf_binding_site_variant'           => 'brown',
    'regulatory_region_variant'         => 'brown',
    'regulatory_region_ablation'        => 'brown',
    'regulatory_region_amplification'   => 'brown',
  },
);

our $MAX_NOT_ORDERED_VARIANTS = 100;
our $MAX_NOT_ORDERED_VARIANTS_DISTANCE = 100;

1;
