# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

package Bio::EnsEMBL::VEP::Constants;

use strict;
use warnings;

use base qw(Exporter);

our $VERSION = 85;

our @EXPORT_OK = qw(
  @FLAG_FIELDS
  %FIELD_DESCRIPTIONS
);

# contains an ordered map between command line flags and output columns
our @FLAG_FIELDS = (

  # general
  { flag => 'individual',      fields => ['IND','ZYG'] },
  { flag => 'allele_number',   fields => ['ALLELE_NUM'] },
  { flag => 'user',            fields => ['IMPACT','DISTANCE','STRAND','FLAGS'] },
  { flag => 'flag_pick',       fields => ['PICK'] },
  { flag => 'flag_pick_allele',fields => ['PICK'] },
  { flag => 'variant_class',   fields => ['VARIANT_CLASS']},
  { flag => 'minimal',         fields => ['MINIMISED']},

  # gene-related
  { flag => 'symbol',          fields => ['SYMBOL','SYMBOL_SOURCE','HGNC_ID'] },
  { flag => 'biotype',         fields => ['BIOTYPE'] },
  { flag => 'canonical',       fields => ['CANONICAL'] },
  { flag => 'tsl',             fields => ['TSL']},
  { flag => 'appris',          fields => ['APPRIS']},
  { flag => 'ccds',            fields => ['CCDS'] },
  { flag => 'protein',         fields => ['ENSP'] },
  { flag => 'uniprot',         fields => ['SWISSPROT', 'TREMBL', 'UNIPARC'] },
  { flag => 'xref_refseq',     fields => ['RefSeq'] },
  { flag => 'refseq',          fields => ['REFSEQ_MATCH'] },
  { flag => 'merged',          fields => ['REFSEQ_MATCH'] },
  { flag => 'gene_phenotype',  fields => ['GENE_PHENO'] },

  # non-synonymous predictions
  { flag => 'sift',            fields => ['SIFT'] },
  { flag => 'polyphen',        fields => ['PolyPhen'] },

  # transcript/protein stuff
  { flag => 'numbers',         fields => ['EXON','INTRON'] },
  { flag => 'domains',         fields => ['DOMAINS'] },
  { flag => 'hgvs',            fields => ['HGVSc','HGVSp','HGVS_OFFSET'] },

  # frequency stuff
  { flag => 'gmaf',            fields => ['GMAF'] },
  { flag => 'maf_1kg',         fields => ['AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF'] },
  { flag => 'maf_esp',         fields => ['AA_MAF','EA_MAF'] },
  { flag => 'maf_exac',        fields => ['ExAC_MAF','ExAC_Adj_MAF','ExAC_AFR_MAF','ExAC_AMR_MAF','ExAC_EAS_MAF','ExAC_FIN_MAF','ExAC_NFE_MAF','ExAC_OTH_MAF','ExAC_SAS_MAF'] },
  { flag => 'check_frequency', fields => ['FREQS'] },

  # misc variation stuff
  { flag => 'check_existing',  fields => ['CLIN_SIG','SOMATIC','PHENO'] },
  { flag => 'pubmed',          fields => ['PUBMED'] },
  { flag => 'check_svs',       fields => ['SV'] },

  # regulatory
  { flag => 'regulatory',      fields => ['MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE'] },
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
  'TSL'                => 'Transcript support level',
  'APPRIS'             => 'Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods',
  'CCDS'               => 'Indicates if transcript is a CCDS transcript',
  'SYMBOL'             => 'Gene symbol (e.g. HGNC)',
  'SYMBOL_SOURCE'      => 'Source of gene symbol',
  'SOURCE'             => 'Source of transcript in merged gene set',
  'HGNC_ID'            => 'Stable identifer of HGNC gene symbol',
  'ENSP'               => 'Protein identifer',
  'FLAGS'              => 'Transcript quality flags',
  'SWISSPROT'          => 'UniProtKB/Swiss-Prot accession',
  'TREMBL'             => 'UniProtKB/TrEMBL accession',
  'UNIPARC'            => 'UniParc accession',
  'HGVSc'              => 'HGVS coding sequence name',
  'HGVSp'              => 'HGVS protein sequence name',
  'SIFT'               => 'SIFT prediction and/or score',
  'PolyPhen'           => 'PolyPhen prediction and/or score',
  'EXON'               => 'Exon number(s) / total',
  'INTRON'             => 'Intron number(s) / total',
  'DOMAINS'            => 'The source and identifer of any overlapping protein domains',
  'MOTIF_NAME'         => 'The source and identifier of a transcription factor binding profile (TFBP) aligned at this position',
  'MOTIF_POS'          => 'The relative position of the variation in the aligned TFBP',
  'HIGH_INF_POS'       => 'A flag indicating if the variant falls in a high information position of the TFBP',
  'MOTIF_SCORE_CHANGE' => 'The difference in motif score of the reference and variant sequences for the TFBP',
  'CELL_TYPE'          => 'List of cell types and classifications for regulatory feature',
  'IND'                => 'Individual name',
  'ZYG'                => 'Zygosity of individual genotype at this locus',
  'SV'                 => 'IDs of overlapping structural variants',
  'FREQS'              => 'Frequencies of overlapping variants used in filtering',
  'GMAF'               => 'Minor allele and frequency of existing variant in 1000 Genomes combined population',
  'AFR_MAF'            => 'Frequency of existing variant in 1000 Genomes combined African population',
  'AMR_MAF'            => 'Frequency of existing variant in 1000 Genomes combined American population',
  'ASN_MAF'            => 'Frequency of existing variant in 1000 Genomes combined Asian population',
  'EAS_MAF'            => 'Frequency of existing variant in 1000 Genomes combined East Asian population',
  'EUR_MAF'            => 'Frequency of existing variant in 1000 Genomes combined European population',
  'SAS_MAF'            => 'Frequency of existing variant in 1000 Genomes combined South Asian population',
  'AA_MAF'             => 'Frequency of existing variant in NHLBI-ESP African American population',
  'EA_MAF'             => 'Frequency of existing variant in NHLBI-ESP European American population',
  'ExAC_MAF',          => 'Frequency of existing variant in ExAC combined population',
  'ExAC_Adj_MAF',      => 'Adjusted frequency of existing variant in ExAC combined population',
  'ExAC_AFR_MAF',      => 'Frequency of existing variant in ExAC African/American population',
  'ExAC_AMR_MAF',      => 'Frequency of existing variant in ExAC American population',
  'ExAC_EAS_MAF',      => 'Frequency of existing variant in ExAC East Asian population',
  'ExAC_FIN_MAF',      => 'Frequency of existing variant in ExAC Finnish population',
  'ExAC_NFE_MAF',      => 'Frequency of existing variant in ExAC Non-Finnish European population',
  'ExAC_OTH_MAF',      => 'Frequency of existing variant in ExAC combined other combined populations',
  'ExAC_SAS_MAF',      => 'Frequency of existing variant in ExAC South Asian population',
  'DISTANCE'           => 'Shortest distance from variant to transcript',
  'CLIN_SIG'           => 'ClinVar clinical significance of the dbSNP variant',
  'BIOTYPE'            => 'Biotype of transcript or regulatory feature',
  'PUBMED'             => 'Pubmed ID(s) of publications that cite existing variant',
  'ALLELE_NUM'         => 'Allele number from input; 0 is reference, 1 is first alternate etc',
  'STRAND'             => 'Strand of the feature (1/-1)',
  'PICK'               => 'Indicates if this consequence has been picked as the most severe',
  'SOMATIC'            => 'Somatic status of existing variant',
  'REFSEQ_MATCH'       => 'RefSeq transcript match status',
  'VARIANT_CLASS'      => 'SO variant class',
  'PHENO'              => 'Indicates if existing variant(s) is associated with a phenotype, disease or trait',
  'GENE_PHENO'         => 'Indicates if gene is associated with a phenotype, disease or trait',
  'MINIMISED'          => 'Alleles in this variant have been converted to minimal representation before consequence calculation',
  'HGVS_OFFSET'        => 'Indicates by how many bases the HGVS notations for this variant have been shifted',
);

1;
