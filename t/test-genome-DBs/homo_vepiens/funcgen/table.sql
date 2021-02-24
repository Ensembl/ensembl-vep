CREATE TABLE `alignment` (
  `alignment_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `name` varchar(255) NOT NULL,
  `bam_file_id` int(11) DEFAULT NULL,
  `bigwig_file_id` int(11) DEFAULT NULL,
  `experiment_id` int(15) unsigned DEFAULT NULL,
  `has_duplicates` tinyint(1) DEFAULT NULL,
  `is_control` tinyint(1) DEFAULT NULL,
  `source_alignment_id` int(22) unsigned DEFAULT NULL,
  `deduplicated_alignment_id` int(28) unsigned DEFAULT NULL,
  `to_gender` enum('male','female','hermaphrodite','mixed','unknown') DEFAULT NULL,
  `is_complete` tinyint(1) DEFAULT NULL,
  PRIMARY KEY (`alignment_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `alignment_qc_flagstats` (
  `alignment_qc_flagstats_id` int(28) unsigned NOT NULL AUTO_INCREMENT,
  `alignment_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  `category` varchar(100) NOT NULL,
  `qc_passed_reads` int(10) unsigned DEFAULT NULL,
  `qc_failed_reads` int(10) unsigned DEFAULT NULL,
  `path` varchar(512) NOT NULL,
  `bam_file` varchar(512) NOT NULL,
  PRIMARY KEY (`alignment_qc_flagstats_id`),
  UNIQUE KEY `name_exp_idx` (`alignment_qc_flagstats_id`,`category`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `alignment_read_file` (
  `alignment_read_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `alignment_id` int(10) unsigned NOT NULL,
  `read_file_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`alignment_read_file_id`,`alignment_id`),
  UNIQUE KEY `rset_table_idname_idx` (`alignment_id`,`read_file_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `analysis` (
  `analysis_id` smallint(5) unsigned NOT NULL AUTO_INCREMENT,
  `created` datetime DEFAULT NULL,
  `logic_name` varchar(100) NOT NULL,
  `db` varchar(120) DEFAULT NULL,
  `db_version` varchar(40) DEFAULT NULL,
  `db_file` varchar(120) DEFAULT NULL,
  `program` varchar(80) DEFAULT NULL,
  `program_version` varchar(40) DEFAULT NULL,
  `program_file` varchar(80) DEFAULT NULL,
  `parameters` text,
  `module` varchar(80) DEFAULT NULL,
  `module_version` varchar(40) DEFAULT NULL,
  `gff_source` varchar(40) DEFAULT NULL,
  `gff_feature` varchar(40) DEFAULT NULL,
  PRIMARY KEY (`analysis_id`),
  UNIQUE KEY `logic_name_idx` (`logic_name`)
) ENGINE=MyISAM AUTO_INCREMENT=81 DEFAULT CHARSET=latin1;

CREATE TABLE `analysis_description` (
  `analysis_id` smallint(5) unsigned NOT NULL,
  `description` text,
  `display_label` varchar(255) NOT NULL,
  `displayable` tinyint(1) NOT NULL DEFAULT '1',
  `web_data` text,
  UNIQUE KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `array` (
  `array_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(40) DEFAULT NULL,
  `format` varchar(20) DEFAULT NULL,
  `vendor` varchar(40) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `type` varchar(20) DEFAULT NULL,
  `class` varchar(20) DEFAULT NULL,
  `is_probeset_array` tinyint(1) NOT NULL DEFAULT '0',
  `is_linked_array` tinyint(1) NOT NULL DEFAULT '0',
  `has_sense_interrogation` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`array_id`),
  UNIQUE KEY `vendor_name_idx` (`vendor`,`name`),
  UNIQUE KEY `class_name_idx` (`class`,`name`)
) ENGINE=MyISAM AUTO_INCREMENT=63 DEFAULT CHARSET=latin1;

CREATE TABLE `array_chip` (
  `array_chip_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `design_id` varchar(100) DEFAULT NULL,
  `array_id` int(10) unsigned NOT NULL,
  `name` varchar(100) DEFAULT NULL,
  PRIMARY KEY (`array_chip_id`),
  UNIQUE KEY `array_design_idx` (`array_id`,`design_id`)
) ENGINE=MyISAM AUTO_INCREMENT=103 DEFAULT CHARSET=latin1;

CREATE TABLE `associated_feature_type` (
  `table_id` int(10) unsigned NOT NULL,
  `table_name` enum('annotated_feature','external_feature','regulatory_feature','feature_type') NOT NULL,
  `feature_type_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`table_id`,`table_name`,`feature_type_id`),
  KEY `feature_type_index` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `associated_group` (
  `associated_group_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `description` varchar(128) DEFAULT NULL,
  PRIMARY KEY (`associated_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `associated_motif_feature` (
  `annotated_feature_id` int(10) unsigned NOT NULL,
  `motif_feature_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`annotated_feature_id`,`motif_feature_id`),
  KEY `motif_feature_idx` (`motif_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `associated_xref` (
  `associated_xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `object_xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `source_xref_id` int(10) unsigned DEFAULT NULL,
  `condition_type` varchar(128) DEFAULT NULL,
  `associated_group_id` int(10) unsigned DEFAULT NULL,
  `rank` int(10) unsigned DEFAULT '0',
  PRIMARY KEY (`associated_xref_id`),
  UNIQUE KEY `object_associated_source_type_idx` (`object_xref_id`,`xref_id`,`source_xref_id`,`condition_type`,`associated_group_id`),
  KEY `associated_source_idx` (`source_xref_id`),
  KEY `associated_object_idx` (`object_xref_id`),
  KEY `associated_idx` (`xref_id`),
  KEY `associated_group_idx` (`associated_group_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `binding_matrix` (
  `binding_matrix_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(200) NOT NULL,
  `threshold` double DEFAULT NULL,
  `source` varchar(20) NOT NULL,
  `stable_id` varchar(128) NOT NULL,
  PRIMARY KEY (`binding_matrix_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=403 DEFAULT CHARSET=latin1;

CREATE TABLE `binding_matrix_frequencies` (
  `binding_matrix_frequencies_id` int(11) NOT NULL AUTO_INCREMENT,
  `binding_matrix_id` int(11) unsigned NOT NULL,
  `position` int(11) unsigned NOT NULL,
  `nucleotide` enum('A','C','G','T') NOT NULL,
  `frequency` int(10) unsigned NOT NULL,
  PRIMARY KEY (`binding_matrix_frequencies_id`),
  UNIQUE KEY `unique_constraint_idx` (`binding_matrix_id`,`position`,`nucleotide`),
  KEY `binding_matrix_id_idx` (`binding_matrix_id`)
) ENGINE=MyISAM AUTO_INCREMENT=26689 DEFAULT CHARSET=latin1;

CREATE TABLE `binding_matrix_transcription_factor_complex` (
  `binding_matrix_transcription_factor_complex_id` int(11) NOT NULL AUTO_INCREMENT,
  `binding_matrix_id` int(11) unsigned NOT NULL,
  `transcription_factor_complex_id` int(11) NOT NULL,
  PRIMARY KEY (`binding_matrix_transcription_factor_complex_id`),
  UNIQUE KEY `binding_matrix_id_transcription_factor_complex_id_idx` (`binding_matrix_id`,`transcription_factor_complex_id`),
  KEY `binding_matrix_id_idx` (`binding_matrix_id`),
  KEY `transcription_factor_complex_id_idx` (`transcription_factor_complex_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `chance` (
  `chance_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `signal_alignment_id` int(10) unsigned DEFAULT NULL,
  `control_alignment_id` int(10) unsigned DEFAULT NULL,
  `analysis_id` smallint(10) unsigned DEFAULT NULL,
  `p` double DEFAULT NULL,
  `q` double DEFAULT NULL,
  `divergence` double DEFAULT NULL,
  `z_score` double DEFAULT NULL,
  `percent_genome_enriched` double DEFAULT NULL,
  `input_scaling_factor` double DEFAULT NULL,
  `differential_percentage_enrichment` double DEFAULT NULL,
  `control_enrichment_stronger_than_chip_at_bin` double DEFAULT NULL,
  `first_nonzero_bin_at` double DEFAULT NULL,
  `pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome` double DEFAULT NULL,
  `path` varchar(512) DEFAULT NULL,
  `run_failed` tinyint(1) DEFAULT '0',
  `error_message` text,
  PRIMARY KEY (`chance_id`),
  UNIQUE KEY `signal_control_alignment_unique` (`signal_alignment_id`,`control_alignment_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `data_file` (
  `data_file_id` int(11) NOT NULL AUTO_INCREMENT,
  `table_id` int(10) unsigned NOT NULL,
  `table_name` varchar(32) NOT NULL,
  `path` varchar(255) NOT NULL,
  `file_type` enum('BAM','BAMCOV','BIGBED','BIGWIG','VCF','CRAM','DIR') NOT NULL DEFAULT 'BAM',
  `md5sum` varchar(45) DEFAULT NULL,
  PRIMARY KEY (`data_file_id`),
  UNIQUE KEY `table_id_name_path_idx` (`table_id`,`table_name`,`path`),
  UNIQUE KEY `data_file_id` (`data_file_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `epigenome` (
  `epigenome_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(120) NOT NULL,
  `short_name` varchar(120) NOT NULL,
  `description` mediumtext,
  `production_name` varchar(120) DEFAULT NULL,
  `gender` enum('male','female','hermaphrodite','mixed','unknown') DEFAULT 'unknown',
  `search_terms` mediumtext,
  `full_name` mediumtext,
  PRIMARY KEY (`epigenome_id`),
  UNIQUE KEY `name_idx` (`name`),
  UNIQUE KEY `short_name_idx` (`short_name`)
) ENGINE=MyISAM AUTO_INCREMENT=273 DEFAULT CHARSET=latin1;

CREATE TABLE `execution_plan` (
  `execution_plan_id` int(18) unsigned NOT NULL AUTO_INCREMENT,
  `time` bigint(20) DEFAULT NULL,
  `experiment_id` int(16) unsigned NOT NULL,
  `execution_plan` longtext NOT NULL,
  PRIMARY KEY (`execution_plan_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `experiment` (
  `experiment_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) DEFAULT NULL,
  `experimental_group_id` smallint(6) unsigned DEFAULT NULL,
  `control_id` int(10) unsigned DEFAULT NULL,
  `is_control` tinyint(3) unsigned DEFAULT '0',
  `feature_type_id` int(10) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  `archive_id` varchar(60) DEFAULT NULL,
  PRIMARY KEY (`experiment_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `experimental_group_idx` (`experimental_group_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `epigenome_idx` (`epigenome_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `experimental_group` (
  `experimental_group_id` smallint(6) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(40) NOT NULL,
  `production_name` varchar(255) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  `url` varchar(255) DEFAULT NULL,
  `is_project` tinyint(1) DEFAULT '0',
  PRIMARY KEY (`experimental_group_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM AUTO_INCREMENT=6 DEFAULT CHARSET=latin1;

CREATE TABLE `external_db` (
  `external_db_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `db_name` varchar(100) NOT NULL,
  `db_release` varchar(255) DEFAULT NULL,
  `status` enum('KNOWNXREF','KNOWN','XREF','PRED','ORTH','PSEUDO') NOT NULL,
  `dbprimary_acc_linkable` tinyint(1) NOT NULL DEFAULT '1',
  `priority` int(11) NOT NULL,
  `db_display_name` varchar(255) DEFAULT NULL,
  `type` enum('ARRAY','ALT_TRANS','ALT_GENE','MISC','LIT','PRIMARY_DB_SYNONYM','ENSEMBL') DEFAULT NULL,
  `secondary_db_name` varchar(255) DEFAULT NULL,
  `secondary_db_table` varchar(255) DEFAULT NULL,
  `description` text,
  PRIMARY KEY (`external_db_id`),
  UNIQUE KEY `db_name_release_idx` (`db_name`,`db_release`)
) ENGINE=MyISAM AUTO_INCREMENT=69 DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=80;

CREATE TABLE `external_feature` (
  `external_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `display_label` varchar(60) DEFAULT NULL,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `feature_set_id` int(10) unsigned NOT NULL,
  `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL,
  PRIMARY KEY (`external_feature_id`),
  UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `feature_set_idx` (`feature_set_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=80;

CREATE TABLE `external_feature_file` (
  `external_feature_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(100) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`external_feature_file_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `epigenome_idx` (`epigenome_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM AUTO_INCREMENT=48 DEFAULT CHARSET=latin1;

CREATE TABLE `external_synonym` (
  `xref_id` int(10) unsigned NOT NULL,
  `synonym` varchar(100) NOT NULL,
  PRIMARY KEY (`xref_id`,`synonym`),
  KEY `name_index` (`synonym`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=20;

CREATE TABLE `fastqc` (
  `fastqc_id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `read_file_id` int(14) unsigned NOT NULL,
  `basic_statistics` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `per_base_sequence_quality` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `per_tile_sequence_quality` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `per_sequence_quality_scores` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `per_base_sequence_content` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `per_sequence_gc_content` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `per_base_n_content` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `sequence_length_distribution` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `sequence_duplication_levels` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `overrepresented_sequences` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `adapter_content` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `kmer_content` enum('PASS','WARN','FAIL') DEFAULT NULL,
  `run_failed` tinyint(1) DEFAULT '0',
  `error_message` text,
  PRIMARY KEY (`fastqc_id`),
  UNIQUE KEY `read_file_id_unique` (`read_file_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `feature_set` (
  `feature_set_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `feature_type_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `name` varchar(100) DEFAULT NULL,
  `type` enum('annotated','regulatory','external','segmentation','mirna_target') DEFAULT NULL,
  `description` varchar(80) DEFAULT NULL,
  `display_label` varchar(80) DEFAULT NULL,
  PRIMARY KEY (`feature_set_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `feature_type_idx` (`feature_type_id`)
) ENGINE=MyISAM AUTO_INCREMENT=609 DEFAULT CHARSET=latin1;

CREATE TABLE `feature_type` (
  `feature_type_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(40) NOT NULL,
  `class` enum('Insulator','DNA','Regulatory Feature','Histone','RNA','Polymerase','Transcription Factor','Transcription Factor Complex','Regulatory Motif','Enhancer','Expression','Pseudo','Open Chromatin','Search Region','Association Locus','Segmentation State','DNA Modification','Transcription Start Site') DEFAULT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `so_accession` varchar(64) DEFAULT NULL,
  `so_term` varchar(255) DEFAULT NULL,
  `production_name` varchar(120) DEFAULT NULL,
  PRIMARY KEY (`feature_type_id`),
  UNIQUE KEY `name_class_analysis_idx` (`name`,`class`,`analysis_id`),
  KEY `so_accession_idx` (`so_accession`)
) ENGINE=MyISAM AUTO_INCREMENT=180054 DEFAULT CHARSET=latin1;

CREATE TABLE `frip` (
  `frip_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `peak_calling_id` int(18) unsigned NOT NULL,
  `frip` double DEFAULT NULL,
  `total_reads` int(14) DEFAULT NULL,
  PRIMARY KEY (`frip_id`),
  UNIQUE KEY `peak_calling_id_unique` (`peak_calling_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `identity_xref` (
  `object_xref_id` int(10) unsigned NOT NULL,
  `xref_identity` int(5) DEFAULT NULL,
  `ensembl_identity` int(5) DEFAULT NULL,
  `xref_start` int(11) DEFAULT NULL,
  `xref_end` int(11) DEFAULT NULL,
  `ensembl_start` int(11) DEFAULT NULL,
  `ensembl_end` int(11) DEFAULT NULL,
  `cigar_line` text,
  `score` double DEFAULT NULL,
  `evalue` double DEFAULT NULL,
  PRIMARY KEY (`object_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `idr` (
  `idr_id` int(9) unsigned NOT NULL AUTO_INCREMENT,
  `experiment_id` int(15) unsigned NOT NULL,
  `max_peaks` int(11) unsigned DEFAULT NULL,
  `type` enum('on biological replicates','on technical replicates','no_idr') NOT NULL,
  `failed_idr_pairs` text,
  PRIMARY KEY (`idr_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `meta` (
  `meta_id` int(10) NOT NULL AUTO_INCREMENT,
  `species_id` int(10) unsigned DEFAULT '1',
  `meta_key` varchar(46) NOT NULL,
  `meta_value` varchar(950) NOT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `species_key_value_idx` (`species_id`,`meta_key`,`meta_value`),
  KEY `species_value_idx` (`species_id`,`meta_value`)
) ENGINE=MyISAM AUTO_INCREMENT=777 DEFAULT CHARSET=latin1;

CREATE TABLE `meta_coord` (
  `table_name` varchar(40) NOT NULL,
  `coord_system_id` int(10) unsigned NOT NULL,
  `max_length` int(11) DEFAULT NULL,
  UNIQUE KEY `table_name` (`table_name`,`coord_system_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `mirna_target_feature` (
  `mirna_target_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `accession` varchar(60) DEFAULT NULL,
  `display_label` varchar(60) DEFAULT NULL,
  `evidence` varchar(60) DEFAULT NULL,
  `method` varchar(60) DEFAULT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `supporting_information` varchar(100) DEFAULT NULL,
  `analysis_id` smallint(10) unsigned DEFAULT NULL,
  `gene_stable_id` varchar(128) DEFAULT NULL,
  PRIMARY KEY (`mirna_target_feature_id`),
  UNIQUE KEY `unique_idx` (`accession`,`gene_stable_id`,`seq_region_start`,`seq_region_end`,`evidence`,`method`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000;

CREATE TABLE `motif_feature` (
  `motif_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `binding_matrix_id` int(10) unsigned NOT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `score` double DEFAULT NULL,
  `stable_id` varchar(18) DEFAULT NULL,
  PRIMARY KEY (`motif_feature_id`),
  UNIQUE KEY `unique_idx` (`binding_matrix_id`,`seq_region_id`,`seq_region_start`,`seq_region_strand`),
  UNIQUE KEY `stable_id_idx` (`stable_id`),
  KEY `seq_region_idx` (`seq_region_id`,`seq_region_start`),
  KEY `binding_matrix_idx` (`binding_matrix_id`)
) ENGINE=MyISAM AUTO_INCREMENT=73368 DEFAULT CHARSET=latin1;

CREATE TABLE `motif_feature_peak` (
  `motif_feature_peak_id` int(11) NOT NULL AUTO_INCREMENT,
  `motif_feature_id` int(11) unsigned NOT NULL,
  `peak_id` int(11) unsigned NOT NULL,
  PRIMARY KEY (`motif_feature_peak_id`),
  UNIQUE KEY `motif_feature_idx` (`motif_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `motif_feature_regulatory_feature` (
  `motif_feature_regulatory_feature_id` int(11) NOT NULL AUTO_INCREMENT,
  `motif_feature_id` int(11) unsigned NOT NULL,
  `regulatory_feature_id` int(11) unsigned NOT NULL,
  `epigenome_id` int(11) unsigned DEFAULT NULL,
  `has_matching_Peak` tinyint(3) unsigned DEFAULT '0',
  PRIMARY KEY (`motif_feature_regulatory_feature_id`),
  UNIQUE KEY `mf_rf_ep_idx` (`motif_feature_id`,`regulatory_feature_id`,`epigenome_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `object_xref` (
  `object_xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ensembl_id` int(10) unsigned NOT NULL,
  `ensembl_object_type` enum('Epigenome','Experiment','RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType','MirnaTargetFeature','ProbeSet','Probe','ProbeFeature','ReadFile') NOT NULL,
  `xref_id` int(10) unsigned NOT NULL,
  `linkage_annotation` varchar(255) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  PRIMARY KEY (`object_xref_id`),
  UNIQUE KEY `xref_idx` (`xref_id`,`ensembl_object_type`,`ensembl_id`,`analysis_id`),
  KEY `analysis_idx` (`analysis_id`),
  KEY `ensembl_idx` (`ensembl_object_type`,`ensembl_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=40;

CREATE TABLE `ontology_xref` (
  `object_xref_id` int(10) unsigned NOT NULL DEFAULT '0',
  `source_xref_id` int(10) unsigned DEFAULT NULL,
  `linkage_type` enum('IC','IDA','IEA','IEP','IGI','IMP','IPI','ISS','NAS','ND','TAS','NR','RCA') NOT NULL,
  UNIQUE KEY `object_xref_id_2` (`object_xref_id`,`source_xref_id`,`linkage_type`),
  KEY `object_xref_id` (`object_xref_id`),
  KEY `source_xref_id` (`source_xref_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `peak` (
  `peak_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `score` double DEFAULT NULL,
  `peak_calling_id` int(10) unsigned NOT NULL,
  `summit` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`peak_id`),
  UNIQUE KEY `seq_region_feature_set_idx` (`seq_region_id`,`seq_region_start`,`peak_calling_id`),
  KEY `feature_set_idx` (`peak_calling_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=100000000 AVG_ROW_LENGTH=39;

CREATE TABLE `peak_calling` (
  `peak_calling_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(300) NOT NULL,
  `display_label` varchar(300) NOT NULL,
  `feature_type_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `signal_alignment_id` int(23) unsigned DEFAULT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  `experiment_id` int(10) unsigned DEFAULT NULL,
  `run_failed` tinyint(1) DEFAULT '0',
  `error_message` text,
  `control_alignment_id` int(23) unsigned DEFAULT NULL,
  `used_for_regulatory_build` tinyint(1) DEFAULT '1',
  PRIMARY KEY (`peak_calling_id`),
  UNIQUE KEY `peak_calling_id_idx` (`peak_calling_id`),
  UNIQUE KEY `peak_calling_name_unique` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `peak_calling_statistic` (
  `peak_calling_statistic_id` int(28) unsigned NOT NULL AUTO_INCREMENT,
  `peak_calling_id` int(18) unsigned DEFAULT NULL,
  `epigenome_id` int(15) unsigned DEFAULT NULL,
  `feature_type_id` int(18) unsigned DEFAULT NULL,
  `statistic` varchar(255) NOT NULL,
  `value` float unsigned DEFAULT NULL,
  PRIMARY KEY (`peak_calling_statistic_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `phantom_peak` (
  `phantom_peak_id` int(17) unsigned NOT NULL AUTO_INCREMENT,
  `analysis_id` smallint(5) unsigned DEFAULT NULL,
  `alignment_id` int(15) unsigned NOT NULL,
  `num_reads` int(12) unsigned DEFAULT NULL,
  `est_frag_len` double DEFAULT NULL,
  `est_frag_len_2` double DEFAULT NULL,
  `est_frag_len_3` double DEFAULT NULL,
  `corr_est_frag_len` double DEFAULT NULL,
  `corr_est_frag_len_2` double DEFAULT NULL,
  `corr_est_frag_len_3` double DEFAULT NULL,
  `phantom_peak` int(17) unsigned DEFAULT NULL,
  `corr_phantom_peak` double DEFAULT NULL,
  `argmin_corr` int(14) DEFAULT NULL,
  `min_corr` double DEFAULT NULL,
  `nsc` double DEFAULT NULL,
  `rsc` double DEFAULT NULL,
  `quality_tag` int(14) DEFAULT NULL,
  `run_failed` tinyint(1) DEFAULT '0',
  `error_message` text,
  PRIMARY KEY (`phantom_peak_id`),
  UNIQUE KEY `alignment_id_unique` (`alignment_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `probe` (
  `probe_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_set_id` int(10) unsigned DEFAULT NULL,
  `name` varchar(100) NOT NULL,
  `length` smallint(6) unsigned NOT NULL,
  `array_chip_id` int(10) unsigned NOT NULL,
  `class` varchar(20) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  `probe_seq_id` int(10) DEFAULT NULL,
  PRIMARY KEY (`probe_id`,`name`,`array_chip_id`),
  UNIQUE KEY `probe_idx` (`probe_id`),
  KEY `probe_set_idx` (`probe_set_id`),
  KEY `array_chip_idx` (`array_chip_id`),
  KEY `name_idx` (`name`),
  KEY `probe_seq_idx` (`probe_seq_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `probe_feature` (
  `probe_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_start` int(10) NOT NULL,
  `seq_region_end` int(10) NOT NULL,
  `seq_region_strand` tinyint(4) NOT NULL,
  `probe_id` int(10) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `mismatches` tinyint(4) NOT NULL,
  `cigar_line` varchar(50) DEFAULT NULL,
  `hit_id` varchar(255) DEFAULT NULL,
  `source` enum('genomic','transcript') DEFAULT NULL,
  PRIMARY KEY (`probe_feature_id`),
  KEY `probe_idx` (`probe_id`),
  KEY `seq_region_probe_probe_feature_idx` (`seq_region_id`,`seq_region_start`,`seq_region_end`,`probe_id`,`probe_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `probe_feature_transcript` (
  `probe_feature_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_feature_id` int(10) unsigned DEFAULT NULL,
  `stable_id` varchar(128) DEFAULT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_feature_transcript_id`),
  KEY `probe_feature_transcript_id_idx` (`probe_feature_transcript_id`),
  KEY `probe_feature_id_idx` (`probe_feature_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `probe_mapping` (
  `probe_mapping_id` int(22) unsigned NOT NULL AUTO_INCREMENT,
  `assembly` varchar(255) DEFAULT NULL,
  `gene_build_version` varchar(255) DEFAULT NULL,
  `five_prime_utr` int(22) unsigned DEFAULT NULL,
  `three_prime_utr` int(22) unsigned DEFAULT NULL,
  `sample_probe_id` int(22) unsigned DEFAULT NULL,
  `sample_probe_set_id` int(22) unsigned DEFAULT NULL,
  `release_version` varchar(255) DEFAULT NULL,
  `release_date` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_mapping_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `probe_mapping_statistic` (
  `probe_mapping_statistic_id` int(29) unsigned NOT NULL AUTO_INCREMENT,
  `array_id` int(11) unsigned DEFAULT NULL,
  `statistic` varchar(255) NOT NULL,
  `value` double unsigned DEFAULT NULL,
  PRIMARY KEY (`probe_mapping_statistic_id`)
) ENGINE=InnoDB AUTO_INCREMENT=34 DEFAULT CHARSET=latin1;

CREATE TABLE `probe_seq` (
  `probe_seq_id` int(10) NOT NULL AUTO_INCREMENT,
  `sequence` text NOT NULL,
  `sequence_upper` text NOT NULL,
  `sequence_upper_sha1` char(40) NOT NULL,
  PRIMARY KEY (`probe_seq_id`),
  UNIQUE KEY `sequence_upper_sha1` (`sequence_upper_sha1`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `probe_set` (
  `probe_set_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(100) NOT NULL,
  `size` smallint(6) unsigned NOT NULL,
  `family` varchar(20) DEFAULT NULL,
  `array_chip_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`probe_set_id`),
  KEY `name` (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `probe_set_transcript` (
  `probe_set_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_set_id` int(10) unsigned NOT NULL,
  `stable_id` varchar(128) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_set_transcript_id`),
  KEY `probe_set_transcript_id_idx` (`probe_set_transcript_id`),
  KEY `probe_set_transcript_stable_id_idx` (`stable_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `probe_transcript` (
  `probe_transcript_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `probe_id` int(10) unsigned NOT NULL,
  `stable_id` varchar(128) NOT NULL,
  `description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`probe_transcript_id`),
  KEY `probe_transcript_id` (`probe_transcript_id`),
  KEY `probe_transcript_stable_id_idx` (`stable_id`),
  KEY `probe_transcript_idx` (`probe_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `read_file` (
  `read_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(300) NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `is_paired_end` tinyint(1) DEFAULT NULL,
  `file_size` bigint(20) DEFAULT NULL,
  `number_of_reads` bigint(20) DEFAULT NULL,
  `read_length` int(10) DEFAULT NULL,
  `md5sum` varchar(45) DEFAULT NULL,
  `file` text,
  `notes` text,
  PRIMARY KEY (`read_file_id`),
  UNIQUE KEY `read_file_id_idx` (`read_file_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `read_file_experimental_configuration` (
  `read_file_experimental_configuration_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `read_file_id` int(10) unsigned DEFAULT NULL,
  `experiment_id` int(10) unsigned NOT NULL,
  `biological_replicate` tinyint(3) unsigned NOT NULL DEFAULT '1',
  `technical_replicate` tinyint(3) unsigned NOT NULL DEFAULT '1',
  `paired_end_tag` int(11) DEFAULT NULL,
  `multiple` int(11) DEFAULT '1',
  PRIMARY KEY (`read_file_experimental_configuration_id`),
  UNIQUE KEY `name_exp_idx` (`experiment_id`,`biological_replicate`,`technical_replicate`,`paired_end_tag`,`multiple`),
  KEY `experiment_idx` (`experiment_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `regulatory_activity` (
  `regulatory_activity_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `regulatory_feature_id` int(10) unsigned DEFAULT NULL,
  `activity` enum('INACTIVE','REPRESSED','POISED','ACTIVE','NA') NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`regulatory_activity_id`),
  UNIQUE KEY `uniqueness_constraint_idx` (`epigenome_id`,`regulatory_feature_id`),
  KEY `regulatory_feature_idx` (`regulatory_feature_id`)
) ENGINE=MyISAM AUTO_INCREMENT=7800 DEFAULT CHARSET=latin1;

CREATE TABLE `regulatory_build` (
  `regulatory_build_id` int(4) unsigned NOT NULL AUTO_INCREMENT,
  `name` text,
  `release_version` int(11) DEFAULT NULL,
  `description` text,
  `version` varchar(50) DEFAULT NULL,
  `initial_release_date` varchar(50) DEFAULT NULL,
  `last_annotation_update` varchar(50) DEFAULT NULL,
  `feature_type_id` int(4) unsigned NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `is_current` tinyint(1) NOT NULL DEFAULT '0',
  `sample_regulatory_feature_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`regulatory_build_id`)
) ENGINE=MyISAM AUTO_INCREMENT=2 DEFAULT CHARSET=latin1;

CREATE TABLE `regulatory_build_epigenome` (
  `regulatory_build_epigenome_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `regulatory_build_id` int(10) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`regulatory_build_epigenome_id`)
) ENGINE=MyISAM AUTO_INCREMENT=39 DEFAULT CHARSET=latin1;

CREATE TABLE `regulatory_build_statistic` (
  `regulatory_build_statistic_id` int(30) unsigned NOT NULL AUTO_INCREMENT,
  `regulatory_build_id` int(22) unsigned DEFAULT NULL,
  `statistic` varchar(255) DEFAULT NULL,
  `value` float unsigned DEFAULT NULL,
  PRIMARY KEY (`regulatory_build_statistic_id`),
  UNIQUE KEY `stats_uniq` (`statistic`,`regulatory_build_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `regulatory_evidence` (
  `regulatory_activity_id` int(10) unsigned NOT NULL,
  `attribute_feature_id` int(10) unsigned NOT NULL,
  `attribute_feature_table` enum('annotated','motif') NOT NULL DEFAULT 'annotated',
  PRIMARY KEY (`regulatory_activity_id`,`attribute_feature_table`,`attribute_feature_id`),
  KEY `attribute_feature_idx` (`attribute_feature_id`,`attribute_feature_table`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `regulatory_feature` (
  `regulatory_feature_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `seq_region_id` int(10) unsigned NOT NULL,
  `seq_region_strand` tinyint(1) NOT NULL,
  `seq_region_start` int(10) unsigned NOT NULL,
  `seq_region_end` int(10) unsigned NOT NULL,
  `stable_id` varchar(18) DEFAULT NULL,
  `bound_start_length` mediumint(3) unsigned NOT NULL,
  `bound_end_length` mediumint(3) unsigned NOT NULL,
  `epigenome_count` smallint(6) DEFAULT NULL,
  `regulatory_build_id` int(10) unsigned DEFAULT NULL,
  PRIMARY KEY (`regulatory_feature_id`),
  UNIQUE KEY `uniqueness_constraint_idx` (`feature_type_id`,`seq_region_id`,`seq_region_strand`,`seq_region_start`,`seq_region_end`,`stable_id`,`bound_start_length`,`bound_end_length`,`regulatory_build_id`),
  KEY `feature_type_idx` (`feature_type_id`),
  KEY `stable_id_idx` (`stable_id`)
) ENGINE=MyISAM AUTO_INCREMENT=201 DEFAULT CHARSET=latin1;

CREATE TABLE `segmentation` (
  `segmentation_id` int(18) unsigned NOT NULL AUTO_INCREMENT,
  `regulatory_build_id` int(22) unsigned DEFAULT NULL,
  `name` varchar(255) NOT NULL,
  `superclass` varchar(255) NOT NULL,
  `class` varchar(255) NOT NULL,
  PRIMARY KEY (`segmentation_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `segmentation_cell_tables` (
  `superclass` varchar(255) NOT NULL,
  `class` varchar(255) NOT NULL,
  `segmentation_id` int(18) unsigned NOT NULL,
  `epigenome_id` int(16) unsigned NOT NULL,
  `feature_type_id` int(18) unsigned NOT NULL,
  `signal_alignment_id` int(23) unsigned NOT NULL,
  `control_alignment_id` int(23) unsigned DEFAULT NULL
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `segmentation_file` (
  `segmentation_file_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `regulatory_build_id` int(4) unsigned DEFAULT NULL,
  `name` varchar(100) DEFAULT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `epigenome_id` int(10) unsigned DEFAULT NULL,
  `segmentation_id` int(18) unsigned DEFAULT NULL,
  PRIMARY KEY (`segmentation_file_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `epigenome_idx` (`epigenome_id`),
  KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM AUTO_INCREMENT=39 DEFAULT CHARSET=latin1;

CREATE TABLE `segmentation_state_assignment` (
  `segmentation_state_assignment_id` int(35) unsigned NOT NULL AUTO_INCREMENT,
  `state` int(8) NOT NULL,
  `segmentation` varchar(255) NOT NULL,
  `assignment` varchar(255) NOT NULL,
  PRIMARY KEY (`segmentation_state_assignment_id`),
  UNIQUE KEY `state` (`state`,`segmentation`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `segmentation_state_emission` (
  `segmentation_state_emission_id` int(27) unsigned NOT NULL AUTO_INCREMENT,
  `segmentation` varchar(255) NOT NULL,
  `state` int(7) DEFAULT NULL,
  `CTCF` double DEFAULT NULL,
  `DNase1` double DEFAULT NULL,
  `H3K27ac` double DEFAULT NULL,
  `H3K27me3` double DEFAULT NULL,
  `H3K36me3` double DEFAULT NULL,
  `H3K4me1` double DEFAULT NULL,
  `H3K4me2` double DEFAULT NULL,
  `H3K4me3` double DEFAULT NULL,
  `H3K9ac` double DEFAULT NULL,
  `H3K9me3` double DEFAULT NULL,
  PRIMARY KEY (`segmentation_state_emission_id`),
  UNIQUE KEY `state` (`state`,`segmentation`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `segmentation_statistic` (
  `segmentation_statistic_id` int(30) unsigned NOT NULL AUTO_INCREMENT,
  `segmentation_id` int(18) unsigned DEFAULT NULL,
  `state` int(8) unsigned DEFAULT NULL,
  `epigenome_id` int(22) unsigned DEFAULT NULL,
  `label` varchar(255) DEFAULT NULL,
  `statistic` varchar(255) NOT NULL,
  `value` float unsigned DEFAULT NULL,
  PRIMARY KEY (`segmentation_statistic_id`),
  UNIQUE KEY `stats_uniq` (`statistic`,`segmentation_id`,`epigenome_id`,`label`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

CREATE TABLE `transcription_factor` (
  `transcription_factor_id` int(11) NOT NULL AUTO_INCREMENT,
  `name` varchar(120) NOT NULL,
  `feature_type_id` int(10) unsigned DEFAULT NULL,
  `gene_stable_id` varchar(128) DEFAULT NULL,
  PRIMARY KEY (`transcription_factor_id`),
  UNIQUE KEY `name_idx` (`name`),
  KEY `feature_type_id_idx` (`feature_type_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `transcription_factor_complex` (
  `transcription_factor_complex_id` int(11) NOT NULL AUTO_INCREMENT,
  `production_name` varchar(120) NOT NULL,
  `display_name` varchar(120) NOT NULL,
  PRIMARY KEY (`transcription_factor_complex_id`),
  UNIQUE KEY `production_name_idx` (`production_name`),
  UNIQUE KEY `display_name_idx` (`display_name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `transcription_factor_complex_composition` (
  `transcription_factor_complex_composition_id` int(11) NOT NULL AUTO_INCREMENT,
  `transcription_factor_complex_id` int(11) NOT NULL,
  `transcription_factor_id` int(11) NOT NULL,
  PRIMARY KEY (`transcription_factor_complex_composition_id`),
  UNIQUE KEY `tfc_id_tf_id_idx` (`transcription_factor_complex_id`,`transcription_factor_id`),
  KEY `transcription_factor_complex_id_idx` (`transcription_factor_complex_id`),
  KEY `transcription_factor_id_idx` (`transcription_factor_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `unmapped_object` (
  `unmapped_object_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `type` enum('xref','probe2transcript','array_mapping') NOT NULL,
  `analysis_id` smallint(5) unsigned NOT NULL,
  `external_db_id` int(10) unsigned DEFAULT NULL,
  `identifier` varchar(255) NOT NULL,
  `unmapped_reason_id` int(10) unsigned NOT NULL,
  `query_score` double DEFAULT NULL,
  `target_score` double DEFAULT NULL,
  `ensembl_id` int(10) unsigned DEFAULT '0',
  `ensembl_object_type` enum('RegulatoryFeature','ExternalFeature','AnnotatedFeature','FeatureType','Probe','ProbeSet','ProbeFeature') NOT NULL,
  `parent` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`unmapped_object_id`),
  UNIQUE KEY `unique_unmapped_obj_idx` (`ensembl_id`,`ensembl_object_type`,`identifier`,`unmapped_reason_id`,`parent`,`external_db_id`),
  KEY `anal_exdb_idx` (`analysis_id`,`external_db_id`),
  KEY `id_idx` (`identifier`(50)),
  KEY `ext_db_identifier_idx` (`external_db_id`,`identifier`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `unmapped_reason` (
  `unmapped_reason_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `summary_description` varchar(255) DEFAULT NULL,
  `full_description` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`unmapped_reason_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

CREATE TABLE `xref` (
  `xref_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `external_db_id` int(10) unsigned DEFAULT NULL,
  `dbprimary_acc` varchar(512) NOT NULL,
  `display_label` varchar(512) NOT NULL,
  `version` varchar(10) DEFAULT NULL,
  `description` text,
  `info_type` enum('NONE','PROJECTION','MISC','DEPENDENT','DIRECT','SEQUENCE_MATCH','INFERRED_PAIR','PROBE','UNMAPPED','COORDINATE_OVERLAP','CHECKSUM') NOT NULL DEFAULT 'NONE',
  `info_text` varchar(255) NOT NULL DEFAULT '',
  PRIMARY KEY (`xref_id`),
  UNIQUE KEY `id_index` (`dbprimary_acc`,`external_db_id`,`info_type`,`info_text`,`version`),
  KEY `display_index` (`display_label`),
  KEY `info_type_idx` (`info_type`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 AVG_ROW_LENGTH=100;

