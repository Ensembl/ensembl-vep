CREATE TABLE `alt_id` (
  `alt_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `term_id` int(10) unsigned NOT NULL,
  `accession` varchar(64) NOT NULL,
  PRIMARY KEY (`alt_id`),
  UNIQUE KEY `term_alt_idx` (`term_id`,`alt_id`),
  KEY `accession_idx` (`accession`(50))
) ENGINE=MyISAM  ;

CREATE TABLE `aux_GO_Cross_product_review_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_goslim_aspergillus_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_goslim_candida_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_goslim_generic_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_goslim_metagenomics_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_goslim_pir_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_goslim_plant_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_goslim_pombe_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_goslim_yeast_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_gosubset_prok_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_high_level_annotation_qc_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_mf_needs_review_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_GO_virus_checked_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_SO_DBVAR_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_SO_SOFA_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `aux_SO_biosapiens_map` (
  `term_id` int(10) unsigned NOT NULL,
  `subset_term_id` int(10) unsigned NOT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  UNIQUE KEY `map_idx` (`term_id`,`subset_term_id`)
) ENGINE=MyISAM ;

CREATE TABLE `closure` (
  `closure_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `child_term_id` int(10) unsigned NOT NULL,
  `parent_term_id` int(10) unsigned NOT NULL,
  `subparent_term_id` int(10) unsigned DEFAULT NULL,
  `distance` tinyint(3) unsigned NOT NULL,
  `ontology_id` int(10) unsigned NOT NULL,
  `confident_relationship` tinyint(1) NOT NULL DEFAULT '0',
  PRIMARY KEY (`closure_id`),
  UNIQUE KEY `child_parent_idx` (`child_term_id`,`parent_term_id`,`subparent_term_id`,`ontology_id`),
  KEY `parent_subparent_idx` (`parent_term_id`,`subparent_term_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `meta` (
  `meta_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `meta_key` varchar(64) NOT NULL,
  `meta_value` varchar(128) DEFAULT NULL,
  `species_id` int(1) unsigned DEFAULT NULL,
  PRIMARY KEY (`meta_id`),
  UNIQUE KEY `key_value_idx` (`meta_key`,`meta_value`)
) ENGINE=MyISAM  ;

CREATE TABLE `ontology` (
  `ontology_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  `namespace` varchar(64) NOT NULL,
  `data_version` varchar(64) DEFAULT NULL,
  PRIMARY KEY (`ontology_id`),
  UNIQUE KEY `name_namespace_idx` (`name`,`namespace`)
) ENGINE=MyISAM  ;

CREATE TABLE `relation` (
  `relation_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `child_term_id` int(10) unsigned NOT NULL,
  `parent_term_id` int(10) unsigned NOT NULL,
  `relation_type_id` int(10) unsigned NOT NULL,
  `intersection_of` tinyint(3) unsigned NOT NULL DEFAULT '0',
  `ontology_id` int(10) unsigned NOT NULL,
  PRIMARY KEY (`relation_id`),
  UNIQUE KEY `child_parent_idx` (`child_term_id`,`parent_term_id`,`relation_type_id`,`intersection_of`,`ontology_id`),
  KEY `parent_idx` (`parent_term_id`)
) ENGINE=MyISAM  ;

CREATE TABLE `relation_type` (
  `relation_type_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  PRIMARY KEY (`relation_type_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM  ;

CREATE TABLE `subset` (
  `subset_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(64) NOT NULL,
  `definition` varchar(128) NOT NULL,
  PRIMARY KEY (`subset_id`),
  UNIQUE KEY `name_idx` (`name`)
) ENGINE=MyISAM  ;

CREATE TABLE `synonym` (
  `synonym_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `term_id` int(10) unsigned NOT NULL,
  `name` mediumtext COLLATE utf8_swedish_ci NOT NULL,
  `type` enum('EXACT','BROAD','NARROW','RELATED') COLLATE utf8_swedish_ci DEFAULT NULL,
  `dbxref` varchar(256) COLLATE utf8_swedish_ci NOT NULL,
  PRIMARY KEY (`synonym_id`),
  UNIQUE KEY `term_synonym_idx` (`term_id`,`synonym_id`),
  KEY `name_idx` (`name`(50))
) ENGINE=MyISAM  DEFAULT CHARSET=utf8 COLLATE=utf8_swedish_ci;

CREATE TABLE `term` (
  `term_id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `ontology_id` int(10) unsigned NOT NULL,
  `subsets` text,
  `accession` varchar(64) NOT NULL,
  `name` varchar(255) NOT NULL,
  `definition` text,
  `is_root` int(11) DEFAULT NULL,
  `is_obsolete` int(11) DEFAULT NULL,
  PRIMARY KEY (`term_id`),
  UNIQUE KEY `accession_idx` (`accession`),
  UNIQUE KEY `ontology_acc_idx` (`ontology_id`,`accession`),
  KEY `name_idx` (`name`)
) ENGINE=MyISAM  ;

