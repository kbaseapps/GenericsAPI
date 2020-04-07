/*
@author jjeffryes
*/
module KBaseMatrices{
    /*
      The workspace ID for a Genome data object.
      @id ws KBaseGenomes.Genome
    */
    typedef string ws_genome_id;

    /*
     * Reference to a handle ID
     * @id handle
     */
    typedef string handle_ref;

    /*
      The workspace ID for a Genome data object.
      @id ws
    */
    typedef string ws_ref;

    /*
      The workspace ID for a A data object
      @id ws KBaseExperiments.AttributeMapping
    */
    typedef string ws_attributemapping_id;

    /*
      A simple 2D matrix of floating point numbers with labels/ids for rows and
      columns.  The matrix is stored as a list of lists, with the outer list
      containing rows, and the inner lists containing values for each column of
      that row.  Row/Col ids should be unique.

      row_ids - unique ids for rows.
      col_ids - unique ids for columns.
      values - two dimensional array indexed as: values[row][col]
      @metadata ws length(row_ids) as n_rows
      @metadata ws length(col_ids) as n_cols
    */
    typedef structure {
      list<string> row_ids;
      list<string> col_ids;
      list<list<float>> values;
    } FloatMatrix2D;

    /*
      The workspace id for a single end or paired end reads object
      @id ws KBaseMatrices.DifferentialExpressionMatrix KBaseFeatureValues.DifferentialExpressionMatrix
    */
    typedef string differential_expression_matrix_ref;

    /*
      A wrapper around a FloatMatrix2D designed for simple matrices of Expression
      data.  Rows map to features, and columns map to conditions.  The data type
      includes some information about normalization factors and contains
      mappings from row ids to features and col ids to conditions.

      KBaseMatrices Fields:
      description - short optional description of the dataset
      scale - raw, ln, log2, log10
      col_normalization - mean_center, median_center, mode_center, zscore
      row_normalization - mean_center, median_center, mode_center, zscore
      col_mapping - map from col_id to an id in the col_condition_set
      row_mapping - map from row_id to a id in the row_condition_set
      col_attributemapping_ref - a reference to a AttributeMapping that relates to the columns
      row_attributemapping_ref - a reference to a AttributeMapping that relates to the rows
      attributes - a mapping of additional information pertaining to the object
      search_attributes - a mapping of object information used by search
      data - contains values for (feature,condition) pairs, where
          features correspond to rows and conditions are columns
          (ie data.values[feature][condition])

      Additional Fields:
      genome_ref - a reference to the aligned genome
      feature_mapping - map from row_id to feature id in the genome
      diff_expr_matrix_ref - added to connect filtered expression matrix to differential expression matrix
          used for filtering

      Validation:
      @unique data.row_ids
      @unique data.col_ids
      @conditionally_required row_attributemapping_ref row_mapping
      @conditionally_required col_attributemapping_ref col_mapping
      @contains data.row_ids row_mapping
      @contains data.col_ids col_mapping
      @contains values(row_mapping) row_attributemapping_ref:instances
      @contains values(col_mapping) col_attributemapping_ref:instances
      @contains data.row_ids genome_ref:features.[*].id genome_ref:mrnas.[*].id genome_ref:cdss.[*].id genome_ref:non_codeing_features.[*].id
      @contains values(feature_mapping) genome_ref:features.[*].id genome_ref:mrnas.[*].id genome_ref:cdss.[*].id genome_ref:non_codeing_features.[*].id

      @optional description row_normalization col_normalization
      @optional col_mapping row_mapping col_attributemapping_ref row_attributemapping_ref
      @optional attributes search_attributes genome_ref feature_mapping diff_expr_matrix_ref

      @metadata ws scale
      @metadata ws row_normalization
      @metadata ws col_normalization
      @metadata ws genome_ref as genome
      @metadata ws col_attributemapping_ref as col_attribute_mapping
      @metadata ws row_attributemapping_ref as row_attribute_mapping
      @metadata ws length(data.row_ids) as feature_count
      @metadata ws length(data.col_ids) as condition_count
    */
    typedef structure {
      string description;
      string scale;
      string row_normalization;
      string col_normalization;
      mapping<string, string> col_mapping;
      ws_attributemapping_id col_attributemapping_ref;
      mapping<string, string> row_mapping;
      ws_attributemapping_id row_attributemapping_ref;
      mapping<string, string> attributes;
      list<string> search_attributes;
      ws_genome_id genome_ref;
      mapping<string, string> feature_mapping;
      differential_expression_matrix_ref diff_expr_matrix_ref;
      FloatMatrix2D data;
    } ExpressionMatrix;

    /*
      A wrapper around a FloatMatrix2D designed for simple matrices of Differential
      Expression data.  Rows map to features, and columns map to conditions.  The
      data type includes some information about normalization factors and contains
      mappings from row ids to features and col ids to conditions.

      KBaseMatrices Fields:
      description - short optional description of the dataset
      scale - raw, ln, log2, log10
      col_normalization - mean_center, median_center, mode_center, zscore
      row_normalization - mean_center, median_center, mode_center, zscore
      col_mapping - map from col_id to an id in the col_condition_set
      row_mapping - map from row_id to a id in the row_condition_set
      col_attributemapping_ref - a reference to a AttributeMapping that relates to the columns
      row_attributemapping_ref - a reference to a AttributeMapping that relates to the rows
      attributes - a mapping of additional information pertaining to the object
      search_attributes - a mapping of object information used by search

      data - contains values for (feature,condition) pairs, where
          features correspond to rows and conditions are columns
          (ie data.values[feature][condition])

      Additional Fields:
      genome_ref - a reference to the aligned genome
      feature_mapping - map from row_id to feature id in the genome

      Validation:
      @unique data.row_ids
      @unique data.col_ids
      @conditionally_required row_attributemapping_ref row_mapping
      @conditionally_required col_attributemapping_ref col_mapping
      @contains data.row_ids row_mapping
      @contains data.col_ids col_mapping
      @contains values(row_mapping) row_attributemapping_ref:instances
      @contains values(col_mapping) col_attributemapping_ref:instances
      @contains data.row_ids genome_ref:features.[*].id genome_ref:mrnas.[*].id genome_ref:cdss.[*].id genome_ref:non_codeing_features.[*].id
      @contains values(feature_mapping) genome_ref:features.[*].id genome_ref:mrnas.[*].id genome_ref:cdss.[*].id genome_ref:non_codeing_features.[*].id

      @optional description row_normalization col_normalization
      @optional col_mapping row_mapping col_attributemapping_ref row_attributemapping_ref
      @optional attributes search_attributes genome_ref feature_mapping

      @metadata ws scale
      @metadata ws row_normalization
      @metadata ws col_normalization
      @metadata ws genome_ref as genome
      @metadata ws col_attributemapping_ref as col_attribute_mapping
      @metadata ws row_attributemapping_ref as row_attribute_mapping
      @metadata ws length(data.row_ids) as feature_count
      @metadata ws length(data.col_ids) as condition_count
    */
    typedef structure {
      string description;
      string scale;
      string row_normalization;
      string col_normalization;
      mapping<string, string> col_mapping;
      ws_attributemapping_id col_attributemapping_ref;
      mapping<string, string> row_mapping;
      ws_attributemapping_id row_attributemapping_ref;
      mapping<string, string> attributes;
      list<string> search_attributes;
      ws_genome_id genome_ref;
      mapping<string, string> feature_mapping;
      FloatMatrix2D data;
    } DifferentialExpressionMatrix;

    /*
      A wrapper around a FloatMatrix2D designed for simple matrices of Fitness data
      for gene/feature knockouts.  Generally fitness is measured as growth rate
      for the knockout strain relative to wildtype.

      KBaseMatrices Fields:
      description - short optional description of the dataset
      scale - raw, ln, log2, log10
      col_normalization - mean_center, median_center, mode_center, zscore
      row_normalization - mean_center, median_center, mode_center, zscore
      col_mapping - map from col_id to an id in the col_condition_set
      row_mapping - map from row_id to a id in the row_condition_set
      col_attributemapping_ref - a reference to a AttributeMapping that relates to the columns
      row_attributemapping_ref - a reference to a AttributeMapping that relates to the rows
      attributes - a mapping of additional information pertaining to the object
      search_attributes - a mapping of object information used by search

      data - contains values for (feature,condition) pairs, where
          features correspond to rows and conditions are columns
          (ie data.values[feature][condition])

      Additional Fields:
      genome_ref - a reference to the aligned genome
      feature_mapping - map from row_id to a set feature ids in the genome

      Validation:
      @unique data.row_ids
      @unique data.col_ids
      @conditionally_required row_attributemapping_ref row_mapping
      @conditionally_required col_attributemapping_ref col_mapping
      @contains data.row_ids row_mapping
      @contains data.col_ids col_mapping
      @contains values(row_mapping) row_attributemapping_ref:instances
      @contains values(col_mapping) col_attributemapping_ref:instances
      @contains data.row_ids genome_ref:features.[*].id genome_ref:mrnas.[*].id genome_ref:cdss.[*].id genome_ref:non_codeing_features.[*].id
      @contains values(feature_mapping) genome_ref:features.[*].id genome_ref:mrnas.[*].id genome_ref:cdss.[*].id genome_ref:non_codeing_features.[*].id

      @optional description row_normalization col_normalization
      @optional col_mapping row_mapping col_attributemapping_ref row_attributemapping_ref
      @optional attributes search_attributes genome_ref feature_mapping

      @metadata ws scale
      @metadata ws row_normalization
      @metadata ws col_normalization
      @metadata ws genome_ref as genome
      @metadata ws col_attributemapping_ref as col_attribute_mapping
      @metadata ws row_attributemapping_ref as row_attribute_mapping
      @metadata ws length(data.row_ids) as feature_count
      @metadata ws length(data.col_ids) as condition_count
    */
    typedef structure {
      string description;
      string scale;
      string row_normalization;
      string col_normalization;
      mapping<string, string> col_mapping;
      ws_attributemapping_id col_attributemapping_ref;
      mapping<string, string> row_mapping;
      ws_attributemapping_id row_attributemapping_ref;
      mapping<string, string> attributes;
      list<string> search_attributes;
      ws_genome_id genome_ref;
      mapping<string, list<string>> feature_mapping;
      FloatMatrix2D data;
    } FitnessMatrix;

    /*
      A wrapper around a FloatMatrix2D designed for matrices of data assigned to individual reactions.

      The columns represent experimental conditions while the rows correspond to reactions from a single
      metabolic reconstruction or from a biochemistry object.

      KBaseMatrices Fields:
      description - short optional description of the dataset
      scale - raw, ln, log2, log10

      data - contains values for (reaction,condition) pairs, where
             reactions correspond to rows and conditions are columns
             (ie data.values[reaction][condition])

      Additional Fields:
      fbamodel_ref - a reference to a FBAModel object
      biochemistry_ref - a reference to a Biochemistry object
      expression_ref - a reference to a ExpressionMatrix object (from which reaction values can be derived)
      fba_ref - a reference to a FBA object (from which reaction fluxes can be derived)

      Validation:
      @unique data.row_ids
      @unique data.col_ids

      @optional description fbamodel_ref biochemistry_ref expression_ref fba_refs

      @metadata ws scale
      @metadata ws length(data.row_ids) as reaction_count
      @metadata ws length(data.col_ids) as condition_count
    */
    typedef structure {
      string description;
      string scale;
      list<ws_ref> fba_refs;
      ws_ref fbamodel_ref;
      ws_ref expression_ref;
      ws_ref biochemistry_ref;
      FloatMatrix2D data;
    } ReactionMatrix;

    /*
      A wrapper around a FloatMatrix2D designed for matrices of chemical concentration data. The
      columns represent experimental conditions while the rows correspond to individual
      identified metabolites

      KBaseMatrices Fields:
      description - short optional description of the dataset
      scale - raw, ln, log2, log10
      col_normalization - mean_center, median_center, mode_center, zscore
      row_normalization - mean_center, median_center, mode_center, zscore
      col_mapping - map from col_id to an id in the col_condition_set
      row_mapping - map from row_id to a id in the row_condition_set
      col_attributemapping_ref - a reference to a AttributeMapping that relates to the columns
      row_attributemapping_ref - a reference to a AttributeMapping that relates to the rows
      attributes - a mapping of additional information pertaining to the object
      search_attributes - a mapping of object information used by search

      data - contains values for (compound,condition) pairs, where
             compounds correspond to rows and conditions are columns
             (ie data.values[compound][condition])

      Additional Fields:
      biochemistry_ref - a reference to a biochemistry object
      biochemistry_mapping - map from row_id to a set compound ids in a biochemistry object

      Validation:
      @unique data.row_ids
      @unique data.col_ids
      @conditionally_required row_attributemapping_ref row_mapping
      @conditionally_required col_attributemapping_ref col_mapping
      @contains data.row_ids row_mapping
      @contains data.col_ids col_mapping
      @contains values(row_mapping) row_attributemapping_ref:instances
      @contains values(col_mapping) col_attributemapping_ref:instances
      @contains values(biochemistry_mapping) biochemistry_ref:compounds.[*].id

      @optional description row_normalization col_normalization
      @optional col_mapping row_mapping col_attributemapping_ref sample_set_ref
      @optional attributes search_attributes biochemistry_mapping unit type

      @metadata ws scale
      @metadata ws unit
      @metadata ws type
      @metadata ws row_normalization
      @metadata ws col_normalization
      @metadata ws col_attributemapping_ref as col_attribute_mapping
      @metadata ws row_attributemapping_ref as row_attribute_mapping
      @metadata ws sample_set_ref as sample_set
      @metadata ws length(data.row_ids) as compound_count
      @metadata ws length(data.col_ids) as sample_count
    */
    typedef structure {
      string description;
      string scale;
      string row_normalization;
      string col_normalization;
      mapping<string, string> col_mapping;
      ws_attributemapping_id col_attributemapping_ref;
      mapping<string, string> row_mapping;
      ws_attributemapping_id row_attributemapping_ref;
      mapping<string, string> attributes;
      list<string> search_attributes;
      ws_ref biochemistry_ref;
      mapping<string, list<string>> biochemistry_mapping;
      FloatMatrix2D data;
      ws_ref sample_set_ref;
      string unit;
      string type;
    } ChemicalAbundanceMatrix;

    /*
      A wrapper around a FloatMatrix2D designed for matrices of amplicon data. The
      columns represent experimental conditions while the rows correspond to individual
      amplicons.

      KBaseMatrices Fields:
      description - short optional description of the dataset
      scale - raw, ln, log2, log10
      col_normalization - mean_center, median_center, mode_center, zscore
      row_normalization - mean_center, median_center, mode_center, zscore
      col_mapping - map from col_id to an id in the col_condition_set
      row_mapping - map from row_id to a id in the row_condition_set
      col_attributemapping_ref - a reference to a AttributeMapping that relates to the columns
      row_attributemapping_ref - a reference to a AttributeMapping that relates to the rows
      attributes - a mapping of additional information pertaining to the object
      search_attributes - a mapping of object information used by search

      data - contains values for (amplicons,condition) pairs, where
             amplicons correspond to rows and conditions are columns
             (ie data.values[amplicons][condition])

      Additional Fields:
      reads_set_ref - a reference to the set of reads libraries that produced this table
      sequence_mapping - map from row_id to the representative sequence for that row

      Validation:
      @unique data.row_ids
      @unique data.col_ids
      @conditionally_required row_attributemapping_ref row_mapping
      @conditionally_required col_attributemapping_ref col_mapping
      @contains data.row_ids row_mapping
      @contains data.col_ids col_mapping
      @contains values(row_mapping) row_attributemapping_ref:instances
      @contains values(col_mapping) col_attributemapping_ref:instances

      @optional description row_normalization col_normalization
      @optional col_mapping row_mapping col_attributemapping_ref row_attributemapping_ref sample_set_ref
      @optional attributes search_attributes sequence_mapping reads_set_ref amplicon_set_ref
      @optional extraction_kit amplicon_type target_gene_region forward_primer_sequence
      @optional reverse_primer_sequence sequencing_platform sequencing_run sequencing_kit
      @optional sequencing_quality_filter_cutoff clustering_cutoff clustering_method sequencing_file_handle

      @metadata ws scale
      @metadata ws row_normalization
      @metadata ws col_normalization
      @metadata ws col_attributemapping_ref as col_attribute_mapping
      @metadata ws row_attributemapping_ref as row_attribute_mapping
      @metadata ws length(data.row_ids) as amplicon_count
      @metadata ws length(data.col_ids) as condition_count
      @metadata ws description
      @metadata ws extraction_kit
      @metadata ws amplicon_type
      @metadata ws target_gene_region
      @metadata ws forward_primer_sequence
      @metadata ws reverse_primer_sequence
      @metadata ws sequencing_platform
      @metadata ws sequencing_run
      @metadata ws sequencing_kit
      @metadata ws sequencing_quality_filter_cutoff
      @metadata ws clustering_cutoff
      @metadata ws clustering_method

    */
    typedef structure {
      string description;
      string scale;
      string row_normalization;
      string col_normalization;
      mapping<string, string> col_mapping;
      ws_attributemapping_id col_attributemapping_ref;
      mapping<string, string> row_mapping;
      ws_attributemapping_id row_attributemapping_ref;
      mapping<string, string> attributes;
      list<string> search_attributes;
      ws_ref reads_set_ref;
      mapping<string, string> sequence_mapping;
      FloatMatrix2D data;
      ws_ref amplicon_set_ref;
      ws_ref sample_set_ref;
      string extraction_kit;
      string amplicon_type;
      string target_gene_region;
      string forward_primer_sequence;
      string reverse_primer_sequence;
      string sequencing_platform;
      string sequencing_run;
      string sequencing_kit;
      string sequencing_quality_filter_cutoff;
      float clustering_cutoff;
      string clustering_method;
      handle_ref sequencing_file_handle;
    } AmpliconMatrix;
    /*
      A wrapper around a FloatMatrix2D designed for matrices of trait data for use in population
      studies. The columns represent genotypes while the rows correspond to traits.
      KBaseMatrices Fields:
      description - short optional description of the dataset
      scale - raw, ln, log2, log10
      col_normalization - mean_center, median_center, mode_center, zscore
      row_normalization - mean_center, median_center, mode_center, zscore
      col_mapping - map from col_id to an id in the col_condition_set
      row_mapping - map from row_id to a id in the row_condition_set
      col_attributemapping_ref - a reference to a AttributeMapping that relates to the columns
      row_attributemapping_ref - a reference to a AttributeMapping that relates to the rows
      attributes - a mapping of additional information pertaining to the object
      search_attributes - a mapping of object information used by search

      data - contains values for (genotype,trait) pairs, where
             traits correspond to the rows and genotypes are columns
             (ie data.values[amplicons][condition])

      Additional Fields:

      Validation:
      @unique data.row_ids
      @unique data.col_ids
      @conditionally_required row_attributemapping_ref row_mapping
      @conditionally_required col_attributemapping_ref col_mapping
      @contains data.row_ids row_mapping
      @contains data.col_ids col_mapping
      @contains values(row_mapping) row_attributemapping_ref:instances
      @contains values(col_mapping) col_attributemapping_ref:instances
      @contains set(trait_id,trait_description) row_attributemapping_ref:attributes.[*].attribute
      @contains set(individual_id,family_id,paternal_id,maternal_id,sex) col_attributemapping_ref:attributes.[*].attribute

      @optional description row_normalization col_normalization
      @optional col_mapping row_mapping
      @optional attributes search_attributes

      @metadata ws scale
      @metadata ws row_normalization
      @metadata ws col_normalization
      @metadata ws col_attributemapping_ref as col_attribute_mapping
      @metadata ws row_attributemapping_ref as row_attribute_mapping
      @metadata ws length(data.row_ids) as genotype_count
      @metadata ws length(data.col_ids) as trait_count
    */
    typedef structure {
      string description;
      string scale;
      string row_normalization;
      string col_normalization;
      mapping<string, string> col_mapping;
      ws_attributemapping_id col_attributemapping_ref;
      mapping<string, string> row_mapping;
      ws_attributemapping_id row_attributemapping_ref;
      mapping<string, string> attributes;
      list<string> search_attributes;
      FloatMatrix2D data;
    } TraitMatrix;
};
