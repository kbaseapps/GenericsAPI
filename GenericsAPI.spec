/*
A KBase module: GenericsAPI
*/
#include <KBaseExperiments.spec>

module GenericsAPI {
  /* A boolean - 0 for false, 1 for true.
  */
  typedef int boolean;

  /* An X/Y/Z style reference
  */
  typedef string obj_ref;

  /* workspace name of the object */
  typedef string workspace_name;

  /* Input of the fetch_data function
    obj_ref: generics object reference

    Optional arguments:
    generics_module: the generics data module to be retrieved from
                    e.g. for an given data type like below:
                    typedef structure {
                      FloatMatrix2D data;
                      condition_set_ref condition_set_ref;
                    } SomeGenericsMatrix;
                    generics_module should be
                    {'data': 'FloatMatrix2D',
                     'condition_set_ref': 'condition_set_ref'}
  */
  typedef structure {
    obj_ref obj_ref;
    mapping<string, string> generics_module;
  } FetchDataParams;

  /* Ouput of the fetch_data function
    data_matrix: a pandas dataframe in json format
  */
  typedef structure {
    string data_matrix;
  } FetchDataReturn;

  /* fetch_data: fetch generics data as pandas dataframe for a generics data object*/
  funcdef fetch_data(FetchDataParams params) returns(FetchDataReturn returnVal) authentication required;

  /* Input of the export_matrix function
    obj_ref: generics object reference

    Optional arguments:
    generics_module: select the generics data to be retrieved from
                        e.g. for an given data type like below:
                        typedef structure {
                          FloatMatrix2D data;
                          condition_set_ref condition_set_ref;
                        } SomeGenericsMatrix;
                        and only 'FloatMatrix2D' is needed
                        generics_module should be
                        {'data': FloatMatrix2D'}
  */
  typedef structure {
      obj_ref obj_ref;
      mapping<string, string> generics_module;
  } ExportParams;

  typedef structure {
      string shock_id;
  } ExportOutput;

  funcdef export_matrix (ExportParams params) returns (ExportOutput returnVal) authentication required;

  /* Input of the validate_data function
    obj_type: obj type e.g.: 'KBaseMatrices.ExpressionMatrix-1.1'
    data: data to be validated
  */
  typedef structure {
      string obj_type;
      mapping<string, string> data;
  } ValidateParams;

  typedef structure {
      boolean validated;
      mapping<string, string> failed_constraint;
  } ValidateOutput;

  /* validate_data: validate data*/
  funcdef validate_data (ValidateParams params) returns (ValidateOutput returnVal) authentication required;


  /* Input of the import_matrix_from_excel function
    obj_type: a type in KBaseMatrices
    input_shock_id: file shock id
    input_file_path: absolute file path
    input_staging_file_path: staging area file path
    matrix_name: matrix object name
    description: optional, a description of the matrix
    workspace_name: workspace name matrix object to be saved to

    optional:
    col_attributemapping_ref: column AttributeMapping reference
    row_attributemapping_ref: row AttributeMapping reference
    genome_ref: genome reference
    diff_expr_matrix_ref: DifferentialExpressionMatrix reference
    biochemistry_ref: (for ChemicalAbundanceMatrix)
    reads_set_ref: (raw data for AmpliconMatrix)
    sample_set_ref: SampleSet object reference

  */
  typedef structure {
      string obj_type;
      string input_shock_id;
      string input_file_path;
      string input_staging_file_path;
      string matrix_name;
      string amplicon_set_name;
      string scale;
      string description;
      workspace_name workspace_name;

      obj_ref genome_ref;
      obj_ref col_attributemapping_ref;
      obj_ref row_attributemapping_ref;
      obj_ref diff_expr_matrix_ref;
      obj_ref biochemistry_ref;
      obj_ref reads_set_ref;
      obj_ref sample_set_ref;

      string unit;
      string type;
  } ImportMatrixParams;

  typedef structure {
      string report_name;
      string report_ref;
      obj_ref matrix_obj_ref;
  } ImportMatrixOutput;

  /* import_matrix_from_excel: import matrix object from excel*/
  funcdef import_matrix_from_excel (ImportMatrixParams params) returns (ImportMatrixOutput returnVal) authentication required;

  typedef structure {
      string obj_type;
      string taxonomic_abundance_tsv;
      string taxonomic_fasta;
      string input_local_file;
      string matrix_name;
      string amplicon_set_name;
      string scale;
      string description;
      workspace_name workspace_name;

      obj_ref genome_ref;
      obj_ref col_attributemapping_ref;
      obj_ref row_attributemapping_ref;
      obj_ref diff_expr_matrix_ref;
      obj_ref biochemistry_ref;
      obj_ref reads_set_ref;
      obj_ref sample_set_ref;

      list<string> metadata_keys;
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
  } ImportOTUParams;

  /* import_matrix_from_biom: import matrix object from BIOM file format*/
  funcdef import_matrix_from_biom (ImportOTUParams params) returns (ImportMatrixOutput returnVal) authentication required;

  /* Input of the import_matrix_from_excel function
    obj_type: saving object data type
    obj_name: saving object name
    data: data to be saved
    workspace_name: workspace name matrix object to be saved to
  */
  typedef structure {
      string obj_type;
      string obj_name;
      mapping<string, string> data;
      workspace_name workspace_name;
  } SaveObjectParams;

  typedef structure {
      obj_ref obj_ref;
  } SaveObjectOutput;

  /* save_object: validate data constraints and save matrix object*/
  funcdef save_object (SaveObjectParams params) returns (SaveObjectOutput returnVal) authentication required;

  /* Input of the search_matrix function
    matrix_obj_ref: object reference of a matrix
    workspace_name: workspace name objects to be saved to
  */
  typedef structure {
      obj_ref matrix_obj_ref;
      workspace_name workspace_name;
  } MatrixSelectorParams;

  typedef structure {
      string report_name;
      string report_ref;
  } MatrixSelectorOutput;

  /* search_matrix: generate a HTML report that allows users to select feature ids*/
  funcdef search_matrix (MatrixSelectorParams params) returns (MatrixSelectorOutput returnVal) authentication required;

  /* Input of the filter_matrix function
    matrix_obj_ref: object reference of a matrix
    workspace_name: workspace name objects to be saved to
    filter_ids: string of column or row ids that result matrix contains
    filtered_matrix_name: name of newly created filtered matrix object
  */
  typedef structure {
      obj_ref matrix_obj_ref;
      workspace_name workspace_name;
      string filtered_matrix_name;
      string remove_ids;
      string dimension;
  } MatrixFilterParams;

  typedef structure {
      string report_name;
      string report_ref;
      list<obj_ref> matrix_obj_refs;
  } MatrixFilterOutput;

  /* filter_matrix: create sub-matrix based on input filter_ids*/
  funcdef filter_matrix (MatrixFilterParams params) returns (MatrixFilterOutput returnVal) authentication required;

  /* Input of the standardize_matrix function
    input_matrix_ref: object reference of a matrix
    workspace_name: workspace name objects to be saved to
    with_mean: center data before scaling
    with_std: scale data to unit variance
    new_matrix_name: name of newly created matrix object
  */
  typedef structure {
      obj_ref input_matrix_ref;
      workspace_name workspace_name;
      boolean with_mean;
      boolean with_std;
      string dimension;
      string new_matrix_name;
  } StandardizeMatrixParams;

  typedef structure {
      string report_name;
      string report_ref;
      obj_ref new_matrix_obj_ref;
  } StandardizeMatrixOutput;

  /* standardize_matrix: standardize a matrix*/
  funcdef standardize_matrix (StandardizeMatrixParams params) returns (StandardizeMatrixOutput returnVal) authentication required;


  typedef structure {
      obj_ref input_matrix_ref;
      workspace_name workspace_name;
      int workspace_id;
      string new_matrix_name;
      mapping<string, string> abundance_filtering_params;
      mapping<string, string> standardization_params;
      mapping<string, string> ratio_transformation_params;
      boolean perform_relative_abundance;
  } TransformMatrixParams;

  typedef structure {
      string report_name;
      string report_ref;
      obj_ref new_matrix_obj_ref;
  } TransformMatrixOutput;

  funcdef transform_matrix (TransformMatrixParams params) returns (TransformMatrixOutput returnVal) authentication required;

  typedef structure {
      obj_ref input_matrix_ref;
      int workspace_id;
      string new_matrix_name;
      int seed_number;
  } RarefyMatrixParams;

  typedef structure {
      string report_name;
      string report_ref;
      obj_ref new_matrix_obj_ref;
  } RarefyMatrixOutput;

  funcdef perform_rarefy (RarefyMatrixParams params) returns (RarefyMatrixOutput returnVal) authentication required;

  typedef structure {
      obj_ref input_matrix_ref;
      obj_ref attribute_mapping_obj_ref;
      int workspace_id;
      string dist_metric;
      string dimension;
      string grouping;
      int permutations;
      boolean perform_anosim;
      boolean perform_permanova;
      boolean perform_permdisp;
  } VariableStatsParams;

  typedef structure {
      string report_name;
      string report_ref;
  } VariableStatsOutput;

  funcdef perform_variable_stats_matrix (VariableStatsParams params) returns (VariableStatsOutput returnVal) authentication required;


  typedef structure {
      list<obj_ref> input_matrix_refs;
      int workspace_id;
      string dist_metric;
      string dimension;
      string correlation_method;
      int permutations;
      string alternative_hypothesis;
  } MantelTestParams;

  typedef structure {
      string report_name;
      string report_ref;
  } MantelTestOutput;

  funcdef perform_mantel_test (MantelTestParams params) returns (MantelTestOutput returnVal) authentication required;


  /* ATTRIBUTE MAPPING */

    /*
        input_shock_id and input_file_path - alternative input params,
    */
    typedef structure {
        string input_shock_id;
        string input_file_path;
        string output_ws_id;
        string output_obj_name;
    } FileToAttributeMappingParams;

    typedef structure {
        obj_ref attribute_mapping_ref;
    } FileToAttributeMappingOutput;

    funcdef file_to_attribute_mapping(FileToAttributeMappingParams params)
        returns (FileToAttributeMappingOutput result) authentication required;

    funcdef file_to_fbamodel_attribute_mapping(FileToAttributeMappingParams params)
        returns (FileToAttributeMappingOutput result) authentication required;

    typedef structure {
        string staging_file_subdir_path;
        string dimension;
        obj_ref input_matrix_ref;
        string workspace_name;
        string output_am_obj_name;
        string output_matrix_obj_name;
    } UpdateMatrixAMParams;

    typedef structure {
        string report_name;
        string report_ref;
        obj_ref new_matrix_obj_ref;
        obj_ref new_attribute_mapping_ref;
    } UpdateMatrixAMOutput;

    funcdef update_matrix_attribute_mapping(UpdateMatrixAMParams params)
        returns (UpdateMatrixAMOutput returnVal) authentication required;

    typedef structure {
        obj_ref input_ref;
        string destination_dir;
    } AttributeMappingToTsvFileParams;

    typedef structure {
        string file_path;
    } AttributeMappingToTsvFileOutput;

    funcdef attribute_mapping_to_tsv_file(AttributeMappingToTsvFileParams params)
        returns (AttributeMappingToTsvFileOutput result) authentication required;

    typedef structure {
        obj_ref input_ref;
    } ExportObjectParams;

    funcdef export_attribute_mapping_tsv(ExportObjectParams params)
        returns (ExportOutput result) authentication required;

    funcdef export_attribute_mapping_excel(ExportObjectParams params)
        returns (ExportOutput result) authentication required;

    funcdef export_cluster_set_excel(ExportObjectParams params)
        returns (ExportOutput result) authentication required;

    funcdef export_corr_matrix_excel(ExportObjectParams params)
        returns (ExportOutput result) authentication required;

    funcdef export_pca_matrix_excel(ExportObjectParams params)
        returns (ExportOutput result) authentication required;

    funcdef export_amplicon_set_tsv(ExportObjectParams params)
        returns (ExportOutput result) authentication required;


  /* Input of the compute_correlation_matrix function
    input_obj_ref: object reference of a matrix
    workspace_name: workspace name objects to be saved to
    corr_matrix_name: correlation matrix object name
    dimension: compute correlation on column or row, one of ['col', 'row']
    method: correlation method, one of ['pearson', 'kendall', 'spearman']
    plot_corr_matrix: plot correlation matrix in report, default False
    plot_scatter_matrix: plot scatter matrix in report, default False
    compute_significance: also compute Significance in addition to correlation matrix
  */
  typedef structure {
      obj_ref input_obj_ref;
      workspace_name workspace_name;
      string corr_matrix_name;
      string dimension;
      string method;
      boolean plot_corr_matrix;
      boolean plot_scatter_matrix;
      boolean compute_significance;
  } CompCorrParams;

  typedef structure {
      string report_name;
      string report_ref;
      obj_ref corr_matrix_obj_ref;
  } CompCorrOutput;

  /* compute_correlation_matrix: create sub-matrix based on input filter_ids*/
  funcdef compute_correlation_matrix (CompCorrParams params) returns (CompCorrOutput returnVal) authentication required;

  /* Input of the compute_correlation_across_matrices function
    matrix_ref_1: object reference of a matrix
    matrix_ref_2: object reference of a matrix
    workspace_name: workspace name objects to be saved to
    corr_matrix_name: correlation matrix object name
    dimension: compute correlation on column or row, one of ['col', 'row']
    method: correlation method, one of ['pearson', 'kendall', 'spearman']
    plot_corr_matrix: plot correlation matrix in report, default False
    compute_significance: also compute Significance in addition to correlation matrix
  */
  typedef structure {
      obj_ref matrix_ref_1;
      obj_ref matrix_ref_2;
      workspace_name workspace_name;
      string corr_matrix_name;
      string dimension;
      string method;
      boolean plot_corr_matrix;
      boolean compute_significance;
      float corr_threshold;
  } CompCorrMetriceParams;

  /* compute_correlation_across_matrices: compute correlation matrix across matrices*/
  funcdef compute_correlation_across_matrices (CompCorrMetriceParams params) returns (CompCorrOutput returnVal) authentication required;

  /* Input of the build_network function
    corr_matrix_ref: CorrelationMatrix object
    workspace_name: workspace name objects to be saved to
    network_obj_name: Network object name
    filter_on_threshold: Dictory holder that holds filter on thredshold params
    params in filter_on_threshold:
      coefficient_threshold: correlation coefficient threshold (select pairs with greater correlation coefficient)
  */
  typedef structure {
      obj_ref corr_matrix_ref;
      workspace_name workspace_name;
      string network_obj_name;
      mapping<string, string> filter_on_threshold;
  } BuildNetworkParams;

  typedef structure {
      string report_name;
      string report_ref;
      obj_ref network_obj_ref;
  } BuildNetworkOutput;

  /* build_network: filter correlation matrix and build network*/
  funcdef build_network (BuildNetworkParams params) returns (BuildNetworkOutput returnVal) authentication required;

  /* Input of the run_pca function
    input_obj_ref: object reference of a matrix
    workspace_name: the name of the workspace
    pca_matrix_name: name of PCA (KBaseExperiments.PCAMatrix) object
    dimension: compute PCA on column or row, one of ['col', 'row']
    n_components - number of components (default 2)
    attribute_mapping_obj_ref - associated attribute_mapping_obj_ref
    scale_size_by - used for PCA plot to scale data size
    color_marker_by - used for PCA plot to group data
  */
  typedef structure {
    obj_ref input_obj_ref;
    string workspace_name;
    string pca_matrix_name;
    string dimension;
    int n_components;
    obj_ref attribute_mapping_obj_ref;
    mapping<string, string> scale_size_by;
    mapping<string, string> color_marker_by;
  } PCAParams;

  /* Ouput of the run_pca function
    pca_ref: PCA object reference (as KBaseExperiments.PCAMatrix data type)
    report_name: report name generated by KBaseReport
    report_ref: report reference generated by KBaseReport
  */
  typedef structure {
    obj_ref pca_ref;
    string report_name;
    string report_ref;
  } PCAOutput;

  /* run_pca: PCA analysis on matrix*/
  funcdef run_pca (PCAParams params) returns (PCAOutput returnVal) authentication required;

  typedef structure {
    obj_ref input_matrix_ref;
    string workspace_name;
    boolean with_attribute_info;
  } ViewMatrixParams;

  typedef structure {
    string report_name;
    string report_ref;
  } ViewMatrixOutput;

  /* view_matrix: generate a report for matrix viewer*/
  funcdef view_matrix (ViewMatrixParams params) returns (ViewMatrixOutput returnVal) authentication required;

  typedef structure {
    string workspace_name;
    int workspace_id;
    obj_ref sample_set_ref;
    mapping<string, int> chemical_data_included;
    mapping<string, int> chemical_ids_included;
  } ChemAbunTempParams;

  funcdef build_chemical_abundance_template(ChemAbunTempParams params) returns (ViewMatrixOutput returnVal) authentication required;


};
