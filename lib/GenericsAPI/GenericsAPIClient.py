# -*- coding: utf-8 -*-
############################################################
#
# Autogenerated by the KBase type compiler -
# any changes made here will be overwritten
#
############################################################

from __future__ import print_function
# the following is a hack to get the baseclient to import whether we're in a
# package or not. This makes pep8 unhappy hence the annotations.
try:
    # baseclient and this client are in a package
    from .baseclient import BaseClient as _BaseClient  # @UnusedImport
except:
    # no they aren't
    from baseclient import BaseClient as _BaseClient  # @Reimport


class GenericsAPI(object):

    def __init__(
            self, url=None, timeout=30 * 60, user_id=None,
            password=None, token=None, ignore_authrc=False,
            trust_all_ssl_certificates=False,
            auth_svc='https://kbase.us/services/authorization/Sessions/Login'):
        if url is None:
            raise ValueError('A url is required')
        self._service_ver = None
        self._client = _BaseClient(
            url, timeout=timeout, user_id=user_id, password=password,
            token=token, ignore_authrc=ignore_authrc,
            trust_all_ssl_certificates=trust_all_ssl_certificates,
            auth_svc=auth_svc)

    def fetch_data(self, params, context=None):
        """
        fetch_data: fetch generics data as pandas dataframe for a generics data object
        :param params: instance of type "FetchDataParams" (Input of the
           fetch_data function obj_ref: generics object reference Optional
           arguments: generics_module: the generics data module to be
           retrieved from e.g. for an given data type like below: typedef
           structure { FloatMatrix2D data; condition_set_ref
           condition_set_ref; } SomeGenericsMatrix; generics_module should be
           {'data': 'FloatMatrix2D', 'condition_set_ref':
           'condition_set_ref'}) -> structure: parameter "obj_ref" of type
           "obj_ref" (An X/Y/Z style reference), parameter "generics_module"
           of mapping from String to String
        :returns: instance of type "FetchDataReturn" (Ouput of the fetch_data
           function data_matrix: a pandas dataframe in json format) ->
           structure: parameter "data_matrix" of String
        """
        return self._client.call_method(
            'GenericsAPI.fetch_data',
            [params], self._service_ver, context)

    def export_matrix(self, params, context=None):
        """
        :param params: instance of type "ExportParams" (Input of the
           export_matrix function obj_ref: generics object reference Optional
           arguments: generics_module: select the generics data to be
           retrieved from e.g. for an given data type like below: typedef
           structure { FloatMatrix2D data; condition_set_ref
           condition_set_ref; } SomeGenericsMatrix; and only 'FloatMatrix2D'
           is needed generics_module should be {'data': FloatMatrix2D'}) ->
           structure: parameter "obj_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "generics_module" of mapping from String to
           String
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        return self._client.call_method(
            'GenericsAPI.export_matrix',
            [params], self._service_ver, context)

    def validate_data(self, params, context=None):
        """
        validate_data: validate data
        :param params: instance of type "ValidateParams" (Input of the
           validate_data function obj_type: obj type e.g.:
           'KBaseMatrices.ExpressionMatrix-1.1' data: data to be validated)
           -> structure: parameter "obj_type" of String, parameter "data" of
           mapping from String to String
        :returns: instance of type "ValidateOutput" -> structure: parameter
           "validated" of type "boolean" (A boolean - 0 for false, 1 for
           true. @range (0, 1)), parameter "failed_constraint" of mapping
           from String to String
        """
        return self._client.call_method(
            'GenericsAPI.validate_data',
            [params], self._service_ver, context)

    def import_matrix_from_excel(self, params, context=None):
        """
        import_matrix_from_excel: import matrix object from excel
        :param params: instance of type "ImportMatrixParams" (Input of the
           import_matrix_from_excel function obj_type: a type in
           KBaseMatrices input_shock_id: file shock id input_file_path:
           absolute file path input_staging_file_path: staging area file path
           matrix_name: matrix object name description: optional, a
           description of the matrix workspace_name: workspace name matrix
           object to be saved to optional: col_attributemapping_ref: column
           AttributeMapping reference row_attributemapping_ref: row
           AttributeMapping reference genome_ref: genome reference
           diff_expr_matrix_ref: DifferentialExpressionMatrix reference) ->
           structure: parameter "obj_type" of String, parameter
           "input_shock_id" of String, parameter "input_file_path" of String,
           parameter "input_staging_file_path" of String, parameter
           "matrix_name" of String, parameter "scale" of String, parameter
           "description" of String, parameter "workspace_name" of type
           "workspace_name" (workspace name of the object), parameter
           "genome_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "col_attributemapping_ref" of type "obj_ref" (An X/Y/Z
           style reference), parameter "row_attributemapping_ref" of type
           "obj_ref" (An X/Y/Z style reference), parameter
           "diff_expr_matrix_ref" of type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "ImportMatrixOutput" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "matrix_obj_ref" of type "obj_ref" (An X/Y/Z
           style reference)
        """
        return self._client.call_method(
            'GenericsAPI.import_matrix_from_excel',
            [params], self._service_ver, context)

    def save_object(self, params, context=None):
        """
        save_object: validate data constraints and save matrix object
        :param params: instance of type "SaveObjectParams" (Input of the
           import_matrix_from_excel function obj_type: saving object data
           type obj_name: saving object name data: data to be saved
           workspace_name: workspace name matrix object to be saved to) ->
           structure: parameter "obj_type" of String, parameter "obj_name" of
           String, parameter "data" of mapping from String to String,
           parameter "workspace_name" of type "workspace_name" (workspace
           name of the object)
        :returns: instance of type "SaveObjectOutput" -> structure: parameter
           "obj_ref" of type "obj_ref" (An X/Y/Z style reference)
        """
        return self._client.call_method(
            'GenericsAPI.save_object',
            [params], self._service_ver, context)

    def search_matrix(self, params, context=None):
        """
        search_matrix: generate a HTML report that allows users to select feature ids
        :param params: instance of type "MatrixSelectorParams" (Input of the
           search_matrix function matrix_obj_ref: object reference of a
           matrix workspace_name: workspace name objects to be saved to) ->
           structure: parameter "matrix_obj_ref" of type "obj_ref" (An X/Y/Z
           style reference), parameter "workspace_name" of type
           "workspace_name" (workspace name of the object)
        :returns: instance of type "MatrixSelectorOutput" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        return self._client.call_method(
            'GenericsAPI.search_matrix',
            [params], self._service_ver, context)

    def filter_matrix(self, params, context=None):
        """
        filter_matrix: create sub-matrix based on input filter_ids
        :param params: instance of type "MatrixFilterParams" (Input of the
           filter_matrix function matrix_obj_ref: object reference of a
           matrix workspace_name: workspace name objects to be saved to
           filter_ids: string of column or row ids that result matrix
           contains filtered_matrix_name: name of newly created filtered
           matrix object) -> structure: parameter "matrix_obj_ref" of type
           "obj_ref" (An X/Y/Z style reference), parameter "workspace_name"
           of type "workspace_name" (workspace name of the object), parameter
           "filter_ids" of String, parameter "filtered_matrix_name" of String
        :returns: instance of type "MatrixFilterOutput" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "matrix_obj_refs" of list of type "obj_ref" (An
           X/Y/Z style reference)
        """
        return self._client.call_method(
            'GenericsAPI.filter_matrix',
            [params], self._service_ver, context)

    def file_to_attribute_mapping(self, params, context=None):
        """
        :param params: instance of type "FileToAttributeMappingParams"
           (input_shock_id and input_file_path - alternative input params,)
           -> structure: parameter "input_shock_id" of String, parameter
           "input_file_path" of String, parameter "output_ws_id" of String,
           parameter "output_obj_name" of String
        :returns: instance of type "FileToAttributeMappingOutput" ->
           structure: parameter "attribute_mapping_ref" of type "obj_ref" (An
           X/Y/Z style reference)
        """
        return self._client.call_method(
            'GenericsAPI.file_to_attribute_mapping',
            [params], self._service_ver, context)

    def attribute_mapping_to_tsv_file(self, params, context=None):
        """
        :param params: instance of type "AttributeMappingToTsvFileParams" ->
           structure: parameter "input_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "destination_dir" of String
        :returns: instance of type "AttributeMappingToTsvFileOutput" ->
           structure: parameter "file_path" of String
        """
        return self._client.call_method(
            'GenericsAPI.attribute_mapping_to_tsv_file',
            [params], self._service_ver, context)

    def export_attribute_mapping_tsv(self, params, context=None):
        """
        :param params: instance of type "ExportObjectParams" -> structure:
           parameter "input_ref" of type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        return self._client.call_method(
            'GenericsAPI.export_attribute_mapping_tsv',
            [params], self._service_ver, context)

    def export_attribute_mapping_excel(self, params, context=None):
        """
        :param params: instance of type "ExportObjectParams" -> structure:
           parameter "input_ref" of type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        return self._client.call_method(
            'GenericsAPI.export_attribute_mapping_excel',
            [params], self._service_ver, context)

    def export_cluster_set_excel(self, params, context=None):
        """
        :param params: instance of type "ExportObjectParams" -> structure:
           parameter "input_ref" of type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        return self._client.call_method(
            'GenericsAPI.export_cluster_set_excel',
            [params], self._service_ver, context)

    def export_corr_matrix_excel(self, params, context=None):
        """
        :param params: instance of type "ExportObjectParams" -> structure:
           parameter "input_ref" of type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "ExportOutput" -> structure: parameter
           "shock_id" of String
        """
        return self._client.call_method(
            'GenericsAPI.export_corr_matrix_excel',
            [params], self._service_ver, context)

    def compute_correlation_matrix(self, params, context=None):
        """
        compute_correlation_matrix: create sub-matrix based on input filter_ids
        :param params: instance of type "CompCorrParams" (Input of the
           filter_matrix function input_obj_ref: object reference of a matrix
           workspace_name: workspace name objects to be saved to
           corr_matrix_name: correlation matrix object name dimension:
           compute correlation on column or row, one of ['col', 'row']
           method: correlation method, one of ['pearson', 'kendall',
           'spearman'] plot_corr_matrix: plot correlation matrix in report,
           default False plot_scatter_matrix: plot scatter matrix in report,
           default False compute_significance: also compute Significance in
           addition to correlation matrix) -> structure: parameter
           "input_obj_ref" of type "obj_ref" (An X/Y/Z style reference),
           parameter "workspace_name" of type "workspace_name" (workspace
           name of the object), parameter "corr_matrix_name" of String,
           parameter "dimension" of String, parameter "method" of String,
           parameter "plot_corr_matrix" of type "boolean" (A boolean - 0 for
           false, 1 for true. @range (0, 1)), parameter "plot_scatter_matrix"
           of type "boolean" (A boolean - 0 for false, 1 for true. @range (0,
           1)), parameter "compute_significance" of type "boolean" (A boolean
           - 0 for false, 1 for true. @range (0, 1))
        :returns: instance of type "CompCorrOutput" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String,
           parameter "corr_matrix_obj_ref" of type "obj_ref" (An X/Y/Z style
           reference)
        """
        return self._client.call_method(
            'GenericsAPI.compute_correlation_matrix',
            [params], self._service_ver, context)

    def build_network(self, params, context=None):
        """
        build_network: filter correlation matrix and build network
        :param params: instance of type "BuildNetworkParams" (Input of the
           build_network function corr_matrix_ref: CorrelationMatrix object
           workspace_name: workspace name objects to be saved to
           network_obj_name: Network object name filter_on_threshold: Dictory
           holder that holds filter on thredshold params params in
           filter_on_threshold: coefficient_threshold: correlation
           coefficient threshold (select pairs with greater correlation
           coefficient)) -> structure: parameter "corr_matrix_ref" of type
           "obj_ref" (An X/Y/Z style reference), parameter "workspace_name"
           of type "workspace_name" (workspace name of the object), parameter
           "network_obj_name" of String, parameter "filter_on_threshold" of
           mapping from String to String
        :returns: instance of type "BuildNetworkOutput" -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "network_obj_ref" of type "obj_ref" (An X/Y/Z
           style reference)
        """
        return self._client.call_method(
            'GenericsAPI.build_network',
            [params], self._service_ver, context)

    def run_pca(self, params, context=None):
        """
        run_pca: PCA analysis on matrix
        :param params: instance of type "PCAParams" (Input of the run_pca
           function input_obj_ref: object reference of a matrix
           workspace_name: the name of the workspace pca_matrix_name: name of
           PCA (KBaseExperiments.PCAMatrix) object dimension: compute PCA on
           column or row, one of ['col', 'row'] n_components - number of
           components (default 2) attribute_mapping_obj_ref - associated
           attribute_mapping_obj_ref customize_instance_group - customer and
           select which instance group to plot scale_size_by - used for PCA
           plot to scale data size) -> structure: parameter "input_obj_ref"
           of type "obj_ref" (An X/Y/Z style reference), parameter
           "workspace_name" of String, parameter "pca_matrix_name" of String,
           parameter "dimension" of String, parameter "n_components" of Long,
           parameter "attribute_mapping_obj_ref" of type "obj_ref" (An X/Y/Z
           style reference), parameter "customize_instance_group" of list of
           mapping from String to String, parameter "scale_size_by" of
           mapping from String to String
        :returns: instance of type "PCAOutput" (Ouput of the run_pca function
           pca_ref: PCA object reference (as KBaseExperiments.PCAMatrix data
           type) report_name: report name generated by KBaseReport
           report_ref: report reference generated by KBaseReport) ->
           structure: parameter "pca_ref" of type "obj_ref" (An X/Y/Z style
           reference), parameter "report_name" of String, parameter
           "report_ref" of String
        """
        return self._client.call_method(
            'GenericsAPI.run_pca',
            [params], self._service_ver, context)

    def status(self, context=None):
        return self._client.call_method('GenericsAPI.status',
                                        [], self._service_ver, context)
