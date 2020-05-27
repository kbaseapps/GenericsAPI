import collections
import errno
import logging
import os
import re
import shutil
import uuid
import time
import traceback
import sys

import pandas as pd
from openpyxl import load_workbook
from xlrd.biffh import XLRDError
from sklearn import preprocessing
from skbio.stats.composition import ilr, clr
from skbio import DistanceMatrix
from skbio.stats.distance import anosim, permanova, permdisp
import scipy.spatial.distance as dist


from installed_clients.DataFileUtilClient import DataFileUtil
from GenericsAPI.Utils.AttributeUtils import AttributesUtil
from GenericsAPI.Utils.SampleServiceUtil import SampleServiceUtil
from GenericsAPI.Utils.DataUtil import DataUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.fba_toolsClient import fba_tools
from installed_clients.kb_GenericsReportClient import kb_GenericsReport

TYPE_ATTRIBUTES = {'description', 'scale', 'row_normalization', 'col_normalization'}
SCALE_TYPES = {'raw', 'ln', 'log2', 'log10'}


class MatrixUtil:

    def _validate_import_matrix_from_excel_params(self, params):
        """
        _validate_import_matrix_from_excel_params:
            validates params passed to import_matrix_from_excel method
        """
        logging.info('start validating import_matrix_from_excel params')

        # check for required parameters
        for p in ['obj_type', 'matrix_name', 'workspace_name', 'scale']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        obj_type = params.get('obj_type')
        if obj_type not in self.matrix_types:
            raise ValueError('Unknown matrix object type: {}'.format(obj_type))

        scale = params.get('scale')
        if scale not in SCALE_TYPES:
            raise ValueError('Unknown scale type: {}'.format(scale))

        if params.get('input_file_path'):
            file_path = params.get('input_file_path')
        elif params.get('input_shock_id'):
            file_path = self.dfu.shock_to_file(
                {'shock_id': params['input_shock_id'],
                 'file_path': self.scratch}).get('file_path')
        elif params.get('input_staging_file_path'):
            file_path = self.dfu.download_staging_file(
                        {'staging_file_subdir_path': params.get('input_staging_file_path')}
                        ).get('copy_file_path')
        else:
            error_msg = "Must supply either a input_shock_id or input_file_path "
            error_msg += "or input_staging_file_path"
            raise ValueError(error_msg)

        refs = {k: v for k, v in params.items() if "_ref" in k}

        return (obj_type, file_path, params.get('workspace_name'),
                params.get('matrix_name'), refs, scale)

    def _upload_to_shock(self, file_path):
        """
        _upload_to_shock: upload target file to shock using DataFileUtil
        """
        logging.info('Start uploading file to shock: {}'.format(file_path))

        file_to_shock_params = {
            'file_path': file_path,
            'pack': 'zip'
        }
        shock_id = self.dfu.file_to_shock(file_to_shock_params).get('shock_id')

        return shock_id

    @staticmethod
    def _mkdir_p(path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    @staticmethod
    def _find_between(s, start, end):
        """
        _find_between: find string in between start and end
        """

        return re.search('{}(.*){}'.format(start, end), s).group(1)

    @staticmethod
    def _write_mapping_sheet(file_path, sheet_name, mapping, index):
        """
        _write_mapping_sheet: write mapping to sheet
        """
        df_dict = collections.OrderedDict()

        df_dict[index[0]] = []
        df_dict[index[1]] = []

        for key, value in mapping.items():
            df_dict.get(index[0]).append(key)
            df_dict.get(index[1]).append(value)

        df = pd.DataFrame.from_dict(df_dict)

        with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
            writer.book = load_workbook(file_path)
            df.to_excel(writer, sheet_name=sheet_name)

    def _generate_tab_content(self, index_page, viewer_name):
        tab_content = ''

        if index_page:
            tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
            tab_content += '\n<iframe height="1200px" width="100%" '
            tab_content += 'src="{}" '.format(index_page)
            tab_content += 'style="border:none;"></iframe>'
            tab_content += '\n</div>\n'
        else:
            tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
            tab_content += '''\n<p style="color:red;" >'''
            tab_content += '''Matrix is too large to be displayed.</p>\n'''
            tab_content += '\n</div>\n'

        return tab_content

    def _generate_variable_stat_tab_content(self, res, viewer_name):
        tab_content = ''

        tab_content += '''\n<div id="{}" class="tabcontent">\n'''.format(viewer_name)
        tab_content += '''<table>\n'''
        for key, value in res.items():
            tab_content += '''<tr>\n'''
            tab_content += '''<td>{}</td>\n'''.format(key)
            tab_content += '''<td>{}</td>\n'''.format(value)
            tab_content += '''</tr>\n'''
        tab_content += '''</table>\n'''
        tab_content += '\n</div>\n'

        return tab_content

    def _generate_variable_stats_visualization_content(self, anosim_res,
                                                       permanova_res, permdisp_res):
        tab_def_content = ''
        tab_content = ''

        first_tab_token = False

        if anosim_res is not None:
            viewer_name = 'anosim_res'

            first_tab_token = True
            tab_def_content += '''\n<div class="tab">\n'''
            tab_def_content += '''\n<button class="tablinks" '''
            tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
            tab_def_content += ''' id="defaultOpen"'''
            tab_def_content += '''>Analysis of Similarities</button>\n'''

            tab_content += self._generate_variable_stat_tab_content(anosim_res, viewer_name)

        if permanova_res is not None:

            viewer_name = 'permanova_res'

            if first_tab_token:
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += '''>Permutational Multivariate Analysis of Variance</button>\n'''
            else:
                first_tab_token = True
                tab_def_content += '''\n<div class="tab">\n'''
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += ''' id="defaultOpen"'''
                tab_def_content += '''>Permutational Multivariate Analysis of Variance</button>\n'''

            tab_content += self._generate_variable_stat_tab_content(permanova_res, viewer_name)

        if permdisp_res is not None:
            viewer_name = 'permdisp_res'

            if first_tab_token:
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += '''>Homogeneity Multivariate Analysis of Variance</button>\n'''
            else:
                first_tab_token = True
                tab_def_content += '''\n<div class="tab">\n'''
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += ''' id="defaultOpen"'''
                tab_def_content += '''>Homogeneity Multivariate Analysis of Variance</button>\n'''

            tab_content += self._generate_variable_stat_tab_content(permdisp_res, viewer_name)

        tab_def_content += '\n</div>\n'
        return tab_def_content + tab_content

    def _generate_trans_visualization_content(self, output_directory, original_matrix_dir,
                                              filtered_matrix_dir, relative_abundance_matrix_dir,
                                              standardize_matrix_dir,
                                              ratio_transformed_matrix_dir):
        tab_def_content = ''
        tab_content = ''

        tab_def_content += '''
        <div class="tab">
            <button class="tablinks" onclick="openTab(event, 'OriginalMatrixViewer')" id="defaultOpen">Original Matrix</button>
        '''
        original_matrix_report_files = os.listdir(original_matrix_dir)
        original_matrix_index_page = None
        for original_matrix_report_file in original_matrix_report_files:
            if original_matrix_report_file.endswith('.html'):
                original_matrix_index_page = original_matrix_report_file

            shutil.copy2(os.path.join(original_matrix_dir, original_matrix_report_file),
                         output_directory)
        tab_content += self._generate_tab_content(original_matrix_index_page,
                                                  'OriginalMatrixViewer')

        if filtered_matrix_dir is not None:
            viewer_name = 'AbundanceFilteredMatrixViewer'
            tab_def_content += '''\n<button class="tablinks" '''
            tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
            tab_def_content += '''>Filtered Abundance Matrix</button>\n'''
            filtered_matrix_report_files = os.listdir(filtered_matrix_dir)
            filtered_matrix_index_page = None
            for filtered_matrix_report_file in filtered_matrix_report_files:
                if filtered_matrix_report_file.endswith('.html'):
                    filtered_matrix_index_page = filtered_matrix_report_file

                shutil.copy2(os.path.join(filtered_matrix_dir, filtered_matrix_report_file),
                             output_directory)
            tab_content += self._generate_tab_content(filtered_matrix_index_page, viewer_name)

        if relative_abundance_matrix_dir is not None:
            viewer_name = 'RelativeAbundanceMatrixViewer'
            tab_def_content += '''\n<button class="tablinks" '''
            tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
            tab_def_content += '''>Relative Abundance Matrix</button>\n'''
            relative_abundance_matrix_report_files = os.listdir(relative_abundance_matrix_dir)
            relative_abundance_matrix_index_page = None
            for relative_abundance_matrix_report_file in relative_abundance_matrix_report_files:
                if relative_abundance_matrix_report_file.endswith('.html'):
                    relative_abundance_matrix_index_page = relative_abundance_matrix_report_file

                shutil.copy2(os.path.join(relative_abundance_matrix_dir,
                                          relative_abundance_matrix_report_file),
                             output_directory)
            tab_content += self._generate_tab_content(relative_abundance_matrix_index_page,
                                                      viewer_name)

        if standardize_matrix_dir is not None:
            viewer_name = 'StandardizeMatrixViewer'
            tab_def_content += '''\n<button class="tablinks" '''
            tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
            tab_def_content += '''>Standardized Matrix</button>\n'''
            standardize_matrix_report_files = os.listdir(standardize_matrix_dir)
            standardize_matrix_index_page = None
            for standardize_matrix_report_file in standardize_matrix_report_files:
                if standardize_matrix_report_file.endswith('.html'):
                    standardize_matrix_index_page = standardize_matrix_report_file

                shutil.copy2(os.path.join(standardize_matrix_dir, standardize_matrix_report_file),
                             output_directory)
            tab_content += self._generate_tab_content(standardize_matrix_index_page, viewer_name)

        if ratio_transformed_matrix_dir is not None:
            viewer_name = 'RatioTransformedMatrixViewer'
            tab_def_content += '''\n<button class="tablinks" '''
            tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
            tab_def_content += '''>Log Ratio Transformed Matrix</button>\n'''
            ratio_transformed_matrix_report_files = os.listdir(ratio_transformed_matrix_dir)
            ratio_transformed_matrix_index_page = None
            for ratio_transformed_matrix_report_file in ratio_transformed_matrix_report_files:
                if ratio_transformed_matrix_report_file.endswith('.html'):
                    ratio_transformed_matrix_index_page = ratio_transformed_matrix_report_file

                shutil.copy2(os.path.join(ratio_transformed_matrix_dir,
                                          ratio_transformed_matrix_report_file),
                             output_directory)
            tab_content += self._generate_tab_content(ratio_transformed_matrix_index_page,
                                                      viewer_name)

        tab_def_content += '\n</div>\n'
        return tab_def_content + tab_content

    def _generate_visualization_content(self, output_directory, heatmap_dir):
        tab_def_content = ''
        tab_content = ''

        tab_def_content += '''
        <div class="tab">
            <button class="tablinks" onclick="openTab(event, 'MatrixViewer')" id="defaultOpen">Matrix Heatmap</button>
        </div>
        '''

        heatmap_report_files = os.listdir(heatmap_dir)

        heatmap_index_page = None
        for heatmap_report_file in heatmap_report_files:
            if heatmap_report_file.endswith('.html'):
                heatmap_index_page = heatmap_report_file

            shutil.copy2(os.path.join(heatmap_dir, heatmap_report_file),
                         output_directory)

        if heatmap_index_page:
            tab_content += '''\n<div id="MatrixViewer" class="tabcontent">'''
            tab_content += '\n<iframe height="1200px" width="100%" '
            tab_content += 'src="{}" '.format(heatmap_index_page)
            tab_content += 'style="border:none;"></iframe>'
            tab_content += '\n</div>\n'
        else:
            tab_content += '''\n<div id="MatrixViewer" class="tabcontent">'''
            tab_content += '''\n<p style="color:red;" >'''
            tab_content += '''Heatmap is too large to be displayed.</p>\n'''
            tab_content += '\n</div>\n'

        return tab_def_content + tab_content

    def _generate_variable_stats_html_report(self, anosim_res, permanova_res, permdisp_res):
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'variable_stats_viewer_report.html')

        visualization_content = self._generate_variable_stats_visualization_content(anosim_res,
                                                                                    permanova_res,
                                                                                    permdisp_res)

        table_style_content = '''
                                table {
                                  font-family: arial, sans-serif;
                                  border-collapse: collapse;
                                  width: 66%;
                                }

                                td, th {
                                  border: 1px solid #dddddd;
                                  text-align: left;
                                  padding: 8px;
                                }

                                tr:nth-child(even) {
                                  background-color: #dddddd;
                                }

                                </style>'''

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'templates', 'matrix_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Visualization_Content</p>',
                                                          visualization_content)
                report_template = report_template.replace('</style>',
                                                          table_style_content)
                result_file.write(report_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for Compute Correlation App'
                            })
        return html_report

    def _generate_transform_html_report(self, original_matrix_dir, filtered_matrix_dir,
                                        relative_abundance_matrix_dir,
                                        standardize_matrix_dir,
                                        ratio_transformed_matrix_dir):
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'transform_matrix_viewer_report.html')

        visualization_content = self._generate_trans_visualization_content(
                                                                    output_directory,
                                                                    original_matrix_dir,
                                                                    filtered_matrix_dir,
                                                                    relative_abundance_matrix_dir,
                                                                    standardize_matrix_dir,
                                                                    ratio_transformed_matrix_dir)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'templates', 'matrix_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Visualization_Content</p>',
                                                          visualization_content)
                result_file.write(report_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for Compute Correlation App'
                            })
        return html_report

    def _generate_heatmap_html_report(self, heatmap_dir):

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'matrix_viewer_report.html')

        visualization_content = self._generate_visualization_content(output_directory,
                                                                     heatmap_dir)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'templates', 'matrix_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Visualization_Content</p>',
                                                          visualization_content)
                result_file.write(report_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for Compute Correlation App'
                            })
        return html_report

    def _generate_transform_report(self, new_matrix_obj_ref, workspace_name, original_matrix_df,
                                   filtered_df, relative_abundance_df,
                                   standardize_df, ratio_transformed_df):
        objects_created = [{'ref': new_matrix_obj_ref, 'description': 'Transformed Matrix'}]

        data_tsv_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(data_tsv_directory)
        logging.info('Start generating matrix tsv files in {}'.format(data_tsv_directory))
        original_matrix_tsv_path = os.path.join(data_tsv_directory,
                                                'original_matrix_{}.tsv'.format(
                                                                                str(uuid.uuid4())))
        original_matrix_df.to_csv(original_matrix_tsv_path)
        original_matrix_dir = self.report_util.build_heatmap_html({
                                            'tsv_file_path': original_matrix_tsv_path})['html_dir']
        filtered_matrix_dir = None
        if filtered_df is not None:
            filtered_matrix_tsv_path = os.path.join(data_tsv_directory,
                                                    'filtered_matrix_{}.tsv'.format(
                                                                                str(uuid.uuid4())))
            filtered_df.to_csv(filtered_matrix_tsv_path)
            filtered_matrix_dir = self.report_util.build_heatmap_html({
                                            'tsv_file_path': filtered_matrix_tsv_path})['html_dir']

        relative_abundance_matrix_dir = None
        if relative_abundance_df is not None:
            relative_abundance_matrix_tsv_path = os.path.join(
                                                    data_tsv_directory,
                                                    'relative_abundance_matrix_{}.tsv'.format(
                                                                                str(uuid.uuid4())))
            relative_abundance_df.to_csv(relative_abundance_matrix_tsv_path)
            relative_abundance_matrix_dir = self.report_util.build_heatmap_html({
                                'tsv_file_path': relative_abundance_matrix_tsv_path})['html_dir']

        standardize_matrix_dir = None
        if standardize_df is not None:
            standardize_matrix_tsv_path = os.path.join(data_tsv_directory,
                                                       'standardize_matrix_{}.tsv'.format(
                                                                                str(uuid.uuid4())))
            standardize_df.to_csv(standardize_matrix_tsv_path)
            standardize_matrix_dir = self.report_util.build_heatmap_html({
                                        'tsv_file_path': standardize_matrix_tsv_path})['html_dir']

        ratio_transformed_matrix_dir = None
        if ratio_transformed_df is not None:
            ratio_transformed_matrix_tsv_path = os.path.join(
                                                    data_tsv_directory,
                                                    'ratio_transformed_matrix_{}.tsv'.format(
                                                                                str(uuid.uuid4())))
            ratio_transformed_df.to_csv(ratio_transformed_matrix_tsv_path)
            ratio_transformed_matrix_dir = self.report_util.build_heatmap_html({
                                'tsv_file_path': ratio_transformed_matrix_tsv_path})['html_dir']

        output_html_files = self._generate_transform_html_report(original_matrix_dir,
                                                                 filtered_matrix_dir,
                                                                 relative_abundance_matrix_dir,
                                                                 standardize_matrix_dir,
                                                                 ratio_transformed_matrix_dir)

        report_params = {'message': '',
                         'objects_created': objects_created,
                         'workspace_name': workspace_name,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 1300,
                         'report_object_name': 'transform_matrix_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_variable_stats_report(self, workspace_id,
                                        anosim_res, permanova_res, permdisp_res):

        output_html_files = self._generate_variable_stats_html_report(anosim_res,
                                                                      permanova_res,
                                                                      permdisp_res)

        report_params = {'message': '',
                         'workspace_id': workspace_id,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 600,
                         'report_object_name': 'variable_stats_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_report(self, matrix_obj_ref, workspace_name, new_row_attr_ref=None,
                         new_col_attr_ref=None, data=None):
        """
        _generate_report: generate summary report
        """

        objects_created = [{'ref': matrix_obj_ref, 'description': 'Imported Matrix'}]

        if new_row_attr_ref:
            objects_created.append({'ref': new_row_attr_ref,
                                    'description': 'Imported Row Attribute Mapping'})

        if new_col_attr_ref:
            objects_created.append({'ref': new_col_attr_ref,
                                    'description': 'Imported Column Attribute Mapping'})

        if data:
            data_df = pd.DataFrame(data['values'], index=data['row_ids'], columns=data['col_ids'])
            result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
            self._mkdir_p(result_directory)
            tsv_file_path = os.path.join(result_directory, 'heatmap_data_{}.tsv'.format(
                                                                                str(uuid.uuid4())))
            data_df.to_csv(tsv_file_path)
            heatmap_dir = self.report_util.build_heatmap_html({
                                                    'tsv_file_path': tsv_file_path})['html_dir']

            output_html_files = self._generate_heatmap_html_report(heatmap_dir)

            report_params = {'message': '',
                             'objects_created': objects_created,
                             'workspace_name': workspace_name,
                             'html_links': output_html_files,
                             'direct_html_link_index': 0,
                             'html_window_height': 1300,
                             'report_object_name': 'import_matrix_from_excel_' + str(uuid.uuid4())}

        else:
            report_params = {'message': '',
                             'objects_created': objects_created,
                             'workspace_name': workspace_name,
                             'report_object_name': 'import_matrix_from_excel_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    @staticmethod
    def _process_mapping_sheet(file_path, sheet_name):
        """
        _process_mapping: process mapping sheet
        """

        try:
            df = pd.read_excel(file_path, sheet_name=sheet_name, dtype='str')
        except XLRDError:
            return dict()
        else:
            mapping = {value[0]: value[1] for value in df.values.tolist()}

        return mapping

    def _process_attribute_mapping_sheet(self, file_path, sheet_name, matrix_name, workspace_id):
        """
        _process_attribute_mapping_sheet: process attribute_mapping sheet
        """

        try:
            df = pd.read_excel(file_path, sheet_name=sheet_name)
        except XLRDError:
            return ''
        else:
            obj_name = f'{matrix_name}_{sheet_name}'
            result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
            self._mkdir_p(result_directory)
            file_path = os.path.join(result_directory, '{}.xlsx'.format(obj_name))
            df.to_excel(file_path)
            import_attribute_mapping_params = {
                'output_obj_name': obj_name,
                'output_ws_id': workspace_id,
                'input_file_path': file_path
            }

            ref = self.attr_util.file_to_attribute_mapping(import_attribute_mapping_params)

            return ref.get('attribute_mapping_ref')

    @staticmethod
    def _file_to_df(file_path):
        logging.info('start parsing file content to data frame')

        try:
            df = pd.read_excel(file_path, sheet_name='data', index_col=0)

        except XLRDError:
            try:
                df = pd.read_excel(file_path, index_col=0)
                logging.warning('WARNING: A sheet named "data" was not found in the attached file,'
                                ' proceeding with the first sheet as the data sheet.')

            except XLRDError:

                try:
                    reader = pd.read_csv(file_path, sep=None, iterator=True)
                    inferred_sep = reader._engine.data.dialect.delimiter
                    df = pd.read_csv(file_path, sep=inferred_sep, index_col=0)
                except Exception:
                    raise ValueError('Cannot parse file. Please provide valide tsv, excel or csv file')

        df.index = df.index.astype('str')
        df.columns = df.columns.astype('str')
        # fill NA with "None" so that they are properly represented as nulls in the KBase Object
        df = df.where((pd.notnull(df)), None)

        return df

    def _file_to_chem_abun_data(self, file_path, refs, matrix_name, workspace_id):
        logging.info('Start reading and converting excel file data')
        data = refs

        df = self._file_to_df(file_path)

        metadata_df = None

        rename_map = {'Aggregate M/Z': 'aggregate_mz',
                      'Compound Name': 'name',
                      'Predicted Formula': 'formula',
                      'Predicted Structure (smiles)': 'smiles',
                      'Predicted Structure (inchi)': 'inchi',
                      'Predicted Structure (inchi-key)': 'inchikey',
                      'Theoretical Mass': 'mass',
                      'Retention Time': 'retention_time',
                      'Polarity': 'polarity',
                      'KEGG': 'kegg',
                      'ChemBi': 'chembi',
                      'ModelSEED': 'modelseed',
                      # 'Theoretical M/Z': 'theoretical_mz',
                      # 'Reference Standard RT (seconds)': 'reference_rt',
                      }
        df.rename(columns=rename_map, inplace=True)

        metadata_keys = rename_map.values()

        shared_metadata_keys = list(set(metadata_keys) & set(df.columns))
        if shared_metadata_keys:
            metadata_df = df[shared_metadata_keys]
            if set(metadata_df.all(skipna=False).tolist()) == {None}:
                raise ValueError('All of metadata fields are None')
            df.drop(columns=shared_metadata_keys, inplace=True)
        else:
            err_msg = 'Please provide at least one of below metadata fields:\n{}'.format(
                                                                        list(rename_map.keys()))
            raise ValueError(err_msg)

        try:
            df = df.astype(float)
        except ValueError:
            err_msg = 'Found some non-float values. Matrix contains only numeric values\n'
            err_msg += 'Please list any non-numeric column names in  Metadata Keys field'
            raise ValueError(err_msg)
        df.fillna(0, inplace=True)
        matrix_data = {'row_ids': df.index.tolist(),
                       'col_ids': df.columns.tolist(),
                       'values': df.values.tolist()}

        data.update({'data': matrix_data})

        data.update(self._get_axis_attributes('col', matrix_data, refs, file_path, matrix_name,
                                              workspace_id))
        data.update(self._get_axis_attributes('row', matrix_data, refs, file_path, matrix_name,
                                              workspace_id, metadata_df=metadata_df))

        return data

    def _file_to_data(self, file_path, refs, matrix_name, workspace_id):
        logging.info('Start reading and converting excel file data')
        data = refs

        df = self._file_to_df(file_path)

        matrix_data = {'row_ids': df.index.tolist(),
                       'col_ids': df.columns.tolist(),
                       'values': df.values.tolist()}

        data.update({'data': matrix_data})
        data.update(self._get_axis_attributes('col', matrix_data, refs, file_path, matrix_name,
                                              workspace_id))
        data.update(self._get_axis_attributes('row', matrix_data, refs, file_path, matrix_name,
                                              workspace_id))

        # processing metadata
        metadata = self._process_mapping_sheet(file_path, 'metadata')
        data['attributes'] = {}
        data['search_attributes'] = []
        for k, v in metadata.items():
            k = k.strip()
            v = v.strip()
            if k in TYPE_ATTRIBUTES:
                data[k] = v
            else:
                data['attributes'][k] = v
                data['search_attributes'].append(" | ".join((k, v)))

        return data

    def _sample_set_to_attribute_mapping(self, axis_ids, sample_set_ref, obj_name, ws_id):

        am_data = self.sampleservice_util.sample_set_to_attribute_mapping(sample_set_ref)

        unmatched_ids = set(axis_ids) - set(am_data['instances'].keys())
        if unmatched_ids:
            name = "Column"
            raise ValueError(f"The following {name} IDs from the uploaded matrix do not match "
                             f"the supplied {name} attribute mapping: {', '.join(unmatched_ids)}"
                             f"\nPlease verify the input data or upload an excel file with a"
                             f"{name} mapping tab.")

        logging.info('start saving AttributeMapping object: {}'.format(obj_name))
        info = self.dfu.save_objects({
            "id": ws_id,
            "objects": [{
                "type": "KBaseExperiments.AttributeMapping",
                "data": am_data,
                "name": obj_name
            }]
        })[0]

        return f'{info[6]}/{info[0]}/{info[4]}'

    def _meta_df_to_attribute_mapping(self, axis_ids, metadata_df, obj_name, ws_id):
        data = {'ontology_mapping_method': "TSV file", 'instances': {}}
        attribute_keys = metadata_df.columns.tolist()
        data['attributes'] = [{'attribute': key, 'source': 'upload'} for key in attribute_keys]

        for axis_id in axis_ids:
            data['instances'][axis_id] = [str(i) for i in metadata_df.loc[axis_id].tolist()]

        logging.info('start saving AttributeMapping object: {}'.format(obj_name))
        info = self.dfu.save_objects({
            "id": ws_id,
            "objects": [{
                "type": "KBaseExperiments.AttributeMapping",
                "data": data,
                "name": obj_name
            }]
        })[0]

        return f'{info[6]}/{info[0]}/{info[4]}'

    def _get_axis_attributes(self, axis, matrix_data, refs, file_path, matrix_name, workspace_id,
                             metadata_df=None):
        """Get the row/col_attributemapping and mapping of ids, validating as needed"""
        # Parameter specified mappings should take precedence over tabs in excel so only process
        # if attributemapping_ref is missing:
        attr_data = {}
        axis_ids = matrix_data[f'{axis}_ids']

        attributemapping_ref = None

        if refs.get('sample_set_ref') and axis == 'col':
            name = matrix_name + "_{}_attributes".format(axis)
            attributemapping_ref = self._sample_set_to_attribute_mapping(
                axis_ids, refs.get('sample_set_ref'), name, workspace_id)
        elif refs.get(f'{axis}_attributemapping_ref'):
            attributemapping_ref = refs[f'{axis}_attributemapping_ref']
        elif metadata_df is not None:
            name = matrix_name + "_{}_attributes".format(axis)
            attributemapping_ref = self._meta_df_to_attribute_mapping(
                axis_ids, metadata_df, name, workspace_id)
        else:
            attributemapping_ref = self._process_attribute_mapping_sheet(
                file_path, f'{axis}_attribute_mapping', matrix_name, workspace_id)

        if attributemapping_ref:
            attr_data[f'{axis}_attributemapping_ref'] = attributemapping_ref

        # col/row_mappings may not be supplied
        id_mapping = self._process_mapping_sheet(file_path, f'{axis}_mapping')
        if id_mapping:
            attr_data[f'{axis}_mapping'] = id_mapping
        # if no mapping, axis ids must match the attribute mapping
        elif attributemapping_ref:
            am_data = self.dfu.get_objects(
                {'object_refs': [attributemapping_ref]}
            )['data'][0]['data']
            unmatched_ids = set(axis_ids) - set(am_data['instances'].keys())
            if unmatched_ids:
                name = "Column" if axis == 'col' else "Row"
                raise ValueError(f"The following {name} IDs from the uploaded matrix do not match "
                                 f"the supplied {name} attribute mapping: {', '.join(unmatched_ids)}"
                                 f"\nPlease verify the input data or upload an excel file with a"
                                 f"{name} mapping tab.")
            else:
                # just gen the IDs in this matrix
                attr_data[f'{axis}_mapping'] = {x: x for x in axis_ids}

        return attr_data

    @staticmethod
    def _build_header_str(attribute_names):

        header_str = ''
        width = 100.0/len(attribute_names)

        header_str += '<tr class="header">'
        header_str += '<th style="width:{0:.2f}%;">Feature ID</th>'.format(width)

        for attribute_name in attribute_names:
            header_str += '<th style="width:{0:.2f}%;"'.format(width)
            header_str += '>{}</th>'.format(attribute_name)
        header_str += '</tr>'

        return header_str

    def _build_html_str(self, row_mapping, attributemapping_data, row_ids):

        logging.info('Start building html replacement')

        attribute_names = [attributes.get('attribute') for attributes in attributemapping_data.get('attributes')]

        header_str = self._build_header_str(attribute_names)

        table_str = ''

        instances = attributemapping_data.get('instances')

        for feature_id, attribute_id in row_mapping.items():
            if feature_id in row_ids:
                feature_instances = instances.get(attribute_id)

                table_str += '<tr>'
                table_str += '<td>{}</td>'.format(feature_id)

                for feature_instance in feature_instances:
                    table_str += '<td>{}</td>'.format(feature_instance)
                table_str += '</tr>'

        return header_str, table_str

    def _generate_search_html_report(self, header_str, table_str):

        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'search.html')

        shutil.copy2(os.path.join(os.path.dirname(__file__), 'templates', 'kbase_icon.png'),
                     output_directory)
        shutil.copy2(os.path.join(os.path.dirname(__file__), 'templates', 'search_icon.png'),
                     output_directory)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'templates', 'search_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('//HEADER_STR', header_str)
                report_template = report_template.replace('//TABLE_STR', table_str)
                result_file.write(report_template)

        report_shock_id = self.dfu.file_to_shock({'file_path': output_directory,
                                                  'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for Search Matrix App'})

        return html_report

    def _generate_search_report(self, header_str, table_str, workspace_name):
        logging.info('Start creating report')

        output_html_files = self._generate_search_html_report(header_str, table_str)

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 366,
                         'report_object_name': 'kb_matrix_filter_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    @staticmethod
    def _filter_value_data(value_data, remove_ids, dimension):
        """Filters a value matrix based on column or row ids"""
        def _norm_id(_id):
            return _id.replace(" ", "_")

        val_df = pd.DataFrame(value_data['values'], index=value_data['row_ids'],
                              columns=value_data['col_ids'], dtype='object')

        if dimension == 'row':
            filtered_df = val_df.drop(remove_ids, axis=0, errors='ignore')
            filtered_df = filtered_df.drop([_norm_id(x) for x in remove_ids], axis=0, errors='ignore')
        elif dimension == 'col':
            filtered_df = val_df.drop(remove_ids, axis=1, errors='ignore')
            filtered_df = filtered_df.drop([_norm_id(x) for x in remove_ids], axis=1, errors='ignore')
        else:
            raise ValueError('Unexpected dimension: {}'.format(dimension))

        filtered_value_data = {
            "values": filtered_df.values.tolist(),
            "col_ids": list(filtered_df.columns),
            "row_ids": list(filtered_df.index),
        }

        return filtered_value_data

    def _standardize_df(self, df, dimension='col', with_mean=True, with_std=True):

        logging.info("Standardizing matrix data")

        if dimension == 'col':
            df = df.T

        df.fillna(0, inplace=True)

        x_train = df.values

        scaler = preprocessing.StandardScaler(with_mean=with_mean, with_std=with_std).fit(x_train)

        standardized_values = scaler.transform(x_train)

        standardize_df = pd.DataFrame(index=df.index, columns=df.columns, data=standardized_values)

        if dimension == 'col':
            standardize_df = standardize_df.T

        return standardize_df

    def _ratio_trans_df(self, df, method='clr', dimension='col'):

        logging.info("Performaing log ratio transformation matrix data")

        if dimension == 'col':
            df = df.T

        df.fillna(0, inplace=True)

        if method == 'clr':
            ratio_trans = clr(df)
        elif method == 'ilr':
            ratio_trans = ilr(df)
        else:
            raise ValueError('Unexpected ratio transformation method')

        ratio_transformed_df = pd.DataFrame(index=df.index, columns=df.columns, data=ratio_trans)

        if dimension == 'col':
            ratio_transformed_df = ratio_transformed_df.T

        return ratio_transformed_df

    def _remove_all_zero(self, df):

        logging.info("Removing all zero rows")
        row_check = (df != 0).any(axis=1)
        removed_row_ids = list(row_check[row_check == False].index)
        df = df.loc[row_check]

        logging.info("Removing all zero columns")
        col_check = (df != 0).any(axis=0)
        removed_col_ids = list(col_check[col_check == False].index)
        df = df.loc[:, col_check]

        return df, removed_row_ids, removed_col_ids

    def _filtering_matrix(self, df, row_threshold=0, columns_threshold=0):

        logging.info("Removing rows with values all below {}".format(row_threshold))
        row_check = (df > row_threshold).any(axis=1)
        removed_row_ids = list(row_check[row_check == False].index)
        df = df.loc[row_check]

        logging.info("Removing columns with values all below {}".format(row_threshold))
        col_check = (df > columns_threshold).any(axis=0)
        removed_col_ids = list(col_check[col_check == False].index)
        df = df.loc[:, col_check]

        return df, removed_row_ids, removed_col_ids

    def _relative_abundance(self, df, dimension='col'):

        logging.info("Creating relative abundance matrix")

        if dimension == 'col':
            df = df.T

        df.fillna(0, inplace=True)
        values = df.values

        rel_values = list()
        for value in values:
            total = value.sum()
            rel_values.append([v/float(total) for v in value])

        relative_abundance_df = pd.DataFrame(index=df.index, columns=df.columns, data=rel_values)

        if dimension == 'col':
            relative_abundance_df = relative_abundance_df.T

        return relative_abundance_df

    def _create_distance_matrix(self, df, dist_metric='euclidean', dimension='col'):
        '''
        dist_metric: The distance metric to use. Default set to 'euclidean'.
                     The distance function can be
                     ["braycurtis", "canberra", "chebyshev", "cityblock", "correlation", "cosine",
                      "dice", "euclidean", "hamming", "jaccard", "kulsinski", "matching",
                      "rogerstanimoto", "russellrao", "sokalmichener", "sokalsneath",
                      "sqeuclidean", "yule"]
        '''

        # calculate distance matrix
        logging.info('start calculating distance matrix')

        if dimension == 'col':
            df = df.T

        df.fillna(0, inplace=True)
        values = df.values
        labels = df.index.tolist()
        Y = dist.pdist(values, metric=dist_metric)

        dist_matrix = dist.squareform(Y)

        dm = DistanceMatrix(dist_matrix, labels)

        return dm

    def _run_anosim(self, dm, grouping, permutations):

        anosim_res = anosim(dm, grouping, permutations=permutations)

        return dict(anosim_res)

    def _run_permanova(self, dm, grouping, permutations):

        permanova_res = permanova(dm, grouping, permutations=permutations)

        return dict(permanova_res)

    def _run_permdisp(self, dm, grouping, permutations):

        permdisp_res = permdisp(dm, grouping, permutations=permutations)

        return dict(permdisp_res)

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.dfu = DataFileUtil(self.callback_url)
        self.fba_tools = fba_tools(self.callback_url)
        self.report_util = kb_GenericsReport(self.callback_url)
        self.data_util = DataUtil(config)
        self.attr_util = AttributesUtil(config)
        self.sampleservice_util = SampleServiceUtil(config)
        self.matrix_types = [x.split(".")[1].split('-')[0]
                             for x in self.data_util.list_generic_types()]

    def standardize_matrix(self, params):
        """
        standardize a matrix
        """

        input_matrix_ref = params.get('input_matrix_ref')
        workspace_name = params.get('workspace_name')
        new_matrix_name = params.get('new_matrix_name')
        with_mean = params.get('with_mean', 1)
        with_std = params.get('with_std', 1)
        dimension = params.get('dimension', 'col')

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
        input_matrix_info = input_matrix_obj['info']
        input_matrix_name = input_matrix_info[1]
        input_matrix_data = input_matrix_obj['data']

        if not new_matrix_name:
            current_time = time.localtime()
            new_matrix_name = input_matrix_name + time.strftime('_%H_%M_%S_%Y_%m_%d', current_time)

        data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
        df = pd.read_json(data_matrix)
        standardize_df = self._standardize_df(df, dimension=dimension,
                                              with_mean=with_mean, with_std=with_std)

        new_matrix_data = {'row_ids': df.index.tolist(),
                           'col_ids': df.columns.tolist(),
                           'values': standardize_df.values.tolist()}

        input_matrix_data['data'] = new_matrix_data

        logging.info("Saving new standardized matrix object")
        info = self.dfu.save_objects({
            "id": workspace_id,
            "objects": [{
                "type": input_matrix_info[2],
                "data": input_matrix_data,
                "name": new_matrix_name
            }]
        })[0]

        new_matrix_obj_ref = "%s/%s/%s" % (info[6], info[0], info[4])

        objects_created = [{'ref': new_matrix_obj_ref, 'description': 'Standardized Matrix'}]

        report_params = {'message': '',
                         'objects_created': objects_created,
                         'workspace_name': workspace_name,
                         'report_object_name': 'import_matrix_from_biom_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        return {'new_matrix_obj_ref': new_matrix_obj_ref,
                'report_name': output['name'], 'report_ref': output['ref']}

    def perform_variable_stats_matrix(self, params):

        input_matrix_ref = params.get('input_matrix_ref')
        workspace_id = params.get('workspace_id')
        dimension = params.get('dimension', 'col')
        dist_metric = params.get('dist_metric', 'euclidean')
        permutations = params.get('permutations', 0)
        grouping = params.get('grouping')
        perform_anosim = params.get('perform_anosim', True)
        perform_permanova = params.get('perform_permanova', True)
        perform_permdisp = params.get('perform_permdisp', True)

        if not bool(perform_anosim or perform_permanova or perform_permdisp):
            raise ValueError('Please select at least one algorithm to perform')

        input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
        input_matrix_data = input_matrix_obj['data']

        data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
        df = pd.read_json(data_matrix)

        dm = self._create_distance_matrix(df, dist_metric=dist_metric, dimension=dimension)

        if dimension not in ['col', 'row']:
            raise ValueError('Please use "col" or "row" for input dimension')

        am_ref = input_matrix_data.get('{}_attributemapping_ref'.format(dimension))

        if not am_ref:
            raise ValueError('Missing {} attribute mapping from original matrix'.format(dimension))

        am_data = self.dfu.get_objects({'object_refs': [am_ref]})['data'][0]['data']
        attribute_names = [am.get('attribute') for am in am_data.get('attributes')]

        if grouping not in attribute_names:
            raise ValueError('Cannot find {} in {} attribute mapping'.format(grouping, dimension))

        attri_pos = attribute_names.index(grouping)

        instances = am_data.get('instances')

        grouping_names = [instance[attri_pos] for instance in instances.values()]

        anosim_res = None
        if perform_anosim:
            anosim_res = self._run_anosim(dm, grouping_names, permutations)

        permanova_res = None
        if perform_permanova:
            permanova_res = self._run_permanova(dm, grouping_names, permutations)

        permdisp_res = None
        if perform_permdisp:
            permdisp_res = self._run_permdisp(dm, grouping_names, permutations)

        report_output = self._generate_variable_stats_report(workspace_id,
                                                             anosim_res,
                                                             permanova_res,
                                                             permdisp_res)

        return report_output

    def transform_matrix(self, params):
        """
        transform a matrix
        """

        input_matrix_ref = params.get('input_matrix_ref')
        workspace_name = params.get('workspace_name')
        new_matrix_name = params.get('new_matrix_name')

        abundance_filtering_params = params.get('abundance_filtering_params')
        perform_relative_abundance = params.get('perform_relative_abundance', False)
        standardization_params = params.get('standardization_params')
        ratio_transformation_params = params.get('ratio_transformation_params')

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
        input_matrix_info = input_matrix_obj['info']
        input_matrix_name = input_matrix_info[1]
        input_matrix_data = input_matrix_obj['data']

        if not new_matrix_name:
            current_time = time.localtime()
            new_matrix_name = input_matrix_name + time.strftime('_%H_%M_%S_%Y_%m_%d', current_time)

        data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
        df = pd.read_json(data_matrix)
        original_matrix_df = df.copy(deep=True)

        filtered_df = None
        if abundance_filtering_params is not None:
            (df, removed_row_ids, removed_col_ids) = self._filtering_matrix(
                                        df,
                                        row_threshold=abundance_filtering_params.get(
                                                        'abundance_filtering_row_threshold', 0),
                                        columns_threshold=abundance_filtering_params.get(
                                                    'abundance_filtering_columns_threshold', 0))
            filtered_df = df.copy(deep=True)

        relative_abundance_df = None
        if perform_relative_abundance:

            df = self._relative_abundance(df, dimension='col')
            relative_abundance_df = df.copy(deep=True)

        standardize_df = None
        if standardization_params is not None:
            df = self._standardize_df(df,
                                      dimension=standardization_params.get(
                                                            'standardization_dimension', 'col'),
                                      with_mean=standardization_params.get(
                                                            'standardization_with_mean', True),
                                      with_std=standardization_params.get(
                                                            'standardization_with_std', True))
            standardize_df = df.copy(deep=True)

        ratio_transformed_df = None
        if ratio_transformation_params is not None:
            df = self._ratio_trans_df(df,
                                      method=ratio_transformation_params.get(
                                                        'ratio_transformation_method', 'clr'),
                                      dimension=ratio_transformation_params.get(
                                                        'ratio_transformation_dimension', 'col'))
            ratio_transformed_df = df.copy(deep=True)

        new_matrix_data = {'row_ids': df.index.tolist(),
                           'col_ids': df.columns.tolist(),
                           'values': df.values.tolist()}

        input_matrix_data['data'] = new_matrix_data

        logging.info("Saving new transformed matrix object")
        info = self.dfu.save_objects({
            "id": workspace_id,
            "objects": [{
                "type": input_matrix_info[2],
                "data": input_matrix_data,
                "name": new_matrix_name
            }]
        })[0]

        new_matrix_obj_ref = "%s/%s/%s" % (info[6], info[0], info[4])
        returnVal = {'new_matrix_obj_ref': new_matrix_obj_ref}

        report_output = self._generate_transform_report(new_matrix_obj_ref, workspace_name,
                                                        original_matrix_df,
                                                        filtered_df, relative_abundance_df,
                                                        standardize_df, ratio_transformed_df)

        returnVal.update(report_output)

        return returnVal

    def filter_matrix(self, params):
        """
        filter_matrix: create sub-matrix based on input feature_ids

        arguments:
        matrix_obj_ref: object reference of a matrix
        workspace_name: workspace name
        feature_ids: string of feature ids that result matrix contains
        filtered_matrix_name: name of newly created filtered matrix object
        """

        matrix_obj_ref = params.get('matrix_obj_ref')
        workspace_name = params.get('workspace_name')
        remove_ids = params.get('remove_ids')
        dimension = params.get('dimension')
        filtered_matrix_name = params.get('filtered_matrix_name')

        matrix_source = self.dfu.get_objects(
            {"object_refs": [matrix_obj_ref]})['data'][0]
        matrix_info = matrix_source.get('info')
        matrix_data = matrix_source.get('data')

        matrix_type = self._find_between(matrix_info[2], '\.', '\-')

        value_data = matrix_data.get('data')
        remove_ids = [x.strip() for x in remove_ids.split(',')]
        filtered_value_data = self._filter_value_data(value_data, remove_ids, dimension)

        # if the matrix has changed shape, update the mappings
        if len(filtered_value_data['row_ids']) < len(matrix_data['data']['row_ids']):
            if matrix_data.get('row_mapping'):
                matrix_data['row_mapping'] = {k: matrix_data['row_mapping'][k]
                                              for k in filtered_value_data['row_ids']}
            if matrix_data.get('feature_mapping'):
                matrix_data['feature_mapping'] = {k: matrix_data['feature_mapping'][k]
                                                  for k in filtered_value_data['row_ids']}

        if len(filtered_value_data['col_ids']) < len(matrix_data['data']['col_ids']):
            if matrix_data.get('col_mapping'):
                matrix_data['col_mapping'] = {k: matrix_data['col_mapping'][k]
                                              for k in filtered_value_data['col_ids']}
        matrix_data['data'] = filtered_value_data

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        filtered_matrix_obj_ref = self.data_util.save_object({
                                                'obj_type': 'KBaseMatrices.{}'.format(matrix_type),
                                                'obj_name': filtered_matrix_name,
                                                'data': matrix_data,
                                                'workspace_name': workspace_id})['obj_ref']

        returnVal = {'matrix_obj_refs': [filtered_matrix_obj_ref]}

        report_output = self._generate_report(filtered_matrix_obj_ref, workspace_name,
                                              data=filtered_value_data)

        returnVal.update(report_output)

        return returnVal

    def search_matrix(self, params):
        """
        search_matrix: generate a HTML report that allows users to select feature ids

        arguments:
        matrix_obj_ref: object reference of a matrix
        workspace_name: workspace name
        """

        matrix_obj_ref = params.get('matrix_obj_ref')
        workspace_name = params.get('workspace_name')

        matrix_source = self.dfu.get_objects(
            {"object_refs": [matrix_obj_ref]})['data'][0]
        matrix_data = matrix_source.get('data')

        row_mapping = matrix_data.get('row_mapping')
        row_attributemapping_ref = matrix_data.get('row_attributemapping_ref')

        row_ids = matrix_data['data']['row_ids']

        if not (row_mapping and row_attributemapping_ref):
            raise ValueError('Matrix obejct is missing either row_mapping or row_attributemapping_ref')

        attributemapping_data = self.dfu.get_objects(
                                    {"object_refs": [row_attributemapping_ref]})['data'][0]['data']

        header_str, table_str = self._build_html_str(row_mapping, attributemapping_data, row_ids)

        returnVal = self._generate_search_report(header_str, table_str, workspace_name)

        return returnVal

    def import_matrix_from_excel(self, params):
        """
        import_matrix_from_excel: import matrix object from excel

        arguments:
        obj_type: one of ExpressionMatrix, FitnessMatrix, DifferentialExpressionMatrix
        matrix_name: matrix object name
        workspace_name: workspace name matrix object to be saved to
        input_shock_id: file shock id
        or
        input_file_path: absolute file path
        or
        input_staging_file_path: staging area file path

        optional arguments:
        col_attributemapping_ref: column AttributeMapping reference
        row_attributemapping_ref: row AttributeMapping reference
        genome_ref: genome reference
        matrix_obj_ref: Matrix reference
        """

        (obj_type, file_path, workspace_name,
         matrix_name, refs, scale) = self._validate_import_matrix_from_excel_params(params)

        if not isinstance(workspace_name, int):
            workspace_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            workspace_id = workspace_name

        if obj_type in ['ChemicalAbundanceMatrix', 'MetaboliteMatrix']:
            data = self._file_to_chem_abun_data(file_path, refs, matrix_name, workspace_id)
        else:
            data = self._file_to_data(file_path, refs, matrix_name, workspace_id)

        data['scale'] = scale
        for key in ['description', 'unit', 'type']:
            if params.get(key):
                data[key] = params[key]

        new_row_attr_ref = None
        if not params.get('row_attributemapping_ref'):
            new_row_attr_ref = data.get('row_attributemapping_ref')

        new_col_attr_ref = None
        if not params.get('col_attributemapping_ref'):
            new_col_attr_ref = data.get('col_attributemapping_ref')

        matrix_obj_ref = self.data_util.save_object({
                                                'obj_type': 'KBaseMatrices.{}'.format(obj_type),
                                                'obj_name': matrix_name,
                                                'data': data,
                                                'workspace_name': workspace_id})['obj_ref']

        try:
            logging.info('Start trying to look up ModelSeed ID')
            if obj_type in ['ChemicalAbundanceMatrix', 'MetaboliteMatrix']:
                ret = self.fba_tools.lookup_modelseed_ids(
                                                {'workspace': workspace_name,
                                                 'chemical_abundance_matrix_id': matrix_name,
                                                 'chemical_abundance_matrix_out_id': matrix_name})

                matrix_obj_ref = ret.get('new_chemical_abundance_matrix_ref')

                matrix_data = self.dfu.get_objects(
                                    {"object_refs": [matrix_obj_ref]})['data'][0]['data']

                if not params.get('row_attributemapping_ref'):
                    new_row_attr_ref = matrix_data.get('row_attributemapping_ref')

                if not params.get('col_attributemapping_ref'):
                    new_col_attr_ref = matrix_data.get('col_attributemapping_ref')
        except Exception:
            logging.info('Failed looking up ModelSeed ID')
            logging.warning('failed to run run_model_characterization')
            logging.warning(traceback.format_exc())
            logging.warning(sys.exc_info()[2])

        returnVal = {'matrix_obj_ref': matrix_obj_ref}

        report_output = self._generate_report(matrix_obj_ref, workspace_name,
                                              new_row_attr_ref=new_row_attr_ref,
                                              new_col_attr_ref=new_col_attr_ref,
                                              data=data['data'])

        returnVal.update(report_output)

        return returnVal

    def export_matrix(self, params):
        """
        export_matrix: univeral downloader for matrix data object

        arguments:
        obj_ref: generics object reference

        optional arguments:
        generics_module: select the generics data to be retrieved from
                        e.g. for an given data type like below:
                        typedef structure {
                          FloatMatrix2D data;
                          condition_set_ref condition_set_ref;
                        } SomeGenericsMatrix;
                        and only data is needed
                        generics_module should be
                        {'data': 'FloatMatrix2D'}
        """
        logging.info('Start exporting matrix')

        if 'input_ref' in params:
            params['obj_ref'] = params.pop('input_ref')

        obj_source = self.dfu.get_objects(
            {"object_refs": [params.get('obj_ref')]})['data'][0]
        obj_data = obj_source.get('data')

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)
        file_path = os.path.join(result_directory, '{}.xlsx'.format(obj_source.get('info')[1]))

        data_matrix = self.data_util.fetch_data(params).get('data_matrix')
        df = pd.read_json(data_matrix)

        df.to_excel(file_path, sheet_name='data')

        if obj_data.get('col_mapping'):
            self._write_mapping_sheet(file_path, 'col_mapping',
                                      obj_data.get('col_mapping'), ['col_name', 'instance_name'])
            obj_data.pop('col_mapping')

        if obj_data.get('row_mapping'):
            self._write_mapping_sheet(file_path, 'row_mapping',
                                      obj_data.get('row_mapping'), ['row_name', 'instance_name'])
            obj_data.pop('row_mapping')

        try:
            obj_data.pop('data')
        except KeyError:
            logging.warning('Missing key [data]')

        obj_data.update(obj_data.get('attributes', {}))  # flatten for printing
        self._write_mapping_sheet(file_path, 'metadata', obj_data, ['name', 'value'])

        shock_id = self._upload_to_shock(file_path)

        return {'shock_id': shock_id}
