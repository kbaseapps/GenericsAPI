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
import numpy as np
from openpyxl import load_workbook
from xlrd.biffh import XLRDError
from sklearn import preprocessing
from skbio.stats.composition import ilr, clr
from skbio import DistanceMatrix
from skbio.stats.distance import anosim, permanova, permdisp, pwmantel
import scipy.spatial.distance as dist
import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter
import plotly.graph_objects as go
from plotly.offline import plot
import plotly.express as px

from installed_clients.DataFileUtilClient import DataFileUtil
from GenericsAPI.Utils.AttributeUtils import AttributesUtil
from GenericsAPI.Utils.SampleServiceUtil import SampleServiceUtil
from GenericsAPI.Utils.DataUtil import DataUtil
import GenericsAPI.Utils.MatrixValidation as vd
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
            tab_content += '\n<iframe height="900px" width="100%" '
            tab_content += 'src="{}" '.format(index_page)
            tab_content += 'style="border:none;"></iframe>'
            tab_content += '\n</div>\n'
        else:
            tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
            tab_content += '''\n<p style="color:red;" >'''
            tab_content += '''Matrix is too large to be displayed.</p>\n'''
            tab_content += '\n</div>\n'

        return tab_content

    def _generate_simper_tab_content(self, res, viewer_name):
        tab_content = ''

        tab_content += '''\n<div id="{}" class="tabcontent">\n'''.format(viewer_name)
        html = '''<pre class="tab">''' + str(res).replace("\n", "<br>") + "</pre>"
        tab_content += html
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

    def _generate_mantel_test_visualization_content(self, pwmantel_res):
        tab_def_content = ''
        tab_content = ''

        viewer_name = 'pwmantel_res'

        tab_def_content += '''\n<div class="tab">\n'''
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += ''' id="defaultOpen"'''
        tab_def_content += '''>Mantel Test</button>\n'''
        tab_def_content += '\n</div>\n'

        tab_content += '''\n<div id="{}" class="tabcontent">\n'''.format(viewer_name)
        tab_content += '''<table>\n'''

        # add table headers
        tab_content += '''<tr>\n'''
        tab_content += '''<th>Distance Matrix 1</th>\n'''
        tab_content += '''<th>Distance Matrix 2</th>\n'''
        for col in pwmantel_res.columns:
            tab_content += '''<th>{}</th>\n'''.format(col)
        tab_content += '''</tr>\n'''

        # add table contents
        for idx, values in enumerate(pwmantel_res.values):
            tab_content += '''<tr>\n'''
            tab_content += '''<td>{}</td>\n'''.format(pwmantel_res.index[idx][0])
            tab_content += '''<td>{}</td>\n'''.format(pwmantel_res.index[idx][1])
            values[0] = round(values[0], 4)
            for value in values:
                tab_content += '''<td>{}</td>\n'''.format(value)
            tab_content += '''</tr>\n'''

        tab_content += '''</table>\n'''
        tab_content += '\n</div>\n'

        return tab_def_content + tab_content

    def _generate_simper_plot(self, species_stats, grouping_names):

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating plotly simper plot in {}'.format(output_directory))
        self._mkdir_p(output_directory)
        simper_plot_path = os.path.join(output_directory, 'SimperPlot.html')

        species = list(species_stats.keys())

        plot_data = list()
        for grouping_name in set(grouping_names):
            y_values = list()
            y_error = list()
            for species_name in species:
                species_data = species_stats[species_name]
                y_values.append(species_data[grouping_name][0])
                y_error.append(species_data[grouping_name][1])

            plot_data.append(go.Bar(name=str(grouping_name), x=species, y=y_values,
                                    error_y=dict(type='data', array=y_error)))

        fig = go.Figure(data=plot_data)
        fig.update_layout(barmode='group',
                          xaxis=dict(title='species'),
                          yaxis=dict(title='average abundance count'))
        plot(fig, filename=simper_plot_path)

        return simper_plot_path

    def _generate_simper_plot_content(self, viewer_name, species_stats, grouping_names,
                                      output_directory):

        simper_plot_path = self._generate_simper_plot(species_stats, grouping_names)

        simper_plot_name = 'SimperPlot.html'
        shutil.copy2(simper_plot_path,
                     os.path.join(output_directory, simper_plot_name))

        tab_content = ''

        tab_content += '''\n<div id="{}" class="tabcontent">\n'''.format(viewer_name)

        tab_content += '<iframe height="500px" width="100%" '
        tab_content += 'src="{}" '.format(simper_plot_name)
        tab_content += 'style="border:none;"></iframe>\n<p></p>\n'

        tab_content += '\n</div>\n'

        return tab_content

    def _generate_simper_visualization_content(self, simper_ret, simper_sum,
                                               species_stats, grouping_names, output_directory):
        tab_def_content = ''
        tab_content = ''

        viewer_name = 'simper_plot'
        tab_def_content += '''\n<div class="tab">\n'''
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += ''' id="defaultOpen"'''
        tab_def_content += '''>Most Influential Species Bar Plot</button>\n'''

        tab_content += self._generate_simper_plot_content(viewer_name, species_stats,
                                                          grouping_names, output_directory)

        viewer_name = 'simper_ret'
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += '''>Most Influential Species Info</button>\n'''

        tab_content += self._generate_simper_tab_content(simper_ret, viewer_name)

        viewer_name = 'simper_sum'
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += '''>Similarity Percentage Summary</button>\n'''

        tab_content += self._generate_simper_tab_content(simper_sum, viewer_name)

        tab_def_content += '\n</div>\n'
        return tab_def_content + tab_content

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
                # first_tab_token = True
                tab_def_content += '''\n<div class="tab">\n'''
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += ''' id="defaultOpen"'''
                tab_def_content += '''>Homogeneity Multivariate Analysis of Variance</button>\n'''

            tab_content += self._generate_variable_stat_tab_content(permdisp_res, viewer_name)

        tab_def_content += '\n</div>\n'
        return tab_def_content + tab_content

    def _generate_rarefy_visualization_content(self, output_directory,
                                               rarefied_matrix_dir, rarecurve_image,
                                               obs_vs_rare_image, random_rare_df):
        tab_def_content = ''
        tab_content = ''

        row_data_summary = random_rare_df.T.describe().to_string()
        col_data_summary = random_rare_df.describe().to_string()

        tab_def_content = ''
        tab_content = ''

        viewer_name = 'data_summary'
        tab_def_content += '''\n<div class="tab">\n'''
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += ''' id="defaultOpen"'''
        tab_def_content += '''>Rarefied Matrix Statistics</button>\n'''

        tab_content += '''\n<div id="{}" class="tabcontent" style="overflow:auto">'''.format(viewer_name)
        tab_content += '''\n<h5>Rarefied Matrix Size: {} x {}</h5>'''.format(
                                                                    len(random_rare_df.index),
                                                                    len(random_rare_df.columns))
        tab_content += '''\n<h5>Row Aggregating Statistics</h5>'''
        html = '''\n<pre class="tab">''' + str(row_data_summary).replace("\n", "<br>") + "</pre>"
        tab_content += html
        tab_content += '''\n<br>'''
        tab_content += '''\n<hr style="height:2px;border-width:0;color:gray;background-color:gray">'''
        tab_content += '''\n<br>'''
        tab_content += '''\n<h5>Column Aggregating Statistics</h5>'''
        html = '''\n<pre class="tab">''' + str(col_data_summary).replace("\n", "<br>") + "</pre>"
        tab_content += html
        tab_content += '\n</div>\n'

        if False and len(random_rare_df.columns) <= 200:
            viewer_name = 'MatrixLinearPlotViewer'
            tab_def_content += '''\n<button class="tablinks" '''
            tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
            tab_def_content += '''>Matrix Linear Plot</button>\n'''

            linear_plot_page = self._generate_linear_plot(random_rare_df, output_directory)

            tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
            tab_content += '\n<iframe height="900px" width="100%" '
            tab_content += 'src="{}" '.format(linear_plot_page)
            tab_content += 'style="border:none;"></iframe>'
            tab_content += '\n</div>\n'

        viewer_name = 'RarefiedMatrixViewer'
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += '''>Rarefied Matrix Heatmap</button>\n'''
        rarefied_matrix_report_files = os.listdir(rarefied_matrix_dir)
        rarefied_matrix_index_page = None
        for rarefied_matrix_report_file in rarefied_matrix_report_files:
            if rarefied_matrix_report_file.endswith('.html'):
                rarefied_matrix_index_page = rarefied_matrix_report_file

            shutil.copy2(os.path.join(rarefied_matrix_dir, rarefied_matrix_report_file),
                         output_directory)
        tab_content += self._generate_tab_content(rarefied_matrix_index_page, viewer_name)

        rarecurve_image_name = os.path.basename(rarecurve_image)
        shutil.copy2(rarecurve_image,
                     os.path.join(output_directory, rarecurve_image_name))

        obs_vs_rare_image_name = os.path.basename(obs_vs_rare_image)
        shutil.copy2(obs_vs_rare_image,
                     os.path.join(output_directory, obs_vs_rare_image_name))

        viewer_name = 'RarecurvePlot'
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += '''>Rarecurve Plot</button>\n'''

        tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
        tab_content += '''\n<img src="{}" '''.format(rarecurve_image_name)
        tab_content += '''alt="rarecurve" width="600" height="600">\n'''
        tab_content += '''<br>\n<br>\n'''
        tab_content += '''\n<img src="{}" '''.format(obs_vs_rare_image_name)
        tab_content += '''alt="rarecurve" width="600" height="600">\n'''
        tab_content += '\n</div>\n'

        tab_def_content += '\n</div>\n'
        return tab_def_content + tab_content

    def _generate_trans_visualization_content(self, output_directory,
                                              operations, heatmap_html_dir_l,
                                              transformed_matrix_df, variable_specific):
        row_data_summary = transformed_matrix_df.T.describe().to_string()
        col_data_summary = transformed_matrix_df.describe().to_string()

        tab_def_content = ''
        tab_content = ''

        op_2_name = {
            'abundance_filtering': 'Filtered',
            'standardization': 'Standardized',
            'ratio_transformation': 'Log Ratio Transformed',
            'relative_abundance': 'Relative Abundance',
            'logit': 'Logit',
            'sqrt': 'Square Root',
            'log': 'Log',
        }

        ## Start tabs ##
        tab_def_content += '''\n<div class="tab">\n'''

        ## Operations tabs ##
        for i, (op, heatmap_html_dir) in enumerate(zip(operations, heatmap_html_dir_l)):
            viewer_name = 'op%s_%s' % (i, op)
            tab_def_content += '''\n<button class="tablinks" '''
            tab_def_content += '''onclick="openTab(event, '%s')"''' % viewer_name
            tab_def_content += '''>%d. %s</button>\n''' % (i+1, op_2_name[op])

            flnms = os.listdir(heatmap_html_dir)
            heatmap_html_flnm = None
            for flnm in flnms:
                if flnm.endswith('.html'):
                    heatmap_html_flnm = flnm

                shutil.copy2(os.path.join(heatmap_html_dir, flnm), output_directory)
            tab_content += self._generate_tab_content(heatmap_html_flnm, viewer_name)

        ## Transformed matrix statistics tab ##
        viewer_name = 'data_summary'
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += ''' id="defaultOpen"'''
        if variable_specific:
            tab_def_content += '''>Transformed Selected Variables Statistics</button>\n'''
        else:
            tab_def_content += '''>Transformed Matrix Statistics</button>\n'''
        tab_content += '''\n<div id="{}" class="tabcontent" style="overflow:auto">'''.format(
                                                                                    viewer_name)
        if variable_specific:
            tab_content += '''\n<h5>Transformed Selected Variables Size: {} x {}</h5>'''.format(
                                                                len(transformed_matrix_df.index),
                                                                len(transformed_matrix_df.columns))
        else:
            tab_content += '''\n<h5>Transformed Matrix Size: {} x {}</h5>'''.format(
                                                                len(transformed_matrix_df.index),
                                                                len(transformed_matrix_df.columns))
        tab_content += '''\n<h5>Row Aggregating Statistics</h5>'''
        html = '''\n<pre class="tab">''' + str(row_data_summary).replace("\n", "<br>") + "</pre>"
        tab_content += html
        tab_content += '''\n<br>'''
        tab_content += '''\n<hr style="height:2px;border-width:0;color:gray;background-color:gray">'''
        tab_content += '''\n<br>'''
        tab_content += '''\n<h5>Column Aggregating Statistics</h5>'''
        html = '''\n<pre class="tab">''' + str(col_data_summary).replace("\n", "<br>") + "</pre>"
        tab_content += html
        tab_content += '\n</div>\n'

        tab_def_content += '\n</div>\n'
        return tab_def_content + tab_content

    def _generate_linear_plot(self, data_df, output_directory, row_name='abundance',
                              top_percent=100):
        linear_plot_path = 'linear_plot.html'

        sum_order = data_df.sum(axis=1).sort_values(ascending=False).index
        data_df = data_df.reindex(sum_order)

        top_index = data_df.index[:int(data_df.index.size * top_percent / 100)]
        data_df = data_df.loc[top_index]

        links = data_df.stack().reset_index()

        col_names = links.columns
        links.rename(columns={col_names[0]: row_name,
                              col_names[1]: 'samples',
                              col_names[2]: 'value'},
                     inplace=True)
        fig = px.line(links, x=row_name, y='value', color='samples')

        plot(fig, filename=os.path.join(output_directory, linear_plot_path))

        return linear_plot_path

    def _generate_visualization_content(self, output_directory, heatmap_dir, data_df,
                                        top_heatmap_dir, top_percent):

        row_data_summary = data_df.T.describe().to_string()
        col_data_summary = data_df.describe().to_string()

        tab_def_content = ''
        tab_content = ''

        viewer_name = 'data_summary'
        tab_def_content += '''\n<div class="tab">\n'''
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += ''' id="defaultOpen"'''
        tab_def_content += '''>Matrix Statistics</button>\n'''

        tab_content += '''\n<div id="{}" class="tabcontent" style="overflow:auto">'''.format(
                                                                                    viewer_name)
        tab_content += '''\n<h5>Matrix Size: {} x {}</h5>'''.format(len(data_df.index),
                                                                    len(data_df.columns))
        tab_content += '''\n<h5>Row Aggregating Statistics</h5>'''
        html = '''\n<pre class="tab">''' + str(row_data_summary).replace("\n", "<br>") + "</pre>"
        tab_content += html
        tab_content += '''\n<br>'''
        tab_content += '''\n<hr style="height:2px;border-width:0;color:gray;background-color:gray">'''
        tab_content += '''\n<br>'''
        tab_content += '''\n<h5>Column Aggregating Statistics</h5>'''
        html = '''\n<pre class="tab">''' + str(col_data_summary).replace("\n", "<br>") + "</pre>"
        tab_content += html
        tab_content += '\n</div>\n'

        if top_heatmap_dir:
            viewer_name = 'TopHeatmapViewer'
            tab_def_content += '''\n<button class="tablinks" '''
            tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
            tab_def_content += '''>Top {} Percent Heatmap</button>\n'''.format(top_percent)

            heatmap_report_files = os.listdir(top_heatmap_dir)

            heatmap_index_page = None
            for heatmap_report_file in heatmap_report_files:
                if heatmap_report_file.endswith('.html'):
                    heatmap_index_page = heatmap_report_file

                shutil.copy2(os.path.join(top_heatmap_dir, heatmap_report_file),
                             output_directory)

            if heatmap_index_page:
                tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
                msg = 'Top {} percent of matrix sorted by sum of abundance values.'.format(top_percent)
                tab_content += '''<p style="color:red;" >{}</p>'''.format(msg)

                tab_content += '\n<iframe height="900px" width="100%" '
                tab_content += 'src="{}" '.format(heatmap_index_page)
                tab_content += 'style="border:none;"></iframe>'
                tab_content += '\n</div>\n'
            else:
                tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
                tab_content += '''\n<p style="color:red;" >'''
                tab_content += '''Heatmap is too large to be displayed.</p>\n'''
                tab_content += '\n</div>\n'

        if False and len(data_df.columns) <= 200:
            if top_heatmap_dir:
                viewer_name = 'MatrixLinearPlotViewer'
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += '''>Top {} Percent Linear Plot</button>\n'''.format(top_percent)

                linear_plot_page = self._generate_linear_plot(data_df, output_directory,
                                                              top_percent=top_percent)

                tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
                msg = 'Top {} percent of matrix sorted by sum of abundance values.'.format(top_percent)
                tab_content += '''<p style="color:red;" >{}</p>'''.format(msg)
                tab_content += '\n<iframe height="900px" width="100%" '
                tab_content += 'src="{}" '.format(linear_plot_page)
                tab_content += 'style="border:none;"></iframe>'
                tab_content += '\n</div>\n'
            else:
                viewer_name = 'MatrixLinearPlotViewer'
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += '''>Matrix Linear Plot</button>\n'''

                linear_plot_page = self._generate_linear_plot(data_df, output_directory)

                tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
                tab_content += '\n<iframe height="900px" width="100%" '
                tab_content += 'src="{}" '.format(linear_plot_page)
                tab_content += 'style="border:none;"></iframe>'
                tab_content += '\n</div>\n'

        viewer_name = 'MatrixHeatmapViewer'
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += '''>Matrix Heatmap</button>\n'''

        heatmap_report_files = os.listdir(heatmap_dir)

        heatmap_index_page = None
        for heatmap_report_file in heatmap_report_files:
            if heatmap_report_file.endswith('.html'):
                heatmap_index_page = heatmap_report_file

            shutil.copy2(os.path.join(heatmap_dir, heatmap_report_file),
                         output_directory)

        if heatmap_index_page:
            tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
            tab_content += '\n<iframe height="900px" width="100%" '
            tab_content += 'src="{}" '.format(heatmap_index_page)
            tab_content += 'style="border:none;"></iframe>'
            tab_content += '\n</div>\n'
        else:
            tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
            tab_content += '''\n<p style="color:red;" >'''
            tab_content += '''Heatmap is too large to be displayed.</p>\n'''
            tab_content += '\n</div>\n'

        tab_def_content += '\n</div>\n'
        return tab_def_content + tab_content

    def _generate_mantel_test_html_report(self, pwmantel_res):
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'mantel_test_viewer_report.html')

        visualization_content = self._generate_mantel_test_visualization_content(pwmantel_res)

        table_style_content = '''
                                table {
                                  font-family: arial, sans-serif;
                                  border-collapse: collapse;
                                  width: 100%;
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
                            'description': 'HTML summary report for Mantel Test App'
                            })
        return html_report

    def _generate_simper_html_report(self, simper_ret, simper_sum, species_stats, grouping_names):
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'simper_viewer_report.html')

        visualization_content = self._generate_simper_visualization_content(simper_ret,
                                                                            simper_sum,
                                                                            species_stats,
                                                                            grouping_names,
                                                                            output_directory)

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
                            'description': 'HTML summary report for Simper App'
                            })
        return html_report

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
                            'description': 'HTML summary report for Variable Stats App'
                            })
        return html_report

    def _generate_rarefy_html_report(self, rarefied_matrix_dir,
                                     rarecurve_image, obs_vs_rare_image, random_rare_df):
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'rarefy_matrix_viewer_report.html')

        visualization_content = self._generate_rarefy_visualization_content(
                                                                    output_directory,
                                                                    rarefied_matrix_dir,
                                                                    rarecurve_image,
                                                                    obs_vs_rare_image,
                                                                    random_rare_df)

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
                            'description': 'HTML summary report for Transform Matrix App'
                            })
        return html_report

    def _generate_transform_html_report(self, operations, heatmap_html_dir_l,
                                        transformed_matrix_df, variable_specific):
        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'transform_matrix_viewer_report.html')

        visualization_content = self._generate_trans_visualization_content(
                                                                    output_directory,
                                                                    operations,
                                                                    heatmap_html_dir_l,
                                                                    transformed_matrix_df,
                                                                    variable_specific)

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
                            'description': 'HTML summary report for Transform Matrix App'
                            })
        return html_report

    def _generate_heatmap_html_report(self, data):

        logging.info('Start generating heatmap report page')

        data_df = pd.DataFrame(data['values'], index=data['row_ids'], columns=data['col_ids'])
        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)
        tsv_file_path = os.path.join(result_directory, 'heatmap_data_{}.tsv'.format(
                                                                            str(uuid.uuid4())))
        data_df.to_csv(tsv_file_path)
        heatmap_dir = self.report_util.build_heatmap_html({
                                                    'tsv_file_path': tsv_file_path,
                                                    'cluster_data': True})['html_dir']

        top_heatmap_dir = None
        top_percent = 100
        if len(data_df.index) > 500:
            display_count = 200  # roughly count for display items
            top_percent = min(int(display_count / len(data_df.index) * 100), 100)
            top_percent = max(top_percent, 1)
            top_heatmap_dir = self.report_util.build_heatmap_html({
                                                        'tsv_file_path': tsv_file_path,
                                                        'sort_by_sum': True,
                                                        'top_percent': top_percent})['html_dir']

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating html report in {}'.format(output_directory))

        html_report = list()

        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'matrix_viewer_report.html')

        visualization_content = self._generate_visualization_content(output_directory,
                                                                     heatmap_dir,
                                                                     data_df,
                                                                     top_heatmap_dir,
                                                                     top_percent)

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
                            'description': 'HTML summary report for Import Matrix App'
                            })
        return html_report

    def _generate_rarefy_report(self, new_matrix_obj_ref, workspace_id,
                                random_rare_df, rarecurve_image, obs_vs_rare_image,
                                warnings):

        objects_created = [{'ref': new_matrix_obj_ref, 'description': 'Randomly Rarefied Matrix'}]

        data_tsv_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(data_tsv_directory)
        logging.info('Start generating matrix tsv files in {}'.format(data_tsv_directory))
        rarefied_matrix_tsv_path = os.path.join(data_tsv_directory,
                                                'rarefied_matrix_{}.tsv'.format(
                                                                                str(uuid.uuid4())))
        random_rare_df.to_csv(rarefied_matrix_tsv_path)
        rarefied_matrix_dir = self.report_util.build_heatmap_html({
                                            'tsv_file_path': rarefied_matrix_tsv_path,
                                            'cluster_data': True})['html_dir']

        output_html_files = self._generate_rarefy_html_report(rarefied_matrix_dir,
                                                              rarecurve_image,
                                                              obs_vs_rare_image,
                                                              random_rare_df)

        report_params = {'message': '',
                         'objects_created': objects_created,
                         'workspace_id': workspace_id,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 660,
                         'report_object_name': 'rarefy_matrix_' + str(uuid.uuid4()),
                         'warnings': warnings}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_transform_report(self, new_matrix_obj_ref, workspace_id,
                                   operations, df_results, variable_specific=False):
        objects_created = [{'ref': new_matrix_obj_ref, 'description': 'Transformed Matrix'}]

        data_tsv_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(data_tsv_directory)

        heatmap_html_dir_l = []
        for i, (op, df) in enumerate(zip(operations, df_results)):
            tsv_path = os.path.join(data_tsv_directory, 'op%d_%s.tsv' %(i, op))
            df.to_csv(tsv_path)
            heatmap_html_dir = self.report_util.build_heatmap_html({
                'tsv_file_path': tsv_path,
                'cluster_data': True
            })['html_dir']
            heatmap_html_dir_l.append(heatmap_html_dir)

        output_html_files = self._generate_transform_html_report(operations, heatmap_html_dir_l,
                                                                 df_results[-1],
                                                                 variable_specific)

        report_params = {'message': '',
                         'objects_created': objects_created,
                         'workspace_id': workspace_id,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 660,
                         'report_object_name': 'transform_matrix_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_mantel_test_report(self, workspace_id, pwmantel_res):
        output_html_files = self._generate_mantel_test_html_report(pwmantel_res)

        report_params = {'message': '',
                         'workspace_id': workspace_id,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 300,
                         'report_object_name': 'mantel_test_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_simper_report(self, workspace_id, simper_ret, simper_sum,
                                species_stats, grouping_names):
        output_html_files = self._generate_simper_html_report(simper_ret, simper_sum,
                                                              species_stats, grouping_names)

        report_params = {'message': '',
                         'workspace_id': workspace_id,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 450,
                         'report_object_name': 'simper_' + str(uuid.uuid4())}

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
                         'html_window_height': 450,
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
            output_html_files = self._generate_heatmap_html_report(data)

            report_params = {'message': '',
                             'objects_created': objects_created,
                             'workspace_name': workspace_name,
                             'html_links': output_html_files,
                             'direct_html_link_index': 0,
                             'html_window_height': 660,
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
            df = pd.read_excel(file_path, sheet_name=sheet_name, index_col=0)
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

    @staticmethod
    def _check_df_col_inclusive(df, col_name, valid_values):
        # check if given column contains all values in valid_values
        unmatched_type = set(df[col_name]) - valid_values
        if unmatched_type:
            err_msg = 'Found unsupported {}: {}\n'.format(' '.join(col_name.split('_')),
                                                          unmatched_type)
            err_msg += 'Please use one of {} as {}'.format(valid_values,
                                                           ' '.join(col_name.split('_')))
            raise ValueError(err_msg)

    @staticmethod
    def _check_df_col_non_empty(df, col_name):
        # check if any column cell is empty(nan)
        if df[col_name].isna().any():
            empty_idx = list(df.loc[df[col_name].isna()].index)
            raise ValueError('Missing [{}] value for index: {}'.format(col_name, empty_idx))

    @staticmethod
    def _check_chem_ids(df):
        # check chemical abundance has at least one of database id
        id_fields = {'mass', 'formula', 'inchikey', 'inchi', 'smiles', 'compound_name',
                     'kegg', 'chembi', 'modelseed'}

        common_ids = list(df.columns & id_fields)

        ids_df = df.loc[:, common_ids]
        missing_ids_idx = list(ids_df.loc[ids_df.isnull().all(axis=1)].index)

        if missing_ids_idx:
            err_msg = 'Missing compound identification for {}\n'.format(missing_ids_idx)
            err_msg += 'Please provide at least one of {}'.format(id_fields)
            raise ValueError(err_msg)

    def _check_chem_abun_metadata(self, metadata_df):
        logging.info('Start checking metadata fields for Chemical Abundance Matrix')

        metadata_df.replace(r'^\s+$', np.nan, regex=True, inplace=True)

        self._check_chem_ids(metadata_df)

        str_cols = ['chemical_type', 'measurement_type', 'units', 'unit_medium']
        for str_col in str_cols:
            metadata_df[str_col] = metadata_df[str_col].apply(lambda s: s.lower()
                                                              if type(s) == str else s)

        valid_chem_types = {'specific', 'aggregate'}
        self._check_df_col_inclusive(metadata_df, 'chemical_type', valid_chem_types)

        specific_abun = metadata_df.loc[metadata_df['chemical_type'] == 'specific']
        aggregate_abun = metadata_df.loc[metadata_df['chemical_type'] == 'aggregate']

        if not specific_abun.index.empty:
            logging.info('Start examing specific chemical abundances')

            valid_measurement_types = {'unknown', 'fticr', 'orbitrap', 'quadrapole'}
            self._check_df_col_inclusive(specific_abun, 'measurement_type', valid_measurement_types)

            valid_unit_medium = {'soil', 'solvent', 'water'}
            self._check_df_col_inclusive(specific_abun, 'unit_medium', valid_unit_medium)

            self._check_df_col_non_empty(specific_abun, 'units')

        if not aggregate_abun.index.empty:
            logging.info('Start examing aggregate chemical abundances')
            pass

    def _file_to_chem_abun_data(self, file_path, refs, matrix_name, workspace_id):
        logging.info('Start reading and converting excel file data')
        data = refs

        df = self._file_to_df(file_path)

        metadata_df = None

        rename_map = {'Aggregate M/Z': 'aggregate_mz',
                      'Compound Name': 'compound_name',
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
                      'Chemical Type': 'chemical_type',
                      'Measurement Type': 'measurement_type',
                      'Units': 'units',
                      'Unit Medium': 'unit_medium',
                      'Chemical Ontology Class': 'chemical_ontology_class',
                      'Measured Identification Level': 'measured_identification_level',
                      'Chomotagraphy Type': 'chomotagraphy_type',
                      'Chemical Class': 'chemical_class',
                      'Protocol': 'protocol',
                      'Identifier': 'identifier'
                      }
        df.rename(columns=rename_map, inplace=True)

        metadata_keys = rename_map.values()

        shared_metadata_keys = list(set(metadata_keys) & set(df.columns))
        if shared_metadata_keys:
            metadata_df = df[shared_metadata_keys]
            if set(metadata_df.all(skipna=False).tolist()) == {None}:
                raise ValueError('All of metadata fields are None')
            df.drop(columns=shared_metadata_keys, inplace=True)
            self._check_chem_abun_metadata(metadata_df)
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

        if dimension == 'row':
            df = df.T

        df.fillna(0, inplace=True)

        x_train = df.values

        scaler = preprocessing.StandardScaler(with_mean=with_mean, with_std=with_std).fit(x_train)

        standardized_values = scaler.transform(x_train)

        standardize_df = pd.DataFrame(index=df.index, columns=df.columns, data=standardized_values)

        if dimension == 'row':
            standardize_df = standardize_df.T

        standardize_df.fillna(0, inplace=True)
        standardize_df.replace(np.inf, 2 ** 32, inplace=True)
        standardize_df.replace(-np.inf, -2 ** 32, inplace=True)

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

        ratio_transformed_df.fillna(0, inplace=True)
        ratio_transformed_df.replace(np.inf, 2 ** 32, inplace=True)
        ratio_transformed_df.replace(-np.inf, -2 ** 32, inplace=True)

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

    def _filtering_matrix(self, df, row_threshold=0, columns_threshold=0,
                          row_sum_threshold=10000, columns_sum_threshold=10000):

        logging.info("Removing rows with values all below {}".format(row_threshold))
        row_check = (df > row_threshold).any(axis=1)
        removed_row_ids = list(row_check[row_check == False].index)
        logging.info("Removed rows: {}".format(removed_row_ids))
        df = df.loc[row_check]

        logging.info("Removing columns with values all below {}".format(columns_threshold))
        col_check = (df > columns_threshold).any(axis=0)
        removed_col_ids = list(col_check[col_check == False].index)
        logging.info("Removed columns: {}".format(removed_col_ids))
        df = df.loc[:, col_check]

        logging.info("Removing rows with sum below {}".format(row_sum_threshold))
        row_check = df.sum(axis=1) > row_sum_threshold
        additional_removed_row_ids = list(row_check[row_check == False].index)
        removed_row_ids += additional_removed_row_ids
        logging.info("Removed rows: {}".format(additional_removed_row_ids))
        df = df.loc[row_check]

        logging.info("Removing columns with sum below {}".format(columns_sum_threshold))
        col_check = df.sum(axis=0) > columns_sum_threshold
        additional_removed_col_ids = list(col_check[col_check == False].index)
        removed_col_ids += additional_removed_col_ids
        logging.info("Removed columns: {}".format(additional_removed_col_ids))
        df = df.loc[:, col_check]

        return df

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

        relative_abundance_df.fillna(0, inplace=True)
        relative_abundance_df.replace(np.inf, 2 ** 32, inplace=True)
        relative_abundance_df.replace(-np.inf, -2 ** 32, inplace=True)

        return relative_abundance_df

    @staticmethod
    def _logit(df: pd.DataFrame):
        # entries are all in range (0,1), exclusively
        vd.assert_in_range(df, rng=(0, 1), inclusive=(False, False), opname='logit')

        f = np.vectorize(lambda p: np.log(p/(1-p)))
        df = pd.DataFrame(f(df.values), index=df.index, columns=df.columns)

        return df

    @staticmethod
    def _sqrt(df):
        # entries are nonnegative
        vd.assert_is_nonnegative(df, opname='sqrt')

        return pd.DataFrame(np.sqrt(df.values), index=df.index, columns=df.columns)

    @staticmethod
    def _log(df, base, a):
        '''
        log(a+x)
        '''
        m = df.values + a

        # entries are nonnegative
        # TODO allow 0? gives -np.inf
        vd.assert_is_nonnegative(m, opname='log')

        m = np.log(m) / np.log(base)

        return pd.DataFrame(m, index=df.index, columns=df.columns)

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
        logging.info('start performing anosim')

        anosim_res = anosim(dm, grouping, permutations=permutations)

        return dict(anosim_res)

    def _run_permanova(self, dm, grouping, permutations):
        logging.info('start performing permanova')

        permanova_res = permanova(dm, grouping, permutations=permutations)

        return dict(permanova_res)

    def _run_permdisp(self, dm, grouping, permutations):
        logging.info('start performing permdisp')

        permdisp_res = permdisp(dm, grouping, permutations=permutations)

        return dict(permdisp_res)

    def _run_mantel_tests(self, dms, labels, permutations=0, correlation_method='pearson',
                          alternative_hypothesis='two-sided'):

        logging.info('start performing mantel test')

        pwmantel_res = pwmantel(dms, labels=labels, permutations=permutations,
                                method=correlation_method, alternative=alternative_hypothesis)

        return pwmantel_res

    def _compute_target_cols(self, df, simper_ret, grouping_names):
        target_cols = [col for col in df.columns if col in str(simper_ret)]
        target_cols = list(set(target_cols))

        try:
            max_target_col_len = 18
            if len(target_cols) > max_target_col_len:
                # choose first few most influential species from each condition pair
                comp_group_len = len(simper_ret)
                num_choosen_col = max(max_target_col_len//comp_group_len, 1)
                target_cols = list()
                for comp_group in simper_ret:
                    species_pos = list(comp_group.names).index('species')
                    ord_pos = list(comp_group.names).index('ord')

                    species = list(comp_group[species_pos])
                    ord_list = list(comp_group[ord_pos])

                    target_species_pos = [i - 1 for i in ord_list[:num_choosen_col]]

                    for p in target_species_pos:
                        target_cols.append(species[p])

                target_cols = list(set(target_cols))

        except Exception:
            warning_msg = 'got unexpected error fetching most influential species'
            logging.warning(warning_msg)

        return target_cols

    def _generate_species_stats(self, df, simper_ret, grouping_names):
        logging.info('start calculating species stats')

        target_cols = self._compute_target_cols(df, simper_ret, grouping_names)

        species_stats = dict()
        for target_col in target_cols:
            logging.info('start calculating {} stats'.format(target_col))
            dist_grouping_names = set(grouping_names)
            average_abun = dict()

            abun_values = df.loc[:, target_col]
            for dist_grouping_name in dist_grouping_names:
                grouping_name_pos = [index for index, value in enumerate(grouping_names)
                                     if value == dist_grouping_name]
                filtered_abun_values = []
                for pos in grouping_name_pos:
                    filtered_abun_values.append(abun_values[pos])

                mean = 0
                std = 0
                try:
                    mean = round(np.mean(filtered_abun_values), 2)
                    std = round(np.std(filtered_abun_values), 2)
                except Exception:
                    warning_msg = 'got unexpected error calculating mean/std abundance value\n'
                    warning_msg += 'grouping_name_pos: {}\n'.format(grouping_name_pos)
                    warning_msg += 'abundance_values: {}\n'.format(abun_values)
                    warning_msg += 'returning 0 as mean/std abundance value\n'
                    logging.warning(warning_msg)
                average_abun[dist_grouping_name] = [mean, std]

            species_stats[target_col] = average_abun

        return species_stats

    def _sync_attribute_mapping(self, matrix_data, removed_ids, new_attri_mapping_name, dimension,
                                workspace_id):

        attri_mapping_ref = matrix_data.get('{}_attributemapping_ref'.format(dimension))

        if attri_mapping_ref:
            logging.info('Start removing {} from {} attribute mapping object'.format(removed_ids,
                                                                                     dimension))
            am_data = self.dfu.get_objects({"object_refs": [attri_mapping_ref]})['data'][0]['data']
            instances = am_data.get('instances', {})
            for removed_id in removed_ids:
                instances.pop(removed_id, None)

            # save new attribute mapping
            info = self.dfu.save_objects({"id": workspace_id,
                                          "objects": [{
                                                "type": 'KBaseExperiments.AttributeMapping',
                                                "data": am_data,
                                                "name": new_attri_mapping_name
                                            }]})[0]

            new_attri_mapping_ref = "%s/%s/%s" % (info[6], info[0], info[4])

            matrix_data['{}_attributemapping_ref'.format(dimension)] = new_attri_mapping_ref

            mapping = matrix_data.get('{}_mapping'.format(dimension))
            if mapping:
                for remove_id in removed_ids:
                    mapping.pop(remove_id, None)

        return matrix_data

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
                         'report_object_name': 'standardize_matrix_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        return {'new_matrix_obj_ref': new_matrix_obj_ref,
                'report_name': output['name'], 'report_ref': output['ref']}

    def perform_simper(self, params):
        logging.info('Start performing SIMPER with {}'.format(params))
        input_matrix_ref = params.get('input_matrix_ref')
        workspace_id = params.get('workspace_id')
        grouping = params.get('grouping')
        dimension = params.get('dimension', 'col')
        permutations = int(params.get('permutations', 0))

        if dimension not in ['col', 'row']:
            raise ValueError('Please use "col" or "row" for input dimension')

        input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
        input_matrix_info = input_matrix_obj['info']
        input_matrix_data = input_matrix_obj['data']

        matrix_type = input_matrix_info[2]

        if 'KBaseMatrices' in matrix_type:
            am_ref = input_matrix_data.get('{}_attributemapping_ref'.format(dimension))
            if not am_ref:
                raise ValueError('Missing {} attribute mapping from original matrix'.format(dimension))
        elif 'KBaseProfile' in matrix_type:
            profile_category = input_matrix_data.get('profile_category')
            if profile_category == 'community' and dimension == 'row':
                raise ValueError('Please choose column dimension for community profile')
            if profile_category == 'organism' and dimension == 'col':
                raise ValueError('Please choose row dimension for organism profile')
            am_ref = input_matrix_data.get('{}_attributemapping_ref'.format(dimension))
            if not am_ref:
                raise ValueError('Missing {} attribute mapping from functional profile'.format(dimension))
        else:
            raise ValueError('Unsupported data type: {}'.format(matrix_type))

        data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
        df = pd.read_json(data_matrix)

        am_ref = '{};{}'.format(input_matrix_ref, am_ref)
        am_data = self.dfu.get_objects({'object_refs': [am_ref]})['data'][0]['data']
        attribute_names = [am.get('attribute') for am in am_data.get('attributes')]

        if grouping not in attribute_names:
            raise ValueError('Cannot find {} in {} attribute mapping'.format(grouping, dimension))

        attri_pos = attribute_names.index(grouping)
        instances = am_data.get('instances')
        if dimension == 'col':
            items = df.columns
        else:
            items = df.index
        grouping_names = list()
        for item in items:
            instance = instances.get(item)
            if not instance:
                raise ValueError('Cannot find instance for {} in attribute mapping'.format(item))
            attri = instance[attri_pos]
            grouping_names.append(attri)

        logging.info('Fetched {} for {} from attributes'.format(grouping_names, grouping))

        if dimension == 'col':
            df = df.T

        df.fillna(0, inplace=True)

        vegan = rpackages.importr('vegan')
        numpy2ri.activate()

        with localconverter(ro.default_converter + pandas2ri.converter):
            simper_ret = vegan.simper(df, grouping_names, permutations=permutations)
            simper_sum = vegan.summary_simper(simper_ret)

        species_stats = self._generate_species_stats(df, simper_ret, grouping_names)

        report_output = self._generate_simper_report(workspace_id, simper_ret, simper_sum,
                                                     species_stats, grouping_names)

        return report_output

    def perform_rarefy(self, params):
        logging.info('Start performing rarefying matrix with {}'.format(params))
        warnings = []
        input_matrix_ref = params.get('input_matrix_ref')
        workspace_id = params.get('workspace_id')
        new_matrix_name = params.get('new_matrix_name')

        seed_number = params.get('seed_number', 'do_not_seed')
        subsample_size = params.get('subsample_size')
        dimension = params.get('dimension', 'col')

        bootstrap = params.get('bootstrap')
        if bootstrap is not None:
            num_rare_reps = bootstrap['num_rare_reps']
            central_tendency = bootstrap['central_tendency']

        input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
        input_matrix_info = input_matrix_obj['info']
        input_matrix_name = input_matrix_info[1]
        input_matrix_data = input_matrix_obj['data']

        for key, obj_data in input_matrix_data.items():
            if key.endswith('_ref'):
                subobj_ref = input_matrix_data[key]
                input_matrix_data[key] = '{};{}'.format(input_matrix_ref, subobj_ref)
                logging.info('updated {} to {}'.format(key, input_matrix_data[key]))

        for dim in ['row', 'col']:
            attribute_mapping = input_matrix_data.get('{}_mapping'.format(dim))
            attributemapping_ref = input_matrix_data.get('{}_attributemapping_ref'.format(dim))
            if not attribute_mapping and attributemapping_ref:
                am_data = self.dfu.get_objects({'object_refs': [attributemapping_ref]})['data'][0]['data']
                attribute_mapping = {x: x for x in am_data['instances'].keys()}
                input_matrix_data['{}_mapping'.format(dim)] = attribute_mapping

        if not new_matrix_name:
            current_time = time.localtime()
            new_matrix_name = input_matrix_name + time.strftime('_%H_%M_%S_%Y_%m_%d', current_time)

        data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
        df = pd.read_json(data_matrix)
        # original_matrix_df = df.copy(deep=True)

        if dimension == 'col':
            df = df.T

        df.fillna(0, inplace=True)

        run_seed = (not seed_number == 'do_not_seed')

        # determining subsample size
        raremax = int(min(df.sum(axis=1)))  # least sample size
        if subsample_size is None:  # default behavior: use least sample size
            subsample_size = raremax
        else:  # user-specified behavior, find any samples too small
            unrarefied = df.index[df.sum(axis=1) < subsample_size].tolist()
            if len(unrarefied) > 0:
                msg = (
                    'At subsampling size %d, samples %s are too small and will not be rarefied. '
                    'Smallest sample size is %d'
                    % (subsample_size, str(unrarefied), raremax)
                )
                warnings.append(msg)
                logging.info(msg)
        logging.info('Using subsample size %d' % subsample_size)

        vegan = rpackages.importr('vegan')
        numpy2ri.activate()

        # generating rarefied matrix
        logging.info('Start executing rrarefy(s)')
        if run_seed:
            ro.r('set.seed({})'.format(seed_number))
        if bootstrap is None:
            with localconverter(ro.default_converter + pandas2ri.converter):
                random_rare = vegan.rrarefy(df, subsample_size)
        else:
            random_rare_l = []
            for rep in range(num_rare_reps):
                with localconverter(ro.default_converter + pandas2ri.converter):
                    random_rare = vegan.rrarefy(df, subsample_size)  # returns np.ndarray
                random_rare_l.append(random_rare)
            if central_tendency == 'mean':
                random_rare = sum(random_rare_l) / num_rare_reps
            elif central_tendency == 'median':
                random_rare = np.median(random_rare_l, axis=0)
            else:
                raise NotImplementedError('Unknown value for `central_tendency`')

        random_rare_df = pd.DataFrame(random_rare, index=df.index, columns=df.columns)

        if dimension == 'col':
            random_rare_df = random_rare_df.T

        # generating plots
        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)

        logging.info('Start generating rarecurve plot')
        rarecurve_image = os.path.join(result_directory, 'rarecurve.jpg')
        ro.r("jpeg('{}')".format(rarecurve_image))
        if run_seed:
            ro.r('set.seed({})'.format(seed_number))
        with localconverter(ro.default_converter + pandas2ri.converter):
            vegan.rarecurve(df, sample=subsample_size, step=20, col="blue", cex=0.6)
        ro.r('dev.off()')

        logging.info('Start generating expected species richness vs raw abundance plot')
        with localconverter(ro.default_converter + pandas2ri.converter):
            Srare = vegan.rarefy(df, subsample_size)
        specnumber = ro.r['specnumber']
        with localconverter(ro.default_converter + pandas2ri.converter):
            S = specnumber(df)
        obs_vs_rare_image = os.path.join(result_directory, 'obs_vs_rare.jpg')
        ro.r("jpeg('{}')".format(obs_vs_rare_image))
        plot = ro.r['plot']
        plot(S, Srare, xlab="Observed No. of Species", ylab="Rarefied No. of Species")
        ro.r('dev.off()')

        new_matrix_data = {'row_ids': random_rare_df.index.tolist(),
                           'col_ids': random_rare_df.columns.tolist(),
                           'values': random_rare_df.values.tolist()}

        input_matrix_data['data'] = new_matrix_data

        logging.info("Saving new rarefy matrix object")

        new_matrix_obj_ref = self.data_util.save_object({
                                                'obj_type': input_matrix_info[2],
                                                'obj_name': new_matrix_name,
                                                'data': input_matrix_data,
                                                'workspace_id': workspace_id})['obj_ref']

        returnVal = {'new_matrix_obj_ref': new_matrix_obj_ref}

        report_output = self._generate_rarefy_report(new_matrix_obj_ref, workspace_id,
                                                     random_rare_df,
                                                     rarecurve_image, obs_vs_rare_image,
                                                     warnings)

        returnVal.update(report_output)

        return returnVal

    def perform_mantel_test(self, params):

        logging.info('Start performing mantel test with {}'.format(params))

        input_matrix_refs = params.get('input_matrix_refs')
        workspace_id = params.get('workspace_id')
        dimension = params.get('dimension', 'col')
        dist_metric = params.get('dist_metric', 'euclidean')

        correlation_method = params.get('correlation_method', 'pearson')
        permutations = params.get('permutations', 0)
        alternative_hypothesis = params.get('alternative_hypothesis', 'two-sided')

        if dimension not in ['col', 'row']:
            raise ValueError('Please use "col" or "row" for input dimension')

        if len(input_matrix_refs) < 2:
            raise ValueError('Please provide at least 2 matrices to perform mentel test')

        dms = list()
        labels = list()
        for input_matrix_ref in input_matrix_refs:
            input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
            input_matrix_info = input_matrix_obj['info']
            input_matrix_name = input_matrix_info[1]
            labels.append(input_matrix_name)

            data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
            df = pd.read_json(data_matrix)

            dm = self._create_distance_matrix(df, dist_metric=dist_metric, dimension=dimension)
            dms.append(dm)

        pwmantel_res = self._run_mantel_tests(dms, labels, permutations=permutations,
                                              correlation_method=correlation_method,
                                              alternative_hypothesis=alternative_hypothesis)

        report_output = self._generate_mantel_test_report(workspace_id, pwmantel_res)

        return report_output

    def perform_variable_stats_matrix(self, params):

        logging.info('Start performing variable statistics with {}'.format(params))

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

        if dimension not in ['col', 'row']:
            raise ValueError('Please use "col" or "row" for input dimension')

        input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
        input_matrix_info = input_matrix_obj['info']
        input_matrix_data = input_matrix_obj['data']

        matrix_type = input_matrix_info[2]
        if 'KBaseMatrices' in matrix_type:
            am_ref = input_matrix_data.get('{}_attributemapping_ref'.format(dimension))
            if not am_ref:
                raise ValueError('Missing {} attribute mapping from original matrix'.format(dimension))
        elif 'KBaseProfile' in matrix_type:
            profile_category = input_matrix_data.get('profile_category')
            if profile_category == 'community' and dimension == 'row':
                raise ValueError('Please choose column dimension for community profile')
            if profile_category == 'organism' and dimension == 'col':
                raise ValueError('Please choose row dimension for organism profile')
            am_ref = input_matrix_data.get('{}_attributemapping_ref'.format(dimension))
            if not am_ref:
                raise ValueError('Missing {} attribute mapping from functional profile'.format(dimension))
        else:
            raise ValueError('Unsupported data type: {}'.format(matrix_type))

        data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
        df = pd.read_json(data_matrix)
        dm = self._create_distance_matrix(df, dist_metric=dist_metric, dimension=dimension)

        am_ref = '{};{}'.format(input_matrix_ref, am_ref)
        am_data = self.dfu.get_objects({'object_refs': [am_ref]})['data'][0]['data']
        attribute_names = [am.get('attribute') for am in am_data.get('attributes')]

        if grouping not in attribute_names:
            raise ValueError('Cannot find {} in {} attribute mapping'.format(grouping, dimension))

        attri_pos = attribute_names.index(grouping)
        instances = am_data.get('instances')
        if dimension == 'col':
            items = df.columns
        else:
            items = df.index
        grouping_names = list()
        for item in items:
            instance = instances.get(item)
            if not instance:
                raise ValueError('Cannot find instance for {} in attribute mapping'.format(item))
            attri = instance[attri_pos]
            grouping_names.append(attri)
        logging.info('Fetched {} for {} from attributes'.format(grouping_names, grouping))

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

    def transform_matrix_variable_specific(self, params):
        """
        Transform a list of variable from a matrix
        """

        OPS = [
            'relative_abundance',
            'standardization',
            'ratio_transformation',
            'log',
            'sqrt',
            'logit',
        ]

        logging.info('Start performing transformation with {}'.format(params))

        input_matrix_ref = params.get('input_matrix_ref')
        workspace_id = params.get('workspace_id')
        new_matrix_name = params.get('new_matrix_name')

        operations = params.get('operations')
        relative_abundance_params = params.get('perform_relative_abundance', {})
        standardization_params = params.get('standardization_params', {})
        ratio_transformation_params = params.get('ratio_transformation_params', {})
        log_params = params.get('log_params', {})

        dimension = params.get('dimension', 'row')
        variables = params.get('variables', list())

        if dimension not in ['col', 'row']:
            raise ValueError('Please use "col" or "row" for input dimension')

        # validate operations
        MAX_OPS = 15
        if len(operations) > MAX_OPS:
            raise Exception('Maximum allowed number of operations is %d' % MAX_OPS)
        for op in operations:
            if op not in OPS:
                raise Exception('Operation %s not in allowed %s' % (str(op), OPS))

        input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
        input_matrix_info = input_matrix_obj['info']
        input_matrix_name = input_matrix_info[1]
        input_matrix_data = input_matrix_obj['data']

        unmatched_var = set(variables) - set(input_matrix_data['data']['{}_ids'.format(dimension)])
        if unmatched_var:
            raise ValueError('variable [{}] is not contained in {} ids from matrix'.format(
                                                                                    unmatched_var,
                                                                                    dimension))

        for key, obj_data in input_matrix_data.items():
            if key.endswith('_ref'):
                subobj_ref = input_matrix_data[key]
                input_matrix_data[key] = '{};{}'.format(input_matrix_ref, subobj_ref)
                logging.info('updated {} to {}'.format(key, input_matrix_data[key]))

        for dim in ['row', 'col']:
            attribute_mapping = input_matrix_data.get('{}_mapping'.format(dim))
            attributemapping_ref = input_matrix_data.get('{}_attributemapping_ref'.format(dim))
            if not attribute_mapping and attributemapping_ref:
                am_data = self.dfu.get_objects({'object_refs': [attributemapping_ref]})['data'][0]['data']
                attribute_mapping = {x: x for x in am_data['instances'].keys()}
                input_matrix_data['{}_mapping'.format(dim)] = attribute_mapping

        if not new_matrix_name:
            current_time = time.localtime()
            new_matrix_name = input_matrix_name + time.strftime('_%H_%M_%S_%Y_%m_%d', current_time)

        data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
        df = pd.read_json(data_matrix)
        original_df = df.copy(deep=True)
        original_row_ids = original_df.index
        original_col_ids = original_df.columns

        if dimension == 'row':
            selected_df = original_df.loc[list(set(variables))]
        else:
            selected_df = original_df.loc[:, list(set(variables))]

        # iterate over operations
        df_results = []
        for op in operations:

            if op == 'relative_abundance':
                selected_df = self._relative_abundance(
                                    selected_df,
                                    dimension=relative_abundance_params.get(
                                                'relative_abundance_dimension', 'col'))

            elif op == 'standardization':
                selected_df = self._standardize_df(
                                    selected_df,
                                    dimension=standardization_params.get(
                                                            'standardization_dimension', 'col'),
                                    with_mean=standardization_params.get(
                                                            'standardization_with_mean', True),
                                    with_std=standardization_params.get(
                                                            'standardization_with_std', True))

            elif op == 'ratio_transformation':
                selected_df = self._ratio_trans_df(
                                    selected_df,
                                    method=ratio_transformation_params.get(
                                                            'ratio_transformation_method', 'clr'),
                                    dimension=ratio_transformation_params.get(
                                                            'ratio_transformation_dimension', 'col'))

            elif op == 'logit':
                selected_df = self._logit(selected_df)

            elif op == 'sqrt':
                selected_df = self._sqrt(selected_df)

            elif op == 'log':
                selected_df = self._log(selected_df,
                                        base=log_params.get('base', 10),
                                        a=log_params.get('offset', 1))

            else:
                raise NotImplementedError('Unknown op `%s`' % op)

            df_results.append(selected_df.copy(deep=True))

        if dimension == 'row':
            df = selected_df.combine_first(original_df)
        else:
            df = selected_df.T.combine_first(original_df.T).T

        new_matrix_data = {'row_ids': df.index.tolist(),
                           'col_ids': df.columns.tolist(),
                           'values': df.values.tolist()}

        input_matrix_data['data'] = new_matrix_data

        removed_row_ids = set(original_row_ids) - set(df.index)
        removed_col_ids = set(original_col_ids) - set(df.columns)

        if removed_row_ids and input_matrix_data.get('row_attributemapping_ref'):
            input_matrix_data = self._sync_attribute_mapping(input_matrix_data, removed_row_ids,
                                                             new_matrix_name + '_row_attribute_mapping',
                                                             'row', workspace_id)

        if removed_col_ids and input_matrix_data.get('col_attributemapping_ref'):
            input_matrix_data = self._sync_attribute_mapping(input_matrix_data, removed_col_ids,
                                                             new_matrix_name + '_col_attribute_mapping',
                                                             'col', workspace_id)

        logging.info("Saving new transformed matrix object")
        new_matrix_obj_ref = self.data_util.save_object({
                                                'obj_type': input_matrix_info[2],
                                                'obj_name': new_matrix_name,
                                                'data': input_matrix_data,
                                                'workspace_id': workspace_id})['obj_ref']

        returnVal = {'new_matrix_obj_ref': new_matrix_obj_ref}

        report_output = self._generate_transform_report(new_matrix_obj_ref, workspace_id,
                                                        operations, df_results,
                                                        variable_specific=True)

        returnVal.update(report_output)

        return returnVal

    def transform_matrix(self, params):
        """
        Transform a matrix
        """

        OPS = [
            'abundance_filtering',
            'relative_abundance',
            'standardization',
            'ratio_transformation',
            'log',
            'sqrt',
            'logit',
        ]

        logging.info('Start performing transformation with {}'.format(params))

        input_matrix_ref = params.get('input_matrix_ref')
        workspace_id = params.get('workspace_id')
        new_matrix_name = params.get('new_matrix_name')

        operations = params.get('operations')
        abundance_filtering_params = params.get('abundance_filtering_params', {})
        relative_abundance_params = params.get('perform_relative_abundance', {})
        standardization_params = params.get('standardization_params', {})
        ratio_transformation_params = params.get('ratio_transformation_params', {})
        log_params = params.get('log_params', {})

        # validate operations
        MAX_OPS = 15
        if len(operations) > MAX_OPS:
            raise Exception('Maximum allowed number of operations is %d' % MAX_OPS)
        for op in operations:
            if op not in OPS:
                raise Exception('Operation %s not in allowed %s' % (str(op), OPS))

        input_matrix_obj = self.dfu.get_objects({'object_refs': [input_matrix_ref]})['data'][0]
        input_matrix_info = input_matrix_obj['info']
        input_matrix_name = input_matrix_info[1]
        input_matrix_data = input_matrix_obj['data']

        for key, obj_data in input_matrix_data.items():
            if key.endswith('_ref'):
                subobj_ref = input_matrix_data[key]
                input_matrix_data[key] = '{};{}'.format(input_matrix_ref, subobj_ref)
                logging.info('updated {} to {}'.format(key, input_matrix_data[key]))

        for dim in ['row', 'col']:
            attribute_mapping = input_matrix_data.get('{}_mapping'.format(dim))
            attributemapping_ref = input_matrix_data.get('{}_attributemapping_ref'.format(dim))
            if not attribute_mapping and attributemapping_ref:
                am_data = self.dfu.get_objects({'object_refs': [attributemapping_ref]})['data'][0]['data']
                attribute_mapping = {x: x for x in am_data['instances'].keys()}
                input_matrix_data['{}_mapping'.format(dim)] = attribute_mapping

        if not new_matrix_name:
            current_time = time.localtime()
            new_matrix_name = input_matrix_name + time.strftime('_%H_%M_%S_%Y_%m_%d', current_time)

        data_matrix = self.data_util.fetch_data({'obj_ref': input_matrix_ref}).get('data_matrix')
        df = pd.read_json(data_matrix)
        original_row_ids = df.index
        original_col_ids = df.columns

        # iterate over operations
        df_results = []
        for op in operations:

            if op == 'abundance_filtering':
                df = self._filtering_matrix(df,
                                            row_threshold=abundance_filtering_params.get(
                                                            'abundance_filtering_row_threshold', 0),
                                            columns_threshold=abundance_filtering_params.get(
                                                        'abundance_filtering_columns_threshold', 0),
                                            row_sum_threshold=abundance_filtering_params.get(
                                                        'abundance_filtering_row_sum_threshold',
                                                        10000),
                                            columns_sum_threshold=abundance_filtering_params.get(
                                                        'abundance_filtering_columns_sum_threshold',
                                                        10000))

            elif op == 'relative_abundance':
                df = self._relative_abundance(df, dimension=relative_abundance_params.get(
                                                            'relative_abundance_dimension', 'col'))

            elif op == 'standardization':
                df = self._standardize_df(df,
                                          dimension=standardization_params.get(
                                                                'standardization_dimension', 'col'),
                                          with_mean=standardization_params.get(
                                                                'standardization_with_mean', True),
                                          with_std=standardization_params.get(
                                                                'standardization_with_std', True))

            elif op == 'ratio_transformation':
                df = self._ratio_trans_df(df,
                                          method=ratio_transformation_params.get(
                                                            'ratio_transformation_method', 'clr'),
                                          dimension=ratio_transformation_params.get(
                                                            'ratio_transformation_dimension', 'col'))

            elif op == 'logit':
                df = self._logit(df)

            elif op == 'sqrt':
                df = self._sqrt(df)

            elif op == 'log':
                df = self._log(df, base=log_params.get('base', 10), a=log_params.get('offset', 1))

            else:
                raise NotImplementedError('Unknown op `%s`' % op)

            df_results.append(df.copy(deep=True))

        new_matrix_data = {'row_ids': df.index.tolist(),
                           'col_ids': df.columns.tolist(),
                           'values': df.values.tolist()}

        input_matrix_data['data'] = new_matrix_data

        removed_row_ids = set(original_row_ids) - set(df.index)
        removed_col_ids = set(original_col_ids) - set(df.columns)

        if removed_row_ids and input_matrix_data.get('row_attributemapping_ref'):
            input_matrix_data = self._sync_attribute_mapping(input_matrix_data, removed_row_ids,
                                                             new_matrix_name + '_row_attribute_mapping',
                                                             'row', workspace_id)

        if removed_col_ids and input_matrix_data.get('col_attributemapping_ref'):
            input_matrix_data = self._sync_attribute_mapping(input_matrix_data, removed_col_ids,
                                                             new_matrix_name + '_col_attribute_mapping',
                                                             'col', workspace_id)

        logging.info("Saving new transformed matrix object")
        new_matrix_obj_ref = self.data_util.save_object({
                                                'obj_type': input_matrix_info[2],
                                                'obj_name': new_matrix_name,
                                                'data': input_matrix_data,
                                                'workspace_id': workspace_id})['obj_ref']

        returnVal = {'new_matrix_obj_ref': new_matrix_obj_ref}

        report_output = self._generate_transform_report(new_matrix_obj_ref, workspace_id,
                                                        operations, df_results)

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
                                                'workspace_id': workspace_id})['obj_ref']

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
                                                'workspace_id': workspace_id})['obj_ref']

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
