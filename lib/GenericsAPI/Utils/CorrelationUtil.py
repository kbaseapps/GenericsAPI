import errno
import json
import logging
import os
import shutil
import traceback
import uuid
import math

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from matplotlib import pyplot as plt
from plotly.offline import plot
from scipy import stats
from natsort import natsorted
from sklearn import metrics
from sklearn.linear_model import LinearRegression

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.kb_GenericsReportClient import kb_GenericsReport
from GenericsAPI.Utils.DataUtil import DataUtil
from installed_clients.KBaseReportClient import KBaseReport

CORR_METHOD = ['pearson', 'kendall', 'spearman', 'mutual_info']  # correlation method
HIDDEN_SEARCH_THRESHOLD = 1500


class CorrelationUtil:

    def _mkdir_p(self, path):
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

    def _validate_compute_corr_matrix_params(self, params):
        """
        _validate_compute_corr_matrix_params:
            validates params passed to compute_correlation_matrix method
        """

        logging.info('start validating compute_corrrelation_matrix params')

        # check for required parameters
        for p in ['input_obj_ref', 'workspace_name', 'corr_matrix_name']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _validate_compute_correlation_across_matrices_params(self, params):
        """
        _validate_compute_correlation_across_matrices_params:
            validates params passed to compute_correlation_across_matrices method
        """

        logging.info('start validating compute_correlation_across_matrices params')

        # check for required parameters
        for p in ['workspace_name', 'corr_matrix_name', 'matrix_ref_1', 'matrix_ref_2']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _build_top_corr_table(self, matrix_2D, output_directory,
                              original_matrix_ref=[], sig_matrix_2D=None):

        row_ids = matrix_2D.get('row_ids')
        col_ids = matrix_2D.get('col_ids')
        values = matrix_2D.get('values')

        data_df = pd.DataFrame(values, index=row_ids, columns=col_ids)
        data_df.fillna(0, inplace=True)

        sig_df = None
        if sig_matrix_2D is not None:
            sig_df = pd.DataFrame(sig_matrix_2D.get('values'),
                                  index=sig_matrix_2D.get('row_ids'),
                                  columns=sig_matrix_2D.get('col_ids'))
            sig_df.fillna(0, inplace=True)

        columns = list()
        if len(original_matrix_ref) == 1:
            data_df = data_df.mask(np.tril(np.ones(data_df.shape)).astype(np.bool))  # remove duplicate entries
            res = self.dfu.get_objects({'object_refs': [original_matrix_ref[0]]})['data'][0]
            obj_type = res['info'][2]
            matrix_type = obj_type.split('-')[0].split('Matrix')[0].split('.')[-1]
            columns.extend(['{} 1'.format(matrix_type), '{} 2'.format(matrix_type)])
        elif len(original_matrix_ref) == 2:
            for matrix_ref in original_matrix_ref:
                res = self.dfu.get_objects({'object_refs': [matrix_ref]})['data'][0]
                obj_type = res['info'][2]
                matrix_type = obj_type.split('-')[0].split('Matrix')[0].split('.')[-1]
                columns.append(matrix_type)
        else:
            columns = ['Variable 1', 'Variable 2']

        value_col_name = 'Correlation Coefficient'
        columns.append(value_col_name)

        links = data_df.stack().reset_index()
        # links = links[links.iloc[:, 0] != links.iloc[:, 1]]  # remove self-comparison

        links.columns = columns

        # sort by absolute value
        links = links.iloc[(-links[value_col_name].abs()).argsort()].reset_index(drop=True)

        top_corr_limit = 200
        top_corr = links[:top_corr_limit]
        top_corr_size = top_corr.index.size

        if sig_df is not None:
            sig_values = list()
            for i in range(top_corr_size):
                corr_pair = top_corr.iloc[i]
                first_item = corr_pair[0]
                second_item = corr_pair[1]
                sig_value = sig_df.loc[first_item, second_item]
                sig_values.append(sig_value)

            top_corr['P-Value'] = sig_values

        headerColor = 'grey'
        rowEvenColor = 'lightgrey'
        rowOddColor = 'white'

        # color codes from px.colors.diverging.RdBu
        colors = ['rgb(247,247,247)', 'rgb(253,219,199)', 'rgb(244,165,130)', 'rgb(214,96,77)',
                  'rgb(67,147,195)', 'rgb(146,197,222)', 'rgb(209,229,240)']
        interval = 0.3
        corr_color_idx = (top_corr[value_col_name]/interval).apply(int)  # divid coefficient by interval and round to int

        fig = go.Figure(data=[go.Table(
            header=dict(values=list(top_corr.columns),
                        line_color='darkslategray',
                        fill_color=headerColor,
                        align='left',
                        font=dict(color='white', size=12)),
            cells=dict(values=top_corr.T.values,
                       line_color='darkslategray',
                       fill_color=[[rowOddColor, rowEvenColor]*top_corr_limit,
                                   [rowOddColor, rowEvenColor]*top_corr_limit,
                                   np.array(colors)[corr_color_idx],
                                   [rowOddColor, rowEvenColor]*top_corr_limit],
                       align='left',
                       font=dict(color='darkslategray', size=11)))
        ])

        if top_corr_size < top_corr_limit:
            fig_title = 'All {} Coefficient Pairs'.format(top_corr_size)
        else:
            fig_title = 'Top {} Coefficient Pairs'.format(top_corr_limit)
        fig.update_layout(
            width=1200,
            height=2000,
            title=dict(text=fig_title, x=0.5,
                       font=dict(family='Times New Roman', size=30, color='Purple')))

        table_file_name = 'top_corr_table.html'
        table_file_path = os.path.join(output_directory, table_file_name)
        fig.write_html(table_file_path)

        tab_content = ''
        tab_content += '\n<iframe height="1300px" width="100%" '
        tab_content += 'src="{}" '.format(table_file_name)
        tab_content += 'style="border:none;"></iframe>'

        return tab_content

    def _build_top_scatter_plot(self, matrix_2D, output_directory, df1, df2,
                                original_matrix_ref=[], sig_matrix_2D=None):
        row_ids = matrix_2D.get('row_ids')
        col_ids = matrix_2D.get('col_ids')
        values = matrix_2D.get('values')

        data_df = pd.DataFrame(values, index=row_ids, columns=col_ids)
        data_df.fillna(0, inplace=True)

        sig_df = None
        if sig_matrix_2D is not None:
            sig_df = pd.DataFrame(sig_matrix_2D.get('values'),
                                  index=sig_matrix_2D.get('row_ids'),
                                  columns=sig_matrix_2D.get('col_ids'))
            sig_df.fillna(0, inplace=True)

        columns = list()
        if len(original_matrix_ref) == 1:
            data_df = data_df.mask(np.tril(np.ones(data_df.shape)).astype(np.bool))
            res = self.dfu.get_objects({'object_refs': [original_matrix_ref[0]]})['data'][0]
            obj_type = res['info'][2]
            matrix_type = obj_type.split('-')[0].split('Matrix')[0].split('.')[-1]
            columns.extend(['{} 1'.format(matrix_type), '{} 2'.format(matrix_type)])
        elif len(original_matrix_ref) == 2:
            for matrix_ref in original_matrix_ref:
                res = self.dfu.get_objects({'object_refs': [matrix_ref]})['data'][0]
                obj_type = res['info'][2]
                matrix_type = obj_type.split('-')[0].split('Matrix')[0].split('.')[-1]
                columns.append(matrix_type)
        else:
            columns = ['Variable 1', 'Variable 2']

        value_col_name = 'Correlation Coefficient'
        columns.append(value_col_name)

        links = data_df.stack().reset_index()
        # links = links[links.iloc[:, 0] != links.iloc[:, 1]]  # remove self-comparison

        links.columns = columns

        # sort by absolute value
        links = links.iloc[(-links[value_col_name].abs()).argsort()].reset_index(drop=True)

        top_corr_limit = 20
        top_corr = links[:top_corr_limit]
        num_plots = top_corr.index.size

        warnings = ''
        if top_corr_limit < links.index.size:
            warnings += 'Note: Limiting to top {} plots'.format(top_corr_limit)

        num_cols = 3
        num_rows = math.ceil(num_plots/num_cols)

        fig = make_subplots(rows=num_rows, cols=num_cols)
        colors = px.colors.qualitative.Plotly
        colors_size = len(colors)
        for i in range(num_plots):
            corr_pair = top_corr.iloc[i]
            first_item = corr_pair[0]
            second_item = corr_pair[1]
            corr_r = corr_pair[2]

            first_item_matrix_value = list(df1.loc[first_item].values)
            second_item_matrix_value = list(df2.loc[second_item].values)

            sub_fig = go.Scatter(
                x=first_item_matrix_value,
                y=second_item_matrix_value,
                mode='markers',
                showlegend=False,
                opacity=0.65,
                marker_color=colors[i % colors_size],
                marker=dict(size=8,)
            )
            X = df1.loc[first_item].values.reshape(-1, 1)
            model = LinearRegression()
            model.fit(X, df2.loc[second_item].values)
            x_range = np.linspace(X.min(), X.max(), 100)
            y_range = model.predict(x_range.reshape(-1, 1))

            sub_fig_trend = go.Scatter(
                x=x_range,
                y=y_range,
                showlegend=False,
                marker_color=colors[i % colors_size],
            )

            fig.add_trace(sub_fig, row=i//num_cols + 1, col=i % num_cols + 1)
            fig.add_trace(sub_fig_trend, row=i//num_cols + 1, col=i % num_cols + 1)

            anno_text = 'correlation coefficient={}'.format(corr_r)
            if i == 0:
                fig.update_layout({'xaxis': {'title': '{} ({})'.format(first_item[:10], links.columns[0])}})
                fig.update_layout({'yaxis': {'title': '{} ({})'.format(second_item[:10], links.columns[1])}})
                x_start = fig['layout']['xaxis']['domain'][0]
                x_end = fig['layout']['xaxis']['domain'][1]
                fig.add_annotation(
                    text=anno_text,
                    x=(x_end - x_start) / 2 + x_start,
                    y=fig['layout']['yaxis']['domain'][1],
                    xref='paper',
                    yref='paper',
                    xanchor='center',
                    yanchor='top',
                    showarrow=False
                  )
                if sig_df is not None:
                    sig_value = sig_df.loc[first_item, second_item]
                    anno_text = 'p-value={}'.format(sig_value)
                    fig.add_annotation(
                        text=anno_text,
                        x=(x_end - x_start) / 2 + x_start,
                        y=fig['layout']['yaxis']['domain'][1] - 0.007,
                        xref='paper',
                        yref='paper',
                        xanchor='center',
                        yanchor='top',
                        showarrow=False
                      )
            else:
                fig.update_layout({'xaxis{}'.format(i+1): {'title': '{} ({})'.format(first_item[:10], links.columns[0])}})
                fig.update_layout({'yaxis{}'.format(i+1): {'title': '{} ({})'.format(second_item[:10], links.columns[1])}})
                x_start = fig['layout']['xaxis{}'.format(i+1)]['domain'][0]
                x_end = fig['layout']['xaxis{}'.format(i+1)]['domain'][1]
                fig.add_annotation(
                    text=anno_text,
                    x=(x_end - x_start) / 2 + x_start,
                    y=fig['layout']['yaxis{}'.format(i+1)]['domain'][1],
                    xref='paper',
                    yref='paper',
                    xanchor='center',
                    yanchor='top',
                    showarrow=False
                  )
                if sig_df is not None:
                    sig_value = sig_df.loc[first_item, second_item]
                    anno_text = 'p-value={}'.format(sig_value)
                    fig.add_annotation(
                        text=anno_text,
                        x=(x_end - x_start) / 2 + x_start,
                        y=fig['layout']['yaxis{}'.format(i+1)]['domain'][1] - 0.007,
                        xref='paper',
                        yref='paper',
                        xanchor='center',
                        yanchor='top',
                        showarrow=False
                      )

        fig_title = 'Scatter Plot For Top {} Coefficient Pairs'.format(num_plots)
        fig.update_layout(
            width=1200,
            height=2000,
            title=dict(text=fig_title, x=0.5,
                       font=dict(family='Times New Roman', size=30, color='Purple')),
            font=dict(family="Courier New, monospace", size=10, color="RebeccaPurple"))

        plot_file_name = 'top_scatter_plot.html'
        plot_file_path = os.path.join(output_directory, plot_file_name)
        fig.write_html(plot_file_path)

        tab_content = ''
        tab_content += '\n<iframe height="1300px" width="100%" '
        tab_content += 'src="{}" '.format(plot_file_name)
        tab_content += 'style="border:none;"></iframe>'

        return tab_content

    def _build_heatmap_content(self, matrix_2D, output_directory, centered_by=None):
        row_ids = matrix_2D.get('row_ids')
        col_ids = matrix_2D.get('col_ids')
        values = matrix_2D.get('values')

        data_df = pd.DataFrame(values, index=row_ids, columns=col_ids)
        data_df.fillna(0, inplace=True)

        tsv_file_path = os.path.join(output_directory, 'heatmap_data_{}.tsv'.format(
                                                                    str(uuid.uuid4())))

        data_df.to_csv(tsv_file_path)
        heatmap_dir = self.report_util.build_heatmap_html({
                                            'tsv_file_path': tsv_file_path,
                                            'cluster_data': True,
                                            'centered_by': centered_by})['html_dir']

        heatmap_report_files = os.listdir(heatmap_dir)

        heatmap_index_page = None
        for heatmap_report_file in heatmap_report_files:
            if heatmap_report_file.endswith('.html'):
                heatmap_index_page = heatmap_report_file

            shutil.copy2(os.path.join(heatmap_dir, heatmap_report_file),
                         output_directory)
        tab_content = '\n'

        if heatmap_index_page:
            tab_content += '\n<iframe height="1300px" width="100%" '
            tab_content += 'src="{}" '.format(heatmap_index_page)
            tab_content += 'style="border:none;"></iframe>'
        else:
            tab_content += '''\n<p style="color:red;" >'''
            tab_content += '''Heatmap is too large to be displayed.</p>\n'''

        return tab_content

    def _generate_visualization_content(self, output_directory, corr_matrix_obj_ref,
                                        corr_matrix_plot_path, scatter_plot_path, df1, df2):

        """
        <div class="tab">
            <button class="tablinks" onclick="openTab(event, 'CorrelationMatrix')" id="defaultOpen">Correlation Matrix</button>
        </div>

        <div id="CorrelationMatrix" class="tabcontent">
            <p>CorrelationMatrix_Content</p>
        </div>"""

        tab_def_content = ''
        tab_content = ''

        corr_data = self.dfu.get_objects({'object_refs': [corr_matrix_obj_ref]})['data'][0]['data']

        coefficient_data = corr_data.get('coefficient_data')
        significance_data = corr_data.get('significance_data')
        original_matrix_ref = corr_data.get('original_matrix_ref')

        tab_def_content += """
        <div class="tab">
            <button class="tablinks" onclick="openTab(event, 'CorrelationTable')" id="defaultOpen">Coefficient Table</button>
        """
        corr_table_content = self._build_top_corr_table(coefficient_data, output_directory,
                                                        original_matrix_ref=original_matrix_ref,
                                                        sig_matrix_2D=significance_data)
        tab_content += """
        <div id="CorrelationTable" class="tabcontent">{}</div>""".format(corr_table_content)

        tab_def_content += """
        <button class="tablinks" onclick="openTab(event, 'ScatterTopMatrix')">Plot Top Coefficient Pairs</button>
        """
        scatter_top_content = self._build_top_scatter_plot(coefficient_data, output_directory,
                                                           df1, df2,
                                                           original_matrix_ref=original_matrix_ref,
                                                           sig_matrix_2D=significance_data)
        tab_content += """
        <div id="ScatterTopMatrix" class="tabcontent">{}</div>""".format(scatter_top_content)

        tab_def_content += """
        <div class="tab">
            <button class="tablinks" onclick="openTab(event, 'CorrelationMatrix')">Correlation Matrix Heatmap</button>
        """
        corr_heatmap_content = self._build_heatmap_content(coefficient_data, output_directory,
                                                           centered_by=0)
        tab_content += """
        <div id="CorrelationMatrix" class="tabcontent">{}</div>""".format(corr_heatmap_content)

        if significance_data:
            tab_def_content += """
            <button class="tablinks" onclick="openTab(event, 'SignificanceMatrix')">Significance Matrix Heatmap</button>
            """
            sig_heatmap_content = self._build_heatmap_content(significance_data, output_directory)
            tab_content += """
            <div id="SignificanceMatrix" class="tabcontent">{}</div>""".format(sig_heatmap_content)

        if corr_matrix_plot_path:
            tab_def_content += """
            <button class="tablinks" onclick="openTab(event, 'CorrelationMatrixPlot')">Correlation Matrix Heatmap</button>
            """

            tab_content += """
            <div id="CorrelationMatrixPlot" class="tabcontent">
            """
            if corr_matrix_plot_path.endswith('.png'):
                corr_matrix_plot_name = 'CorrelationMatrixPlot.png'
                corr_matrix_plot_display_name = 'Correlation Matrix Plot'

                shutil.copy2(corr_matrix_plot_path,
                             os.path.join(output_directory, corr_matrix_plot_name))

                tab_content += '<div class="gallery">'
                tab_content += '<a target="_blank" href="{}">'.format(corr_matrix_plot_name)
                tab_content += '<img src="{}" '.format(corr_matrix_plot_name)
                tab_content += 'alt="{}" width="600" height="400">'.format(
                                                                    corr_matrix_plot_display_name)
                tab_content += '</a><div class="desc">{}</div></div>'.format(
                                                                corr_matrix_plot_display_name)
            elif corr_matrix_plot_path.endswith('.html'):
                corr_matrix_plot_name = 'CorrelationMatrixPlot.html'

                shutil.copy2(corr_matrix_plot_path,
                             os.path.join(output_directory, corr_matrix_plot_name))

                tab_content += '<iframe height="1300px" width="100%" '
                tab_content += 'src="{}" '.format(corr_matrix_plot_name)
                tab_content += 'style="border:none;"></iframe>\n<p></p>\n'
            else:
                raise ValueError('unexpected correlation matrix plot format:\n{}'.format(
                                                                            corr_matrix_plot_path))

            tab_content += """</div>"""

        if scatter_plot_path:

            tab_def_content += """
            <button class="tablinks" onclick="openTab(event, 'ScatterMatrixPlot')">Scatter Matrix Plot</button>
            """

            tab_content += """
            <div id="ScatterMatrixPlot" class="tabcontent">
            """

            scatter_plot_name = 'ScatterMatrixPlot.png'
            scatter_plot_display_name = 'Scatter Matrix Plot'

            shutil.copy2(scatter_plot_path,
                         os.path.join(output_directory, scatter_plot_name))

            tab_content += '<div class="gallery">'
            tab_content += '<a target="_blank" href="{}">'.format(scatter_plot_name)
            tab_content += '<img src="{}" '.format(scatter_plot_name)
            tab_content += 'alt="{}" width="600" height="400">'.format(
                                                                scatter_plot_display_name)
            tab_content += '</a><div class="desc">{}</div></div>'.format(
                                                                scatter_plot_display_name)

            tab_content += """</div>"""

        tab_def_content += """</div>"""

        return tab_def_content + tab_content

    def _generate_corr_html_report(self, corr_matrix_obj_ref, corr_matrix_plot_path,
                                   scatter_plot_path, df1, df2):

        """
        _generate_corr_html_report: generate html summary report for correlation
        """

        logging.info('Start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'compute_correlation_report.html')

        visualization_content = self._generate_visualization_content(
                                                                     output_directory,
                                                                     corr_matrix_obj_ref,
                                                                     corr_matrix_plot_path,
                                                                     scatter_plot_path, df1, df2)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'templates', 'corr_template.html'),
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

    def _generate_corr_report(self, corr_matrix_obj_ref, workspace_name, corr_matrix_plot_path,
                              scatter_plot_path=None, df1=None, df2=None):
        """
        _generate_report: generate summary report
        """
        logging.info('Start creating report')

        output_html_files = self._generate_corr_html_report(corr_matrix_obj_ref,
                                                            corr_matrix_plot_path,
                                                            scatter_plot_path, df1, df2)

        report_params = {'message': '',
                         'objects_created': [{'ref': corr_matrix_obj_ref,
                                              'description': 'Correlation Matrix'}],
                         'workspace_name': workspace_name,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 1400,
                         'report_object_name': 'compute_correlation_matrix_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _corr_for_matrix(self, input_obj_ref, method, dimension):
        """
        _corr_for_matrix: compute correlation matrix df for KBaseMatrices object
        """
        data_matrix = self.data_util.fetch_data({'obj_ref': input_obj_ref}).get('data_matrix')
        data_df = pd.read_json(data_matrix)
        data_df = data_df.reindex(index=natsorted(data_df.index))
        data_df = data_df.reindex(columns=natsorted(data_df.columns))
        data_df.fillna(0, inplace=True)

        corr_df = self.df_to_corr(data_df, method=method, dimension=dimension)

        return corr_df, data_df

    def _compute_significance(self, data_df, dimension):
        """
        _compute_significance: compute pairwsie significance dataframe
                               two-sided p-value for a hypothesis test
        """

        logging.info('Start computing significance matrix')
        if dimension == 'row':
            data_df = data_df.T

        data_df = data_df.dropna()._get_numeric_data()
        dfcols = pd.DataFrame(columns=data_df.columns)
        sig_df = dfcols.transpose().join(dfcols, how='outer')

        for r in data_df.columns:
            for c in data_df.columns:
                pvalue = stats.linregress(data_df[r], data_df[c])[3]
                sig_df[r][c] = round(pvalue, 4)

        return sig_df

    def _df_to_list(self, df, threshold=None):
        """
        _df_to_list: convert Dataframe to FloatMatrix2D matrix data
        """

        df.fillna(0, inplace=True)

        if threshold:
            drop_cols = list()
            for col in df.columns:
                if all(df[col] < threshold) and all(df[col] > -threshold):
                    drop_cols.append(col)
            df.drop(columns=drop_cols, inplace=True, errors='ignore')

            drop_idx = list()
            for idx in df.index:
                if all(df.loc[idx] < threshold) and all(df.loc[idx] > -threshold):
                    drop_idx.append(idx)
            df.drop(index=drop_idx, inplace=True, errors='ignore')

        matrix_data = {'row_ids': df.index.astype(str).tolist(),
                       'col_ids': df.columns.astype(str).tolist(),
                       'values': df.values.tolist()}

        return matrix_data

    def _save_corr_matrix(self, workspace_name, corr_matrix_name, corr_df, sig_df, method,
                          matrix_ref=None, corr_threshold=None):
        """
        _save_corr_matrix: save KBaseExperiments.CorrelationMatrix object
        """
        logging.info('Start saving CorrelationMatrix')

        if not isinstance(workspace_name, int):
            ws_name_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            ws_name_id = workspace_name

        corr_data = {}

        corr_data.update({'coefficient_data': self._df_to_list(corr_df,
                                                               threshold=corr_threshold)})
        corr_data.update({'correlation_parameters': {'method': method}})
        if matrix_ref:
            corr_data.update({'original_matrix_ref': matrix_ref})

        if sig_df is not None:
            corr_data.update({'significance_data': self._df_to_list(sig_df)})

        obj_type = 'KBaseExperiments.CorrelationMatrix'
        info = self.dfu.save_objects({
            "id": ws_name_id,
            "objects": [{
                "type": obj_type,
                "data": corr_data,
                "name": corr_matrix_name
            }]
        })[0]

        return "%s/%s/%s" % (info[6], info[0], info[4])

    def _Matrix2D_to_df(self, Matrix2D):
        """
        _Matrix2D_to_df: transform a FloatMatrix2D to data frame
        """

        index = Matrix2D.get('row_ids')
        columns = Matrix2D.get('col_ids')
        values = Matrix2D.get('values')

        df = pd.DataFrame(values, index=index, columns=columns)

        return df

    def _corr_to_df(self, corr_matrix_ref):
        """
        retrieve correlation matrix ws object to coefficient_df and significance_df
        """

        corr_data = self.dfu.get_objects({'object_refs': [corr_matrix_ref]})['data'][0]['data']

        coefficient_data = corr_data.get('coefficient_data')
        significance_data = corr_data.get('significance_data')

        coefficient_df = self._Matrix2D_to_df(coefficient_data)

        significance_df = None
        if significance_data:
            significance_df = self._Matrix2D_to_df(significance_data)

        return coefficient_df, significance_df

    def _corr_df_to_excel(self, coefficient_df, significance_df, result_dir, corr_matrix_ref):
        """
        write correlation matrix dfs into excel
        """

        corr_info = self.dfu.get_objects({'object_refs': [corr_matrix_ref]})['data'][0]['info']
        corr_name = corr_info[1]

        file_path = os.path.join(result_dir, corr_name + ".xlsx")

        writer = pd.ExcelWriter(file_path)

        coefficient_df.to_excel(writer, "coefficient_data", index=True)

        if significance_df is not None:
            significance_df.to_excel(writer, "significance_data", index=True)

        writer.close()

    def _fetch_matrix_data(self, matrix_ref):

        logging.info('start fectching matrix data')

        res = self.dfu.get_objects({'object_refs': [matrix_ref]})['data'][0]
        obj_type = res['info'][2]

        if "KBaseMatrices" in obj_type or 'KBaseProfile' in obj_type:
            data_matrix = self.data_util.fetch_data({'obj_ref': matrix_ref}).get('data_matrix')
            data_df = pd.read_json(data_matrix)
            data_df = data_df.reindex(index=natsorted(data_df.index))
            data_df = data_df.reindex(columns=natsorted(data_df.columns))

            return data_df
        else:
            err_msg = 'Ooops! [{}] is not supported.\n'.format(obj_type)
            err_msg += 'Please supply KBaseMatrices or KBaseProfile object'
            raise ValueError(err_msg)

    def _compute_metrices_corr(self, df1, df2, method, compute_significance):

        df1.fillna(0, inplace=True)
        df2.fillna(0, inplace=True)

        col_1 = df1.columns
        col_2 = df2.columns
        idx_1 = df1.index
        idx_2 = df2.index

        common_col = col_1.intersection(col_2)
        logging.info('matrices share [{}] common columns'.format(common_col.size))

        if common_col.empty:
            raise ValueError('Matrices share no common columns')

        logging.info('start trimming original matrix')
        df1 = df1.loc[:][common_col]
        df2 = df2.loc[:][common_col]

        corr_df = pd.DataFrame(index=idx_1, columns=idx_2)
        sig_df = pd.DataFrame(index=idx_1, columns=idx_2)

        logging.info('start calculating correlation matrix')
        logging.info('sizing {} x {}'.format(idx_1.size, idx_2.size))
        counter = 0
        for idx_value in idx_1:
            for col_value in idx_2:

                if counter % 100000 == 0:
                    logging.info('computed {} corr/sig values'.format(counter))

                value_array_1 = df1.loc[idx_value].tolist()
                value_array_2 = df2.loc[col_value].tolist()

                if method == 'pearson':
                    corr_value, p_value = stats.pearsonr(value_array_1, value_array_2)
                elif method == 'spearman':
                    corr_value, p_value = stats.spearmanr(value_array_1, value_array_2)
                elif method == 'kendall':
                    corr_value, p_value = stats.kendalltau(value_array_1, value_array_2)
                elif method == 'mutual_info':
                    corr_value = metrics.adjusted_mutual_info_score(value_array_1, value_array_2)
                    # p_value = stats.linregress(value_array_1, value_array_2)[3]
                    p_value = 0
                else:
                    err_msg = 'Input correlation method [{}] is not available.\n'.format(method)
                    err_msg += 'Please choose one of {}'.format(CORR_METHOD)
                    raise ValueError(err_msg)

                corr_df.at[idx_value, col_value] = round(corr_value, 4)
                if compute_significance:
                    sig_df.at[idx_value, col_value] = round(p_value, 4)

                counter += 1

        if not compute_significance:
            sig_df = None

        return corr_df, sig_df, df1, df2

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.scratch = config['scratch']

        self.data_util = DataUtil(config)
        self.dfu = DataFileUtil(self.callback_url)
        self.report_util = kb_GenericsReport(self.callback_url)

        plt.switch_backend('agg')

    def df_to_corr(self, df, method='pearson', dimension='col'):
        """
        Compute pairwise correlation of dimension (col or row)

        method: one of ['pearson', 'kendall', 'spearman']
        """

        logging.info('Computing correlation matrix')

        if method not in CORR_METHOD:
            err_msg = 'Input correlation method [{}] is not available.\n'.format(method)
            err_msg += 'Please choose one of {}'.format(CORR_METHOD)
            raise ValueError(err_msg)

        if dimension == 'row':
            df = df.T
        elif dimension != 'col':
            err_msg = 'Input dimension [{}] is not available.\n'.format(dimension)
            err_msg += 'Please choose either "col" or "row"'
            raise ValueError(err_msg)

        corr_df = df.corr(method=method).round(4)

        return corr_df

    def plotly_corr_matrix(self, corr_df):
        logging.info('Plotting matrix of correlation')

        result_dir = os.path.join(self.scratch, str(uuid.uuid4()) + '_corr_matrix_plots')
        self._mkdir_p(result_dir)

        try:
            trace = go.Heatmap(z=corr_df.values,
                               x=corr_df.columns,
                               y=corr_df.index)
            data = [trace]
        except Exception:
            err_msg = 'Running plotly_corr_matrix returned an error:\n{}\n'.format(
                                                                    traceback.format_exc())
            raise ValueError(err_msg)
        else:
            corr_matrix_plot_path = os.path.join(result_dir, 'corr_matrix_plots.html')
            logging.info('Saving plot to:\n{}'.format(corr_matrix_plot_path))
            plot(data, filename=corr_matrix_plot_path)

        return corr_matrix_plot_path

    def plot_corr_matrix(self, corr_df):
        """
        plot_corr_matrix: genreate correlation matrix plot
        """
        logging.info('Plotting matrix of correlation')

        result_dir = os.path.join(self.scratch, str(uuid.uuid4()) + '_corr_matrix_plots')
        self._mkdir_p(result_dir)

        try:
            plt.clf()
            matrix_size = corr_df.index.size
            figsize = 10 if matrix_size / 5 < 10 else matrix_size / 5
            fig, ax = plt.subplots(figsize=(figsize, figsize))
            cax = ax.matshow(corr_df)
            plt.xticks(list(range(len(corr_df.columns))), corr_df.columns, rotation='vertical',
                       fontstyle='italic')
            plt.yticks(list(range(len(corr_df.columns))), corr_df.columns, fontstyle='italic')
            plt.colorbar(cax)
        except Exception:
            err_msg = 'Running plot_corr_matrix returned an error:\n{}\n'.format(
                                                                    traceback.format_exc())
            raise ValueError(err_msg)
        else:
            corr_matrix_plot_path = os.path.join(result_dir, 'corr_matrix_plots.png')
            logging.info('Saving plot to:\n{}'.format(corr_matrix_plot_path))
            plt.savefig(corr_matrix_plot_path)

        return corr_matrix_plot_path

    def plot_scatter_matrix(self, df, dimension='col', alpha=0.2, diagonal='kde', figsize=(10, 10)):
        """
        plot_scatter_matrix: generate scatter plot for dimension (col or row)
                             ref: https://pandas.pydata.org/pandas-docs/stable/generated/pandas.plotting.scatter_matrix.html
        """
        logging.info('Plotting matrix of scatter')

        result_dir = os.path.join(self.scratch, str(uuid.uuid4()) + '_scatter_plots')
        self._mkdir_p(result_dir)

        if dimension == 'row':
            df = df.T
        elif dimension != 'col':
            err_msg = 'Input dimension [{}] is not available.\n'.format(dimension)
            err_msg += 'Please choose either "col" or "row"'
            raise ValueError(err_msg)

        try:
            plt.clf()
            sm = pd.plotting.scatter_matrix(df, alpha=alpha, diagonal=diagonal, figsize=figsize)

            # Change label rotation
            [s.xaxis.label.set_rotation(45) for s in sm.reshape(-1)]
            [s.yaxis.label.set_rotation(45) for s in sm.reshape(-1)]

            # # May need to offset label when rotating to prevent overlap of figure
            [s.get_yaxis().set_label_coords(-1.5, 0.5) for s in sm.reshape(-1)]

            # Hide all ticks
            [s.set_xticks(()) for s in sm.reshape(-1)]
            [s.set_yticks(()) for s in sm.reshape(-1)]
        except Exception:
            err_msg = 'Running scatter_matrix returned an error:\n{}\n'.format(
                                                                    traceback.format_exc())
            raise ValueError(err_msg)
        else:
            scatter_plot_path = os.path.join(result_dir, 'scatter_plots.png')
            logging.info('Saving plot to:\n{}'.format(scatter_plot_path))
            plt.savefig(scatter_plot_path)

        return scatter_plot_path

    def compute_correlation_across_matrices(self, params):
        """
        matrix_ref_1: object reference of a matrix
        matrix_ref_2: object reference of a matrix
        workspace_name: workspace name objects to be saved to
        corr_matrix_name: correlation matrix object name
        method: correlation method, one of ['pearson', 'kendall', 'spearman', 'mutual_info']
        plot_corr_matrix: plot correlation matrix in report, default False
        compute_significance: also compute Significance in addition to correlation matrix
        """

        logging.info('--->\nrunning CorrelationUtil.compute_correlation_across_matrices\n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        self._validate_compute_correlation_across_matrices_params(params)

        matrix_ref_1 = params.get('matrix_ref_1')
        matrix_ref_2 = params.get('matrix_ref_2')
        workspace_name = params.get('workspace_name')
        corr_matrix_name = params.get('corr_matrix_name')
        corr_threshold = params.get('corr_threshold')
        dimension = params.get('dimension', 'row')

        method = params.get('method', 'pearson')
        if method not in CORR_METHOD:
            err_msg = 'Input correlation method [{}] is not available.\n'.format(method)
            err_msg += 'Please choose one of {}'.format(CORR_METHOD)
            raise ValueError(err_msg)
        plot_corr_matrix = params.get('plot_corr_matrix', False)
        compute_significance = params.get('compute_significance', True)

        matrix_1_type = self.dfu.get_objects({'object_refs': [matrix_ref_1]})['data'][0]['info'][2]

        # making sure otu_ids are on the column of table
        if "AmpliconMatrix" in matrix_1_type:
            matrix_ref_1, matrix_ref_2 = matrix_ref_2, matrix_ref_1

        df1 = self._fetch_matrix_data(matrix_ref_1)
        df2 = self._fetch_matrix_data(matrix_ref_2)

        if dimension == 'col':
            df1 = df1.T
            df2 = df2.T

        corr_df, sig_df, df1, df2 = self._compute_metrices_corr(df1, df2, method, compute_significance)

        if plot_corr_matrix:
            corr_matrix_plot_path = self.plotly_corr_matrix(corr_df)
        else:
            corr_matrix_plot_path = None

        corr_matrix_obj_ref = self._save_corr_matrix(workspace_name, corr_matrix_name, corr_df,
                                                     sig_df, method,
                                                     matrix_ref=[matrix_ref_1, matrix_ref_2],
                                                     corr_threshold=corr_threshold)

        returnVal = {'corr_matrix_obj_ref': corr_matrix_obj_ref}

        df1.index = df1.index.astype('str')
        df1.columns = df1.columns.astype('str')
        df2.index = df2.index.astype('str')
        df2.columns = df2.columns.astype('str')
        report_output = self._generate_corr_report(corr_matrix_obj_ref, workspace_name,
                                                   corr_matrix_plot_path, df1=df1, df2=df2)

        returnVal.update(report_output)

        return returnVal

    def compute_correlation_matrix(self, params):
        """
        input_obj_ref: object reference of a matrix
        workspace_name: workspace name objects to be saved to
        dimension: compute correlation on column or row, one of ['col', 'row']
        corr_matrix_name: correlation matrix object name
        method: correlation method, one of ['pearson', 'kendall', 'spearman']
        compute_significance: compute pairwise significance value, default False
        plot_corr_matrix: plot correlation matrix in repor, default False
        plot_scatter_matrix: plot scatter matrix in report, default False
        """

        logging.info('--->\nrunning CorrelationUtil.compute_correlation_matrix\n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        self._validate_compute_corr_matrix_params(params)

        input_obj_ref = params.get('input_obj_ref')
        workspace_name = params.get('workspace_name')
        corr_matrix_name = params.get('corr_matrix_name')

        method = params.get('method', 'pearson')
        dimension = params.get('dimension', 'row')
        plot_corr_matrix = params.get('plot_corr_matrix', False)
        plot_scatter_matrix = params.get('plot_scatter_matrix', False)
        compute_significance = params.get('compute_significance', False)

        res = self.dfu.get_objects({'object_refs': [input_obj_ref]})['data'][0]
        obj_type = res['info'][2]

        if "KBaseMatrices" in obj_type or 'KBaseProfile' in obj_type:
            corr_df, data_df = self._corr_for_matrix(input_obj_ref, method, dimension)
            sig_df = None
            if compute_significance:
                sig_df = self._compute_significance(data_df, dimension)
        else:
            err_msg = 'Ooops! [{}] is not supported.\n'.format(obj_type)
            err_msg += 'Please supply KBaseMatrices or KBaseProfile object'
            raise ValueError(err_msg)

        if plot_corr_matrix:
            corr_matrix_plot_path = self.plotly_corr_matrix(corr_df)
        else:
            corr_matrix_plot_path = None

        if plot_scatter_matrix:
            scatter_plot_path = self.plot_scatter_matrix(data_df, dimension=dimension)
        else:
            scatter_plot_path = None

        corr_matrix_obj_ref = self._save_corr_matrix(workspace_name, corr_matrix_name, corr_df,
                                                     sig_df, method, matrix_ref=[input_obj_ref])

        returnVal = {'corr_matrix_obj_ref': corr_matrix_obj_ref}

        data_df.index = data_df.index.astype('str')
        data_df.columns = data_df.columns.astype('str')
        if dimension == 'col':
            data_df = data_df.T
        report_output = self._generate_corr_report(corr_matrix_obj_ref, workspace_name,
                                                   corr_matrix_plot_path,
                                                   scatter_plot_path=scatter_plot_path,
                                                   df1=data_df, df2=data_df)

        returnVal.update(report_output)

        return returnVal

    def export_corr_matrix_excel(self, params):
        """
        export CorrelationMatrix as Excel
        """

        corr_matrix_ref = params.get('input_ref')

        coefficient_df, significance_df = self._corr_to_df(corr_matrix_ref)

        result_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_dir)

        self._corr_df_to_excel(coefficient_df, significance_df, result_dir, corr_matrix_ref)

        package_details = self.dfu.package_for_download({
            'file_path': result_dir,
            'ws_refs': [corr_matrix_ref]
        })

        return {'shock_id': package_details['shock_id']}
