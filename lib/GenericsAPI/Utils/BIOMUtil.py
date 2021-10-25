import errno
import logging
import os
import uuid
import biom
import pandas as pd
from Bio import SeqIO
import shutil

from installed_clients.DataFileUtilClient import DataFileUtil
from GenericsAPI.Utils.AttributeUtils import AttributesUtil
from GenericsAPI.Utils.SampleServiceUtil import SampleServiceUtil
from GenericsAPI.Utils.DataUtil import DataUtil
from GenericsAPI.Utils.MatrixUtil import MatrixUtil
from GenericsAPI.Utils.TaxonUtil import TaxonUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.KBaseSearchEngineClient import KBaseSearchEngine
from installed_clients.kb_GenericsReportClient import kb_GenericsReport

TYPE_ATTRIBUTES = {'description', 'scale', 'row_normalization', 'col_normalization'}
SCALE_TYPES = {'raw', 'ln', 'log2', 'log10'}
DEFAULT_META_KEYS = ["lineage", "score", "taxonomy_source", "species_name",
                     "consensus_sequence"]
TARGET_GENE_SUBFRAGMENT_MAP = {'16S': ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9'],
                               '18S': ['V1', 'V2', 'V3', 'V4', 'V9'],
                               'ITS': ['ITS1', 'ITS2']}
SEQ_INSTRUMENTS_MAP = {'Applied Biosystems': ['AB 310 Genetic Analyzer',
                                              'AB 3130 Genetic Analyzer',
                                              'AB 3130xL Genetic Analyzer',
                                              'AB 3500 Genetic Analyzer',
                                              'AB 3500xL Genetic Analyzer',
                                              'AB 3730 Genetic Analyzer',
                                              'AB 3730xL Genetic Analyzer',
                                              'AB 5500xl Genetic Analyzer',
                                              'AB 5500x-Wl Genetic Analyzer',
                                              'AB SOLiD System',
                                              'AB SOLiD System 2.0',
                                              'AB SOLiD System 3.0',
                                              'AB SOLiD 3 Plus System',
                                              'AB SOLiD 4 System',
                                              'AB SOLiD 4hq System',
                                              'AB SOLiD PI System'],
                       'Roche 454': ['454 GS', '454 GS 20', '454 GS FLX', '454 GS FLX+',
                                     '454 GS FLX Titanium'],
                       'Life Sciences': ['454 GS Junior'],
                       'Illumina': ['Illumina Genome Analyzer',
                                    'Illumina Genome Analyzer II',
                                    'Illumina Genome Analyzer IIx',
                                    'Illumina HiScanSQ',
                                    'Illumina HiSeq 1000',
                                    'Illumina HiSeq 1500',
                                    'Illumina HiSeq 2000',
                                    'Illumina HiSeq 2500',
                                    'Illumina HiSeq 3000',
                                    'Illumina HiSeq 4000',
                                    'Illumina HiSeq X',
                                    'HiSeq X Five',
                                    'HiSeq X Ten',
                                    'Illumina iSeq 100',
                                    'Illumina MiSeq',
                                    'Illumina MiniSeq',
                                    'NextSeq 500',
                                    'NextSeq 550',
                                    'NextSeq 1000',
                                    'NextSeq 2000',
                                    'Illumina NovaSeq 6000'],
                       'ThermoFisher': ['Ion Torrent PGM', 'Ion Torrent Proton',
                                        'Ion Torrent S5 XL', 'Ion Torrent S5'],
                       'Pacific Biosciences': ['PacBio RS', 'PacBio RS II', 'PacBio Sequel',
                                               'PacBio Sequel II'],
                       'Oxford Nanopore': ['MinION', 'GridION', 'PromethION'],
                       'BGI Group': ['BGISEQ-500', 'DNBSEQ-G400', 'DNBSEQ-T7', 'DNBSEQ-G50',
                                     'MGISEQ-2000RS']}


class BiomUtil:

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

    def _process_params(self, params):
        logging.info('start validating import_matrix_from_biom params')

        # check for required parameters
        for p in ['obj_type', 'matrix_name', 'workspace_id', 'scale', 'amplicon_type',
                  'sequencing_technology', 'sequencing_instrument',
                  'target_gene', 'target_subfragment', 'taxon_calling']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

        # check sequencing_technology and sequencing_instrument matching
        # sequencing_technology = params.get('sequencing_technology')
        # sequencing_instrument = params.get('sequencing_instrument')

        # check target_gene and target_subfragment matching
        target_gene = params.get('target_gene')
        target_subfragment = params.get('target_subfragment')

        if target_gene not in TARGET_GENE_SUBFRAGMENT_MAP:
            raise ValueError('Unexpected target gene: {}'.format(target_gene))
        expected_subfragments = TARGET_GENE_SUBFRAGMENT_MAP.get(target_gene)
        if not set(target_subfragment) <= set(expected_subfragments):
            raise ValueError('Please select target subfragments among {} for {}'.format(
                expected_subfragments, target_gene))

        # check taxon_calling
        taxon_calling = params.get('taxon_calling')
        taxon_calling_method = list(set(taxon_calling.get('taxon_calling_method')))
        params['taxon_calling_method'] = taxon_calling_method

        if 'denoising' in taxon_calling_method:
            denoise_method = taxon_calling.get('denoise_method')
            sequence_error_cutoff = taxon_calling.get('sequence_error_cutoff')

            if not (denoise_method and sequence_error_cutoff):
                raise ValueError('Please provide denoise_method and sequence_error_cutoff')

            params['denoise_method'] = denoise_method
            params['sequence_error_cutoff'] = sequence_error_cutoff

        if 'clustering' in taxon_calling_method:
            clustering_method = taxon_calling.get('clustering_method')
            clustering_cutoff = taxon_calling.get('clustering_cutoff')

            if not (clustering_method and clustering_cutoff):
                raise ValueError('Please provide clustering_method and clustering_cutoff')

            params['clustering_method'] = clustering_method
            params['clustering_cutoff'] = clustering_cutoff

        obj_type = params.get('obj_type')
        if obj_type not in self.matrix_types:
            raise ValueError('Unknown matrix object type: {}'.format(obj_type))

        scale = params.get('scale')
        if scale not in SCALE_TYPES:
            raise ValueError('Unknown scale type: {}'.format(scale))

        biom_file = None
        tsv_file = None
        fasta_file = None
        metadata_keys = DEFAULT_META_KEYS

        input_local_file = params.get('input_local_file', False)

        if params.get('taxonomic_abundance_tsv') and params.get('taxonomic_fasta'):
            tsv_file = params.get('taxonomic_abundance_tsv')
            fasta_file = params.get('taxonomic_fasta')

            if not (tsv_file and fasta_file):
                raise ValueError('missing TSV or FASTA file')

            if not input_local_file:
                tsv_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': tsv_file}).get('copy_file_path')

                fasta_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': fasta_file}).get('copy_file_path')

            metadata_keys_str = params.get('metadata_keys')
            if metadata_keys_str:
                metadata_keys += [x.strip() for x in metadata_keys_str.split(',')]
            mode = 'tsv_fasta'
        elif params.get('biom_fasta'):
            biom_fasta = params.get('biom_fasta')
            biom_file = biom_fasta.get('biom_file_biom_fasta')
            fasta_file = biom_fasta.get('fasta_file_biom_fasta')

            if not (biom_file and fasta_file):
                raise ValueError('missing BIOM or FASTA file')

            if not input_local_file:
                biom_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': biom_file}).get('copy_file_path')

                fasta_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': fasta_file}).get('copy_file_path')
            mode = 'biom_fasta'
        elif params.get('tsv_fasta'):
            tsv_fasta = params.get('tsv_fasta')
            tsv_file = tsv_fasta.get('tsv_file_tsv_fasta')
            fasta_file = tsv_fasta.get('fasta_file_tsv_fasta')

            if not (tsv_file and fasta_file):
                raise ValueError('missing TSV or FASTA file')

            if not input_local_file:
                tsv_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': tsv_file}).get('copy_file_path')

                fasta_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': fasta_file}).get('copy_file_path')

            metadata_keys_str = tsv_fasta.get('metadata_keys_tsv_fasta')
            if metadata_keys_str:
                metadata_keys += [x.strip() for x in metadata_keys_str.split(',')]
            mode = 'tsv_fasta'
        else:
            raise ValueError('missing valide file group type in parameters')

        return (biom_file, tsv_file, fasta_file, mode, list(set(metadata_keys)))

    def _validate_fasta_file(self, df, fasta_file):
        logging.info('start validating FASTA file')
        try:
            fastq_dict = SeqIO.index(fasta_file, "fasta")
        except Exception:
            raise ValueError('Cannot parse file. Please provide valide FASTA file')

        matrix_ids = df.index
        file_ids = fastq_dict.keys()

        unmatched_ids = set(matrix_ids) - set(file_ids)

        if unmatched_ids:
            raise ValueError('FASTA file does not have [{}] OTU id'.format(unmatched_ids))

    def _file_to_amplicon_data(self, biom_file, tsv_file, fasta_file, mode, refs, matrix_name,
                               workspace_id, scale, description, metadata_keys=None):

        amplicon_data = refs

        if mode.startswith('biom'):
            logging.info('start parsing BIOM file for matrix data')
            table = biom.load_table(biom_file)
            observation_metadata = table._observation_metadata
            sample_metadata = table._sample_metadata

            matrix_data = {'row_ids': table._observation_ids.tolist(),
                           'col_ids': table._sample_ids.tolist(),
                           'values': table.matrix_data.toarray().tolist()}

            logging.info('start building attribute mapping object')
            amplicon_data.update(self.get_attribute_mapping("row", observation_metadata,
                                                            matrix_data, matrix_name, refs,
                                                            workspace_id))
            amplicon_data.update(self.get_attribute_mapping("col", sample_metadata,
                                                            matrix_data, matrix_name, refs,
                                                            workspace_id))

            amplicon_data['attributes'] = {}
            for k in ('create_date', 'generated_by'):
                val = getattr(table, k)
                if not val:
                    continue
                if isinstance(val, bytes):
                    amplicon_data['attributes'][k] = val.decode('utf-8')
                else:
                    amplicon_data['attributes'][k] = str(val)
        elif mode.startswith('tsv'):
            observation_metadata = None
            sample_metadata = None
            try:
                logging.info('start parsing TSV file for matrix data')
                reader = pd.read_csv(tsv_file, sep=None, iterator=True)
                inferred_sep = reader._engine.data.dialect.delimiter
                df = pd.read_csv(tsv_file, sep=inferred_sep, index_col=0)
            except Exception:
                raise ValueError('Cannot parse file. Please provide valide tsv file')
            else:
                self._validate_fasta_file(df, fasta_file)
                metadata_df = None
                if metadata_keys:
                    shared_metadata_keys = list(set(metadata_keys) & set(df.columns))
                    if mode == 'tsv' and 'consensus_sequence' not in shared_metadata_keys:
                        raise ValueError('TSV file does not include consensus_sequence')
                    if shared_metadata_keys:
                        metadata_df = df[shared_metadata_keys]
                        df.drop(columns=shared_metadata_keys, inplace=True)
                try:
                    df = df.astype(float)
                except ValueError:
                    err_msg = 'Found some non-float values. Matrix contains only numeric values\n'
                    err_msg += 'Please list any non-numeric column names in  Metadata Keys field'
                    raise ValueError(err_msg)
                df.fillna(0, inplace=True)
                df.index = df.index.astype('str')
                df.columns = df.columns.astype('str')
                matrix_data = {'row_ids': df.index.tolist(),
                               'col_ids': df.columns.tolist(),
                               'values': df.values.tolist()}

            logging.info('start building attribute mapping object')
            amplicon_data.update(self.get_attribute_mapping("row", observation_metadata,
                                                            matrix_data, matrix_name, refs,
                                                            workspace_id, metadata_df=metadata_df))
            amplicon_data.update(self.get_attribute_mapping("col", sample_metadata,
                                                            matrix_data, matrix_name, refs,
                                                            workspace_id))

            amplicon_data['attributes'] = {}
        else:
            raise ValueError('error parsing _file_to_amplicon_data, mode: {}'.format(mode))

        amplicon_data.update({'data': matrix_data})

        amplicon_data['search_attributes'] = [f'{k}|{v}' for k, v in amplicon_data['attributes'].items()]

        amplicon_data['scale'] = scale
        if description:
            amplicon_data['description'] = description

        return amplicon_data

    def get_attribute_mapping(self, axis, metadata, matrix_data, matrix_name, refs,  workspace_id,
                              metadata_df=None):
        mapping_data = {}
        axis_ids = matrix_data[f'{axis}_ids']
        if refs.get('sample_set_ref') and axis == 'col':
            name = matrix_name + "_{}_attributes".format(axis)
            mapping_data[f'{axis}_attributemapping_ref'] = self._sample_set_to_attribute_mapping(
                axis_ids, refs.get('sample_set_ref'), name, workspace_id)
            mapping_data[f'{axis}_mapping'] = {x: x for x in axis_ids}
        elif refs.get(f'{axis}_attributemapping_ref'):
            am_data = self.dfu.get_objects(
                {'object_refs': [refs[f'{axis}_attributemapping_ref']]}
            )['data'][0]['data']
            unmatched_ids = set(axis_ids) - set(am_data['instances'].keys())
            if unmatched_ids:
                name = "Column" if axis == 'col' else "Row"
                raise ValueError(f"The following {name} IDs from the uploaded matrix do not match "
                                 f"the supplied {name} attribute mapping: {', '.join(unmatched_ids)}"
                                 f"\nPlease verify the input data or upload an excel file with a"
                                 f"{name} mapping tab.")
            else:
                mapping_data[f'{axis}_mapping'] = {x: x for x in axis_ids}
        elif metadata:
            name = matrix_name + "_{}_attributes".format(axis)
            mapping_data[f'{axis}_attributemapping_ref'] = self._metadata_to_attribute_mapping(
                axis_ids, metadata, name, workspace_id)
            # if coming from biom file, metadata and axis IDs are guaranteed to match
            mapping_data[f'{axis}_mapping'] = {x: x for x in axis_ids}
        elif metadata_df is not None:
            name = matrix_name + "_{}_attributes".format(axis)
            mapping_data[f'{axis}_attributemapping_ref'] = self._meta_df_to_attribute_mapping(
                axis_ids, metadata_df, name, workspace_id)
            mapping_data[f'{axis}_mapping'] = {x: x for x in axis_ids}

        return mapping_data

    def _meta_df_to_attribute_mapping(self, axis_ids, metadata_df, obj_name, ws_id):
        data = {'ontology_mapping_method': "TSV file", 'instances': {}}
        metadata_df = metadata_df.astype(str)
        attribute_keys = metadata_df.columns.tolist()
        data['attributes'] = [{'attribute': key, 'source': 'upload'} for key in attribute_keys]

        if 'taxonomy' in attribute_keys:
            data['attributes'].append({'attribute': 'parsed_user_taxonomy', 'source': 'upload'})

        for axis_id in axis_ids:
            data['instances'][axis_id] = metadata_df.loc[axis_id].tolist()
            if 'taxonomy' in attribute_keys:
                parsed_user_taxonomy = None
                taxonomy_index = attribute_keys.index('taxonomy')
                taxonomy_str = metadata_df.loc[axis_id].tolist()[taxonomy_index]
                parsed_user_taxonomy = self.taxon_util.process_taxonomic_str(taxonomy_str)
                data['instances'][axis_id].append(parsed_user_taxonomy)

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

    def _metadata_to_attribute_mapping(self, instances, metadata, obj_name, ws_id):
        data = {'ontology_mapping_method': "BIOM file", 'instances': {}}
        sample_set = metadata[0:min(len(metadata), 25)]
        metadata_keys = sorted(set((k for m_dict in sample_set for k in m_dict)))
        data['attributes'] = [{'attribute': key, 'source': 'upload'} for key in metadata_keys]
        for inst, meta in zip(instances, metadata):
            data['instances'][inst] = [str(meta[attr]) for attr in metadata_keys]

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

    def _generate_visualization_content(self, output_directory, heatmap_dir, data_df,
                                        top_heatmap_dir, top_percent, display_count):

        row_data_summary = data_df.T.describe().round(2).to_string()
        col_data_summary = data_df.describe().round(2).to_string()

        tab_def_content = ''
        tab_content = ''

        viewer_name = 'data_summary'
        tab_def_content += '''\n<div class="tab">\n'''
        tab_def_content += '''\n<button class="tablinks" '''
        tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
        tab_def_content += ''' id="defaultOpen"'''
        tab_def_content += '''>Matrix Statistics</button>\n'''

        tab_content += '''\n<div id="{}" class="tabcontent" style="overflow:auto">'''.format(viewer_name)
        tab_content += '''\n<h5>Amplicon Matrix Size: {} x {}</h5>'''.format(len(data_df.index),
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
            tab_def_content += '''>Top {}% ({} Rows) Heatmap</button>\n'''.format(
                                                                            round(top_percent, 2),
                                                                            display_count)

            heatmap_report_files = os.listdir(top_heatmap_dir)

            heatmap_index_page = None
            for heatmap_report_file in heatmap_report_files:
                if heatmap_report_file.endswith('.html'):
                    heatmap_index_page = heatmap_report_file

                shutil.copy2(os.path.join(top_heatmap_dir, heatmap_report_file),
                             output_directory)

            if heatmap_index_page:
                tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
                msg = 'Top {} percent of matrix sorted by sum of abundance values.'.format(
                                                                            round(top_percent, 2))
                tab_content += '''<p style="color:red;" >{}</p>'''.format(msg)

                tab_content += '\n<iframe height="1300px" width="100%" '
                tab_content += 'src="{}" '.format(heatmap_index_page)
                tab_content += 'style="border:none;"></iframe>'
                tab_content += '\n</div>\n'
            else:
                tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
                tab_content += '''\n<p style="color:red;" >'''
                tab_content += '''Heatmap is too large to be displayed.</p>\n'''
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
            tab_content += '\n<iframe height="1300px" width="100%" '
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

    def _generate_heatmap_html_report(self, data):

        logging.info('Start generating heatmap report page')

        data_df = pd.DataFrame(data['values'], index=data['row_ids'], columns=data['col_ids'])
        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)
        tsv_file_path = os.path.join(result_directory, 'heatmap_data_{}.tsv'.format(
                                                                            str(uuid.uuid4())))
        data_df.to_csv(tsv_file_path)

        if data_df.index.size < 10000:
            heatmap_dir = self.report_util.build_heatmap_html({
                                                        'tsv_file_path': tsv_file_path,
                                                        'cluster_data': True})['html_dir']
        else:
            logging.info('Original matrix is too large. Skip clustering data in report.')
            heatmap_dir = self.report_util.build_heatmap_html({
                                                        'tsv_file_path': tsv_file_path,
                                                        'cluster_data': False})['html_dir']
        top_heatmap_dir = None
        top_percent = 100
        display_count = 200  # roughly count for display items
        if len(data_df.index) > 1000:
            top_percent = min(display_count / data_df.index.size * 100, 100)
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
                                                                     top_percent,
                                                                     display_count)

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
                            'description': 'HTML summary report for Import Amplicon Matrix App'
                            })
        return html_report

    def _generate_report(self, matrix_obj_ref, new_row_attr_ref,
                         new_col_attr_ref, workspace_id, data=None):
        """
        _generate_report: generate summary report
        """

        objects_created = [{'ref': matrix_obj_ref, 'description': 'Imported Amplicon Matrix'}]

        if new_row_attr_ref:
            objects_created.append({'ref': new_row_attr_ref,
                                    'description': 'Imported Amplicons(Row) Attribute Mapping'})

        if new_col_attr_ref:
            objects_created.append({'ref': new_col_attr_ref,
                                    'description': 'Imported Samples(Column) Attribute Mapping'})

        if data:
            output_html_files = self._generate_heatmap_html_report(data)

            report_params = {'message': '',
                             'objects_created': objects_created,
                             'workspace_id': workspace_id,
                             'html_links': output_html_files,
                             'direct_html_link_index': 0,
                             'html_window_height': 1400,
                             'report_object_name': 'import_matrix_from_biom_' + str(uuid.uuid4())}

        else:
            report_params = {'message': '',
                             'objects_created': objects_created,
                             'workspace_id': workspace_id,
                             'report_object_name': 'import_matrix_from_biom_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url, token=self.token)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.dfu = DataFileUtil(self.callback_url)
        self.report_util = kb_GenericsReport(self.callback_url)
        self.data_util = DataUtil(config)
        self.sampleservice_util = SampleServiceUtil(config)
        self.attr_util = AttributesUtil(config)
        self.matrix_util = MatrixUtil(config)
        self.taxon_util = TaxonUtil(config)
        self.matrix_types = [x.split(".")[1].split('-')[0]
                             for x in self.data_util.list_generic_types()]
        self.taxon_wsname = config['taxon-workspace-name']
        self.kbse = KBaseSearchEngine(config['search-url'])
        self.taxon_cache = dict()

    def fetch_sequence(self, matrix_ref):
        logging.info('start to fetch consensus sequence')

        input_matrix_obj = self.dfu.get_objects({'object_refs': [matrix_ref]})['data'][0]
        input_matrix_info = input_matrix_obj['info']
        matrix_name = input_matrix_info[1]
        matrix_type = input_matrix_info[2]
        matrix_data = input_matrix_obj['data']

        if 'KBaseMatrices.AmpliconMatrix' not in matrix_type:
            raise ValueError('Unexpected data type: {}'.format(matrix_type))

        handle = matrix_data.get('sequencing_file_handle')
        if not handle:
            raise ValueError('Missing sequencing_file_handle from the matrix object')

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        logging.info('Start generating consensus sequence file in {}'.format(output_directory))
        self._mkdir_p(output_directory)

        matrix_fasta_file = self.dfu.shock_to_file({'handle_id': handle,
                                                    'file_path': self.scratch}).get('file_path')

        try:
            logging.info('start parsing FASTA file')
            fastq_dict = SeqIO.index(matrix_fasta_file, "fasta")
        except Exception:
            raise ValueError('Cannot parse file. Please provide valide FASTA file')

        row_ids = matrix_data['data']['row_ids']

        fasta_file_path = os.path.join(output_directory, matrix_name + 'consensus_sequence.fasta')

        with open(fasta_file_path, 'w') as f:
            for row_id in row_ids:
                consensus_sequence = str(fastq_dict.get(row_id).seq)
                f.write('>' + str(row_id) + '\n')
                f.write(consensus_sequence + '\n')

        return fasta_file_path

    def import_matrix_from_biom(self, params):
        """
        arguments:
        obj_type: one of ExpressionMatrix, FitnessMatrix, DifferentialExpressionMatrix
        matrix_name: matrix object name
        workspace_id: workspace id matrix object to be saved to
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

        (biom_file, tsv_file, fasta_file, mode, metadata_keys) = self._process_params(params)

        workspace_id = params.get('workspace_id')
        matrix_name = params.get('matrix_name')
        obj_type = params.get('obj_type')
        scale = params.get('scale')
        description = params.get('description')
        refs = {k: v for k, v in params.items() if "_ref" in k}

        amplicon_data = self._file_to_amplicon_data(biom_file, tsv_file, fasta_file, mode,
                                                    refs, matrix_name,
                                                    workspace_id, scale, description, metadata_keys)

        for key in ['amplicon_type', 'amplification', 'extraction',
                    'target_gene', 'target_subfragment', 'pcr_primers',
                    'library_kit', 'library_layout', 'library_screening_strategy',
                    'sequencing_center', 'sequencing_date',
                    'sequencing_technology', 'sequencing_instrument',
                    'sequencing_quality_filter_cutoff',
                    'read_length_cutoff', 'read_pairing',
                    'barcode_error_rate', 'chimera_detection_and_removal',
                    'taxon_calling_method',
                    'denoise_method', 'sequence_error_cutoff',
                    'clustering_method', 'clustering_cutoff',
                    'sample_set_ref', 'reads_set_ref']:
            if params.get(key):
                amplicon_data[key] = params[key]

        new_row_attr_ref = None
        if not params.get('row_attributemapping_ref'):
            new_row_attr_ref = amplicon_data.get('row_attributemapping_ref')

        new_col_attr_ref = None
        if not params.get('col_attributemapping_ref'):
            new_col_attr_ref = amplicon_data.get('col_attributemapping_ref')

        if fasta_file:
            logging.info('start saving consensus sequence file to shock: {}'.format(fasta_file))
            handle_id = self.dfu.file_to_shock({'file_path': fasta_file,
                                                'make_handle': True})['handle']['hid']
            amplicon_data['sequencing_file_handle'] = handle_id

        logging.info('start saving Matrix object: {}'.format(matrix_name))
        matrix_obj_ref = self.data_util.save_object({
                                                'obj_type': 'KBaseMatrices.{}'.format(obj_type),
                                                'obj_name': matrix_name,
                                                'data': amplicon_data,
                                                'workspace_id': workspace_id})['obj_ref']

        if params.get('sample_set_ref'):
            self.matrix_util._link_matrix_to_samples(matrix_obj_ref, amplicon_data, params['sample_set_ref'])

        returnVal = {'matrix_obj_ref': matrix_obj_ref}

        report_output = self._generate_report(matrix_obj_ref,
                                              new_row_attr_ref, new_col_attr_ref, workspace_id,
                                              data=amplicon_data['data'])

        returnVal.update(report_output)

        return returnVal
