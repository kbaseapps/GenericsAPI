import errno
import logging
import os
import uuid

import biom
import pandas as pd
from Bio import SeqIO
import csv
import shutil
import traceback
import sys
from plotly.offline import plot
import plotly.express as px

from installed_clients.DataFileUtilClient import DataFileUtil
from GenericsAPI.Utils.AttributeUtils import AttributesUtil
from GenericsAPI.Utils.SampleServiceUtil import SampleServiceUtil
from GenericsAPI.Utils.DataUtil import DataUtil
from GenericsAPI.Utils.MatrixUtil import MatrixUtil
from GenericsAPI.Utils.TaxonUtil import TaxonUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.KBaseSearchEngineClient import KBaseSearchEngine
from installed_clients.kb_GenericsReportClient import kb_GenericsReport
from installed_clients.WorkspaceClient import Workspace as workspaceService

TYPE_ATTRIBUTES = {'description', 'scale', 'row_normalization', 'col_normalization'}
SCALE_TYPES = {'raw', 'ln', 'log2', 'log10'}
DEFAULT_META_KEYS = ["lineage", "score", "taxonomy_source", "species_name",
                     "consensus_sequence"]


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
        for p in ['obj_type', 'matrix_name', 'workspace_id', 'scale',
                  'amplicon_type', 'target_gene_region', 'forward_primer_sequence',
                  'reverse_primer_sequence', 'sequencing_platform',
                  'sequencing_quality_filter_cutoff', 'clustering_cutoff', 'clustering_method']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

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
        elif params.get('biom_tsv'):
            biom_tsv = params.get('biom_tsv')
            biom_file = biom_tsv.get('biom_file_biom_tsv')
            tsv_file = biom_tsv.get('tsv_file_biom_tsv')

            if not (biom_file and tsv_file):
                raise ValueError('missing BIOM or TSV file')

            if not input_local_file:
                biom_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': biom_file}).get('copy_file_path')

                tsv_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': tsv_file}).get('copy_file_path')
            mode = 'biom_tsv'
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
        elif params.get('tsv'):
            tsv = params.get('tsv')
            tsv_file = tsv.get('tsv_file_tsv')

            if not tsv_file:
                raise ValueError('missing TSV file')

            if not input_local_file:
                tsv_file = self.dfu.download_staging_file(
                                    {'staging_file_subdir_path': tsv_file}).get('copy_file_path')

            metadata_keys_str = tsv.get('metadata_keys_tsv')
            if metadata_keys_str:
                metadata_keys += [x.strip() for x in metadata_keys_str.split(',')]

            mode = 'tsv'
        else:
            raise ValueError('missing valide file group type in parameters')

        return (biom_file, tsv_file, fasta_file, mode, list(set(metadata_keys)))

    def _retrieve_value(self, biom_metadata_dict, tsv_metadata_df, key, required=False):

        if key in biom_metadata_dict:
            return {k.lower(): v for k, v in biom_metadata_dict.items()}.get(key)
        elif key in tsv_metadata_df:
            return {k.lower(): v for k, v in tsv_metadata_df.items()}.get(key)
        elif required:
            raise ValueError('missing necessary [{}] from file'.format(key))
        else:
            return None

    def _search_taxon(self, scientific_name):
        """
        logic borrowed from: GFU.GenomeInterface
        https://github.com/kbaseapps/GenomeFileUtil/blob/master/lib/GenomeFileUtil/core/GenomeInterface.py#L216
        """
        taxon_id = None

        if scientific_name in self.taxon_cache:
            taxon_id = self.taxon_cache.get(scientific_name)
        else:
            search_params = {
                "object_types": ["taxon"],
                "match_filter": {
                    "lookup_in_keys": {
                        "scientific_name": {"value": scientific_name}},
                    "exclude_subobjects": 1
                },
                "access_filter": {
                    "with_private": 0,
                    "with_public": 1
                },
                "sorting_rules": [{
                    "is_object_property": 0,
                    "property": "timestamp",
                    "ascending": 0
                }]
            }

            try:
                objects = self.kbse.search_objects(search_params)['objects']

                if not objects:
                    search_params['match_filter']['lookup_in_keys'] = {
                        "aliases": {"value": scientific_name}
                    }
                    objects = self.kbse.search_objects(search_params)['objects']
                if objects:
                    taxon_id = objects[0].get('object_name')
            except Exception:
                logging.info('Failed looking up taxon in search engine')
                logging.warning(traceback.format_exc())
                logging.warning(sys.exc_info()[2])
                taxon_id = None

            self.taxon_cache[scientific_name] = taxon_id

        return taxon_id

    def _fetch_taxon_level(self, taxon_char):

        taxon_level_mapping = {'l': 'Life', 'd': 'Domain', 'k': 'Kingdom', 'p': 'Phylum',
                               'c': 'Class', 'o': 'Order', 'f': 'Family', 'g': 'Genus',
                               's': 'Species'}

        return taxon_level_mapping.get(taxon_char.lower(), taxon_char)

    def _fetch_taxonomy(self, datarow, observation_id=None, am_data=None):

        taxonomy = dict()
        lineage = None
        if am_data:
            attributes = am_data['attributes']
            instances = am_data['instances']

            instance = instances.get(observation_id)
            if not instance:
                raise ValueError('Cannot find {} instance from attributes'.format(observation_id))

            for index, attribute in enumerate(attributes):
                if attribute['attribute'] == 'taxonomy':
                    lineage = instance[index]
                    break

            if isinstance(lineage, str):
                delimiter = csv.Sniffer().sniff(lineage).delimiter
                lineage = [x.strip() for x in lineage.split(delimiter)]

            if lineage:
                taxonomy['lineage'] = lineage

            for key in ['score', 'taxonomy_source', 'species_name']:
                for index, attribute in enumerate(attributes):
                    if attribute['attribute'] == key:
                        taxonomy[key] = instance[index]
        else:
            lineage = self._retrieve_value([], datarow, 'taxonomy')

            if isinstance(lineage, str):
                delimiter = csv.Sniffer().sniff(lineage).delimiter
                lineage = [x.strip() for x in lineage.split(delimiter)]

            if lineage:
                taxonomy['lineage'] = lineage

            for key in ['score', 'taxonomy_source', 'species_name']:
                val = self._retrieve_value([], datarow, key)
                if val:
                    taxonomy[key] = val

        # if lineage:
        #     for taxon_str in lineage[::-1]:
        #         taxon_items = taxon_str.split('__')
        #         if len(taxon_items) == 1:
        #             scientific_name = taxon_items[-1]
        #             taxon_level_char = 'None'
        #         else:
        #             scientific_name = taxon_items[-1]
        #             taxon_level_char = taxon_items[0]

        #         if scientific_name:
        #             taxon_id = self._search_taxon(scientific_name)
        #             if taxon_id:
        #                 taxon_ref = f"{self.taxon_wsname}/{taxon_id}"
        #                 taxon_level = self._fetch_taxon_level(taxon_level_char)

        #                 taxonomy.update({'taxon_ref': taxon_ref,
        #                                  'taxon_id': taxon_id,
        #                                  'scientific_name': scientific_name,
        #                                  'taxon_level': taxon_level})
        #                 break

        return taxonomy

    def _retrieve_tsv_amplicon_set_data(self, tsv_file, am_data):
        amplicons = dict()

        try:
            logging.info('start parsing TSV file')
            reader = pd.read_csv(tsv_file, sep=None, iterator=True)
            inferred_sep = reader._engine.data.dialect.delimiter
            df = pd.read_csv(tsv_file, sep=inferred_sep, index_col=0)
        except Exception:
            raise ValueError('Cannot parse file. Please provide valide TSV file')

        if 'consensus_sequence' not in df.columns.tolist():
            raise ValueError('TSV file does not include consensus_sequence')

        logging.info('start processing each row in TSV')
        for observation_id in df.index:
            taxonomy = self._fetch_taxonomy(df.loc[observation_id],
                                            observation_id,
                                            am_data)

            amplicon = {'consensus_sequence': df.loc[observation_id, 'consensus_sequence'],
                        'taxonomy': taxonomy}

            amplicons.update({observation_id: amplicon})

        logging.info('finished parsing TSV file')

        return amplicons

    def _retrieve_tsv_fasta_amplicon_set_data(self, tsv_file, fasta_file, am_data):
        amplicons = dict()
        try:
            logging.info('start parsing FASTA file')
            fastq_dict = SeqIO.index(fasta_file, "fasta")
        except Exception:
            raise ValueError('Cannot parse file. Please provide valide FASTA file')

        try:
            logging.info('start parsing TSV file')
            reader = pd.read_csv(tsv_file, sep=None, iterator=True)
            inferred_sep = reader._engine.data.dialect.delimiter
            df = pd.read_csv(tsv_file, sep=inferred_sep, index_col=0)
        except Exception:
            raise ValueError('Cannot parse file. Please provide valide TSV file')

        logging.info('start processing files')
        for observation_id in df.index:
            if observation_id not in fastq_dict:
                raise ValueError('FASTA file does not have [{}] OTU id'.format(observation_id))

            taxonomy = self._fetch_taxonomy(df.loc[observation_id],
                                            observation_id,
                                            am_data)

            amplicon = {'consensus_sequence': str(fastq_dict.get(observation_id).seq),
                        'taxonomy': taxonomy}
            amplicons.update({observation_id: amplicon})

        logging.info('finished processing files')
        return amplicons

    def _retrieve_biom_fasta_amplicon_set_data(self, biom_file, fasta_file, am_data):
        amplicons = dict()
        try:
            logging.info('start parsing FASTA file')
            fastq_dict = SeqIO.index(fasta_file, "fasta")
        except Exception:
            raise ValueError('Cannot parse file. Please provide valide FASTA file')

        logging.info('start parsing BIOM file')
        table = biom.load_table(biom_file)

        observation_ids = table._observation_ids.tolist()
        observation_metadata = table._observation_metadata

        logging.info('start processing files')
        for index, observation_id in enumerate(observation_ids):
            if observation_id not in fastq_dict:
                raise ValueError('FASTA file does not have [{}] OTU id'.format(observation_id))

            taxonomy = self._fetch_taxonomy(observation_metadata[index],
                                            observation_id,
                                            am_data)

            amplicon = {'consensus_sequence': str(fastq_dict.get(observation_id).seq),
                        'taxonomy': taxonomy}

            amplicons.update({observation_id: amplicon})

        logging.info('finished processing files')
        return amplicons

    def _retrieve_biom_tsv_amplicon_set_data(self, biom_file, tsv_file, am_data):
        amplicons = dict()
        try:
            logging.info('start parsing TSV file')
            reader = pd.read_csv(tsv_file, sep=None, iterator=True)
            inferred_sep = reader._engine.data.dialect.delimiter
            df = pd.read_csv(tsv_file, sep=inferred_sep, index_col=0)
        except Exception:
            raise ValueError('Cannot parse file. Please provide valide tsv file')

        if 'consensus_sequence' not in df.columns.tolist():
            raise ValueError('TSV file does not include consensus_sequence')

        logging.info('start parsing BIOM file')
        table = biom.load_table(biom_file)

        observation_ids = table._observation_ids.tolist()
        observation_metadata = table._observation_metadata

        logging.info('start processing files')
        for index, observation_id in enumerate(observation_ids):
            if observation_id not in df.index:
                raise ValueError('TSV file does not have [{}] OTU id'.format(observation_id))

            taxonomy = self._fetch_taxonomy(df.loc[observation_id],
                                            observation_id,
                                            am_data)

            amplicon = {'consensus_sequence': df.loc[observation_id, 'consensus_sequence'],
                        'taxonomy': taxonomy}

            amplicons.update({observation_id: amplicon})

        logging.info('finished processing files')
        return amplicons

    def _file_to_amplicon_set_data(self, biom_file, tsv_file, fasta_file, mode, refs, description):

        logging.info('start parsing amplicon_set_data')

        amplicon_set_data = dict()

        row_attributemapping_ref = refs.get('row_attributemapping_ref')
        am_data = None
        if row_attributemapping_ref:
            am_data = self.dfu.get_objects({'object_refs': [row_attributemapping_ref]}
                                           )['data'][0]['data']

        if mode == 'biom_tsv':
            amplicons = self._retrieve_biom_tsv_amplicon_set_data(biom_file, tsv_file, am_data)
        elif mode == 'biom_fasta':
            amplicons = self._retrieve_biom_fasta_amplicon_set_data(biom_file, fasta_file, am_data)
        elif mode == 'tsv_fasta':
            amplicons = self._retrieve_tsv_fasta_amplicon_set_data(tsv_file, fasta_file, am_data)
        elif mode == 'tsv':
            amplicons = self._retrieve_tsv_amplicon_set_data(tsv_file, am_data)
        else:
            raise ValueError('error parsing _file_to_amplicon_set_data, mode: {}'.format(mode))

        amplicon_set_data.update({'amplicons': amplicons})

        if 'reads_set_ref' in refs:
            amplicon_set_data['reads_set_ref'] = refs.get('reads_set_ref')

        if description:
            amplicon_set_data['description'] = description

        return amplicon_set_data

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

        if False and len(data_df.columns) <= 200:
            if top_heatmap_dir:
                viewer_name = 'MatrixLinearPlotViewer'
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += '''>Top {} Percent Linear Plot</button>\n'''.format(
                                                                            round(top_percent, 2))

                linear_plot_page = self._generate_linear_plot(data_df, output_directory,
                                                              row_name='OTU',
                                                              top_percent=round(top_percent, 2))

                tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
                msg = 'Top {} percent of matrix sorted by sum of abundance values.'.format(
                                                                            round(top_percent, 2))
                tab_content += '''<p style="color:red;" >{}</p>'''.format(msg)
                tab_content += '\n<iframe height="1300px" width="100%" '
                tab_content += 'src="{}" '.format(linear_plot_page)
                tab_content += 'style="border:none;"></iframe>'
                tab_content += '\n</div>\n'
            else:
                viewer_name = 'MatrixLinearPlotViewer'
                tab_def_content += '''\n<button class="tablinks" '''
                tab_def_content += '''onclick="openTab(event, '{}')"'''.format(viewer_name)
                tab_def_content += '''>Matrix Linear Plot</button>\n'''

                linear_plot_page = self._generate_linear_plot(data_df, output_directory, row_name='OTU')

                tab_content += '''\n<div id="{}" class="tabcontent">'''.format(viewer_name)
                tab_content += '\n<iframe height="1300px" width="100%" '
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

    def _df_to_tsv(self, amplicon_set_df, result_dir, amplicon_set_ref):
        logging.info('writting amplicon set data frame to tsv file')
        amplicon_set_obj = self.dfu.get_objects({'object_refs': [amplicon_set_ref]})['data'][0]
        amplicon_set_info = amplicon_set_obj['info']
        amplicon_set_name = amplicon_set_info[1]

        file_path = os.path.join(result_dir, amplicon_set_name + ".tsv")

        amplicon_set_df.to_csv(file_path, sep='\t', index=True, header=True)

        return file_path

    def _amplicon_set_to_df(self, amplicon_set_ref):
        logging.info('converting amplicon set to data frame')
        am_set_data = self.dfu.get_objects({'object_refs': [amplicon_set_ref]})['data'][0]['data']

        amplicon_matrix_ref = am_set_data.get('amplicon_matrix_ref')
        matrix_data = self.dfu.get_objects({'object_refs': [amplicon_matrix_ref]})['data'][0]['data']
        matrix_value_data = matrix_data.get('data')

        index = matrix_value_data.get('row_ids')
        columns = matrix_value_data.get('col_ids')
        values = matrix_value_data.get('values')

        df = pd.DataFrame(values, index=index, columns=columns)

        amplicons = am_set_data.get('amplicons')

        meta_index = list()

        meta_columns = ['taxonomy', 'taxon_id', 'taxon_ref', 'taxon_level', 'score',
                        'taxonomy_source', 'species_name', 'consensus_sequence']
        meta_values = list()
        for otu_id, amplicon in amplicons.items():
            meta_index.append(otu_id)

            taxonomy_data = amplicon.get('taxonomy')

            taxonomy = taxonomy_data.get('lineage')
            taxon_id = taxonomy_data.get('taxon_id')
            taxon_ref = taxonomy_data.get('taxon_ref')
            taxon_level = taxonomy_data.get('taxon_level')
            score = taxonomy_data.get('score')
            taxonomy_source = taxonomy_data.get('taxonomy_source')
            species_name = taxonomy_data.get('species_name')

            consensus_sequence = amplicon.get('consensus_sequence')

            meta_values.append([taxonomy, taxon_id, taxon_ref, taxon_level, score, taxonomy_source,
                                species_name, consensus_sequence])

        meta_df = pd.DataFrame(meta_values, index=meta_index, columns=meta_columns)

        merged_df = df.merge(meta_df, left_index=True, right_index=True, how='left',
                             validate='one_to_one')

        return merged_df

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
        self.ws = workspaceService(config["workspace-url"], token=self.token)

    def fetch_sequence(self, matrix_ref):
        logging.info('start to fetch consensus sequence')

        input_matrix_obj = self.dfu.get_objects({'object_refs': [matrix_ref]})['data'][0]

        logging.info('using workspace to fetch object')
        input_matrix_obj = self.ws.get_objects2({"objects": [{'ref': matrix_ref}]})['data'][0]
        print('return from ws')
        print(input_matrix_obj)

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

        for key in ['extraction_kit', 'amplicon_type', 'target_gene_region',
                    'forward_primer_sequence', 'reverse_primer_sequence', 'sequencing_platform',
                    'sequencing_run', 'sequencing_kit', 'sequencing_quality_filter_cutoff',
                    'clustering_cutoff', 'clustering_method', 'sample_set_ref']:
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

    def export_amplicon_set_tsv(self, params):
        """
        export AmpliconSet as TSV
        """
        logging.info('start exporting amplicon set object')
        amplicon_set_ref = params.get('input_ref')

        amplicon_set_df = self._amplicon_set_to_df(amplicon_set_ref)

        result_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_dir)

        self._df_to_tsv(amplicon_set_df, result_dir, amplicon_set_ref)

        package_details = self.dfu.package_for_download({
            'file_path': result_dir,
            'ws_refs': [amplicon_set_ref]
        })

        return {'shock_id': package_details['shock_id']}
