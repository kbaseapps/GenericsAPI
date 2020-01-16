
import json
import requests
import uuid

from GenericsAPI.Utils.DataUtil import DataUtil
from installed_clients.DataFileUtilClient import DataFileUtil


class SampleServiceUtil:

    def _fetch_attri_from_meta(self, meta_data):
        attributes = list()

        for key, value in meta_data.items():
            attri = {'attribute': key, 'source': 'SampleService'}
            if 'units' in value:
                attri['unit'] = value.get('units')

            attributes.append(attri)

        return attributes

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.token = config['KB_AUTH_TOKEN']
        self.srv_wiz_url = config['srv-wiz-url']
        self.dfu = DataFileUtil(self.callback_url)
        # self.data_util = DataUtil(config)
        # self.matrix_types = [x.split(".")[1].split('-')[0]
        #                      for x in self.data_util.list_generic_types()]

    def get_sample_service_url(self):

        payload = {
            "method": "ServiceWizard.get_service_status",
            "id": '',
            "params": [{"module_name": "SampleService", "version": "dev"}],  # TODO: change to beta/release
            "version": "1.1"
        }

        sw_resp = requests.post(url=self.srv_wiz_url, data=json.dumps(payload))
        wiz_resp = sw_resp.json()
        if wiz_resp.get('error'):
            raise RuntimeError("ServiceWizard Error - " + str(wiz_resp['error']))

        return wiz_resp['result'][0]['url']

    def get_sample(self, sample_id):

        sample_url = self.get_sample_service_url()
        headers = {"Authorization": self.token}
        params = {
            "id": sample_id,
        }
        payload = {
            "method": "SampleService.get_sample",
            "id": str(uuid.uuid4()),
            "params": [params],
            "version": "1.1"
        }
        resp = requests.post(url=sample_url, headers=headers, data=json.dumps(payload))
        resp_json = resp.json()
        if resp_json.get('error'):
            raise RuntimeError(f"Error from SampleService - {resp_json['error']}")
        sample = resp_json['result'][0]

        return sample

    def save_sample(self, sample):
        """
        copied from sample_uploader
        (https://github.com/kbaseapps/sample_uploader/blob/master/lib/sample_uploader/utils/sample_utils.py#L79)
        """
        sample_url = self.get_sample_service_url()
        headers = {"Authorization": self.token}
        params = {
            "sample": sample,
            "prior_version": None,
        }
        payload = {
            "method": "SampleService.create_sample",
            "id": str(uuid.uuid4()),
            "params": [params],
            "version": "1.1"
        }
        resp = requests.post(url=sample_url, headers=headers, data=json.dumps(payload, default=str))
        if not resp.ok:
            raise RuntimeError(f'Error from SampleService - {resp.text}')
        resp_json = resp.json()
        if resp_json.get('error'):
            raise RuntimeError(f"Error from SampleService - {resp_json['error']}")
        sample_id = resp_json['result'][0]['id']
        return sample_id

    def sample_set_to_attribute_mapping(self, sample_set_ref):

        sample_set = self.dfu.get_objects(
                    {"object_refs": [sample_set_ref]})['data'][0]['data']

        samples = sample_set['samples']

        am_data = {'ontology_mapping_method': "SampleService", 'instances': {}}

        attributes = list()

        sample_datas = []
        for sample in samples:
            sample_id = sample['id']

            sample_data = self.get_sample(sample_id)

            node_tree = sample_data.get('node_tree', [{}])

            for node in node_tree:
                meta_controlled = node.get('meta_controlled')
                meta_user = node.get('meta_user')

                if meta_controlled:
                    attributes_controlled = self._fetch_attri_from_meta(meta_controlled)
                    attributes += attributes_controlled

                if meta_user:
                    attributes_user = self._fetch_attri_from_meta(meta_user)
                    attributes += attributes_user

            sample_datas.append(sample_data)

        attributes = [i for n, i in enumerate(attributes) if i not in attributes[n + 1:]]
        attributes = [{'attribute': 'id', 'source': 'SampleService'},
                      {'attribute': 'type', 'source': 'SampleService'},
                      {'attribute': 'parent', 'source': 'SampleService'}] + attributes

        am_data['attributes'] = attributes

        instances = am_data['instances']

        for sample_data in sample_datas:
            instance = list()
            node_tree = sample_data.get('node_tree', [{}])

            for node in node_tree:

                meta_controlled = node.get('meta_controlled')
                meta_user = node.get('meta_user')

                for attribute in attributes:
                    attri_name = attribute['attribute']

                    if attri_name in ['id', 'type', 'parent']:
                        instance.append(node.get(attri_name))
                    else:
                        if meta_user:
                            instance.append(meta_user.get(attri_name, {}).get('value'))
                        elif meta_controlled:
                            instance.append(meta_controlled.get(attri_name, {}).get('value'))

            instance = [str(i) for i in instance]

            instances[sample_data['name']] = instance

        return am_data
