import json
import glob

json_file_PIP = "../test-results/Ti/element-Ti-model-PIP_hcp_bcc_co-test-phonon_hcp-properties.json"
json_file_NRLTB = "../test-results/Ti/element-Ti-model-NRLTB-test-phonon_hcp-properties.json"

PIP_json = json.load(open(json_file_PIP))
NRLTB_json = json.load(open(json_file_NRLTB))