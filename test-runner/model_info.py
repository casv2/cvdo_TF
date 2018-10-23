import os
import matplotlib.pyplot as plt
import argparse
import glob
import json
from pprint import pprint
import matplotlib.pyplot as plt
import subprocess

cwd = os.getcwd()

parser = argparse.ArgumentParser("model info")
parser.add_argument('element', action='store', type=str, help='element')
parser.add_argument('potentials', nargs='+', help='potentials')
args = parser.parse_args()

print args.potentials

fig = plt.figure(figsize=(12,12))
plt.axis([0, 40, -10, 3])

class deserialise(object):
	def __init__(self,j):
		self.__dict__ = json.loads(j)

i = 4

for pot in args.potentials:
	model_dir = os.path.join('..', 'models', args.element, pot)
	os.chdir(model_dir)
	if pot.startswith("GAP"):
		xml_filename = glob.glob("*.xml")
		print xml_filename
		xml_file = open(xml_filename[0], "r")
		for line in xml_file:
			if line.startswith("  <command_line>"):
				pot_string = pot #+ "\n" 
				plt.text(-5, i, pot_string, fontsize=14, ha='left', wrap=True)
				i += -1.9
				plt.text(-5, i, line, fontsize=12, ha='left', wrap=True)
				i += -0.8

		gap_configs = []
		xyz_filename = glob.glob("*.xml.xyz")
		print xyz_filename, len(xyz_filename)
		xyz_data_filename = glob.glob("*.xyz_data")
		print xyz_data_filename, len(xyz_data_filename)
		if len(xyz_filename) != 0:
			xyz_file = open(xyz_filename[0], "r")
			for line in xyz_file:
				if line.startswith("config_type"):
					gap_configs.append(line.split(" ")[0].split("=")[1])
		elif len(xyz_data_filename) != 0:
			xyz_data_file = open(xyz_data_filename[0], "r")
			print "hello"
			for line in xyz_data_file:
				if line.startswith("    <![CDATA[OMIT_config_type"):
					gap_configs.append(line.split("=")[1].split(" ")[0])
		else:
			xml_file = open(xml_filename[0], "r")
			for line in xml_file:
				if "configtype" in line:
					gap_configs.append(line.split("configtype")[1].split(" ")[0][1:])
				#if "config_type" in line:
				#	print line 
					#gap_configs.append(line.split("configt_ype")[1].split(" ")[0][1:])
				if line.startswith("    <![CDATA[config_type"):
					gap_configs.append(line.split("=")[1].split(" ")[0])

		configs = [config for config in set(gap_configs)]
		num_configs = [gap_configs.count(config) for config in set(gap_configs)]
		config_dict = dict(zip(configs, num_configs))

		for key in config_dict.keys():
			if "." in key:
				del config_dict[key]

		plt.text(-5, i, str(config_dict), fontsize=12, ha='left', wrap=True)
		i += -0.5
	if pot.startswith("PIP"):
		#i = i + 0.8
		i += -0.35
		try:
			info_json = json.load(open("./info.json"))
		except:
			subprocess.call(["python", "model.py"])
			info_json = json.load(open("./info.json"))
		info = deserialise(info_json)
		string = str("config_weights " + str(info.configweights) + "\ndata_weights " + str(info.dataweights) + "\nnum_configs " + str(info.numconfigs) + "\ndb " + str(info.dbpath).split("/")[-1])
		string = string.replace("u\'", "").replace("\'", "")
		plt.text(-5, i, pot, fontsize=14, ha='left', wrap=True)
		i += -2
		plt.text(-5, i, string, fontsize=12, ha='left', wrap=True)
		i += -1
	os.chdir(cwd)

plt.axis('off')
plt.savefig("../test-results/{0}/images/_pot-info.png".format(args.element))
#plt.show()

