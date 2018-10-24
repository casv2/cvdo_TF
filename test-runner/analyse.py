
import os
import glob
import json 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import subprocess
import __builtin__
from pandas.tools.plotting import table
import numpy as np
import itertools

import time
start_time = time.time()

parser = argparse.ArgumentParser(description='run fast test')
parser.add_argument('element', action='store', type=str, help='element')
parser.add_argument('potentials', nargs='+')
#parser.add_argument('fast/slow')?
args = parser.parse_args()

print args.potentials

#cwd = os.getcwd
__builtin__.do_io = True

cwd = os.path.join('..', "test-results", args.element)

images_dir = '../test-results/{0}/images'.format(args.element)
if not os.path.exists(images_dir): # and __builtin__.do_io:
    os.mkdir(images_dir)
time.sleep(1)

image_folder = os.path.join(cwd, "images")
#os.path.join(image_folder, )

def pct_change(df):

	for column in df.columns[2:]:
		print column
		print df.columns
		df[column] = (df[column] - df[df.columns[1]])/df[df.columns[1]] * 100
		df[column] = [str(int(df[column].tolist()[i]))+"%" for i in xrange(0,len(df[column].tolist()))]

	df = df.round(2)
	return df

def run_tests(potential, element):
	if element == "Si":
		tests = ["111_layer_test", "bulk_diamond", "surface-energy-100-unrelaxed", 
	"surface-energy-111-unrelaxed", "fourfold-defect-small-GAP_traj_eval", "vac-path-no-relax-configs", 
	"dimer", "surface-decohesion-100-unrelaxed", "surface-decohesion-111-unrelaxed", "surface-energy-110-unrelaxed","surface-decohesion-110-unrelaxed"]
	#"grain-boundary-GAP-relaxed-configs", "bulk_fcc", "bulk_hcp", "bulk_sh", "bulk_bcc", "bulk_hex_diamond", "bulk_hcp_short_and_fat", "bulk_st12"] # ,"111_layer_test"] #"diinterstitial_GAP_relaxed-64-energy" ] 
	elif element == "W":
		tests = ["vacancy-energy", "bulk_bcc", "surface-decohesion-110-unrelaxed", "surface-energy-100-unrelaxed", "surface-energy-111-unrelaxed"
	"surface-energy-111-unrelaxed", "surface-decohesion-100-unrelaxed", "surface-decohesion-111-unrelaxed"] #"surface-decohesion-110-unrelaxed", "surface-energy-100-unrelaxed", "surface-energy-110-unrelaxed", 
	elif element == "Ti":
		tests = ['burgers_path', "bulk_bcc", "bulk_hcp","surface-decohesion-100-unrelaxed", "surface-decohesion-111-unrelaxed", "surface-decohesion-110-unrelaxed" ]#, 'phonon_hcp']#, 'phonon_bcc']# ['burgers_path']"bulk_bcc", "bulk_hcp", "bulk_omega", "dimer", "force-constant-bcc", "force-constant-hcp", 
		#"surface-decohesion-100-unrelaxed", "surface-decohesion-111-unrelaxed", "surface-decohesion-110-unrelaxed", 
		#"force-constant-omega", "surface-energy-100-unrelaxed", "surface-energy-111-unrelaxed", "hcp_to_bcc"] #"surface-decohesion-110-unrelaxed", "surface-energy-110-unrelaxed", "surface-energy-100-unrelaxed", "surface-energy-110-unrelaxed", "surface-energy-111-unrelaxed", "surface-decohesion-100-unrelaxed", "surface-decohesion-110-unrelaxed", "surface-decohesion-111-unrelaxed", 
	ps = []
	for test in tests:
		r_str = ("python run-model-test-cas.py {0} {1} {2}".format(element, potential, test)).split(" ")
		#print r_str
		p = subprocess.Popen(r_str)
		ps.append(p)
	for p in ps:
		p.wait()

def table_plot(potentials, element):
	data_dict = { '-' : []}
	if element == "Si":
		#potentials.append("CASTEP_ASE")
		tests = ["bulk_diamond", "surface-energy-100-unrelaxed", "surface-energy-111-unrelaxed"]  # "bulk_fcc", "bulk_hcp", "bulk_bcc", "bulk_st12" , "bulk_sh", "bulk_hex_diamond", "bulk_hcp_short_and_fat"] #"surface-energy-110-unrelaxed",
	elif element == "W":
		tests = ["bulk_bcc", "surface-energy-100-unrelaxed"]#, "surface-energy-111-unrelaxed", "surface-energy-100-unrelaxed"] #, "surface-energy-100-unrelaxed", "surface-energy-110-unrelaxed", "surface-energy-111-unrelaxed"
	elif element == "Ti":
		tests = [ "bulk_bcc", "bulk_hcp"] #"surface-decohesion-110-unrelaxed","surface-energy-110-unrelaxed",
	for test in sorted(tests):
		if test.split('_')[0] == "bulk":
			plt.figure(figsize=(7,7))
			print test.split('_')[0]
		for pot in potentials:
			json_file = cwd + "/element-{0}-model-{1}-test-{2}-properties.json".format(element, pot, test)
			print json_file
			model_name = os.path.basename(json_file).split('-')[3]
			print model_name
			if model_name not in potentials:
				continue
			if model_name not in data_dict.keys():
				data_dict[model_name] = []
			try:
				data = json.load(open(json_file))
				for key in data.keys():
					entry = data[key]
					if type(entry) == float:
						data_dict[model_name].append(entry)
						if key not in data_dict["-"]:
							data_dict['-'].append(key)
					elif type(entry) == list:
						V = [entry[i][0] for i in xrange(0,len(entry))]
						E = [entry[i][1] for i in xrange(0,len(entry))]
						print json_file, E, V
						E_shift = [ E[i]-min(E) for i in xrange(0,len(E))]
						plt.plot(V,E_shift, label=model_name)
			except:
				print "continue"
		if test.split('_')[0] == "bulk":
			path = os.path.join(image_folder, test)
			plt.legend()
			plt.xlabel("Volume (A)")
			plt.ylabel("Energy (eV)")
			plt.title(test)
			plt.savefig(path)
	print data_dict

	for key in data_dict.keys():
		print key, len(data_dict[key])

	df = pd.DataFrame.from_dict(data_dict)	

	#columnsTitles=[ "-", "NRLTB", "GAP_soap", "Ti_Env4B_reg_tmp"]
	#df=df.reindex(columns=columnsTitles)

	print df

	df = pct_change(df)
	
	"""
	if element == "Ti":
		df.set_index('-', inplace=True, drop=True)
		data_dict = pd.DataFrame.to_dict(df)
		print data_dict
		force_constant_list = ["force-constant-bcc", "force-constant-hcp"] #model-CASTEP_ASE-test-force-constant-hcp
		for force_constant in force_constant_list:
			#data_dict["-"].append(force_constant)
			for pot in potentials:
				print pot
				av_rel_error = force_constant_plot(force_constant, pot, element)
				print av_rel_error
				data_dict[pot][force_constant] =  av_rel_error
				#av_rel_error = force_constant_plot(pot, element)[0]
				#print pot, av_rel_error
		df = pd.DataFrame.from_dict(data_dict)	
	"""
	print df

	filename = 'out.tex'
	pdffile = 'out.pdf'
	outname = '_table.png'

	template = r'''
	\documentclass[preview]{{standalone}}
	\usepackage{{graphicx}}
	\usepackage{{booktabs}}
	\begin{{document}}
	\begin{{table}}
	\begin{{center}}
	\scalebox{{0.6}}{{
	{}}}
	\end{{center}}
	\end{{table}}
	\end{{document}}
	''' 
	with open(filename, 'wb') as f:
		f.write(template.format(df.to_latex()))

	subprocess.call(['pdflatex', filename])
	subprocess.call(['convert', '-strip', '-density', '300', pdffile, '-quality', '90', outname])

	os.rename("_table.png", image_folder + "/_table.png")
	os.system("rm out.*")


def decohesion_plot(potentials, element):
	plt.figure(figsize=(8,8))
	potentials.append("CASTEP_ASE")
	tests = ["surface-decohesion-100-unrelaxed", "surface-decohesion-110-unrelaxed", "surface-decohesion-111-unrelaxed"]
	for test in tests:
		f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,5))
		plt.title(test)
		for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
			model_name = os.path.basename(json_file).split('-')[3]
			if model_name not in potentials:
				continue
			data = json.load(open(json_file))
			ax1.plot(data[data.keys()[1]], data[data.keys()[0]], label=model_name)
			ax1.set_xlabel("Unrelaxed opening (A)")
			ax1.set_ylabel("Energy")
			ax2.plot(data[data.keys()[1]], data[data.keys()[2]], label=model_name)
			ax2.set_xlabel("Unrelaxed opening (A)")
			ax2.set_ylabel("Stress")
		plt.legend()
		path = os.path.join(image_folder, test)
		plt.savefig(path)	

def fourfold_plot(potentials, element):
	test = "fourfold-defect-small-GAP_traj_eval"
	plt.figure(figsize=(8,8))
	for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
		model_name = os.path.basename(json_file).split('-')[3]
		if model_name not in potentials:
			continue
		data = json.load(open(json_file))
		E = [data[data.keys()[0]][i][1] for i in xrange(0,len(data[data.keys()[0]]))]
		plt.plot(E, label=model_name)
		plt.xlabel("Evaluation")
		plt.ylabel("Energy")
		plt.title(test)
	plt.legend()
	path = os.path.join(image_folder, test)
	plt.savefig(path)

def grain_plot(potentials, element):
	test = "grain-boundary-GAP-relaxed-configs"
	plt.figure(figsize=(8,8))
	for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
		model_name = os.path.basename(json_file).split('-')[3]
		if model_name not in potentials:
			continue
		data = json.load(open(json_file))
		E = [data[data.keys()[0]][i][1] for i in xrange(0,len(data[data.keys()[0]]))]
		plt.plot(E, label=model_name)
		plt.xlabel("Evaluation")
		plt.ylabel("Energy")
		plt.title(test)
	plt.legend()
	path = os.path.join(image_folder, test)
	plt.savefig(path)

def diinterstitial_plot(potentials, element):
	test = "diinterstitial-GAP-relaxed-configs"
	plt.figure(figsize=(8,8))
	for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
		model_name = os.path.basename(json_file).split('-')[3]
		if model_name not in potentials:
			continue
		data = json.load(open(json_file))
		E = [data[data.keys()[0]][i][1] for i in xrange(0,len(data[data.keys()[0]]))]
		#print E
		plt.plot(E, label=model_name)
		plt.xlabel("diinterstitial split.*.xyz")
		plt.ylabel("Energy")
		plt.title(test)
	plt.legend()
	path = os.path.join(image_folder, test)
	plt.savefig(path)

def vacancy_path_plot(potentials, element):
	test = "vac-path-no-relax-configs"
	plt.figure(figsize=(8,8))
	for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
		model_name = os.path.basename(json_file).split('-')[3]
		if model_name not in potentials:
			continue
		data = json.load(open(json_file))
		E = data[data.keys()[0]]
		p = plt.plot(E, label=model_name)
		plt.xlabel("Evaluation")
		plt.ylabel("Energy")
		plt.title(test)
	plt.legend()
	path = os.path.join(image_folder, test)
	plt.savefig(path)

def dimer_plot(potentials, element):
	test = "dimer"
	for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
		model_name = os.path.basename(json_file).split('-')[3]
		if model_name not in potentials:
			continue
		data = json.load(open(json_file))
		rr = data[data.keys()[1]]
		E = data[data.keys()[0]]
		p = plt.plot(rr, E, label=model_name)
		plt.xlabel("Distance (A)")
		plt.ylabel("Energy (eV)")
		plt.title(test)
	plt.legend()
	path = os.path.join(image_folder, test)
	plt.savefig(path)

def force_constant_plot(force_constant, pot, element):
	test = force_constant
	print test
	castepdata = json.load(open(cwd + "/element-{0}-model-CASTEP_ASE-test-{1}-properties.json".format(element, test)))
	data = json.load(open(cwd + "/element-{0}-model-{1}-test-{2}-properties.json".format(element, pot, test)))

	rel_error_list = []

	print np.mean(np.array(castepdata[castepdata.keys()[0]][0]) - np.array(data[data.keys()[0]][0]))

	for j in xrange(0, len(castepdata)):
		rel_error = np.mean([np.sqrt(np.mean((np.array(castepdata[castepdata.keys()[j]][i]) - np.array(data[data.keys()[j]][i])) ** 2)) for i in xrange(0,len(castepdata[castepdata.keys()[0]]))])
		rel_error_list.append(rel_error)

	print np.mean(rel_error_list)
	return np.mean(rel_error_list)

	#c_config1 = castepdata[castepdata.keys()[0]]
	#c_config2 = castepdata[castepdata.keys()[1]]

	#config1 = data[data.keys()[0]]
	#config2 = data[data.keys()[1]]

	#force_error_1 =  np.mean([np.sqrt(np.mean((np.array(c_config1[i]) - np.array(config1[i])) ** 2)) for i in xrange(0,len(c_config1))])
	#force_error_2 = np.mean([np.sqrt(np.mean((np.array(c_config2[i]) - np.array(config2[i])) ** 2)) for i in xrange(0,len(c_config1))])

	#av_force_error = np.mean([force_error_1, force_error_2])
	#print av_force_error


#force_constant_plot("force-constant-omega", args.potentials[2], args.element)
#table_plot(args.potentials, args.element)

def layer_plot(potentials, element):
	test = "111_layer_test"
	plt.figure(figsize=(8,8))
	potentials.append("CASTEP_ASE")
	for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
		model_name = os.path.basename(json_file).split('-')[3]
		if model_name not in potentials:
			continue
		data = json.load(open(json_file))
		E = [data[data.keys()[0]][i][1]/6 for i in xrange(0,len(data[data.keys()[0]]))] #24 atoms!
		E_per_a  = np.array(E) - min(E)
		spacing = [0.02*i for i in xrange(0,len(E))]
		plt.plot(spacing, E_per_a, label=model_name)
		plt.xlabel("Interplanar displacement (Angstrom)")
		plt.ylabel("Energy per atom (eV)")
		plt.title(test)
	plt.legend()
	test = "z111_layer_test"
	path = os.path.join(image_folder, test)
	plt.savefig(path)


def burgers_path(potentials, element):
	test = "burgers_path"
	plt.figure(figsize=(8,8))
	for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
		model_name = os.path.basename(json_file).split('-')[3]
		if model_name not in potentials:
			continue
		data = json.load(open(json_file))
		E = [data[data.keys()[0]][i][1] for i in xrange(0,len(data[data.keys()[0]]))] #24 atoms!
		E_per_a  = np.array(E) - min(E)
		#spacing = [0.02*i for i in xrange(0,len(E))]
		plt.plot(E_per_a, label=model_name)
	plt.xlabel("Transform bcc to hcp")
	plt.ylabel("Energy per atom (eV)")
	plt.title(test)
	plt.legend()
	path = os.path.join(image_folder, test)
	plt.savefig(path)

def pot_info(potentials, element):
	r_str = ["python", "model_info.py"]
	r_str.append(element)
	for pot in potentials:
		r_str.append(pot)
	print r_str
	subprocess.call(r_str)

def phonon(potentials, element, phase):
	test = "phonon_{0}".format(phase)
	plt.figure(figsize=(10,10))
	colours = ["red", "blue", "green", "yellow"]
	linestyles = ["-",":", "-.", "--"]
	i = 0
	for json_file in sorted(glob.glob(cwd + "/element-{0}-model-*-test-{1}-properties.json".format(element, test))):
		model_name = os.path.basename(json_file).split('-')[3]
		print model_name
		if model_name not in potentials:
			continue
		data = json.load(open(json_file))
		omega_kn = np.array(data["omega_kn"])
		point_names = data["point_names"]
		path = data["path"]
		Q = data["Q"]
		q = data["q"]
		n = 0
		for n in range(len(omega_kn[0])):
			omega_n = omega_kn[:, n]
			plt.plot(q, omega_n, color=colours[i], label=model_name if n == 0 else "", linestyle=linestyles[i])
			n += 1
		i += 1
	plt.xticks(Q, point_names)
	plt.xlim(q[0], q[-1])
	plt.grid('on')
	plt.title(test)
	plt.legend()
	path = os.path.join(image_folder, test)
	plt.savefig(path)


for pot in args.potentials:
	run_tests(pot, args.element)


if args.element == "Si":
	pot_info(args.potentials, args.element)
	layer_plot(args.potentials, args.element)
	#dimer_plot(args.potentials, args.element)
	#vacancy_path_plot(args.potentials, args.element)
	#fourfold_plot(args.potentials, args.element)
	table_plot(args.potentials, args.element)
	decohesion_plot(args.potentials, args.element)
elif args.element == "W":
	table_plot(args.potentials, args.element)
	decohesion_plot(args.potentials, args.element)
elif args.element == "Ti":
	pot_info(args.potentials, args.element)
	#phonon(args.potentials, args.element, "bcc")
	#phonon(args.potentials, args.element, "bcc")
	burgers_path(args.potentials, args.element)
	#dimer_plot(args.potentials, args.element)
	table_plot(args.potentials, args.element)
	decohesion_plot(args.potentials, args.element)
#grain_plot(args.potentials, args.element)
#diinterstitial_plot(args.potentials, args.element)

os.chdir(image_folder)
cwd = os.getcwd()

png_files = ""
for png_file in sorted(glob.glob("*.png")):
	png_files += "{0} ".format(png_file)

os.system("convert " + png_files + "output.pdf")
os.rename("output.pdf", "../../../test-runner/{0}-output.pdf".format(args.element))

print("--- %s seconds ---" % (time.time() - start_time))

