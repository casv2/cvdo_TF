import argparse
import time
import json
import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt


start_time = time.time()

#setting colour types
struct_colors = [ "black", "red", "blue", "cyan", "orange", "magenta", "green", "grey", "brown" ]

def pct_change(df):

	for column in df.columns[2:]:
		print column
		print df.columns
		df[column] = (df[column] - df[df.columns[1]])/df[df.columns[1]] * 100
		df[column] = [str(int(df[column].tolist()[i]))+"%" for i in xrange(0,len(df[column].tolist()))]

	df = df.round(2)
	return df

#parsing arguments
parser = argparse.ArgumentParser(description='analyse specific tests')
parser.add_argument('--surface', '-s', action='store_true', help='do surface')
parser.add_argument('--bulk', '-b', action='store_true', help='do E/V bulks')
parser.add_argument('')
args = parser.parse_args()

surface = args.surface
bulk = args.bulk

data_dict = { '-' : []}

for json_file in sorted(glob.glob("*.json")):
	#load data
	json_data = json.load(open(json_file))

	#extract model/file_names
	model_name = json_file.split('-')[2]
	file_name = json_file.split('-')[4]

	#add model_name to data dictionary
	if model_name not in data_dict.keys():
		data_dict[model_name] = []

	if surface and "surface" in json_file:
		data_dict[model_name].append(json_data['Ef'])
		if file_name not in data_dict['-']:
			data_dict['-'].append(file_name)

	if bulk and "bulk" in json_file:
		for data in json_data:
			if type(json_data[data]) == float:
				data_dict[model_name].append(json_data[data])
				if "{0} {1}".format(file_name, data) not in data_dict["-"]:
					data_dict['-'].append("{0} {1}".format(file_name, data))
			if type(json_data[data]) == list:
				E = [json_data[data][i][0] for i in xrange(0,len(json_data[data]))]
				V = [json_data[data][i][1] for i in xrange(0,len(json_data[data]))]
				plt.plot(E,V)


df = pd.DataFrame.from_dict(data_dict)
#print len(data_dict["SW"])
#print len(data_dict["GAP_6"])
#print len(data_dict["-"])
print pct_change(df)
