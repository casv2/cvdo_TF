import json

data = json.load(open("./info.json"))

class deserialise(object):
	def __init__(self,j):
		self.__dict__ = json.loads(j)

p = deserialise(data)

print p.numconfigs