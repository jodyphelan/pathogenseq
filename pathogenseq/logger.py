import json

def log(key,log_file):
	x = json.load(open(log_file))
	x[key] = True
	json.dump(open(log_file,"w"),x)

def checkpoint(key,log_file):
	x = json.load(open(log_file))
	return False if key in x else True
