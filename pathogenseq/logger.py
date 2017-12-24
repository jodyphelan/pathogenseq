import json

def log(key,log_file):
	x = json.load(open(log_file))
	x[key] = True
	json.dump(open(log_file,"w"),x)

def check_log(key,log_file):
	x = json.load(open(log_file))
	return True if key in x else False
