#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from processes import *;


def read_model_from_file(lib, model_file):

	def list_to_str(l, prefix, suffix):
		r = "[";
		first_flag = True;
		for e in l:
			if not first_flag: r = r + ",";
			r = r + prefix + str(e) + suffix;
			first_flag = False;
		r = r + "]";
		return r;

	def dict_to_str(d, prefix, suffix, nested = True):
		r = "{";
		first_flag = True;
		for k in d.keys():
			if not first_flag: r = r + ",";
			r = r + "\"" + k + "\":";
			if nested: r = r + list_to_str(d[k], prefix, suffix);
			else: r = r + str(d[k]);
			first_flag = False;
		r = r + "}";
		return r;

	from distutils.text_file import TextFile;
	from string import split, atof;
	from os import system;

	f = TextFile(filename = model_file);
	model_flag = True;
	processes = {};
	entities = {};
	processes_level = [];
	process_index = 0;
	for line in f.readlines():
		fields = line.split(" = ");

		if (len(fields) == 1) and model_flag:
			# process instance line
			ind_process = fields[0];

			# empty line: end of processes list
			if ind_process == "":
				model_flag = False;
				continue;

			# identify indentation level
			ind_list = ind_process.split("  ");
			(level, process) = (len(ind_list) - 1, ind_list[-1]);

			# identify generic process name and compose process instance ID
			process_index = process_index + 1;
			(process_name, parameters) = process.split("(");
			process_id = ("p_%05d" % process_index) + "_" + process_name;

			# keep the level-wise lists updated
			if level == len(processes_level): processes_level.append([]);
			if level <> 0:
				ancestor = processes_level[level-1][-1];
				processes[ancestor]["subprocesses"].append(process_id);

			# compose a process
			process = {
				"instance_of": process_name,
				"relates": {},
				"subprocesses": [],
				"parameters": {}
			}

			# take care of related entities
			relates = process["relates"];
			for parameter in parameters.split(", "):
				(parameter_id, entity_list) = parameter.split(":");
				relates[parameter_id] = [];
				if entity_list[-1] == ")": entity_list = entity_list[:-1];
				for entity_id in entity_list[1:-1].split(","):
					if entity_id == "": continue;
					relates[parameter_id].append(entity_id);

			processes_level[level].append(process_id);
			processes[process_id] = process;

			continue;

		# end of processes
		if (len(fields) == 2) and model_flag and (fields[0] == "SSE"):
			model_flag = False;
			continue;

		# entity/process constant parameter
		if (len(fields) == 2) and model_flag:
			(cp_name, cp_value) = fields;
			cp_name = cp_name.split("  ")[-1];
			
			check_ecp = cp_name.split(".");
			if len(check_ecp) == 2:
				# entity constant parameter
				(entity, cp_name) = check_ecp;

				if entity not in entities.keys():
					entities[entity] = {
						"parameters": {}
					}

				entities[entity]["parameters"][cp_name] = atof(cp_value);
				continue;

			# process constant parameter
			process["parameters"][cp_name] = atof(cp_value);
			continue;

		# init_state assignment
		if (len(fields) == 2) and (not model_flag) and (fields[0] == "init_state"):
			init_state_def = fields[0] + " = " + fields[1];
			continue;

	f.close();

	for eid in entities.keys():
		e = lib.entity_instance_by_id(eid);
		for cp_name in e.generic_entity.parameters.keys():
			if cp_name in entities[eid]["parameters"].keys():
				cp_value = entities[eid]["parameters"][cp_name];
				e.set_parameter_val(cp_name, cp_value);
			else:
				e.set_parameter_val(cp_name, None);

	ordered_keys = processes.keys()[:];
	ordered_keys.sort();
	ordered_keys.reverse();
	f = open("tmp_mfile.py", "w");
	for pid in ordered_keys:
		p = processes[pid];
		pdef = pid + " = process_instance(";

		pdef = pdef + "lib.generic_process_by_id(\"";
		pdef = pdef + p["instance_of"] + "\"), ";

		pdef = pdef + dict_to_str(p["relates"], "lib.entity_instance_by_id(\"", "\")") + ", ";
		pdef = pdef + list_to_str(p["subprocesses"], "", "") + ", ";
		pdef = pdef + dict_to_str(p["parameters"], "", "", nested = False) + ")";

		f.write(pdef + "\n");
	f.write(init_state_def + "\n");
	f.close();
	execfile("tmp_mfile.py");
	system("rm -f pfile.py");
	return (eval(ordered_keys[-1]), eval("init_state"));



def constants_table(lib, log_file, short_process_names_flag = True):

	def remove_spaces(l):
		while " " in l:
			i = l.index(" ");
			l = l[:i] + l[i + 1:];
		return l;


	from distutils.text_file import TextFile;
	from string import split, atof;
	from os import system;

	perf = ["r2", "reMSE"];

	f = TextFile(filename = log_file);
	cp_table = {};
	perf_table = {};
	model_index = 0;
	model_flag = False;
	prefix = [""];
	for line in f.readlines():
		fields = line.split(" = ");

		if (len(fields) == 1) and not model_flag and (line[:4] == "root"):
			model_flag = True;
			continue;

		if (len(fields) == 1) and model_flag:
			# process instance line
			ind_process = fields[0];

			# identify indentation level
			ind_list = ind_process.split("  ");
			(level, process) = (len(ind_list) - 1, ind_list[-1]);
			process = remove_spaces(process);

			if level > len(prefix):
				prefix.append(process);
			else:
				if level == len(prefix):
					prefix[-1] = process;
				else:
					prefix = prefix[:level];

			# print prefix;
			continue;

		# end of processes
		if (len(fields) == 2) and model_flag and (fields[0] == "SSE"):
			model_flag = False;
			prefix = [""];
			model_index = model_index + 1;
			continue;

		if (len(fields) == 2) and not model_flag and (fields[0] in perf):
			perf_table[(fields[0], model_index - 1)] = fields[1];
			continue;

		# entity/process constant parameter
		if (len(fields) == 2) and model_flag:
			(cp_name, cp_value) = fields;
			cp_name = cp_name.split("  ")[-1];

			if short_process_names_flag:
				if prefix[-1] == "":
					cp_full_name = "";
				else:
					cp_full_name = prefix[-1] + ".";
			else:
				cp_full_name = "";
				for p in prefix:
					if p <> "": cp_full_name = cp_full_name + p + ".";

			if "." in cp_name:
				# entity constant parameter
				(eid, cpid) = split(cp_name, ".");
				cp_range = lib.entity_instance_by_id(eid).generic_entity.parameters[cpid];
			else:
				if prefix[-1] == "": gpid = "root";
				else: (gpid, donotcare) = split(prefix[-1], "(");
				cp_range = lib.generic_process_by_id(gpid).parameters[cp_name];

			cp_full_name = cp_full_name + cp_name;

			cp_full_name = ("%05d" % len(prefix)) + "-" + cp_full_name;

			if cp_full_name not in cp_table.keys(): cp_table[cp_full_name] = {"range": cp_range, "values": {}};
			cp_table[cp_full_name]["values"][model_index] = atof(cp_value);
			continue;
	f.close();

	# header
	ordered_keys = cp_table.keys();
	ordered_keys.sort();

	print "parameter",;
	for cp_full_name in ordered_keys: print cp_full_name[6:],;
	for p in perf: print "perf-" + p,;
	print "";

	print "range_min",;
	for cp_full_name in ordered_keys: print cp_table[cp_full_name]["range"][0],;
	for p in perf: print "-",;
	print "";

	print "range_max",;
	for cp_full_name in ordered_keys: print cp_table[cp_full_name]["range"][1],;
	for p in perf: print "-",;
	print "";

	print "value_min",;
	for cp_full_name in ordered_keys: print min(cp_table[cp_full_name]["values"].values()),;
	for p in perf: print "-",;
	print "";

	print "value_max",;
	for cp_full_name in ordered_keys: print max(cp_table[cp_full_name]["values"].values()),;
	for p in perf: print "-",;
	print "";

	# constant parameter values
	for i in range(model_index):
		print i+1,;
		for cp_full_name in ordered_keys:
			if i in cp_table[cp_full_name]["values"].keys(): print cp_table[cp_full_name]["values"][i],;
			else: print "-",;
		for p in perf: print perf_table[(p, i)],;
		print "";

"""
	# header
	print "parameter range_min range_max value_min value_max",;
	for i in range(model_index): print i+1,;
	print "";

	# constant parameter values
	ordered_keys = cp_table.keys();
	ordered_keys.sort();
	for cp_full_name in ordered_keys:
		print cp_full_name[6:],;
		print cp_table[cp_full_name]["range"][0],;
		print cp_table[cp_full_name]["range"][1],;
		print min(cp_table[cp_full_name]["values"].values()),;
		print max(cp_table[cp_full_name]["values"].values()),;
		for i in range(model_index):
			if i in cp_table[cp_full_name]["values"].keys(): print cp_table[cp_full_name]["values"][i],;
			else: print "-",;
		print "";
	print "";

	# model performance
	for p in perf:
		print "perf-" + p,;
		print "- - - -",;
		for i in range(model_index): print perf_table[(p, i)],;
		print "";
"""
