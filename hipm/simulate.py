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


	from string import split, atof;

	f = open(model_file);
	model_flag = True;
	processes = {};
	entities = {};
	processes_level = [];
	process_index = 0;
	for line in f.readlines():
		fields = line[:-1].split(" = ");

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
		if (len(fields) == 1) and model_flag and (fields[0] == ""):
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
	f = open("pfile.py", "w");
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
	execfile("pfile.py");
	return (eval(ordered_keys[-1]), eval("init_state"));
