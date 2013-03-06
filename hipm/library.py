#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from entities import *;
from processes import *;


def enumerate_sets(elements, cardinality):

	def improve(c):
		i = cardinality - 1;
		if i == -1: return 0;
		while c[i] == n - (cardinality - i): i = i - 1;
		if i == -1: return 0;

		c[i] = c[i] + 1;
		while i < cardinality - 1:
			c[i + 1] = c[i] + 1;
			i = i + 1;
		return 1;

	def c2e(c): return [elements[i] for i in c];

	n = len(elements);

	c = range(cardinality);
	sets = [c2e(c)];
	while improve(c): sets.append(c2e(c));
	return sets;


# library of generic processes and entities

class library:

	def __init__(self, id):

		self.id = id;
		self.generic_processes = [];
		self.generic_entities = [];

	
	def to_graph(self, file_name):

		if file_name[-4:] <> ".dot": file_name = file_name + ".dot";

		f = open(file_name, "w");
		f.write("graph lib_%s {\n" % self.id);
		f.write("\n");
		f.write("\tnode[shape = \"box\", font = \"Helvetica\"];\n");

		type_nodes = [];
		for gp in self.generic_processes:
			f.write("\n");
			f.write("\tnode_%s [label = \"%s\"];\n" % (gp.id, gp.id));
			if gp.type == "": continue;

			if gp.type not in type_nodes: f.write("\tnode_%s [label = \"%s\"];\n" % (gp.type, gp.type));
			type_nodes.append(gp.type);
			f.write("\tnode_%s -- node_%s [style = \"dashed\"];\n" % (gp.type, gp.id));

		for gp in self.generic_processes:
			if len(gp.processes) == 0: continue;
			f.write("\n");
			for (sp_id, relates, optional) in gp.processes:
				f.write("\tnode_%s -- node_%s;\n" % (gp.id, sp_id));

		f.write("}\n");
		f.close();


	def add_generic_process(self, id, type, relates, processes, parameters, oaes, odes):

		for gp in self.generic_processes:
			if (gp.id == id) and (gp.type == type):
				print "duplicate generic process id-type %s-%s." % (id, type);
		self.generic_processes.append(generic_process(self, id, type, relates, processes, parameters, oaes, odes));


	def add_generic_entity(self, id, variables, parameters):

		for ge in self.generic_entities:
			if ge.id == id:
				print "duplicate generic entity id %s." % id;
		nge = generic_entity(self, id, variables, parameters);
		self.generic_entities.append(nge);

		return nge;


	def generic_entity_by_id(self, id):

		for ge in self.generic_entities:
			if ge.id == id: return ge;
		return None;


	def entity_instance_by_id(self, id):

		for ge in self.generic_entities:
			for ei in ge.instances:
				if ei.id == id: return ei;
		return None;


	def entities_list_instances(self, ges_list):

		eis = [];
		for ge in ges_list: eis.extend(ge.instances);
		return eis;


	def entities_list_set_instances(self, ges_list, c_from, c_to):

		eis = self.entities_list_instances(ges_list);
		sets = [];
		for cardinality in range(max(0, c_from), min(c_to, len(eis)) + 1):
			sets.extend(enumerate_sets(eis, cardinality));
		return sets;


	def generic_process_by_id(self, id):

		for gp in self.generic_processes:
			if gp.id == id: return gp;
		return None;


	def generic_processes_by_type(self, type):

		gps = [];
		for gp in self.generic_processes:
			if (gp.id == type) or (gp.type == type):
				gps.append(gp);
		return gps;


	# translate generic process relates_instance into
	# relates_instance for the subprocesses of a given
	# type sp_type using the sp_relates_ids as a key

	# example
	#
	# sp_type = limitation_factor("X", "P")
	# sp_relates_ids = "N", "P" ("N" -> "X", "P" -> "P")
	# relates_instance = {"P": [p1], "N": [n1, n2, n3]}
	# return = {"P": [p1], "X": [n1, n2, n3]}

	def subprocess_type_relates_instance(self, sp_type, sp_relates_ids, relates_instance):

		new_relates_instance = 	{};
		for i in range(len(sp_relates_ids)):
			sp_relates_id = sp_relates_ids[i];
			if sp_relates_id not in relates_instance: continue;

			for sp in self.generic_processes_by_type(sp_type):
				if i >= len(sp.relates): continue;
				r_sp_relates_id = sp.relates[i][0];
				new_relates_instance[r_sp_relates_id] = relates_instance[sp_relates_id];

		return new_relates_instance;


	# return a list of all process instances of the given generic process type
	# parameters:
	#  - type: the given process type
	#  - relates_instance: process relates entity instances

	def process_type_instances(self, type, relates_instance):

		def improve(c, list):

			i = len(c) - 1;
			if i == -1: return 0;

			while c[i] == len(list[i]) - 1:
				c[i] = 0;
				i = i - 1;
				if i == -1: return 0;
			c[i] = c[i] + 1;
			return 1;

		def c2pis_list(c):

			ri = {};
			for i in range(len(ris_list)):
				(key, value) = ris_list[i][c[i]].items()[0];
				ri[key] = value;

			pis_list = [];
			for gp in gp_list: pis_list.extend(gp.instances(ri));
			return pis_list;

		def c2pis(c): return [pis_lists[i][c[i]] for i in range(len(pis_lists))];


		# get all processes of the given type
		# group them by compatibility profile

		gp_groups = [];
		for gp in self.generic_processes_by_type(type):
			ris_list = gp.relates_instances_list(relates_instance);
			if ris_list == None: continue;

			found_flag = 0;
			for [p_ris_list, gp_list] in gp_groups:
				if ris_list <> p_ris_list: continue;
				gp_list.append(gp);
				found_flag = 1;
				break;

			if found_flag == 0: gp_groups.append([ris_list, [gp]]);

		type_pis = [];
		# for each compatibility profile group collect the process instances
		for (ris_list, gp_list) in gp_groups:

			c = [0 for ris in ris_list];
			pis_lists = [c2pis_list(c)];
			while improve(c, ris_list): pis_lists.append(c2pis_list(c));
			while [] in pis_lists: pis_lists.remove([]);

			c = [0 for pis in pis_lists];
			type_pis.append(c2pis(c));
			while improve(c, pis_lists): type_pis.append(c2pis(c));

		return type_pis;


	def simplest_process_type_instances(self, type, relates_instance):

		def improve(c, list):

			i = len(c) - 1;
			if i == -1: return 0;

			while c[i] == len(list[i]) - 1:
				c[i] = 0;
				i = i - 1;
				if i == -1: return 0;
			c[i] = c[i] + 1;
			return 1;

		def c2pis_list(c):

			ri = {};
			for i in range(len(ris_list)):
				(key, value) = ris_list[i][c[i]].items()[0];
				ri[key] = value;

			pis_list = [];
			for gp in gp_list: pis_list.extend(gp.simplest_instances(ri));
			return pis_list;

		def c2pis(c): return [pis_lists[i][c[i]] for i in range(len(pis_lists))];

		# get all processes of the given type
		# group them by compatibility profile

		gp_groups = [];
		for gp in self.generic_processes_by_type(type):
			ris_list = gp.relates_instances_list(relates_instance);
			if ris_list == None: continue;

			found_flag = 0;
			for [p_ris_list, gp_list] in gp_groups:
				if ris_list <> p_ris_list: continue;
				gp_list.append(gp);
				found_flag = 1;
				break;

			if found_flag == 0: gp_groups.append([ris_list, [gp]]);

		type_pis = [];
		# for each compatibility profile group collect the process instances
		for (ris_list, gp_list) in gp_groups:

			c = [0 for ris in ris_list];
			pis_lists = [c2pis_list(c)];
			while improve(c, ris_list): pis_lists.append(c2pis_list(c));
			while [] in pis_lists: pis_lists.remove([]);

			c = [0 for pis in pis_lists];
			type_pis.append(c2pis(c));
			while improve(c, pis_lists): type_pis.append(c2pis(c));
		if type_pis == []: type_pis.append([]);

		return type_pis;


	def simplest_process_type_instances_pl(self, type, relates_instance):

		def improve(c, list):

			i = len(c) - 1;
			if i == -1: return 0;

			while c[i] == len(list[i]) - 1:
				c[i] = 0;
				i = i - 1;
				if i == -1: return 0;
			c[i] = c[i] + 1;
			return 1;

		def c2pis_list(c):

			ri = {};
			for i in range(len(ris_list)):
				(key, value) = ris_list[i][c[i]].items()[0];
				ri[key] = value;

			pis_list = [];
			for gp in gp_list: pis_list.extend(gp.simplest_instances(ri));
			return pis_list;

		# get all processes of the given type
		# group them by compatibility profile

		gp_groups = [];
		for gp in self.generic_processes_by_type(type):
			ris_list = gp.relates_instances_list(relates_instance);
			if ris_list == None: continue;

			found_flag = 0;
			for [p_ris_list, gp_list] in gp_groups:
				if ris_list <> p_ris_list: continue;
				gp_list.append(gp);
				found_flag = 1;
				break;

			if found_flag == 0: gp_groups.append([ris_list, [gp]]);

		type_pis = [];
		# for each compatibility profile group collect the process instances
		for (ris_list, gp_list) in gp_groups:
			c = [0 for ris in ris_list];
			type_pis.extend(c2pis_list(c));
			while improve(c, ris_list): type_pis.extend(c2pis_list(c));
		return type_pis;
