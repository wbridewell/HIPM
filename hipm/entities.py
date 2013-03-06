#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.


# entities and their instantiations


# specific entity or generic entity instance contains:
#
#  . ge - reference to the generic entity
#  . id - entity id
#  . v_vals - measured values of the entity variables
#  . p_vals - values of the entity parameters

class entity_instance:

	def __init__(self, ge, id, variables, parameter_vals):

		self.generic_entity = ge;
		self.id = id;

		self.variables = {};

		if variables == None:
			for v in ge.variables.keys(): self.variable_vals[v] = (None, None, None, (None, None));
		else:
			for v in ge.variables.keys():
				if v in variables.keys():
					(vtype, vdata, vinitrange) = variables[v];

					if isinstance(vdata, str) == 0:
						vdataid = None;
						vinit = vdata;
					else:
						vdataid = vdata;
						vinit = None;

					if vinitrange == None:
						vlower = None;
						vupper = None;
					else:
						(vlower, vupper) = vinitrange;

					self.variables[v] = (vtype, vdataid, vinit, (vlower, vupper));
				else:
					self.variables[v] = (None, None, None, (None, None));

		self.parameter_vals = {};

		if parameter_vals == None:
			for p in ge.parameters.keys(): self.parameter_vals[p] = None;
		else:
			for p in ge.parameters.keys():
				if p in parameter_vals.keys():
					self.parameter_vals[p] = parameter_vals[p];
				else:
					self.parameter_vals[p] = None;

		ge.instances.append(self);

	def set_parameter_val(self, p, v): self.parameter_vals[p] = v;

	def variable_id(self, v): return "%s_%s" % (self.id, v);

	def parameter_id(self, p): return "econst_%s_%s" % (self.id, p);

	def parameter_val(self, p): return self.parameter_vals[p];

	def __str__(self): return self.id;



# generic entity class contains:
#
#  . lib - reference to the library that generic entity belongs to
#  . id - entity identifyer
#  . variables - variables of the entity along with their aggregation functions
#  . parameters - constant parameters of the entity

class generic_entity:

	def __init__(self, lib, id, variables, parameters):

		self.lib = lib;
		self.id = id;

		self.variables = variables;
		self.parameters = parameters;

		self.instances = [];

	def aggregation(self, var): return self.variables[var];

	def __str__(self): return self.id;
