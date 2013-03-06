#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

import sys


# processes and their instantiations


# process_instance class has the following fields:
#
#  . generic_process - pointer to the instantiated generic process
#  . relates - (dictionary) entities involved in the process
#  . process_instances - (list) pointers to sub process instances
#  . parameter_vals - (dictionary) values of the constant parameters
#
# examples of field values
#  . relates = {"P": ["p1", "p2"], N: ["n1"], E: []}
#  . parameter_vals = {"lambda": 0.5}
#
# input parameters are the same as the fields

# replacese all occurrences of the "old" string with the "new" string 
# in "str"
def replace_str(old, new, str):
	returned = 0
	fromIndex = 0
	while returned <> -1:
		# find the next occurrence of the old string
		returned = str.find(old, fromIndex)
		if returned <> -1:
			passLeft = False
			passRight = False
			# check the left side
			if returned - 1 > 0:
				if str[returned-1].isalnum() or str[returned-1] == "_":
					fromIndex = returned + len(old)
					continue
				else:
				       	passLeft = True
			else:
				passLeft = True
			# check the right side
			if returned + len(old) < len(str):
				if str[returned + len(old)].isalnum() or str[returned + len(old)] == "_":
					fromIndex = returned + len(old)
					continue
				else:
					passRight = True
			else:
				passRight = True
			# we have a match
			if passLeft and passRight:
				str = str[0:returned] + new + str[returned+len(old):]
				fromIndex = returned + len(new)
	return str

class process_instance:

	def __init__(self, gp, relates, process_instances, parameter_vals):

		self.generic_process = gp;
		self.relates = relates;
		self.process_instances = process_instances;

		self.parameter_vals = {};
		if parameter_vals == None:
			for p in gp.parameters.keys(): self.parameter_vals[p] = None;
		else:
			for p in gp.parameters.keys():
				if p in parameter_vals.keys():
					self.parameter_vals[p] = parameter_vals[p];
				else:
					self.parameter_vals[p] = None;


	def same_type_structure(self, pi):

		if self.generic_process.type <> pi.generic_process.type: return(0);
		if self.generic_process.type == "": return(0);
		if pi.generic_process.type == "": return(0);
		for r in self.relates.keys():
			if r not in pi.relates.keys(): continue;
			if len(self.relates[r]) <> len(pi.relates[r]): return(0);
			for ei in self.relates[r]:
				if ei not in pi.relates[r]: return(0);
		return(1);
	

	def same_structure(self, pi):
		if self.generic_process <> pi.generic_process: return(0);

		if len(self.relates.keys()) <> len(pi.relates.keys()): return(0);
		for r in self.relates.keys():
			if r not in pi.relates.keys(): return(0);
			if len(self.relates[r]) <> len(pi.relates[r]): return(0);
			for ei in self.relates[r]:
				if ei not in pi.relates[r]: return(0);

		if len(self.process_instances) <> len(pi.process_instances): return(0);
		for spi1 in self.process_instances:
			found_flag = 0;
			for spi2 in pi.process_instances:
				if spi1.same_structure(spi2):
					found_flag = 1;
					break;
			if found_flag == 0: return(0);

		return(1);

	# set parameter values to the ones returned by the fit procedure below

	def _set_parameter_vals_(self, parameters, parameter_vals, sp_flag):

		def parameter_value(pid):
			## xxx: in some cases the structure cannot be parameterized at all
			##      we have to just ignore and move on
			if parameter_vals == {}: return None;
			if parameters == None: return None;

			for i in range(len(parameters)):
				(p, prange) = parameters[i];
				if p == pid: return parameter_vals[i];
			return None;

		if not sp_flag:
			for (ent_id, ei_list) in self.relates.items():
				for ei in ei_list:
					for (p, prange) in ei.generic_entity.parameters.items():
						pid = ei.parameter_id(p);
						ei.parameter_vals[p] = parameter_value(pid);

		for (p, prange) in self.generic_process.parameters.items():
			pid = self.parameter_id(p);
			self.parameter_vals[p] = parameter_value(pid);

		for pi in self.process_instances:
			pi._set_parameter_vals_(parameters, parameter_vals, 1);

	def set_parameter_vals(self, parameters, parameter_vals):
		self._set_parameter_vals_(parameters, parameter_vals, 0);


	def parameter_id(self, p):
		id = "pconst";
		for (ent, ei_list) in self.relates.items():
			id = id + "_" + self.generic_process.id + "_" + ent;
			for ei in ei_list: id = id + "_" + str(ei);
		id = id + "_" + p;
		return id;

	def parameter_val(self, p): return self.parameter_vals[p];


	def number_of_process_instances(self):

		npi = len(self.process_instances);
		for pi in self.process_instances:
			npi = npi + pi.number_of_process_instances();
		return npi;


	def __str__(self): return self._str_(0);

        # print the current process
	# sp_flag==0 for the root process, which indicates that the 
	# values of the entity parameters should be printed.
	def _str_(self, sp_flag):
		from re import compile, sub;

		first_flag = 1;
		
		# the name of the process is equivalent to the gp name
		str = "%s(" % self.generic_process.id;
		# print the entities participating in the process
		for id in self.relates.keys():
			if first_flag == 0: str = str + ", ";
			first_flag = 0;
                        # the role name
			str = str + ("%s:" % id);
			# the entity filling the role
			# either blank or the name of an entity instance
			if len(self.relates[id]) == 0:
				str = str + "{}";
			else:
				str = str + ("{%s" % self.relates[id][0].__str__());
				for ei in self.relates[id][1:]: str = str + "," + ei.__str__();
				str = str + "}";
		str = str + ")";

		# print values of the entity parameters only for the "root" process instance
		if not sp_flag:
			for id in self.relates.keys():
				for ei in self.relates[id]:
					for (p, pval) in ei.parameter_vals.items():
						if pval <> None:
							str = str + ("\n  %s.%s = %g" % (ei.__str__(), p, pval));

               	# print the parameter values for the process
		for (p, pval) in self.parameter_vals.items():
			if pval <> None: str = str + ("\n  %s = %g" % (p, pval));

                # return the current string if there are no subprocesses
		if len(self.process_instances) == 0: return str;
                
                # otherwise, recurse on the subprocesses
		for pi in self.process_instances:
			spi = pi._str_(1);
			nl = compile("\n");
			spi = nl.sub("\n  ", spi);
			str = str + "\n  " + spi;

		return str;

        # print the current process in an IPM-friendly format
        # we assume that this method is only called on the root process
	def _str_ipm_(self, name = None, unique_id = None, sse = None, init_state = None):

		# print the supplied name or the name of the root generic process
		if name == None:
			str = "model %s" % self.generic_process.id
		else:
			str = "model %s" % name

		# if we have a unique id, which is helpful when printing multiple models,
		# make sure to include it as part of the model's name
		if unique_id <> None:
			str = str + "_%s" % unique_id

		# include the SSE score if given
		if sse == None:
			str = str + ""
		else:
			str = str + "{%g}" % sse
		str = str + ";\n"
		
		# variables
		# we assume that the role name is equivalent to the variable name
		# we assume that entities only have one variable (called "value")
		# gather a list of the variables and their types, the observed variables,
		# and the exogenous variables
		# variables = [[name, type, obs, exo], ]
		variables = []
		exvars = []
		obsvars = []
		for id in self.relates.keys():
			if len (self.relates[id]) <> 0:
				cur_ent = self.relates[id][0];
		    	# get the information about the one variable in this entity
				cur_var = cur_ent.variables["value"];
				if cur_var[0] == "exogenous":
					exvars.append(cur_ent)
				# if the there is no data file reference, then it's not observed
				if cur_var[1] <> None:
					obsvars.append(cur_ent)
				variables.append([cur_ent.id, cur_ent.generic_entity.id])
	    
		# now print the variable declarations
		first_flag = 1;
		for v in variables:
			if first_flag <> 0:
				str = str + "variables ";
			else:
				str = str + ",";
			first_flag = 0;
			str = str + ("%s{%s}" % (v[0], v[1]));
		if len(variables) > 0: str = str + ";\n"
	    
		# print the list of observable variables
		first_flag = 1;
		for v in obsvars:
			if first_flag <> 0:
				str = str + "observable ";
			else:
				str = str + ",";
			first_flag = 0;
			str = str + ("%s" % v)
		if len(obsvars) > 0: str = str + ";\n"
                    
		# print the list of exogenous variables
		first_flag = 1;
		for v in exvars:
			if first_flag <> 0:
				str = str + "exogenous ";
			else:
				str = str + ",";
			first_flag = 0;
			str = str + ("%s" % v)
		if len(exvars) > 0: str = str + ";\n"
		
		str = str + "\n"
		
		# the header has been printed, so we just need to print the processes
		# we assume that changed processes will have "_adjusted" at the end of their
		# generic process name, whereas fixed processes will have "_tmp" at the end
		# of theirs.  The names of changed and fixed processes should be okay, but
		# we may need to add a unique identifier for the added processes.
		for pi in self.process_instances:
			pi_name = pi.generic_process.id;
			pi_modifier = ""
			gp = ""
			if pi_name.endswith("_tmp"):
				pi_modifier = "fix"
				pi_name = pi_name[0:len(pi_name)-len("_tmp")]
				gp = "none" # can't get the real generic process name here
			elif pi_name.endswith("_adjusted"):
				pi_modifier = "change"
				pi_name = pi_name[0:len(pi_name)-len("_adjusted")]
				gp = "none" # can't get the real generic process name here
			else:
				pi_modifier = "add"
				gp = pi.generic_process.id
			
			# print process header
			str = str + "process %s{%s,%s};\n" % (pi_name, gp, pi_modifier);
		
			# print the equations
			str = str + "equations\n"
			for lhs in pi.generic_process.oaes.keys():
				# remove the ".value" after entity names
				rhs = pi.generic_process.oaes[lhs]
				lhs = lhs.replace(".value", "")
				rhs = rhs.replace(".value", "")
				# replace the parameter names with the actual values
				for p in pi.generic_process.parameters.keys():
					rhs = replace_str(p, "%f" % pi.parameter_vals[p], rhs)
				# replace the role names with the variable names
				for ei in pi.relates.keys():
					# make sure the names are different and need to be changed,
					# else you get an infinite loop
					if ei <> pi.relates[ei][0].id:
						lhs = replace_str(ei, pi.relates[ei][0],lhs)
						rhs = replace_str(ei, pi.relates[ei][0],rhs)
				str = str + "%s = %s;\n" % (lhs, rhs)
                
			for lhs in pi.generic_process.odes.keys():
				# remove the ".value" after entity names
				rhs = pi.generic_process.odes[lhs]
				rhs = rhs.replace(".value", "")
				lhs = lhs.replace(".value", "")
				# replace the parameter names with the actual values
				for p in pi.generic_process.parameters.keys():
					rhs = replace_str(p, "%f" % pi.parameter_vals[p],rhs)
				# replace the role names with the variable names
				for ei in pi.relates.keys():
					# make sure the names are different and need to be changed,
					# else you get an infinite loop
					if ei <> pi.relates[ei][0].id:
						lhs = replace_str(ei, pi.relates[ei][0].id,lhs)
						rhs = replace_str(ei, pi.relates[ei][0].id,rhs)
				str = str + "d[%s,t,1] = %s;\n" % (lhs, rhs);
			
			str = str + "\n"
            
		return str;
            
	def add_process(self):

		def inlist(pi, pis, type_flag):
			for c_pi in pis:
				if type_flag and c_pi.same_type_structure(pi): return(1);
				if c_pi.same_structure(pi): return(1);
			return(0);

		gp = self.generic_process;
		lib = gp.lib;

		candidates = [];
		for (sp_type, sp_relates_ids, sp_optional) in gp.processes:
			if sp_optional == 0: continue;

			sp_ri = lib.subprocess_type_relates_instance(sp_type, sp_relates_ids, self.relates);
			for spi in lib.simplest_process_type_instances_pl(sp_type, sp_ri):
				if inlist(spi, self.process_instances, 1): continue;
				if inlist(spi, candidates, 0): continue;
				candidates.append(spi);

		pis = [];
		for c in candidates:
			pis.append(process_instance(self.generic_process, self.relates, self.process_instances + [c], None));

		for i in range(len(self.process_instances)):
			spi = self.process_instances[i];
			for nspi in spi.add_process():
				nspis = self.process_instances[:i] + [nspi] + self.process_instances[i+1:];
				pis.append(process_instance(self.generic_process, self.relates, nspis, None));

		return pis;


	def _to_eqns_(self):

		def union(l1, l2):

			u = l1;
			for e in l2:
				if e not in u: u.append(e);
			return u;

		# for each equation (ODE, OAE)
		#  . instantiate LHS
		#  . find out what is the aggregation function for the LHS
		#  . instantiate RHS
		#  . store LHS, aggregation, RHS, TYPE = ODE|OAE

		def instantiate_eqn(lhs, rhs, type):
			from string import split;
			from re import compile, sub;

			# replace ent_id.param_var -> param_var_id in rhs_inst

			def instatiate_ent(rhs_inst, param_var, param_var_id):
				new_rhs_inst = replace_str(param_var, param_var_id, rhs_inst)
				return (new_rhs_inst, new_rhs_inst <> rhs_inst);


			eqns_inst = {};
			parameters = [];
			variables = [];

			(lhs_ent_id, lhs_v) = lhs.split(".");

			for lhs_ei in self.relates[lhs_ent_id]:
				# if the variable on the LHS is exogenous do not instantiate
				(lhs_vtype, lhs_vdataid, lhs_vinit, lhs_vinitrange) = lhs_ei.variables[lhs_v];
				if lhs_vtype == "exogenous": continue;

				lhs_ge = lhs_ei.generic_entity;
				lhs_aggr = lhs_ge.variables[lhs_v];
	
				lhs_inst = lhs_ei.variable_id(lhs_v);
				rhs_inst = rhs;

				# replace the LHS entity with its instance
				# do that for entity parameters ...
				for (p, (lower, upper)) in lhs_ge.parameters.items():
					pid = lhs_ei.parameter_id(p);
					pval = lhs_ei.parameter_val(p);
					(rhs_inst, changed) = instatiate_ent(rhs_inst, lhs_ent_id + "." + p, pid);
					if changed: parameters = union(parameters, [(pid, (lower, upper, pval))]);

				# ... and entity variables
				for (v, vaggr) in lhs_ge.variables.items():
					vid = lhs_ei.variable_id(v);
					(vtype, vdataid, vinit, vinitrange) = lhs_ei.variables[v];
					(rhs_inst, changed) = instatiate_ent(rhs_inst, lhs_ent_id + "." + v, vid);
					if changed: variables = union(variables, [(vid, vdataid, vinit, vinitrange)]);

				# replace other process relates entities with their instances
				for (ent_id, ei_list) in self.relates.items():
					if ent_id == lhs_ent_id: continue;
					if len(ei_list) <> 1: continue;

					ei = ei_list[0];
					ge = ei.generic_entity;

					for (p, (lower, upper)) in ge.parameters.items():
						pid = ei.parameter_id(p);
						pval = ei.parameter_val(p);
						(rhs_inst, changed) = instatiate_ent(rhs_inst, ent_id + "." + p, pid);
						if changed: parameters = union(parameters, [(pid, (lower, upper, pval))]);

					for (v, vaggr) in ge.variables.items():
						vid = "datamin_" + ei.variable_id(v);
						(vtype, vdataid, vinit, vinitrange) = ei.variables[v];
						(rhs_inst, changed) = instatiate_ent(rhs_inst, "datamin(" + ent_id + "." + v + ")", vid);
						if changed: variables = union(variables, [(vid, vdataid, vinit, vinitrange)]);

					for (v, vaggr) in ge.variables.items():
						vid = "datamax_" + ei.variable_id(v);
						(vtype, vdataid, vinit, vinitrange) = ei.variables[v];
						(rhs_inst, changed) = instatiate_ent(rhs_inst, "datamax(" + ent_id + "." + v + ")", vid);
						if changed: variables = union(variables, [(vid, vdataid, vinit, vinitrange)]);

					for (v, vaggr) in ge.variables.items():
						vid = ei.variable_id(v);
						(vtype, vdataid, vinit, vinitrange) = ei.variables[v];
						(rhs_inst, changed) = instatiate_ent(rhs_inst, ent_id + "." + v, vid);
						if changed: variables = union(variables, [(vid, vdataid, vinit, vinitrange)]);

				# and take care about "local" process parameters
				for (p, (lower, upper)) in gp.parameters.items():
					pid = self.parameter_id(p);
					pval = self.parameter_val(p);
					new_rhs_inst = replace_str(p, pid, rhs_inst);
					
					if new_rhs_inst == rhs_inst: continue;

					rhs_inst = new_rhs_inst;
					parameters = union(parameters, [(pid, (lower, upper, pval))]);

				eqns_inst[lhs_inst] = (rhs_inst, type, lhs_aggr, lhs_vdataid, lhs_vinit, lhs_vinitrange);

			return (eqns_inst, parameters, variables);


		# aggregate equqtions that have the same variable on the LHS
		#  this functions deals with the equations for the current process

		def aggregate_eqns(eqns, parameters, variables, eqns_inst, params, vars):

			for (v, (rhs_inst, eq_type, aggr, vdataid, vinit, vinitrange)) in eqns_inst.items():
				if v not in odes.keys():
					eqns[v] = ["(" + rhs_inst + ")", aggr, vdataid, vinit, vinitrange];
				else:
					[curr_rhs_inst, curr_aggr, curr_vdataid, curr_vinit, curr_vinitrange] = eqns[v];
					if aggr <> curr_aggr: print "***", aggr, currr_aggr;
					if aggr == "prod": curr_rhs_inst = curr_rhs_inst + " * (" + rhs_inst + ")";
					if aggr == "sum": curr_rhs_inst = curr_rhs_inst + " + (" + rhs_inst + ")";
					if aggr == "min": curr_rhs_inst = "min(" + rhs_inst + ", " + curr_rhs_inst + ")";
					if aggr == "max": curr_rhs_inst = "max(" + rhs_inst + ", " + curr_rhs_inst + ")";
					eqns[v] = [curr_rhs_inst, curr_aggr, curr_vdataid, curr_vinit, curr_vinitrange];

			parameters = union(parameters, params);
			variables = union(variables, vars);

			return (eqns, parameters, variables);

		# and this function deals with the euqaitons for the subprocesses
		#  Q: why we need another function?
		#  A: by the time i wrote this comment - i forgot

		def aggregate_sp_eqns(eqns, sp_eqns):

			for (v, [rhs_inst, aggr, vdataid, vinit, vinitrange]) in sp_eqns.items():
				if v not in eqns.keys():
					eqns[v] = [rhs_inst, aggr, vdataid, vinit, vinitrange];
				else:
					[curr_rhs_inst, curr_aggr, curr_vdataid, curr_vinit, curr_vinitrange] = eqns[v];
					if aggr <> curr_aggr: print "***", aggr, curr_aggr;
					if aggr == "prod": curr_rhs_inst = curr_rhs_inst + " * (" + rhs_inst + ")";
					if aggr == "sum": curr_rhs_inst = curr_rhs_inst + " + (" + rhs_inst + ")";
					if aggr == "min": curr_rhs_inst = "min(" + rhs_inst + ", " + curr_rhs_inst + ")";
					if aggr == "max": curr_rhs_inst = "max(" + rhs_inst + ", " + curr_rhs_inst + ")";
					eqns[v] = [curr_rhs_inst, curr_aggr, curr_vdataid, curr_vinit, curr_vinitrange];

			return eqns;

		# main _to_eqns_

		gp = self.generic_process;

		odes = {};
		oaes = {};
		parameters = [];
		variables = [];

		for (lhs, rhs) in gp.odes.items():
			(odes_inst, params, vars) = instantiate_eqn(lhs, rhs, "D");
			(odes, parameters, variables) = aggregate_eqns(odes, parameters, variables, odes_inst, params, vars);

		for (lhs, rhs) in gp.oaes.items():
			(oaes_inst, params, vars) = instantiate_eqn(lhs, rhs, "A");
			(oaes, parameters, variables) = aggregate_eqns(oaes, parameters, variables, oaes_inst, params, vars);

		for pi in self.process_instances:
			(sp_odes, sp_oaes, sp_params, sp_vars) = pi._to_eqns_();

			odes = aggregate_sp_eqns(odes, sp_odes);
			oaes = aggregate_sp_eqns(oaes, sp_oaes);

			parameters = union(parameters, sp_params);
			variables = union(variables, sp_vars);

		return (odes, oaes, parameters, variables);


	# takes care that EACH system variable appears on the LHS of AT LEAST ONE equation

	def to_eqns(self):

		(odes, oaes, parameters, variables) = self._to_eqns_();
		for ent in self.generic_process.lib.generic_entities:
			for ei in ent.instances:
				for (v, (vtype, vdataid, vinit, vinitrange)) in ei.variables.items():
					if vtype <> "system": continue;
					if vdataid == None: continue;
					vid = ei.variable_id(v);
					if vid in odes.keys(): continue;
					odes[vid] = ['(double) 0.0 * %s' % vid, 'sum', vdataid, vinit, vinitrange];
		return (odes, oaes, parameters, variables);


	# transforms the process instance to C code for parameter fiting
	# we need to causally order the OAEs before spitting them out as C code
	# the current solution is a bit hacky
	###
	# keep track of calculated vars
	#  sys variables and exogenous variables are automatically added to the list
	# if all vars on the RHS of an equation are already calculated, then output the equation
	# keep going until all algebraic-only variables are calculated
	#
	# do a string search for all the variables in the model and
	# ensure that those matching the RHS are in the OK'd list before printing
	###
	# assume: all data sets have the same variables in the same order
	###
	def to_c_model(self, data_sets, normalize_error, init_state_fit, init_state, n_tf_restarts, n_fs_restarts):
		from re import compile, sub;

		ms_init_dims = "void MS_init_dims(void) {\n\n";
		ms_fill_in = "void MS_fill_in(void) {\n\n";
		ms_params = "void MS_params(void) {\n\n";
		ms_model = "void MS_model(double t) {\n\n";

		# the list of "already calculated variables" for ordering the OAEs
		calc_vars = []
		# a complete list of variables in the model (variables is not complete)
		all_vars = []

		# assumes that all data sets have the same variables in the same order		
		ms_init_dims = ms_init_dims + ("\tdata.nsets = %d;\n" % len(data_sets));
		data_length = 0;
		for i in range(len(data_sets)):
			ds = data_sets[i];
			data_length = data_length + ds.length;
			nvars = len(ds.vars);

			ms_fill_in = ms_fill_in + ("\tdata.set_files[%d] = strdup(\"%s\");\n" % (i, ds.file));
			ms_fill_in = ms_fill_in + ("\tdata.set_tos[%d] = %d;\n" % (i, data_length));
		ms_fill_in = ms_fill_in + "\n";
		
		ms_init_dims = ms_init_dims + ("\tdata.nvars = %d;\n" % len(data_sets[0].vars));
		ms_init_dims = ms_init_dims + ("\tdata.length = %d;\n\n" % data_length);

		(odes, oaes, parameters, variables) = self.to_eqns();

		ms_init_dims = ms_init_dims + ("\tsys_var.n = %d;\n" % len(odes.keys()));
		ms_init_dims = ms_init_dims + ("\tparam.n = %d;\n\n" % len(parameters));

		# set parameters of the fit code
		ms_init_dims = ms_init_dims + ("\topt.normalize_error = %d;\n" % normalize_error);
		ms_init_dims = ms_init_dims + ("\topt.init_state_fit = %d;\n" % init_state_fit);
		ms_init_dims = ms_init_dims + ("\topt.tf_restarts = %d;\n" % n_tf_restarts);
		ms_init_dims = ms_init_dims + ("\topt.fs_restarts = %d;\n\n" % n_fs_restarts);

		# set parameter initial values and bounds
		for i in range(len(parameters)):
			(pid, (lower, upper, pval)) = parameters[i];
			if pval == None: pval = (lower + upper) * 0.5;
			ms_params = ms_params + ("\tparam.vals[%d] = (double) %g;\n" % (i, pval));
			ms_params = ms_params + ("\tparam.bounds[%d] = (double) %g;\n" % (2 * i, lower));
			ms_params = ms_params + ("\tparam.bounds[%d] = (double) %g;\n\n" % (2 * i + 1, upper));

		ms_fill_in = ms_fill_in + ("\tdata.time = data.table[%d];\n" % 0);
		ms_fill_in = ms_fill_in + ("\tdata.sys_var_names[%d] = strdup(\"%s\");\n" % (0, data_sets[0].vars[0]));

		# take care about exogenous variables first
		model_eqns = "";
		
		for (v, vdataid, vinit, vinitrange) in variables:
			if (v in odes.keys()) or (v in oaes.keys()): continue;
			# consider exogenous variables OK
			all_vars.append(v)
			calc_vars.append(v)
			vi = ds.var_index(vdataid);
			if vi == -1:
				#print "Warning: Unobserved exogenous variable %s." % v;
				model_eqns = model_eqns + ("\tdouble %s = (double) 0;\n" % v);
			else:
				if v[:8] == "datamin_":
					model_eqns = model_eqns + ("\tdouble %s = data.min[%d];\n" % (v, vi));
				elif v[:8] == "datamax_":
					model_eqns = model_eqns + ("\tdouble %s = data.max[%d];\n" % (v, vi));
				else: model_eqns = model_eqns + ("\tdouble %s = data.table[%d][i];\n" % (v, vi));

		if model_eqns <> "":
			model_eqns = model_eqns + "\n";
			ms_model = ms_model + "\textern int t_to_index(double t);\n\n";
			ms_model = ms_model + "\tint i = t_to_index(t);\n\n";

		# add the system variables to the OK list
		# add system and oae variables to the list of all variables
		for i in odes.keys():
			all_vars.append(i)
			calc_vars.append(i)
		for i in oaes.keys():
			all_vars.append(i)

		# this may be unnecessary, but for sanity's sake, make sure we don't leave
		# this loop without actually printing all of the OAEs.
		printed_oaes = 0
		while len(calc_vars) < len(all_vars) and printed_oaes < len(oaes.keys()):
			for i in range(len(oaes.keys())):
				(v, [rhs_inst, aggr, vdataid, vinit, vinitrange]) = oaes.items()[i];
				
				# don't print variables that have already been dealt with
				# assumes that variables on the LHS of an algebraic equation are
				# never on the LHS of an ODE (or exogenous).
				if (v in calc_vars): continue
        
				# don't print an equation if all of its dependencies haven't been printed
				print_eq = 1
				for (vname, vtmp1, vtmp2, vtmp3) in variables:
					if rhs_inst.find(vname) >= 0 and not vname in calc_vars: print_eq = -1
				if print_eq < 0: continue

				calc_vars.append(v)
				printed_oaes += 1
				for j in range(len(parameters)):
					rhs_inst = replace_str(parameters[j][0], "param.vals[%d]" % j, rhs_inst);

				model_eqns = model_eqns + ("\tdouble %s = %s;\n" % (v, rhs_inst));
			model_eqns = model_eqns + "\n";
	
		# and finally system variables
		nobserved = 0;
		sys_vars = [];
		obs_sys_vars = [];
		for i in range(len(odes.keys())):
			(v, [rhs_inst, aggr, vdataid, vinit, (vlower, vupper)]) = odes.items()[i];
			sys_vars.append(v);

                        # introduce parameter(s) for fitting the initial value of the system variable
                        is_param_index = len(parameters) + i;

                        # one parameter per training data set
                        for j in range(len(data_sets)):
				# initial parameter value used for fitting ...
				if init_state == None:
					if vinit <> None:
						ms_params = ms_params + ("\tparam.vals[%d] = (double) %g;\n" % (is_param_index, vinit));
				else:
					if v in init_state.keys():
						ms_params = ms_params + ("\tparam.vals[%d] = (double) %g;\n" % (is_param_index, init_state[v][j]));

				# ... and bounds
				is_bound_index = 2 * is_param_index;
				if vlower <> None:
					ms_params = ms_params + ("\tparam.bounds[%d] = (double) %g;\n" % (is_bound_index, vlower));
				if vupper <> None:
					ms_params = ms_params + ("\tparam.bounds[%d] = (double) %g;\n\n" % (is_bound_index  + 1, vupper));

				is_param_index = is_param_index + len(odes.keys());

			vi = ds.var_index(vdataid);

			# data.sys_var_names stores the names of observed system variables
			if vi == -1:
				#print "Warning: Unobserved system variable %s." % v;
				ms_fill_in = ms_fill_in + ("\tdata.sys_vars[%d] = (double *) NULL;\n" % i);
				ms_fill_in = ms_fill_in + ("\tdata.sys_var_names[%d] = strdup(\"%s\");\n" % (i+1, v));
			else:
				obs_sys_vars.append(v);
				ms_fill_in = ms_fill_in + ("\tdata.sys_vars[%d] = data.table[%d];\n" % (i, vi));
				nobserved = nobserved + 1;
				ms_fill_in = ms_fill_in + ("\tdata.sys_var_names[%d] = strdup(\"%s\");\n" % (i+1, vdataid));

			# replace the parameters with their C variable names
			for j in range(len(parameters)):
				rhs_inst = replace_str(parameters[j][0], "param.vals[%d]" % j, rhs_inst)

			ms_model = ms_model + ("\tdouble %s = sys_var.vals[%d];\n" % (v, i));
			model_eqns = model_eqns + ("\tsys_var.dot_vals[%d] = %s;\n" % (i, rhs_inst));

		ms_init_dims = ms_init_dims + ("\tsys_var.nobserved = %d;\n" % nobserved);

		ms_init_dims = ms_init_dims + "}\n";
		ms_params = ms_params + "}\n";
		ms_fill_in = ms_fill_in + "}\n";
		ms_model = ms_model + "\n" + model_eqns + "}\n";

		f = open("ms.c", "w");
		f.write(ms_init_dims + "\n");
		f.write(ms_fill_in + "\n");
		f.write(ms_params + "\n");
		f.write(ms_model);
		f.close();

		return parameters, sys_vars, obs_sys_vars;


	def simulate(self, data_sets, init_state = None, window = 0, output = sys.stdout):

		from distutils.text_file import TextFile;
		from string import split, atof;


		if window == 0:
			self.fit_params(data_sets, init_state = init_state, n_tf_restarts = 0, n_fs_restarts = 0, output = output);
			return;

		first_flag = True;
		for data_set in data_sets:
			for i in range(data_set.length - window):
				ptime = [data_set.time[i + window]];
				if i == 0:
					ptime = ptime + [data_set.time[0]];
					if first_flag:
						ptime = ptime + [-1];
						first_flag = False;

				ds_window = data_set.subset(range(i, i + window + 1), file_suffix = "w");
				ds_window.write_to_file();
				f = open("model.out", "w");
				self.simulate([ds_window], output = f);
				f.close();

				f = TextFile(filename = "model.out");
				simulation_flag = False;
				for line in f.readlines():
					if line[-4:] == "_sim":
						simulation_flag = True;
						if (-1 in ptime): print >> output, line;
						continue;
					if line[:5] == "SSE, ":
						simulation_flag = False;
						continue;
					if simulation_flag:
						fields = line.split(" ");
						if atof(fields[0]) in ptime: print >> output, line;
				f.close();
		# print >> output, "SSE, reMSE, r2, AIC, BIC: not calculated";


	# fit parameters of a model
    #
	#  . data_sets: list of data sets
	#  . normalize_error: error normalization on (1) and off (0)
	#  . init_state_fit: whether to fit the initial state of ODE models (1) or not (0)
	#  . init_state: list of initial values for the initial state
	#                used only for model simulation (n_tf_restarts + n_fs_restarts = 0)
	#  . n_tf_restarts: number of teacher forcing restarts
	#  . n_fs_restarts: number of full simulation restarts

	def fit_params(self, data_sets, normalize_error = 1, init_state_fit = 1, 
		       init_state = None, n_tf_restarts = 128, n_fs_restarts = 8, 
		       output = sys.stdout):
		from os import system, popen;
		from string import split, atof;

		(parameters, sys_vars, obs_sys_vars) = self.to_c_model(data_sets, 
								       normalize_error, 
								       init_state_fit, 
								       init_state, 
								       n_tf_restarts, 
								       n_fs_restarts);

		if parameters == None: return None;
		if system("make > /dev/null") <> 0:
			print "Simulation error for the model:";
			print self;
			return None;

		if n_tf_restarts + n_fs_restarts == 0:
			sse = float("inf")
			nmse = float("inf")
			r2 = 0
			aic = float("inf")
			bic = float("inf")
			for line in popen("./model").readlines(): 
				print >> output, line[:-1];
				fields = line[:-1].split(" ");
				parameter_vals = None;
				if fields[0] == "SSE,": 
					[sse, nmse, r2, aic, bic] = map(atof, [fields[5], fields[6], fields[7],
									       fields[8], fields[9]]);
			return (sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state);

		[sse, nmse, r2, aic, bic] = [None, None, None, None, None];
		for line in popen("./model").readlines():
			fields = line[:-1].split(" ");
			if fields[0] == "PARAMS:": parameter_vals = map(atof, fields[1:]);
			if fields[0] == "SSE,": 
				[sse, nmse, r2, aic, bic] = map(atof, 
								[fields[5], fields[6], fields[7],
								 fields[8], fields[9]]);

		if sse == None: return None;
		
		if init_state_fit == 1:
			init_state = {};
			is_param_index = len(parameters);
			for i in range(len(data_sets)):
				for v in sys_vars:
					if v not in init_state: init_state[v] = [];
					init_state[v].append(parameter_vals[is_param_index])
					is_param_index = is_param_index + 1;

			parameter_vals = parameter_vals[:len(parameters)];
		else:
			init_state = None;

		self.set_parameter_vals(parameters, parameter_vals);
		return (sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state);




class generic_process:

	def __init__(self, lib, id, type, relates, processes, parameters, oaes, odes):

		self.lib = lib;

		self.id = id;
		self.type = type;

		self.relates = relates;
		self.processes = processes;

		self.parameters = parameters;
		self.oaes = oaes;
		self.odes = odes;


	def relates_index_by_id(self, id):

		for i in range(len(self.relates)):
			if self.relates[i][0] == id: return i;
		return -1;


	def relates_aggregation_by_id(self, id_var):
		from string import split;

		[id, var] = id_var.split(".");
		aggregation = [];
		for e in self.relates[self.relates_index_by_id(id)][1]:
			aggregation.append(e.aggregation(var));
		return aggregation;


	# check the compatibility of each relates_instance item with respect
	# to the relates specification of the generic process
	#
	# if incompatibility item is found return None
	# compatible items are returned
	#
	# special match: set vs. scalar <1 to 1> relates specification

	def relates_instances_list(self, relates_instance):

		def sublist(l1, l2):

			for e in l1:
				if e not in l2: return 0;
			return 1;

		ris_list = [];
		for (id, ges, c_from, c_to) in self.relates:

			# relates_instance item is not present - match OK
			if id not in relates_instance: continue;
			
			# get the list of entities for the relates specification
			# this is the set of entities that the process accepts
			geis_list = self.lib.entities_list_instances(ges);
			
			# get the list of entities from relates_instance item
			# this is the set of entities coming into the process
			eis_list = relates_instance[id];

			# we want to pass along those entities that
			# fit the type restrictions of the process, so
			# we need to figure out which entities can be
			# passed along and make sure that the
			# cardinality of that set fits within the
			# specified bounds.

			# get the subset of incoming entities that fit
			# the type restrictions.
			ent_int = filter(lambda i: i in geis_list, eis_list);

			# if too few entities fit the specs, we can't go on
			if len(ent_int) < c_from: return None;

			# length of the relates_instance item OK
			if len(eis_list) <= c_to:
				ris_list.append([{id: ent_int}]);
				continue;

			# in the special case where the upper bound
			# equals one and more than one entity is
			# coming into the process, we will create
			# multiple processes---one for each entering
			# entity.

			# length of the relates_instance is not OK, but the relates specification is scalar
			if c_to == 1:
				if c_from == 1: ris_list.append([{id: [i]} for i in ent_int]);
				else: ris_list.append([{id: [i]} for i in ent_int] + [{id: []}]);
				continue;

			# if none of the above - mismatch
			return None;
		return ris_list;

			# OLD STUFF

			# relates_instance item contains incompatible entities - mismatch
#			if not sublist(eis_list, geis_list): return None;

#			if mixing == 1: print "good match"

			# too few instances in the relates_instance item
#			if len(eis_list) < c_from: return None;

#			if mixing == 1: print "lower bound met"

			# length of the relates_instance item OK
#			if len(eis_list) <= c_to:
#				ris_list.append([{id: relates_instance[id]}]);
#				continue;

#			if mixing == 1: print "upper bound met"

			# this is the special case where we unroll processes
			# length of the relates_instance is not OK, but the relates specification is scalar
#			if c_to == 1:
#				if c_from == 1: ris_list.append([{id: [i]} for i in relates_instance[id]]);
#				else: ris_list.append([{id: [i]} for i in relates_instance[id]] + [{id: []}]);
#				continue;

			# if none of the above - mismatch
#			return None;

#		return ris_list;


	# get a list of all possible instances of a process
	#  . relates_instance - process relates entity instances

	def instances(self, relates_instance):

		def improve(c):

			i = len(c) - 1;
			if i == -1: return 0;

			while c[i] == len(relates_possibilities[self.relates[i][0]]) - 1:
				c[i] = 0;
				i = i - 1;
				if i == -1: return 0;
			c[i] = c[i] + 1;
			return 1;

		def c2gpis(c):

			def subprocess_instances(type, relates_ids, optional):
	
				subp_ri = {};
				for i in range(len(relates_ids)):
					id = relates_ids[i];
					if id not in relates_instance: continue;
					for gp in self.lib.generic_processes_by_type(type):
						if i >= len(gp.relates): continue;
						subp_id = gp.relates[i][0];
						subp_ri[subp_id] = relates_instance[id];

				subpis = self.lib.process_type_instances(type, subp_ri);
				if optional == 1: subpis = [[None]] + subpis;
				return subpis;


			def improve(c, list):

				i = len(c) - 1;
				if i == -1: return 0;

				while c[i] == len(list[i]) - 1:
					c[i] = 0;
					i = i - 1;
					if i == -1: return 0;
				c[i] = c[i] + 1;
				return 1;

			def c2pis(c):
				subp_instances = [];
				for i in range(len(pis_lists)):
					pis_list = pis_lists[i][c[i]];
					if pis_list == [None]: continue;
					subp_instances.extend(pis_list);
				return process_instance(self, relates_instance, subp_instances, None);


			# compile the relate_instance and build process instance
			relates_instance = {};
			for i in range(len(self.relates)):
				id = self.relates[i][0];

				# assumption (!might be wrong for some domains, make it explicit!):
				#  if a process relates two variables of the same type
				#  they should be instantiated differently
				#  used in a chemistry domain to avoid reactions X -> X

				for j in range(i):
					if self.relates[i][1] == self.relates[j][1]:
						idj = self.relates[j][0];
						if (relates_possibilities[id][c[i]] == relates_possibilities[idj][c[j]]):
							return [];

				relates_instance[id] = relates_possibilities[id][c[i]];
			gpi = process_instance(self, relates_instance, [], None);

			# if there are no subprocesses that is all we should do
			if len(self.processes) == 0: return [gpi];

			# otherwise
			# generate all possible subprocess instances

			pis_lists = [];
			for (subp_type, relates_ids, optional) in self.processes:
				pis_list = subprocess_instances(subp_type, relates_ids, optional);
				if len(pis_list) == 0: continue;
				pis_lists.append(pis_list);

			c = [0 for pis_list in pis_lists];
			gpis = [c2pis(c)];
			while improve(c, pis_lists): gpis.append(c2pis(c));

			return gpis;


		relates_possibilities = {};
		for (id, ges, c_from, c_to) in self.relates:
			if id in relates_instance:
				relates_possibilities[id] = [relates_instance[id]];
			else:
				relates_possibilities[id] = self.lib.entities_list_set_instances(ges, c_from, c_to);

			if len(relates_possibilities[id]) == 0: return [];

		c = [0 for e in self.relates];
		gpis = c2gpis(c);
		while improve(c): gpis.extend(c2gpis(c));
		return gpis;


	def translate_relates_instance(self, relates_ids, relates_instance):

		sp_relates_instance = {};
		for i in range(len(self.relates)):
			id = self.relates[i];
			if id not in relates_ids: continue;


	# return a list of simplest instances of a process
	# given the relates_instance

	def simplest_instances(self, relates_instance):

		def improve(c, list):

			i = len(c) - 1;
			if i == -1: return 0;

			while c[i] == len(list[i]) - 1:
				c[i] = 0;
				i = i - 1;
				if i == -1: return 0;
			c[i] = c[i] + 1;
			return 1;

		def c2pis(c):
			spis = [];
			for i in range(len(c)):
				spi_list = spis_candidates[i][c[i]];
				if spi_list == []: continue;
				spis.extend(spi_list);
			return process_instance(self, relates_instance, spis, None);


		lib = self.lib;

		# assumption (!might be wrong for some domains, make it explicit!):
		#  if a process relates two variables of the same type
		#  they should be instantiated differently
		#  used in a chemistry domain to avoid reactions X -> X

		for i in range(len(self.relates)):
			(i_id, i_ges, i_c_from, i_c_to) = self.relates[i];
			for j in range(i):
				(j_id, j_ges, j_cfrom, j_cto) = self.relates[j];
				if i_ges == j_ges:
					if relates_instance[i_id] == relates_instance[j_id]:
						return [];

		# generate a process instance with no subprocesses
		if len(self.processes) == 0: return [process_instance(self, relates_instance, [], None)];

		# for each subprocess collect the list of (simplest) candidate process instances
		spis_candidates = [];
		for (sp_type, sp_relates_ids, sp_optional) in self.processes:
			# if the subprocess is optional the simplest candidate is "process is not present"
			if sp_optional == 1:
				spis_candidates.append([[]]);
				continue;

			# otherwise instantiate the process type
			sp_type_relates_instance = lib.subprocess_type_relates_instance(sp_type, sp_relates_ids, relates_instance);
			spi_candidates = lib.simplest_process_type_instances(sp_type, sp_type_relates_instance);
			spis_candidates.append(spi_candidates);

		pis = [];
		c = [0 for spi_candidates in spis_candidates];
		pis.append(c2pis(c));
		while improve(c, spis_candidates): pis.append(c2pis(c));

		return pis;
