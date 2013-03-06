#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.


from library import library;

def count_models(lib, root_id, maximum_processes = 0):
        from sys import stdout;

        pis_list = lib.process_type_instances(root_id, {});
        n = 0;
        for pis in pis_list:
                for pi in pis:
                        if maximum_processes > 0 and pi.number_of_process_instances() > maximum_processes: continue;
                        n = n + 1;
        print "Total Model Structures: %d" % n;
        return n;

# perform an exhaustive search by visiting each model structure in an arbitrary order
# lib = generic process library
# root_id = root process in the process hierarchy
# data = list of data files for evaluation purposes
# eval_f = the fitness measure for the objective function
# verbosity_level = 1: all models are printed to standard out, else no models are printed
# normalize_error = 1: normalize the error measure used for evaluating model fitness
# init_state_fit = 1: treat the initial conditions of unobserved variables as parameters
# n_tf_restarts = the number of teacher forcing optimization runs
# n_fs_restarts = the number of full simulation optimization runs
# maximum_processes = the upper limit on the size of the model (0 indicates no limit on size)
def exhaustive(lib, root_id, data, eval_f = "sse", verbosity_level = 1, normalize_error = 1, 
	       init_state_fit = 1, n_tf_restarts = 128, n_fs_restarts = 8, maximum_processes = 0):
	from sys import stdout;

	n_models = count_models(lib, root_id, maximum_processes);
	n = 0;

	(opt_sse, opt_nmse, opt_r2, opt_aic, opt_bic) = (None, None, None, None, None);
	pis_list = lib.process_type_instances(root_id, {});
	for pis in pis_list:
		for pi in pis:
			if maximum_processes > 0 and pi.number_of_process_instances() > maximum_processes: continue;
			res = pi.fit_params(data, normalize_error = normalize_error, init_state_fit = init_state_fit, 
					    n_tf_restarts = n_tf_restarts, n_fs_restarts = n_fs_restarts);
			if res == None: continue;
			(sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state) = res;
			n = n + 1;
			if verbosity_level == 1:
				print pi;
				print "";
				print "#Parameters = %d" % len(parameters);
				print "#Processes = %d" % pi.number_of_process_instances();
				print "SSE = %g" % sse;
				print "reMSE = %g" % nmse;
				print "r2 = %g" % r2;
				print "AIC = %g" % aic;
				print "BIC = %g" % bic;
				if init_state <> None: print "init_state = %s" % str(init_state);
				print "";
				print "";
				stdout.flush();

			opt_flag = 0;
			if eval_f == "sse":
				if (opt_sse == None) or (sse < opt_sse):
					opt_flag = 1;
			elif eval_f == "nmse":
				if (opt_nmse == None) or (nmse < opt_nmse):
					opt_flag = 1;
			else:
				if (opt_r2 == None) or (r2 > opt_r2):
					opt_flag = 1;

			if (opt_flag == 1):
				(opt_sse, opt_nmse, opt_r2, opt_aic, opt_bic) = (sse, nmse, r2, aic, bic);
				opt_pi = pi;
				(opt_parameters, opt_parameter_vals, opt_init_state) = (parameters, parameter_vals, init_state);

	opt_pi.set_parameter_vals(opt_parameters, opt_parameter_vals);
	return (opt_pi, opt_sse, opt_nmse, opt_r2, opt_aic, opt_bic, opt_init_state);


# perform a beam search by progressively adding processes to the model structure
# lib = generic process library
# root_id = root process in the process hierarchy
# data = list of data files for evaluation purposes
# beam_width = the number of models to keep on the beam after each search iteration
# eval_f = the fitness measure for the objective function ("r2" or "sse")
# verbosity_level = 1: all models are printed to standard out, else no models are printed
# normalize_error = 1: normalize the error measure used for evaluating model fitness
# init_state_fit = 1: treat the initial conditions of unobserved variables as parameters
# n_tf_restarts = the number of teacher forcing optimization runs
# n_fs_restarts = the number of full simulation optimization runs
# maximum_processes = the upper limit on the size of the model (0 indicates no limit on size)
def beam_search(lib, root_id, data, beam_width = 32, eval_f = "sse", verbosity_level = 1, normalize_error = 1, 
		init_state_fit = 1, n_tf_restarts = 128, n_fs_restarts = 8, maximum_processes = 0):
	from sys import stdout;

	def in_beam(beam):

		for i in range(len(beam)):
			(c_error, c_pi, rest) = beam[i];
			if pi.same_structure(c_pi): return(1);
		return(0);

	def add_to_beam(beam, error):

		candidate = (error, pi, (sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state));
		if len(beam) == 0:
			beam = [candidate];
			return(beam);

		if error > beam[-1][0]:
			beam = beam + [candidate];
		else:
			for i in range(len(beam)):
				(c_error, c_pi, rest) = beam[i];
				if error < c_error: break;
			beam = beam[:i] + [candidate] + beam[i:];
		return(beam[:beam_width]);

	n = 0;

	beam = [];
	root_gp = lib.generic_process_by_id(root_id);
	ri = {};
	for (id, ges, c_from, c_to) in root_gp.relates: ri[id] = lib.entities_list_instances(ges);

	# seed the beam with all non-redundant minimal models
	for pi in root_gp.simplest_instances(ri):
		if in_beam(beam): continue;

		res = pi.fit_params(data, normalize_error = normalize_error, init_state_fit = init_state_fit, 
				    n_tf_restarts = n_tf_restarts, n_fs_restarts = n_fs_restarts);
		n = n + 1;

		if res == None: continue;
		(sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state) = res;
		if eval_f == "sse": beam = add_to_beam(beam, sse);
		elif eval_f == "nmse": beam = add_to_beam(beam, nmse);
		else: beam = add_to_beam(beam, -r2);

		if verbosity_level == 2:
			pi.set_parameter_vals(parameters, parameter_vals);
			print pi;
			print "";
			print "#Parameters = %d" % len(parameters);
			print "#Processes = %d" % pi.number_of_process_instances();
			print "SSE = %g" % sse;
			print "reMSE = %g" % nmse;
			print "r2 = %g" % r2;
			print "AIC = %g" % aic;
			print "BIC = %g" % bic;
			if init_state <> None: print "init_state = %s" % str(init_state);
			print "";
			print "";
			stdout.flush();

	while 1:
		if verbosity_level == 1:
			print "Current beam (%d):" % len(beam);
			for candidate in beam:
				(error, pi, (sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state)) = candidate;
				pi.set_parameter_vals(parameters, parameter_vals);
				print pi;
				print "";
				print "#Parameters = %d" % len(parameters);
				print "#Processes = %d" % pi.number_of_process_instances();
				print "SSE = %g" % sse;
				print "reMSE = %g" % nmse;
				print "r2 = %g" % r2;
				print "AIC = %g" % aic;
				print "BIC = %g" % bic;
				if init_state <> None: print "init_state = %s" % str(init_state);
				print "";
				print "";
				stdout.flush();

		p_beam = [candidate for candidate in beam];
		for candidate in p_beam:
			(error, c_pi, (sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state)) = candidate;
			# Force the number of processes to be less than or equal to the maximum number 
			# of allowable processes---if specified
			if maximum_processes > 0 and c_pi.number_of_process_instances() >= maximum_processes: continue;
			c_pi.set_parameter_vals(parameters, parameter_vals);
			for pi in c_pi.add_process():
				if in_beam(beam): continue;

				res = pi.fit_params(data, normalize_error = normalize_error, init_state_fit = init_state_fit, 
						    n_tf_restarts = n_tf_restarts, n_fs_restarts = n_fs_restarts);
				n = n + 1;

				if res == None: continue;
				(sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state) = res;
				if eval_f == "sse": beam = add_to_beam(beam, sse);
				elif eval_f == "nmse": beam = add_to_beam(beam, nmse);
				else: beam = add_to_beam(beam, -r2);

				if verbosity_level == 2:
					pi.set_parameter_vals(parameters, parameter_vals);
					print pi;
					print "";
					print "#Parameters = %d" % len(parameters);
					print "#Processes = %d" % pi.number_of_process_instances();
					print "SSE = %g" % sse;
					print "reMSE = %g" % nmse;
					print "r2 = %g" % r2;
					print "AIC = %g" % aic;
					print "BIC = %g" % bic;
					if init_state <> None: print "init_state = %s" % str(init_state);
					print "";
					print "";
					stdout.flush();

		if p_beam == beam: break;
	if verbosity_level == 1: print "%d candidate models evaluated\n\n\n" % n;

	(opt_error, opt_pi, (opt_sse, opt_nmse, opt_r2, opt_aic, opt_bic, opt_parameters, opt_parameter_vals, opt_init_state)) = beam[0];
	opt_pi.set_parameter_vals(opt_parameters, opt_parameter_vals);
	return (opt_pi, opt_sse, opt_nmse, opt_r2, opt_aic, opt_bic, opt_init_state);

# perform an exhaustive search by visiting each model structure in an arbitrary order
# lib = generic process library
# base_model = the model from which revision will proceed
# data = list of data files for evaluation purposes
# beam_width = the number of models to keep on the beam after each search iteration
# eval_f = the fitness measure for the objective function ("r2" or "sse")
# verbosity_level = 1: all models are printed to standard out, else no models are printed
# normalize_error = 1: normalize the error measure used for evaluating model fitness
# init_state_fit = 1: treat the initial conditions of unobserved variables as parameters
# n_tf_restarts = the number of teacher forcing optimization runs
# n_fs_restarts = the number of full simulation optimization runs
# maximum_processes = the upper limit on the size of the model (0 indicates no limit on size)
def beam_revision(lib, base_model, data, beam_width = 32, eval_f = "sse", verbosity_level = 0, normalize_error = 1, 
		  init_state_fit = 1, n_tf_restarts = 128, n_fs_restarts = 8, maximum_processes = 0):
	from sys import stdout;

	def in_beam(beam):

		for i in range(len(beam)):
			(c_error, c_pi, rest) = beam[i];
			if pi.same_structure(c_pi): return(1);
		return(0);

	def add_to_beam(beam, error):

		candidate = (error, pi, (sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state));
		if len(beam) == 0:
			beam = [candidate];
			return(beam);

		if error > beam[-1][0]:
			beam = beam + [candidate];
		else:
			for i in range(len(beam)):
				(c_error, c_pi, rest) = beam[i];
				if error < c_error: break;
			beam = beam[:i] + [candidate] + beam[i:];
		return(beam[:beam_width]);

	n = 0;
	beam = [];
	# seed the beam with the original model
	pi = base_model
    	res = pi.fit_params(data, normalize_error = normalize_error, init_state_fit = init_state_fit, 
			    n_tf_restarts = n_tf_restarts, n_fs_restarts = n_fs_restarts);
	n = n + 1;
	if res == None:
		sse = pow(2,31)
		nmse = sse
		aic = 0
		bic = 0
		r2 = 0
		init_state = None
		parameters = pi.to_c_model(data, normalize_error, init_state_fit, init_state, n_tf_restarts, n_fs_restarts)
		parameter_vals = pi.parameter_vals
	else:
		(sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state) = res;
	
	if eval_f == "sse": beam = add_to_beam(beam, sse);
	elif eval_f == "nmse": beam = add_to_beam(beam, nmse);
	else: beam = add_to_beam(beam, -r2);

	while 1:
		if verbosity_level == 1:
			print "Current beam (%d):" % len(beam);
			for candidate in beam:
				(error, pi, (sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state)) = candidate;
				pi.set_parameter_vals(parameters, parameter_vals);
				print pi;
				print "";
				print "#Parameters = %d" % len(parameters);
				print "#Processes = %d" % pi.number_of_process_instances();
				print "SSE = %g" % sse;
				print "reMSE = %g" % nmse;
				print "r2 = %g" % r2;
				print "AIC = %g" % aic;
				print "BIC = %g" % bic;
				if init_state <> None: print "init_state = %s" % str(init_state);
				print "";
				print "";
				stdout.flush();
				
		p_beam = [candidate for candidate in beam];
		for candidate in p_beam:
			(error, c_pi, (sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state)) = candidate;
			# Force the number of processes to be less than or equal to the maximum number 
			# of allowable processes---if specified
			if maximum_processes > 0 and c_pi.number_of_process_instances() >= maximum_processes: continue;
			c_pi.set_parameter_vals(parameters, parameter_vals);
			for pi in c_pi.add_process():
				if in_beam(beam): continue;

				res = pi.fit_params(data, normalize_error = normalize_error, init_state_fit = init_state_fit, 
						    n_tf_restarts = n_tf_restarts, n_fs_restarts = n_fs_restarts);
				n = n + 1;

				#if res == None: continue;
				#(sse, r2, parameters, parameter_vals, init_state) = res;

				if res == None:
					sse = None
					nmse = None
					r2 = None
					aic = None
					bic = None
					parameters = None
					parameter_vals = None
					init_state = None
				else:
					(sse, nmse, r2, aic, bic, parameters, parameter_vals, init_state) = res;
				
				if eval_f == "sse": beam = add_to_beam(beam, sse);
				elif eval_f == "nmse": beam = add_to_beam(beam, nmse);
				else: beam = add_to_beam(beam, -r2);
		if p_beam == beam: break;
	if verbosity_level == 1: print "%d candidate models evaluated\n\n\n" % n;

	return beam

#	(opt_error, opt_pi, (opt_sse, opt_r2, opt_parameters, opt_parameter_vals, opt_init_state)) = beam[0];
#	opt_pi.set_parameter_vals(opt_parameters, opt_parameter_vals);
#	return (opt_pi, opt_sse, opt_r2, opt_init_state);
