#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from ross_entities import *;
from data import data;
from search import *;


ross_first_year = data("./data/rsp_yr1_sat.data");
train_data = [ross_first_year];

no3eis = no3e.instances;
eeis = ee.instances;
peis = pe.instances;

for bw in [4, 8, 16, 32]:
	for fsr in [1, 8]:

		for zeis in [[], [z1]]:
			ze.instances = zeis;

			# for no3eis in [[], [no3]]:
				# no3e.instances = no3eis;

			for feeis in [[], [fe]]:
				fee.instances = feeis;

				for deis in [[], [d1]]:
					de.instances = deis;

					print "---------------------------"
					print "beam width:", bw;
					print "full sim restarts:", fsr;

					print "Z:",;
					for zei in zeis: print zei,;
					print "";

					print "P:",;
					for pei in peis: print pei,;
					print "";

					print "N:",;
					for nei in no3eis + feeis: print nei,;
					print "";

					print "D:",;
					for dei in deis: print dei,;
					print "";

					print "E:",;
					for eei in eeis: print eei,;
					print "";
					print "---------------------------"
					print "";

					(model, sse, nmse, r2, aic, bic, init_state) = beam_search(lib, "root", train_data, beam_width = bw, eval_f = "nmse", normalize_error = 1, verbosity_level = 2, n_fs_restarts = fsr);

					# the , on the next line is important, since in that way we "hide"
				    # the model from the script that transforms models to Aleph examples
					# note: the best model is already reported in the log file, so the
					# example based on it would be redundant/duplicate
					print "Best model:",;
					print model;
					print "SSE = %g, nMSE = %g, r2 = %g, init_state = %s\n" % (sse, nmse, r2, str(init_state));

					print "Best model simulation:";
					model.fit_params(train_data, init_state = init_state, n_tf_restarts = 0, n_fs_restarts = 0);
					print "";
					print "";
