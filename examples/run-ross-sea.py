#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from ross_entities import *;
from data import data;
from search import *;
import time;

ross1 = data("./ross-sea-yr1.data");
#ross2 = data("./ross-sea-yr2.data");

data_files = [ross1]
#data_files = [ross1,ross2]

print "Start: ", time.strftime("%m/%d/%y %H:%M:%S", time.localtime());

(model, sse, r2, init_state) = exhaustive(lib, "root", data_files);
#(model, sse, r2, init_state) = beam_search(lib, "root", data_files, beam_width = 1, n_fs_restarts = 3,init_state_fit=1);

print "Best model:";
print model;
print "SSE = %g, r2 = %g, init_state = %s\n" % (sse, r2, str(init_state));

model.fit_params(data_files, init_state = init_state, n_tf_restarts = 0, n_fs_restarts = 0);

print "End: ", time.strftime("%m/%d/%y %H:%M:%S", time.localtime());
print "";
