#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from ringkobing_entities import *;
from data import data;
from search import *;

fjord = data("./fjord.data");
 
(pi, sse, nmse, r2, aic, bic, init_state) = exhaustive(lib, "root", [fjord], n_tf_restarts = 2, n_fs_restarts = 1);

print "Best Model:";
print pi;
print "SSE = %g, nMSE = %g, r2 = %g, init_state = %s\n" % (sse, nmse, r2, str(init_state));

pi.fit_params([fjord], init_state = init_state, n_tf_restarts = 0, n_fs_restarts = 0);
