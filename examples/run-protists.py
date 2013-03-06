#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from protist_entities import *;
from data import data;
from search import *;

pp = [data("1cs2.txt")];

# model multiple data sets
#pp = [data("1as3.txt"),data("1cs2.txt")];


# this beam width is large enough that it is equivalent to exhaustive search for this modeling task.
# see search.py for a description of the parameters to beam_search() and exhaustive().
(pi, sse, nmse, r2, aic, bic, init_state) = beam_search(lib, "root", pp, beam_width=32, eval_f="nmse");

print "Best model:";
print pi;
print "SSE = %g, nMSE = %g, r2 = %g, AIC = %g, BIC = %g, init_state = %s\n" % (sse, nmse, r2, aic, bic, str(init_state));

## simulate the best model.
## can also call pi.simulate(pp, init_state=init_state)
pi.fit_params(pp, init_state = init_state, n_tf_restarts = 0, n_fs_restarts = 0)
