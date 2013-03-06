#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from population_dynamics_lib import *;

# <entity_id> = entity_instance(<generic_entity_id>, "entity name",
#       { # variables
#           # regular model variable
#           "variable name": ("system",
#                             "data column name" or <initial_value>,
#                             (<minimum_value>, <maximum_value>))
#           # exogenous variable
#           "variable name": ("exogenous",
#                             "data column name",
#                             None)
#       },
#       { # constants
#           "constant name": <constant_value>
#       };

# observed grazer
g1 = entity_instance(ge, 
                     "nasutum", 
		     {"conc": ("system", "nasutum", None)}, 
                     None);

# observed primary producer
p1 = entity_instance(pe, 
                     "aurelia", 
                     {"conc": ("system", "aurelia", None)}, 
                     None);
