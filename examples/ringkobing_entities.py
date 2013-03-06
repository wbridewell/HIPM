#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from fjord_lib import *;

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


h = entity_instance(le, "water_level", {"fjord": ("system", "h", None), 
                                        "sea": ("exogenous", "h_sea", None)}, None);
A = entity_instance(se, "fjord_surface", {"area": ("exogenous", "A", None)}, None);
a = entity_instance(ge, "gate2", {"opened": ("exogenous", "a", None)}, None);
Q_f = entity_instance(fe, "fresh_water_inflow", {"value": ("exogenous", "Q_f", None)}, None);
W = entity_instance(we, "wind", {"velocity": ("exogenous", "W_Vel", None), 
                                 "direction": ("exogenous", "W_Dir", None)}, None);

