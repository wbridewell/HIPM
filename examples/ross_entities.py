#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from aquatic_ecosystem_lib import *;  # import library

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


# observed primary producer
p1 = entity_instance(pe, "phyto",
	{
		"conc": ("system", "P", (1,11)), # ugC/L
		"growth_rate": ("system", 0, (0,1)),
		"growth_lim": ("system", 1, (0,1))
	},
	{
		"max_growth": 0.59,
		"exude_rate": 0.19,
		"death_rate": 0.025,
		"Ek_max": 30,
		"biomin": 0.025,
		"PhotoInhib": 200
	}
);

# unobserved grazer
z1 = entity_instance(ze, "zoo",
	{
		"conc": ("system", 1, (1.5,2.5)),
		"growth_rate": ("system", 0.1, (0, 1))
	},
	{
		"assim_eff": 0.75,
		"death_rate": 0.02,
		"respiration_rate": 0.019,
		"gmax": 0.4,
		"gcap":200
	}
);

# observed nitrate
no31 = entity_instance(no3, "NO3",
	{
		"conc": ("system", "NO3", (31,32)),
		"mixing_rate": ("system", 0, (0,1))
	},
	None
);

# unobserved iron
fe1 = entity_instance(fe, "Fe", 
	{
		"conc": ("system", 0.000429206, (0,0.001)),
		"mixing_rate": ("system", 0, (0,1))
	},
	None
);

# observed/exogenous ENVIRONMENT
e1 = entity_instance(ee, "environment",
	{
		"PUR": ("exogenous", "PUR", None),
		"TH2O": ("exogenous", "SST", None),
		"ice": ("exogenous", "ice", None)
	},
	{
		"beta":0.7
	}
);

# unobserved detritus with initial value from [0.001,1] default 0.1
d1 = entity_instance(de, "detritus", 
	{
		"conc": ("system", 0.1, (0.001, 0.1))
	},
	None
);
