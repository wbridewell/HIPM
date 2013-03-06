#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.


# J. D. Murray's Mathematical Biology.  Inhibition and activation
# processes are commented out to reduce the size of the search space.

from library import *;
from entities import *;
from processes import *;


lib = library("chemical_kinetics");

ce = lib.add_generic_entity("C",
	{"conc": "sum", "pos_flux": "prod", "neg_flux": "prod"},
	{"pos_rate": (0, 10), "neg_rate": (0, 10)}
);

ee = lib.add_generic_entity("E",
	{"conc": "sum"},
	{}
);


lib.add_generic_process("root", "",
	[("C",[ce],2,100), ("E",[ee],0,100)],
	[("reaction", ["C","C","E"], 1)],
	{}, {},
	{
		"C.conc": "C.pos_rate * C.pos_flux - C.neg_rate * C.neg_flux"
	}
);


# C1 -> C2

lib.add_generic_process("irreversible", "reaction",
	[("C1",[ce],1,1), ("C2",[ce],1,1), ("E",[ee],0,1)],
	[("influence", ["C1", "C2", "E"], 0)],
	{"kinetic_order_1": (0,1), "kinetic_order_2": (0,1)},
	{
		"C1.neg_flux": "pow(C1.conc, kinetic_order_1)",
		"C2.pos_flux": "pow(C1.conc, kinetic_order_2)"
	},
	{}
);

# C1 -> C2 inhibited by E

#lib.add_generic_process("inhibition", "influence",
#	[("C1",[ce],1,1), ("C2",[ce],1,1), ("E",[ee],1,1)],
#	[],
#	{
#		"kinetic_order_1": (0,1), "kinetic_order_2": (0,1),
#	},
#	{
#		"C1.neg_flux": "pow(E.conc, -1 * kinetic_order_1)",
#		"C2.pos_flux": "pow(E.conc, -1 * kinetic_order_2)"
#	},
#	{}
#);

# C1 -> C2 activated by E

#lib.add_generic_process("activation", "influence",
#	[("C1",[ce],1,1), ("C2",[ce],1,1), ("E",[ee],1,1)],
#	[],
#	{
#		"kinetic_order_1": (0,1), "kinetic_order_2": (0,1),
#	},
#	{
#		"C1.neg_flux": "pow(E.conc, kinetic_order_1)",
#		"C2.pos_flux": "pow(E.conc, kinetic_order_2)"
#	},
#	{}
#);


# C1 <-> C2

lib.add_generic_process("reversible", "reaction",
	[("C1",[ce],1,1), ("C2",[ce],1,1)],
	[],
	{
		"kinetic_order_11p": (0,1), "kinetic_order_12p": (0,1),
		"kinetic_order_11n": (0,1), "kinetic_order_12n": (0,1),
		"kinetic_order_21p": (0,1), "kinetic_order_22p": (0,1),
		"kinetic_order_21n": (0,1), "kinetic_order_22n": (0,1)
	},
	{
		"C1.pos_flux": "pow(C1.conc, kinetic_order_11p) * pow(C2.conc, kinetic_order_12p)",
		"C1.neg_flux": "pow(C1.conc, kinetic_order_11n) * pow(C2.conc, kinetic_order_12n)",
		"C2.pos_flux": "pow(C1.conc, kinetic_order_21p) * pow(C2.conc, kinetic_order_22p)",
		"C2.neg_flux": "pow(C1.conc, kinetic_order_21n) * pow(C2.conc, kinetic_order_22n)"
	},
	{}
);
