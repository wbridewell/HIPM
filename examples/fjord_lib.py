#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from library import *;
from entities import *;
from processes import *;



lib = library("fjord");


# generic entities: id, variables, constant parameters

le = lib.add_generic_entity("level",
	{"fjord": "sum", "sea": "sum"},
	{}
);

ge = lib.add_generic_entity("gate",
	{"opened": "sum", "influence": "sum"},
	{}
);

se = lib.add_generic_entity("surface",
	{"area": "sum"},
	{}
);

fe = lib.add_generic_entity("flow",
	{"value": "sum"},
	{}
);

we = lib.add_generic_entity("wind",
	{"velocity": "sum", "direction": "sum", "influence": "sum"},
	{}
);



# generic processes: id, type, entities related, list of subprocesses, constant parameters, equations

lib.add_generic_process("root", "",
	[("WL",[le],1,1), ("G",[ge],1,1), ("A",[se],1,1), ("Q_f",[fe],1,1), ("W",[we],1,1)],
	[
		("gate_influence_0", ["G"], 0),
		("gate_influence_1", ["G"], 1),
		("gate_influence_2", ["G"], 1),
		("wind_forcing_0", ["W"], 0),
		("wind_forcing_1v", ["W"], 1),
		("wind_forcing_1d", ["W"], 1),
		("wind_forcing_2v", ["W"], 1),
		("wind_forcing_2d", ["W"], 1),
		("wind_forcing_2vd", ["W"], 1)
	],
	{"wl_0": (-5,5)},
	{},
	{
		"WL.fjord": "G.influence * (WL.sea - WL.fjord + wl_0) / A.area + Q_f.value / A.area + W.influence"
	}
);

lib.add_generic_process("gate_influence_0", "",
	[("G",[ge],1,1)],
	[],
	{"const_gi_0": (-100000, 100000)},
	{
		"G.influence": "const_gi_0"
	},
	{}
);

lib.add_generic_process("gate_influence_1", "",
	[("G",[ge],1,1)],
	[],
	{"const_gi_1": (-100000, 100000)},
	{
		"G.influence": "const_gi_1 * G.opened"
	},
	{}
);

lib.add_generic_process("gate_influence_2", "",
	[("G",[ge],1,1)],
	[],
	{"const_gi_2": (-100000, 100000)},
	{
		"G.influence": "const_gi_2 * G.opened * G.opened"
	},
	{}
);

lib.add_generic_process("gate_influence_3", "",
	[("G",[ge],1,1)],
	[],
	{"const_gi_3": (-100000, 100000)},
	{
		"G.influence": "const_gi_3 * G.opened * G.opened * G.opened"
	},
	{}
);

lib.add_generic_process("gate_influence_4", "",
	[("G",[ge],1,1)],
	[],
	{"const_gi_4": (-100000, 100000)},
	{
		"G.influence": "const_gi_4 * G.opened * G.opened * G.opened * G.opened"
	},
	{}
);

lib.add_generic_process("gate_influence_5", "",
	[("G",[ge],1,1)],
	[],
	{"const_gi_5": (-100000, 100000)},
	{
		"G.influence": "const_gi_5 * G.opened * G.opened * G.opened * G.opened * G.opened"
	},
	{}
);

lib.add_generic_process("wind_forcing_0", "wind_forcing_0",
	[("W",[we],1,1)],
	[],
	{"const_wf_0": (-100000, 100000)},
	{
		"W.influence": "const_wf_0"
	},
	{}
);

lib.add_generic_process("wind_forcing_1v", "wind_forcing_1v",
	[("W",[we],1,1)],
	[],
	{"const_wf_1v": (-100000, 100000)},
	{
		"W.influence": "const_wf_1v * W.velocity"
	},
	{}
);

lib.add_generic_process("wind_forcing_2v", "wind_forcing_2v",
	[("W",[we],1,1)],
	[],
	{"const_wf_2v": (-100000, 100000)},
	{
		"W.influence": "const_wf_2v * W.velocity * W.velocity"
	},
	{}
);

lib.add_generic_process("wind_forcing_1d_plain", "wind_forcing_1d",
	[("W",[we],1,1)],
	[],
	{"const_wf_1d": (-100000, 100000)},
	{
		"W.influence": "const_wf_1d * W.direction"
	},
	{}
);

lib.add_generic_process("wind_forcing_1d_sin", "wind_forcing_1d",
	[("W",[we],1,1)],
	[],
	{"const_wf_1d": (-100000, 100000)},
	{
		"W.influence": "const_wf_1d * sin(W.direction * 3.14159 / 180)"
	},
	{}
);

lib.add_generic_process("wind_forcing_1d_cos", "wind_forcing_1d",
	[("W",[we],1,1)],
	[],
	{"const_wf_1d": (-100000, 100000)},
	{
		"W.influence": "const_wf_1d * cos(W.direction * 3.14159 / 180)"
	},
	{}
);

lib.add_generic_process("wind_forcing_2d_plain", "wind_forcing_2d",
	[("W",[we],1,1)],
	[],
	{"const_wf_2d": (-100000, 100000)},
	{
		"W.influence": "const_wf_2d * W.direction * W.direction"
	},
	{}
);

lib.add_generic_process("wind_forcing_2d_sin", "wind_forcing_2d",
	[("W",[we],1,1)],
	[],
	{"const_wf_2d": (-100000, 100000)},
	{
		"W.influence": "const_wf_2d * sin(W.direction * 3.14159 / 180) * sin(W.direction * 3.14159 / 180)"
	},
	{}
);

lib.add_generic_process("wind_forcing_2d_cos", "wind_forcing_2d",
	[("W",[we],1,1)],
	[],
	{"const_wf_2d": (-100000, 100000)},
	{
		"W.influence": "const_wf_2d * cos(W.direction * 3.14159 / 180) * cos(W.direction * 3.14159 / 180)"
	},
	{}
);

lib.add_generic_process("wind_forcing_2d_sin_cos", "wind_forcing_2d",
	[("W",[we],1,1)],
	[],
	{"const_wf_2d": (-100000, 100000)},
	{
		"W.influence": "const_wf_2d * sin(W.direction * 3.14159 / 180) * cos(W.direction * 3.14159 / 180)"
	},
	{}
);

lib.add_generic_process("wind_forcing_2vd_plain", "wind_forcing_2vd",
	[("W",[we],1,1)],
	[],
	{"const_wf_2vd": (-100000, 100000)},
	{
		"W.influence": "const_wf_2vd * W.velocity * W.direction"
	},
	{}
);

lib.add_generic_process("wind_forcing_2vd_sin", "wind_forcing_2vd",
	[("W",[we],1,1)],
	[],
	{"const_wf_2vd": (-100000, 100000)},
	{
		"W.influence": "const_wf_2vd * W.velocity * sin(W.direction * 3.14159 / 180)"
	},
	{}
);

lib.add_generic_process("wind_forcing_2vd_cos", "wind_forcing_2vd",
	[("W",[we],1,1)],
	[],
	{"const_wf_2vd": (-100000, 100000)},
	{
		"W.influence": "const_wf_2vd * W.velocity * cos(W.direction * 3.14159 / 180)"
	},
	{}
);

