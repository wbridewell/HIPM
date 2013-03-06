#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

from library import *;
from entities import *;
from processes import *;



lib = library("pp_simple");


# generic entities: id, variables, constant parameters

ge = lib.add_generic_entity("G",
	{"conc": "sum"},
	{}
);

pe = lib.add_generic_entity("P",
	{"conc": "sum"},
	{}
);


# generic processes: id, type, entities related, list of subprocesses, constant parameters, equations

lib.add_generic_process("logistic_growth", "growth",
	[("P",[pe],1,1)],
	[],
	{"growth_rate": (0,3), "k": (0,0.5)},
	{},
	{
		"P.conc": "growth_rate * P.conc * (1 - k * P.conc)"
	}
);

lib.add_generic_process("exponential_growth", "growth",
	[("P",[pe],1,1)],
	[],
	{"growth_rate": (0,2)},
	{},
	{
		"P.conc": "growth_rate * P.conc"
	}
);

lib.add_generic_process("holling_1", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"predation_rate": (0,1), "conversion_factor": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc"
	}
);

lib.add_generic_process("holling_2", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"predation_rate": (0,1), "conversion_factor": (0,1), "lambda": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc / (1 + lambda * predation_rate * P.conc)",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc / (1 + lambda * predation_rate * P.conc)"
	}
);

lib.add_generic_process("holling_3", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"predation_rate": (0,1), "conversion_factor": (0,1), "lambda": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc * P.conc / (1 + lambda * predation_rate * P.conc * P.conc)",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc * P.conc / (1 + lambda * predation_rate * P.conc * P.conc)"
	}
);

lib.add_generic_process("ratio_dependent_2", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"predation_rate": (0,1), "conversion_factor": (0,1), "lambda": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc / (G.conc + lambda * predation_rate * P.conc)",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc / (G.conc + lambda * predation_rate * P.conc)"
	}
);

lib.add_generic_process("ratio_dependent_3", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"predation_rate": (0,1), "conversion_factor": (0,1), "lambda": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc * P.conc / (G.conc + lambda * predation_rate * P.conc * P.conc)",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc * P.conc / (G.conc* G.conc + lambda * predation_rate * P.conc * P.conc)"
	}
);

lib.add_generic_process("ivlev", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"delta": (0,1), "predation_rate": (0,1), "conversion_factor": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * (1 - exp(-delta * P.conc))",
		"G.conc": "conversion_factor * predation_rate * G.conc * (1 - exp(-delta * P.conc))"
	}
);

lib.add_generic_process("watts", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"sigma": (0,1), "delta": (0,1), "predation_rate": (0,1), "conversion_factor": (0,1), "lambda": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * (1 - exp(-delta * P.conc / pow(G.conc, sigma)))",
		"G.conc": "conversion_factor * predation_rate * G.conc * (1 - exp(-delta * P.conc / pow(G.conc, sigma)))"
	}
);

lib.add_generic_process("hassell_varley_1", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"sigma": (0,1), "predation_rate": (0,1), "conversion_factor": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc * pow(G.conc, -sigma) ",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc * pow(G.conc,-sigma)"
	}
);

lib.add_generic_process("hassell_varley_2", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"sigma": (0,1), "predation_rate": (0,1), "conversion_factor": (0,1), "lambda": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc / (pow(G.conc, sigma) + lambda * predation_rate * P.conc)",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc / (pow(G.conc,sigma) + lambda * predation_rate * P.conc)"
	}
);

lib.add_generic_process("deangelis_beddington", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"delta": (0,1), "predation_rate": (0,1), "conversion_factor": (0,1), "lambda": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc / (1 + lambda * predation_rate * P.conc + delta * G.conc)",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc / (1 + lambda * predation_rate * P.conc + delta * G.conc)"
	}
);

lib.add_generic_process("crowley_martin", "predation",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[],
	{"delta": (0,1), "predation_rate": (0,1), "conversion_factor": (0,1), "lambda": (0,1)},
	{},
	{
		"P.conc": "-1 * predation_rate * G.conc * P.conc / ((1 + lambda * predation_rate * P.conc) * (1 + delta * G.conc))",
		"G.conc": "conversion_factor * predation_rate * G.conc * P.conc / ((1 + lambda * predation_rate * P.conc) * (1 + delta * G.conc))"
	}
);

lib.add_generic_process("exponential_loss", "loss",
	[("G",[ge],1,1)],
	[],
	{"loss_rate": (0,2)},
	{},
	{
		"G.conc": "-1 * loss_rate * G.conc"
	}
);


lib.add_generic_process("root", "",
	[("P",[pe],1,1), ("G",[ge],1,1)],
	[
		("growth", ["P"], 0),
		("predation", ["P", "G"], 0),
		("loss", ["G"], 0)
	],
	{}, {}, {}
);
