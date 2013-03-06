#!/usr/bin/python

# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

"""
This generic model library supports the construction of an ecosystem model of the Ross Sea.
It is hierarchical in processes, but the entities are flat.
This version is updated and corrected.  It is designed for use with the sensitivity analysis experiments
Last revision: SRB 4/16/07

Revision Record
04/16/07  added 2 limited growth functions (logistic)
04/09/07  saved as ross.py for work with satellite data. changed mixing rate parameters for yr1
04/06/07  change rspd_2.py to rspdb.py; also extend lower bound of mixing rate
02/22/07  ross_27 = rspd
02/08/07  ross_25 = ross_23
02/05/07  Ross_22 is based on Ross_15.  In this run, grazing, zoop resp and death are optional.  
12/11/06  Ross_15 is based on Ross_8 (the most sucessful run to date).  Here I am changing the range of the 
	  Zooplankton assimilation efficency from (0.74,0.76) to (0.05,0.4).  This better matches existing 
	  knowledge regarding assimilation efficencies and encompasses the classic 10% rule.   
11/6/06	  This function is the next step from ross_7.  
	    (1) It uses the old rsp.data
	    (2) I have changed the light control functions to match those in Arrigo et al. 1998
10/27/06  replaced PAR with 0.4 * PAR ~ PUR.  Reset light intensity parameter values.
10/26/06  changed (1- E.ice) to (E.ice) to accomodate new data in "limited_growth"
10/20/06  Change P limitation to basic 
9/26/06   (1) names changed for use in rsp-hipm paper, 
	  (2) nutrient limitation functions expanded
9/19/06   corrected mixing process.
"""

from library import *;
from entities import *;
from processes import *;

lib = library("aquatic_ecosystem");

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# 	GENERIC ENTITIES
# 				id, variables, constant parameters
# -----------------------------------------------------------------------

# --- PHYTOPLANKTON ---
pe = lib.add_generic_entity("P",
	{ # variables
		"conc": "sum", # ug C/L
		"growth_rate": "prod", 
		"growth_lim": "min"
	},
	{ # constant parameters
		"max_growth": (0.4,0.8), # schoemann et al 2005 suggest 0.4-1.5, d-1
		"exude_rate": (0.001,0.2), # d-1
		"death_rate": (0.02,0.04), # d-1
		"Ek_max": (1,100), 				
		"sinking_rate": (0.0001,0.25), 
		"biomin": (0.02,0.04), 
		"PhotoInhib": (200,1500),
#		"chla_C_ratio": (0.023,0.025), # mg Chla to mg C -- data says 24 mg Chla to 1 g C (Schoemann et al. 2005)
	}
);

# --- ZOOPLANKTON ---
ze = lib.add_generic_entity("Z",
	{ # variables
		"conc": "sum", # ug C/L
		"grazing_rate": "prod"
	},
	{ # constant parameters
		"assim_eff": (0.05,0.4), 
		"death_rate": (0.001,0.3), 
		"respiration_rate": (0.01,0.04),
		"sinking_rate": (0.001,0.25), 
		"gmax": (0.3,0.5), 
		"glim": (19,21), 
		"gcap": (199,301),
		"attack_rate": (0.3,10)	
	}
);

# --- NUTRIENTs ---

# nitrate
no3 = lib.add_generic_entity("Nitrate",
	{ # variables
		"conc": "sum",
		"mixing_rate": "sum"
	}, # umol/L
	{ # constant parameters
		"toCratio": (6.6,6.7), # mol:mol ratio -- this represents the phytoplankton C:N ratio (redfield)
		"avg_deep_conc": (31,32)
	} # umol/L
);

# iron
fe = lib.add_generic_entity("Iron",
	{ # variables
		"conc": "sum",
		"mixing_rate": "sum"
	}, # umol/L
	{ # constant parameters
		"toCratio": (3000,450000), # mol:mol as in Arrigo and Tagliabue 2005
		"avg_deep_conc": (0.00035,0.00045) # umol/L as in Arrigo and Tagliabue 2005
	}
);

# --- DETRITUS ---
de = lib.add_generic_entity("D",
	{ # variables
		"conc": "sum" # umol C/L
	},
	{ # constant parameters
		"remin_rate": (0.03,0.04), # d-1
		"sinking_rate": (0.00001,0.1) # d-1
	}
);

# --- ENVIRONMENT ---
ee = lib.add_generic_entity("E",
	{ # variables
		"TH2O": "sum", # degrees C
		"PUR": "sum", # uE m-2 s-1
		"ice": "sum" # proportion.  0<ice<1
	},
	{ # constant parameters
		"beta": (0.001,1), # 
	}
);

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# 	GENERIC PROCESSES:
# 					id, type, entities related, list of subprocesses, 
#								constant parameters, equations
# -----------------------------------------------------------------------

# --- GROWTH --- 

lib.add_generic_process("growth", "growth",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("N",[no3,fe],1,100),
		("D",[de],1,1),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
		("limited_growth", ["P","N","E"], 0), 
		("exudation", ["P"], 0), 
		("nutrient_uptake",["P","N"], 0)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"P.conc": "P.growth_rate * P.conc"
	}
);

lib.add_generic_process("exudation", "exudation",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"P.conc": "-1 * P.exude_rate * P.growth_rate * P.conc"
	}
);

lib.add_generic_process("nutrient_uptake", "nutrient_uptake",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("N",[no3,fe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"N.conc": "-1 * 1/( N.toCratio * 12.0107) * P.growth_rate * P.conc" # 12.0107 is the mol weight of C (g C/mol) or (ugC/umol)
	}
);

lib.add_generic_process("limited_growth", "limited_growth",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("N",[no3,fe],1,100),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
		("light_lim",["P","E"],0),
		("nutrient_lim",["P","N"],0)
	],
	{ # constant parameters
	},
	{ # algebraic equations
		"P.growth_rate": "(1.0-E.ice) * P.max_growth * exp(0.06933 * E.TH2O) * P.growth_lim" #(P.growth_lim is a minimum of light and nutrient limitations)
	},
	{ # differential equations
	}
);

lib.add_generic_process("Pearl_Verhulst_logistic", "limited_growth",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("N",[no3,fe],1,100),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"K": (20,10000) # Flag: this range is set to accomodate the ross sea phytoplankton - but may need to be adjusted.
	},
	{ # algebraic equations
		"P.growth_rate": "P.max_growth * (1.0 - P.growth_lim/K)"
	},
	{ # differential equations
	}
); 
	
lib.add_generic_process("Gompertz_logistic", "limited_growth",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("N",[no3,fe],1,100),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"K": (1,10000) # Flag: this range is set to accomodate the ross sea phytoplankton - but may need to be adjusted.
	},
	{ # algebraic equations
		"P.growth_rate": "P.max_growth * (log(K) - log(P.conc))" # referenced in Rosenzweig 1971
	},
	{ # differential equations
	}
); 


# ------ P.growth_lim --
'''
Multiple factors (and formulations of factors) might limit growth.
In this library nutrient and light limitations are combined into P.growth_lim using a minimum function
so that only one operates at a time (i.e., they are substitutable).  The disadvantage
of this encoding is that it will not be possible to determine which factor is operating at a given time.  Temperature
is a multiplicative control factor encoded in the P.growth_rate equation, and in the present library we do not consider
alternative temperature effect functions.  
'''

# --light lim --

lib.add_generic_process("arrigoetal1998", "light_lim",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"a": (5,15)
	},
	{ # algebraic equations
		"P.growth_lim": "(1.0 - exp(-E.PUR / (P.Ek_max / (1.0 + a * exp(E.PUR * exp(1.089 - 2.12 * log10(P.Ek_max)))))))"
	},
	{ # differential equations
	}
);

lib.add_generic_process("arrigoetal1998_w_photoinhibition", "light_lim",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"a": (5,15),
	},
	{ # algebraic equations
		"P.growth_lim": "(1.0 - exp(-E.PUR / (P.Ek_max / (1.0 + a * exp(E.PUR * exp(1.089 - 2.12 * log10(P.Ek_max))))))) * exp(-1.0 * E.PUR /P.PhotoInhib)"
	},
	{ # differential equations
	}
);

# -- nutrient lim -- 	

lib.add_generic_process("monod_lim", "nutrient_lim",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("N",[no3,fe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"k":(0.000001,0.001) # this is based on Kmu_max estimates for NO3 in Schoemann et al. 2005; based on Tagliabue and Arrigo (2005) this range should work for Fe too.
	},
	{ # algebraic equations
		"P.growth_lim": "N.conc / (N.conc + k)"
	},
	{ # differential equations
	}
);

''' this process has been excluded until we find a rational for the equations (see Weigerts literature)
lib.add_generic_process(
	"ratio_lim", "nutrient_lim",
	[("P",[pe],1,1), ("N",[no3,fe],1,1)],
	[],
	{"k":(0.000001,1)},
	{"P.growth_lim": "N.conc / (N.conc + k * P.conc)"},
	{}
	);
'''

lib.add_generic_process("monod_2nd", "nutrient_lim",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("N",[no3,fe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"k": (0.000001,0.001)
	},
	{ # algebraic equations
		"P.growth_lim": "(N.conc * N.conc) / (N.conc * N.conc + k)"
	},
	{ # differential equations
	}
);

lib.add_generic_process("nut_lim_exp", "nutrient_lim",
	[ # parameters (id, generic_entities, min, max)
		("P",[pe],1,1),
		("N",[no3,fe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"k": (0.000001,1) # recall that the parameter in this equation has the inverse interpretation from the monod 1/2 saturation constant.
	},
	{ # algebraic equations
		"P.growth_lim": "1.0 -exp(-1.0 * k * N.conc)"
	},
	{ # differential equations
	}
);


# --- DEATH ---

lib.add_generic_process("death_exp", "death",
	[ # parameters (id, generic_entities, min, max)
		("S",[pe,ze],1,1),
		("D",[de],1,1),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"S.conc": "-1.0 * S.death_rate * S.conc",
		"D.conc": "(1.0 - E.beta) * S.death_rate * S.conc"
	}
);

'''
Second order exponential death (see Steee and Henderson 1992, Edwards 2001).  
This is often used to model both natural mortality as well as loss to grazers of a higher trophic 
level not explicitly modeled.  The could of course be generalized. This is sometimes called a "closure term".
'''	
lib.add_generic_process("death_exp2", "death",
	[ # parameters (id, generic_entities, min, max)
		("S",[pe,ze],1,1),
		("D",[de],1,1),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"S.conc": "-1 * S.death_rate * (S.conc * S.conc)",
		"D.conc": "(1.0 - E.beta) * S.death_rate * (S.conc * S.conc)"
	}
); 


# --- REMINERALIZATION ---

lib.add_generic_process("remineralization", "",
	[ # parameters (id, generic_entities, min, max)
		("D",[de],1,1),
		("N",[fe,no3],1,3)
	],
	[ # subprocesses (id, parameters, optional_flag)
		("nutrient_remineralization",["D","N"],0)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"D.conc": "-1.0 * D.remin_rate * D.conc"
	}
);

lib.add_generic_process("nutrient_remineralization", "",
	[ # parameters (id, generic_entities, min, max)
		("D", [de], 1,1),
		("N", [fe], 1, 1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"N.conc": "(1.0 / (N.toCratio * 12.0107)) * D.remin_rate * D.conc" 
			# 12.0107 is the mol weight of C (g C/mol); the 1000 is to convert from mol to mmol
	}
);


# --- RESPIRATION ---

lib.add_generic_process("respiration", "",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"Z.conc":"-1 * Z.respiration_rate * Z.conc"
	}
);


# --- SINKING ---

lib.add_generic_process("sinking", "",
	[ # parameters (id, generic_entities, min, max)
		("V",[pe,ze,de],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"V.conc": "-1.0 * V.sinking_rate * V.conc"
	}
);
# should sinking be a different function?  Edwards 2001 uses a simple exponential.

	
# --- GRAZING ---

lib.add_generic_process("lotka_volterra", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.gmax * P.conc"
	},
	{ # differential equations
	}
); # this is a Holling Type I functional response.  

lib.add_generic_process("generalized_gause", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"alpha": (0.0001,0.9999) # if alpha = 0 then the function reduces, if alpha = 1 then it equals the lotka_volterra function
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.gmax * pow(P.conc, alpha)"
	},
	{ # differential equations
	}
); # reported in Rosenzweig 1971, a generalization of lotka-volterra

lib.add_generic_process("monod", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
		"Z.grazing_rate": "max(0.0, Z.gmax * P.conc / (Z.gcap + P.conc))" # active predators
	},
	{ # differential equations
	}
);  # this comes from Arrigo et al. 1998 (CIAO model)

lib.add_generic_process("monod_mod", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
	},
	{ # algebraic equations
		"Z.grazing_rate": "max(0.0, (Z.gmax * (P.conc - P.biomin - Z.glim) / (Z.gcap +  (P.conc - P.biomin - Z.glim))))"
	},
	{ # differential equations
	}
);		# this comes from Arrigo et al. 1998 (CIAO model)

lib.add_generic_process("holling_type_2", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"h": (1,5)		
			# this is the extra handling time: L.Gross says - Handling time is defined as the time spent pursuing, subduing, and consuming each prey item plus the time spent preparing to search for the next prey item (including effects of satiation) :: SRB - I don't know if this is appropriate
			# if h < 1 then the max of the response function exceeds 1. 
	},
	{ # algebraic equations
		"Z.grazing_rate": "(Z.attack_rate * P.conc) / (1 + Z.attack_rate * h * P.conc)"
	},
	{ # differential equations
	}
);		# this is Holling's Disk equation

lib.add_generic_process("holling_type_3", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"h": (1,5)		
			# this is the extra handling time: L.Gross says - Handling time is defined as the time spent pursuing, subduing, and consuming each prey item plus the time spent preparing to search for the next prey item (including effects of satiation) :: SRB - I don't know if this is appropriate
			# if h < 1 then the max of the response function exceeds 1. 
	},
	{ # algebraic equations
		"Z.grazing_rate": "(Z.attack_rate * P.conc*P.conc) / (1 + Z.attack_rate * h * P.conc*P.conc)"
	},
	{ # differential equations
	}
);

lib.add_generic_process("ivlev", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"delta": (0.001,1)	
			# delta could be greater than 1.
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.gmax * (1.0 - exp(-1.0 * delta * P.conc))"
	},
	{ # differential equations
	}
); # the Ivlev also appears in Rozensweig 1971 -- with out a reference to Ivlev.


# ---

lib.add_generic_process("ratio_dependent_2", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"h": (1,5)
			# this is the extra handling time: L.Gross says - Handling time is defined as the time spent pursuing, subduing, and consuming each prey item plus the time spent preparing to search for the next prey item (including effects of satiation) :: SRB - I don't know if this is appropriate
			# if h < 1 then the max of the response function exceeds 1.
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.attack_rate * P.conc / (Z.conc + Z.attack_rate * h * P.conc)"
	},
	{ # differential equations
	}
);

lib.add_generic_process("ratio_dependent_3", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"h": (1,5)		
			# this is the extra handling time: L.Gross says - Handling time is defined as the time spent pursuing, subduing, and consuming each prey item plus the time spent preparing to search for the next prey item (including effects of satiation) :: SRB - I don't know if this is appropriate
			# if h < 1 then the max of the response function exceeds 1.
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.attack_rate * (P.conc * P.conc) / (Z.conc * Z.conc + Z.attack_rate * h * P.conc * P.conc)"
	},
	{ # differential equations
	}
);


lib.add_generic_process("watts", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"delta": (0.01,0.5),			# higher value increases grazing rate
		"m": (0.0001,0.9999)	    # higher value retards grazing rate
			# range taken from the the alpha parameter of the generalized_gause process, since watts is a generalization of ivlev analogous to the generalized_gause generalization of lotka_volterra
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.gmax * (1.0 - exp(-1.0 * delta * P.conc / pow(Z.conc,m)))"
			# the assumption here is that Z interfere with their own grazing; thus, more Z implies less grazing.
	},
	{ # differential equations
	}
);

lib.add_generic_process("hassell_varley_1", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"sigma": (1,100)	# higher value retards growth -- SRB not sure of appropriate range
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.gmax * P.conc * pow(Z.conc, -1.0 * sigma)"
	},
	{ # differential equations
	}
);

lib.add_generic_process("hassell_varley_2", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"sigma": (0.0001,0.9999),
		"h": (1,5)		
			# this is the extra handling time: L.Gross says - Handling time is defined as the time spent pursuing, subduing, and consuming each prey item plus the time spent preparing to search for the next prey item (including effects of satiation) :: SRB - I don't know if this is appropriate
			# if h < 1 then the max of the response function exceeds 1.
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.attack_rate * P.conc / (pow(Z.conc,sigma) + Z.attack_rate * h* P.conc)"
	},
	{ # differential equations
	}
);

lib.add_generic_process("deangelis_beddington", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"delta": (0.0,1.0),
		"h": (1,5)		
			# this is the extra handling time :: SRB - I don't know if this is appropriate
			# if h < 1 then the max of the response function exceeds 1.
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.attack_rate * P.conc / (1.0 + Z.attack_rate * h * P.conc + delta * Z.conc)"
	},
	{ # differential equations
	}
);

lib.add_generic_process("crowley_martin", "graze_rate",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"delta": (0.0,1.0),
		"h": (1,5)		
			# this is the extra handling time :: SRB - I don't know if this is appropriate
			# if h < 1 then the max of the response function exceeds 1.
	},
	{ # algebraic equations
		"Z.grazing_rate": "Z.attack_rate * P.conc / ((1.0 + Z.attack_rate * h * P.conc) * (1.0 + delta * Z.conc))"
	},
	{ # differential equations
	}
);

lib.add_generic_process("grazing", "grazing",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],0,1),
		("D",[de],0,1),
		("E",[ee],0,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
		("graze_rate",["Z","P"],0)
	],
	{ # constant parameters
	},
	{# algebraic equations
	},
	{ # differential equations
		"Z.conc": "Z.assim_eff * Z.grazing_rate * Z.conc",
		"P.conc": "-1 * Z.grazing_rate * Z.conc",
		"D.conc": "(1-E.beta) * (1-Z.assim_eff) * Z.grazing_rate * Z.conc"
	}
);


# --- Nutrient Mixing ----------------------------------------------------------------
# this process represents an input of nutrients (nitrate) due to mixing or upwelling.  

lib.add_generic_process("nutrient_mixing", "",
	[ # parameters (id, generic_entities, min, max)
		("N",[no3,fe],1,1),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
		("mixing_rate",["N","E"],0)
	],
	{ # constant parameters
	},
	{ # algebraic equations
	},
	{ # differential equations
		"N.conc": "(N.avg_deep_conc - N.conc) * N.mixing_rate"
	}
);
# notice that mathematically and practically this function can be negative becuase the nutrient gradient 
# can be negative; this is not realistic, but will be ok if it happens.  If it does occur, this process will be
# a loss term rather than a gain term -- which makes sense. srb(9/19/06)

lib.add_generic_process("linear_temp_control", "mixing_rate",
	[ # parameters (id, generic_entities, min, max)
		("N",[no3,fe],1,1),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
	],
	{ # constant parameters
		"max_mixing_rate": (0.000001,1) # FLAG: The lower bound was extended as per KRA
	},
	{ # algebraic equations
		"N.mixing_rate": "max_mixing_rate * (datamax(E.TH2O) - E.TH2O) / (datamax(E.TH2O) - datamin(E.TH2O))" # the numerator of the function is the temperature gradient; the denominator scales the function to be between 0 and 1
	},
	{ # differential equations
	}
);


# --- ROOT ---

lib.add_generic_process("root", "",
	[ # parameters (id, generic_entities, min, max)
		("Z",[ze],1,1),
		("P",[pe],1,2),
		("N",[no3,fe],2,2),
		("D",[de],1,1),
		("E",[ee],1,1)
	],
	[ # subprocesses (id, parameters, optional_flag)
		("growth",["P","N","D","E"], 0),
		("death",["P","D","E"],0),
		("death",["Z","D","E"],0),
		("grazing",["Z","P","D","E"],0),
		("remineralization",["D","N"],0),
		("respiration",["Z"],0),
		("sinking",["P"],0),
		("sinking",["D"],0),
		("nutrient_mixing",["N","E"],0), # this process represesnts an input to nutrients from mixing with below the mixed layer depth
	],
	{ # constant parameters
	},
	{# algebraic equations
	},
	{ # differential equations
	}
);
