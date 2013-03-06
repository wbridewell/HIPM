/* Copyright (c) 2008, Institute for the Study of Learning and Expertise
 * All rights reserved.
 * For details, see the LICENSE file.
 */

/* #define __ALLOW_MISSING_VALUES */

/* this is a silly experimental marker, so caveat */
#define MV_MARKER -100000

/*#include <malloc.h>*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <cvode.h>
#include <cvspgmr.h>
#include <nvector_serial.h>

/* on 64bit machines: fortran_int should be just int */
#define fortran_int long 
/* #define fortran_int int */

/**
 * ASSUMPTIONS:
 *
 * 1. the model simulation should be sampled at the same rate as the observed data
 * 2. 
 */

/* global structures and variables */

/* 1. data */

struct {

  /* multiple data sets */
  int nsets;                 /* number of data files */
  char **set_files;          /* [0..(data.nsets-1)]: array of data file names */
  int *set_tos;              /* [0..(data.nsets-1)]: array of the number of 
				measurements included in the data files up to 
				the current index. For example, if we have three 
				data files of length 10, 20, and 30, the values 
				in this array should be: 10, 10+20=30, and 
				10+20+30=60.                                  */

  /* data dimensions */
  int nvars;                 /* number of observed variables */
  int length;                /* number of measurements in all the data files */
  char **sys_var_names;      /* [0..(sys_var.n)]: array of variable names (index, then system variables) */

  double **table;            /* data table */
  
  /* time and system variables */
  double *time;              /* pointer to the time variable in the data 
				(typically equals data.table[0], which means 
				that the time variable is the first column in 
				the data file).                              */
  double **sys_vars;         /* [0..(sys_var.n-1)]: array of pointers to system 
				variables in the data file. For unobserved 
				variables set to NULL.                       */

  double *sys_vars_disps;    /* dispersions for error normalization */
  /* min and max used for data aggregation function */
  double *min, *max;
  
  int current_set;           /* current data set (used for simulation) */

} data;


/* 2. simulation buffers and cvode stuff */

struct {

  double *buffer_1;
  double *buffer_2;

  fortran_int success_flag;

  /* CVode stuff */
  void *cvode_mem;

  N_Vector y0;
  realtype t0;

  N_Vector abstol;
  realtype reltol;

} sim;


/* 3. system variables */

struct {

  int n;                     /* the number of system variables */     
  int nobserved;             /* the number of observed system variables */

  double *vals;              /* current values */
  double *dot_vals;          /* time derivatives of the current values */

} sys_var;


/* 4. parameters */

struct {

  fortran_int n;             /* number of parameters in the ODEs */ 

  double *vals;              /* [0..(param.n-1)]: initial values of the parameters
				(used for the first restart of the optimization 
				code - in all the other restarts, inital values 
				are randomly generated). */
  double *optimal_vals;      /* [0..(param.n-1)]: storage for the best parameters
				found at any time during parameter estimation. */

  double *bounds;            /* [0..(2*param.n-1)]: lower and upper bounds of the 
				parameter values. param.bounds[2*i] is the lower 
				and param.bounds[2*i+1] is the upper bound for 
				the i-th parameter. */
} param;


/* 5. optimization parameters and alg-717 stuff */

struct {

  int init_state_fit;        /* fit initial values of the system variables? */
  int normalize_error;       /* normalize SSE with dispersion? */
  int tf_restarts;           /* number of TF optimization restarts */
  int fs_restarts;           /* number of FS optimization restarts */
  
  /* alg-717 stuff */
  fortran_int n;
  fortran_int p;
  
  int liv;
  fortran_int *iv;
  
  int lv;
  double *v;
  
} opt;

/* controls whether a CVODE simulation prints the values of all variables as it goes along */
int print_full_sim = 0;

/* alternative combining schemes */
double min(double x, double y) {
  return x < y ? x : y;
}

double max(double x, double y) {
  return x > y ? x : y;
}

/* data aggregation */

double datamin(double *x, int n) {

  double min = x[0];
  int i;

  for (i = 1; i < n; ++i) if (x[i] < min) min = x[i];
  return(min);
}

double datamax(double *x, int n) {

  double max = x[0];
  int i;

  for (i = 1; i < n; ++i) if (x[i] > max) max = x[i];
  return(max);
}

/* initialization and allocation procedures */

/* 1. intialize structure dimensions */

/* model/data specific dimensions */
/* defined in model.c */

extern void MS_init_dims(void);

/* fixed part */

void init_dims(void) {

  /* default optimization parameters */
  opt.init_state_fit = 1;
  opt.normalize_error = 1;
  opt.tf_restarts = 64;
  opt.fs_restarts = 4;

  MS_init_dims();

  opt.p = param.n + data.nsets * sys_var.n;
  opt.n = sys_var.nobserved * data.length;

  opt.liv =  5 * opt.p + 82;
  opt.lv = opt.n * (opt.p + 6) + 4 * opt.p + opt.p * (2 * opt.p + 20) + 105;
}

/* 2. allocate memory for global structures and variables */

void mem_alloc(void) {

  int i;

  data.set_files = (char **) malloc(data.nsets * sizeof(char *));
  data.set_tos = (int *) malloc(data.nsets * sizeof(int));

  data.table = (double **) malloc(data.nvars * sizeof(double *));
  data.min = (double *) malloc(data.nvars * sizeof(double));
  data.max = (double *) malloc(data.nvars * sizeof(double));

  for (i = 0; i < data.nvars; ++i)
    data.table[i] = (double *) malloc(data.length * sizeof(double));
  data.sys_vars = (double **) malloc(sys_var.n * sizeof(double *));
  data.sys_vars_disps = (double *) malloc(sys_var.n * sizeof(double));
  /* index at 0, all others +1 of their position in sys_var.vals */
  data.sys_var_names = (char **) malloc((sys_var.n+1) * sizeof(char *));

  sim.buffer_1 = (double *) malloc(opt.n * sizeof(double));
  sim.buffer_2 = (double *) malloc(opt.n * sizeof(double));

  sys_var.vals = (double *) malloc(sys_var.n * sizeof(double));
  sys_var.dot_vals = (double *) malloc(sys_var.n * sizeof(double));

  param.vals = (double *) malloc(opt.p * sizeof(double));
  param.optimal_vals = (double *) malloc(opt.p * sizeof(double));
  param.bounds = (double *) malloc(2 * opt.p * sizeof(double));

  opt.iv = (fortran_int *) calloc(opt.liv, sizeof(fortran_int));
  opt.v = (double *) calloc(opt.lv, sizeof(double));
}


/* 3. fill in the values */

/* data specific values */
/* defined in model.c */

extern void MS_fill_in(void);

/* parameters initial values and bounds */
/* defined in model.c */

extern void MS_params(void);

/* fixed part */

void fill_in(void) {

  extern double max_double_vector(double *x, int n);
  extern double min_double_vector(double *x, int n);
  extern double disp(double *x, int n);

  char *fn, buffer[256];
  FILE *f;
  int set, set_from, set_to, i, j;

  /* set to 0 to let cvode write to standard error */
  int suppress_cvode = 1;

  FILE *errfp;
  extern void cvode_model(realtype t, N_Vector cvode_y, N_Vector cvode_dy_dt, void *data);


  MS_fill_in();

  /* read data sets from files */
  for (set = i = 0; set < data.nsets; ++set) {
    if ((f = fopen(fn = data.set_files[set], "r")) == NULL) {
      printf("Initialization error: Can not open file %s\n", fn);
      exit(1);
    }

    /* skip variable names */
    for (j = 0; j < data.nvars; ++j) {
      fscanf(f, " %s", buffer);
    }
    
    for (set_to = data.set_tos[set]; i < set_to; ++i)
      for (j = 0; j < data.nvars; ++j) {
        fscanf(f, " %lg", data.table[j] + i);
      }
    fclose(f);
  }

  for (j = 0; j < data.nvars; ++j) {
    if (data.table[j] == NULL) continue;
    data.min[j] = datamin(data.table[j], data.length);
    data.max[j] = datamax(data.table[j], data.length);
  }

  /* calculate dispersions for error normalization */
  for (j = 0; j < sys_var.n; ++j) {
    if (data.sys_vars[j] == NULL) continue;
    data.sys_vars_disps[j] = disp(data.sys_vars[j], data.length);
  }

  /* cvode initialization */
  if ((sim.cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON)) == NULL) {
    printf("Library error: CVodeCreate failed.\n");
    exit(2);
  }

  /* cvode output to /dev/null */
  if (suppress_cvode) {
    errfp = fopen("/dev/null", "w");
    if (CVodeSetErrFile(sim.cvode_mem, errfp) != 0) {
      printf("Library error: CVodeSetErrFile failed.\n");
      exit(2);
    }
  }

  /* cvode malloc */
  sim.reltol = 1.0e-4;
  sim.abstol = N_VNew_Serial(sys_var.n);
  for (j = 0; j < sys_var.n; ++j) NV_Ith_S(sim.abstol, j) = 1e-8;
  sim.t0 = data.time[0];
  sim.y0 = N_VNew_Serial(sys_var.n);
  if (CVodeMalloc(sim.cvode_mem, cvode_model, sim.t0, sim.y0, CV_SV, &(sim.reltol), sim.abstol) != 0) {
    printf("Library error: CVodeMalloc failed.\n");
    exit(2);
  }

  /* cvode linear solver */
  if (CVSpgmr(sim.cvode_mem, PREC_NONE, 0) != 0) {
    printf("Library error: CVDense failed.\n");
    exit(2);
  }

  /* parameters for fitting the initial values of the state vars */
  for (set = 0, i = param.n; set < data.nsets; ++set) {
    set_from = (set == 0) ? 0 : data.set_tos[set - 1] + 1;
    for (j = 0; j < sys_var.n; ++j, ++i) {
      if (data.sys_vars[j] == NULL) continue;
      param.vals[i] = data.sys_vars[j][set_from];

      /* establish the bounds for the observed variables */

      /* HACK code for *single* data set where we fit the unobserved, but not the observed variables */
      /* param.bounds[2*i] = data.sys_vars[j][set_from] - data.sys_vars[j][set_from] * 0.001;
	 param.bounds[2*i+1] = data.sys_vars[j][set_from] + data.sys_vars[j][set_from] * 0.001;*/

      /* original code for observed variable bounds */
            param.bounds[2*i] = min_double_vector(data.sys_vars[j], data.length);
            param.bounds[2*i+1] = max_double_vector(data.sys_vars[j], data.length);
    }
  }
  MS_params();
}


/* the model itself (model specific :-) */
/* defined in model.c */

extern void MS_model(double t);

/* note that this function is not model specific */

void cvode_model(realtype t, N_Vector y, N_Vector dy_dt, void *data) {

  int j;

  for (j = 0; j < sys_var.n; ++j) sys_var.vals[j] = NV_Ith_S(y, j);
  MS_model((double) t);
  for (j = 0; j < sys_var.n; ++j) NV_Ith_S(dy_dt, j) = sys_var.dot_vals[j];
}



/* FROM HERE ON: everything is model independent */


/* 1. misc. utility functions */

/* for given t find the index in the current data set */
int t_to_index(double t) {
 
  /* prepare indices into the time array */
  int set_from = (data.current_set == 0) ? 0 : data.set_tos[data.current_set - 1] + 1;
  int set_to = data.set_tos[data.current_set];
  int i;

  /* if t is larger than the last recorded time, return t */
  if (t > data.time[set_to - 1]) return(set_to - 1);

  /* return the index of the time closest to t that does not go past t */
  for (i = set_from; i < set_to; ++i)
    if (data.time[i] >= t) break;
  return((data.time[i] == t) ? i : i - 1);
}


/* print current values of the parameters */

void print_double_vector(char *head, double *x, int n) {

  int i;

  printf("%s:", head);
  for (i = 0; i < n; ++i) printf(" %.16g", x[i]);
  printf("\n");
}


/* calculate min, max, dispersion */

double min_double_vector(double *x, int n) {

  int i;
  double min;

  for (i = 1, min = x[0]; i < n; ++i) {
#ifdef __ALLOW_MISSING_VALUES
    if (x[i] == MV_MARKER) continue;
#endif
    if (x[i] < min) min = x[i];
  }
  return(min);
}

double max_double_vector(double *x, int n) {

  int i;
  double max;

  for (i = 1, max = x[0]; i < n; ++i)
    if (x[i] > max) max = x[i];
  return(max);
}

double disp(double *x, int n) {

  int i;
#ifdef __ALLOW_MISSING_VALUES
  int length = 0;
#else
  int length = n;
#endif
  double sx, sxx, var;

  for (i = 0, sx = sxx = 0.0; i < n; ++i) {
#ifdef __ALLOW_MISSING_VALUES
    if (x[i]  == MV_MARKER) continue;
    else length++;
#endif
    sx += x[i];
    sxx += x[i] * x[i];
  }

  return(sxx - sx * sx / n);
}

/* Fill the parameters with random values that fall between the two bounds */
void randomize_params(void) {
  int i;
  double lower, upper;

  for (i = 0; i < param.n; ++i) {
    lower = param.bounds[2*i];
    upper = param.bounds[2*i + 1];
    
    param.vals[i] = lower + (upper - lower) * ((double) rand()) / ((double) RAND_MAX);
  }
}

/* estimators of the model error (or goodness of fit) */
/* to allow us to deal with missing values we assume that */
/* x[] contains the observed values and y[] containst the simulated */
/* if x[i] == MV_MARKER we'll consider it missing for the purposes of */
/* calculating the SSE. */
double SSE_double_vector(double *x, double *y, int n) {

  int i;
  double sse, d;

  for (i = 0, sse = 0.0; i < n; ++i) {
#ifdef __ALLOW_MISSING_VALUES
    if (x[i] == MV_MARKER) continue;
#endif
    d = x[i] - y[i];
    sse += d * d;
  }
  return(sse);
}

double SSE(double *sim, int normalize_flag) {

  int j, jo;
  double sse, sse_j;


  for (j = jo = 0, sse = 0.0; j < sys_var.n; ++j) {
    if (data.sys_vars[j] == NULL) continue;
    sse_j = SSE_double_vector(data.sys_vars[j], sim + (jo * data.length), data.length);
    sse += sse_j / ((normalize_flag == 0) ? 1.0 : data.sys_vars_disps[j]);
    ++jo;
  }
  return(sse);
}

/* due to the data.length division, values returned from RMSE, AIC, and 
 * BIC will be incorrect when some values are missing.  We need to keep 
 * a count of the number of non-missing values for each variable for 
 * this to work correctly 
 */
/* Root Mean Squared Error */
double RMSE(double *sim, int normalize_flag) {
  return(sqrt(SSE(sim, normalize_flag) / data.length));
}


/* Mean Squared Error */
double MSE(double *sim, int normalize_flag) {
  return(SSE(sim, normalize_flag) / data.length);
}

double inequality_coefficient(double *sim, int normalize_flag) {

  /* I think this equation is incorrect */
  /* sqrt(MSE) / (sqrt(sum(observed^2/n)) + sqrt(sum(predicted^2/n))) */
  return 0.0;
}

/* Akaike Information Criterion */
double AIC(double *sim, int normalize_flag) {
  /* assumes that the residuals are normally distributed */
  /* note: check the assumption for the likelihood function */
  return 2 * param.n + data.length * log(MSE(sim, normalize_flag));
}


/* Bayesian Information Criterion */
double BIC(double *sim, int normalize_flag) {
  /* note: check the assumption for the likelihood function */
  return data.length * log(MSE(sim, normalize_flag)) + param.n * log(data.length);
}

/* Coefficient of Determination */
double r2(double *sim) {

  int i, j, jo;
#ifdef __ALLOW_MISSING_VALUES
  int n_non_missing = 0;
#endif
  double r2, r2_j, sumx, sumxx, ssxx, sumy, sumyy, ssyy, sumxy, ssxy, *y, *yhat;

  sumxy = sumx = sumy = sumxx = sumyy = 0.0; 
  for (j = jo = 0, r2 = 0.0; j < sys_var.n; ++j) {
    if ((y = data.sys_vars[j]) == NULL) continue;
    yhat = sim + (jo * data.length);
    for (i = 0; i < data.length; ++i) {
#ifdef __ALLOW_MISSING_VALUES
      if (y[i] == MV_MARKER) continue;
      n_non_missing += 1;
#endif
      sumx += y[i];
      sumxx += y[i] * y[i];

      sumy += yhat[i];
      sumyy += yhat[i] * yhat[i];

      sumxy += y[i] * yhat[i];
    }
    ++jo;
  }
#ifdef __ALLOW_MISSING_VALUES
  if (n_non_missing == 0) {
    r2 = 0.0;
  } else {
    ssxx = sumxx - sumx * sumx / n_non_missing;
    ssyy = sumyy - sumy * sumy / n_non_missing;
    ssxy = sumxy - sumx * sumy / n_non_missing;
    r2 = (ssxy * ssxy) / (ssxx * ssyy);
  }
#else    
  ssxx = sumxx - sumx * sumx / (data.length * sys_var.nobserved);
  ssyy = sumyy - sumy * sumy / (data.length * sys_var.nobserved);
  ssxy = sumxy - sumx * sumy / (data.length * sys_var.nobserved);
  r2 = (ssxy * ssxy) / (ssxx * ssyy);
#endif

  return r2;
}


/* print observed vs. simulated table */
/* LJUPCO: variable names and unobservables? */

void print_obs_vs_sim(void) {

  int i, j, jo;
  /* print the names of the variables */
  /* first print time and assume that it is the first variable name */
  /* then print the remaining variables */
  printf ("%s", data.sys_var_names[0]);
  for (i = 0; i < sys_var.n; ++i) {
    if (data.sys_vars[i] == NULL) continue;
    printf (" %s %s_sim", data.sys_var_names[i+1], data.sys_var_names[i+1]);
  }
  printf("\n");
  /* print the values of the variables */
  for (i = 0; i < data.length; ++i) {
    printf("%g", data.time[i]);
    for (j = jo = 0; j < sys_var.n; ++j) {
      if (data.sys_vars[j] == NULL) continue;
#ifdef __ALLOW_MISSING_VALUES
      if (data.sys_vars[j][i] == MV_MARKER)
        printf(" - %g", sim.buffer_2[jo * data.length + i]);
      else
        printf(" %g %g", data.sys_vars[j][i], sim.buffer_2[jo * data.length + i]);
#else
      printf(" %g %g", data.sys_vars[j][i], sim.buffer_2[jo * data.length + i]);
#endif
      ++jo;
    }
    printf("\n");
  }
}



/* 2. model simulation code */

/* teacher forcing simulation - implementation with integrals */
/* if you *really* want to know what are the parameters - see ALG-717 paper */

void sim_model_TF(fortran_int *n, fortran_int *p, double *x, fortran_int *nf, 
		  fortran_int *need, double *r, double *rp, long *ui, 
		  double *ur, int (*uf)()) {

  double integral;
  int i, j, set_to;


  /* calculate the derivatives/slopes and store them in simulation buffer */
  for (i = data.current_set = 0; i < data.length; ++i) {
    if (i > data.set_tos[data.current_set]) ++data.current_set;
    for (j = 0; j < sys_var.n; ++j) sys_var.vals[j] = data.sys_vars[j][i];
    MS_model(data.time[i]);
    for (j = 0; j < sys_var.n; ++j){
      if (isnan(sys_var.vals[j]) || isinf(sys_var.vals[j])) {
	/* in case of nan or inf: report error to the optimization code */
	*nf = 0;
	return;
      }
      sim.buffer_1[j * data.length + i] = sys_var.dot_vals[j];
    }
  }

  /* calculate the integral - teacher forcing simulation and store it in r */
  for (j = 0; j < sys_var.n; ++j) {
    for (data.current_set = i = 0; data.current_set < data.nsets; ++data.current_set) {
      set_to = data.set_tos[data.current_set];
      r[j * data.length + i] = integral = data.sys_vars[j][i];
      for (++i; i < set_to; ++i) {
	/* trapezian formula */
	integral += ((sim.buffer_1[j * data.length + i] + sim.buffer_1[j * data.length + i - 1]) 
		     * (data.time[i] - data.time[i - 1]) * 0.5);
	r[j * data.length + i] = integral;
      }
    }
  }
  /* report OK to the optimization code */
  *nf = 1;
}

void print_cvode_error(int cv_ret) {
	printf("Simulation error: ");
	
    switch(cv_ret) {
	case CV_MEM_NULL:
	  printf("CVode memory was null.\n");
	  break;
	case CV_NO_MALLOC:
	  printf("CVode memory was not allocated.\n");
	  break;
	case CV_ILL_INPUT:
	  printf("Illegal CVode input.\n");
	  break;
	case CV_TOO_MUCH_WORK:
	  printf("Too many internal CVode steps.\n");
	  break;
	case CV_TOO_MUCH_ACC:
	  printf("CVode accuracy demands too strict.\n");
	  break;
	case CV_ERR_FAILURE:
	  printf("Too many CVode errors.\n");
	  break;
	case CV_CONV_FAILURE:
	  printf("CVode cannot converge on a solution.\n");
	  break;
	case CV_LSETUP_FAIL:
	  printf("CVode's linear solver failed to initialize.\n");
	  break;
	case CV_LSOLVE_FAIL:
	  printf("CVode's linear solver failed.\n");
	  break;
	}
}

/* full simulation of the model - using CVODE */
/* if you *really* want to know what are the parameters - see ALG-717 paper */
void sim_model_CVODE(fortran_int *n, fortran_int *p, double *x, fortran_int *nf, 
		     fortran_int *need, double *r, double *rp, long *ui, 
		     double *ur, int (*uf)()) {

  int i, j, jo, set_to, cv_ret;
  realtype t;

  for (data.current_set = i = 0; data.current_set < data.nsets; ++data.current_set) {
    set_to = data.set_tos[data.current_set];

    /* re-initialize CVode */
    sim.t0 = data.time[i];
    for (j = jo = 0; j < sys_var.n; ++j) {
      if (opt.init_state_fit == 0 && data.sys_vars[j] != NULL)
	/* use the data file to initialize the observed variables */
	NV_Ith_S(sim.y0, j) = data.sys_vars[j][i];
      else
	/* use fits to the initial values, or pre-specified initial values */
	NV_Ith_S(sim.y0, j) = param.vals[param.n - (data.nsets - data.current_set) * sys_var.n + j];

      if (data.sys_vars[j] == NULL) continue;
      r[jo * data.length + i] = NV_Ith_S(sim.y0, j);
      ++jo;
    }
    if (CVodeReInit(sim.cvode_mem, cvode_model, sim.t0, sim.y0, CV_SV, &(sim.reltol), sim.abstol) != 0) {
      printf("Library error: CVodeReInit failed.\n");
      *nf = 0;
      return;
    }

    if (print_full_sim) {
      printf ("%s", data.sys_var_names[0]);
      
      for (j = 0; j < sys_var.n; ++j) {
	printf (" %s", data.sys_var_names[j+1]);
      }
      printf("\n");
    }


    /* run CVODE simulation step by step and store results */
    for (++i; i < set_to; ++i) {
      cv_ret = CVode(sim.cvode_mem, (realtype) data.time[i], sim.y0, &t, CV_NORMAL);
      /* if necessary, print a formal error message */
      if (cv_ret != CV_SUCCESS && 
	  cv_ret != CV_ROOT_RETURN && 
	  cv_ret != CV_TSTOP_RETURN) {
	print_cvode_error(cv_ret);
	*nf = 0;
	return;
      }
      
      if (print_full_sim) { printf("%g", data.time[i]); }
      for (j = jo = 0; j < sys_var.n; ++j) {
	if (print_full_sim) { printf(" %g", NV_Ith_S(sim.y0, j)); }

	if (data.sys_vars[j] == NULL) continue;
	r[jo * data.length + i] = NV_Ith_S(sim.y0, j);
	++jo;
      }
      if (print_full_sim) { printf("\n"); }
    }
  }
  /* report OK to the optimization code */
  *nf = 1;
}

/*
 * evaluate model simulation (stored in r) as required by ALG-717:
 *   need[0]: store model 0.5 * SSE in *f
 *   need[1]: store derivatives of the error for each point in r and rp
 *            since error = 0.5 * ((y - yhat) / norm) ^ 2
 *                  d error / d y = (y - yhat) / norm
 *                  d2 error / d y2 = 1 / norm
 *   other parameters: see ALG-717 paper
 */

void eval_model(fortran_int *need, double *f, fortran_int *n, fortran_int *nf, 
		double *xn, double *r, double *rp, fortran_int *ui, 
		double *ur, double *w) {

  int i, j, jo;
  double *y, *yhat;
  double norm;

  if (need[0] == 1) {
    *f = 0.5 * SSE(r, opt.normalize_error);
    if (isnan(*f) || isinf(*f)) {
      /* report error to the optimization code */
      *nf = 0;
      return;
    }
  } else {
    for (j = jo = 0; j < sys_var.n; ++j) {
      if ((y = data.sys_vars[j]) == NULL) continue;
      yhat = r + (jo * data.length);
      norm = (opt.normalize_error == 0) ? 1.0 : data.sys_vars_disps[j];
      for (i = 0; i < data.length; ++i) {
#ifdef __ALLOW_MISSING_VALUES
       if (y[i] == MV_MARKER) 
          r[jo * data.length + i] = 0.0;
        else
#endif 
	  r[jo * data.length + i] = (yhat[i] - y[i]) / norm;
	rp[jo * data.length + i] = (double) 1.0 / norm;
      }
      ++jo;
    }
  }
  /* report OK to the optimization code */
  *nf = 1;
}

/* multistart optimization procedure */
void run_fit(void (*sim_model)(), int nrestarts) {

  extern void dglfb_();
  int r, i, sentinel;
  double sse, min_sse = 0.0;


  for (r = sentinel = 0; (r < nrestarts) && (sentinel < 10 * nrestarts); ++r, ++sentinel) {
    /* initialize random parameter values: in first iteration use default values */
    if (r != 0) {
      randomize_params();
    }

    /* try out the parameter values and reinitialize if simulation is not possible */
    sim_model(&(opt.n), &(param.n), param.vals, &(sim.success_flag), NULL, sim.buffer_2, NULL, NULL, NULL, NULL);
    if (sim.success_flag == 0) {
      /* since we're not randomizing the initial set of parameters, move past them upon failure */
      if (r != 0) --r;
      continue;
    }

    /* restore default parameters for runing ALG-717 */
    for (i = 0; i < opt.liv; ++i) opt.iv[i] = (fortran_int) 0;
    for (i = 0; i < opt.lv; ++i) opt.v[i] = (double) 0.0;

    /* run ALG-717 */
    dglfb_(&(opt.n), &(param.n), &(param.n), param.vals, param.bounds,
           eval_model, NULL, NULL, opt.iv, &(opt.liv), &(opt.lv), opt.v,
           sim_model, NULL, NULL, NULL);

    /* simulate the model with the obtained parameter values */
    sim_model(&(opt.n), &(param.n), param.vals, &(sim.success_flag), NULL, sim.buffer_2, NULL, NULL, NULL, NULL);
    if (sim.success_flag == 0) continue;
    sse = SSE(sim.buffer_2, 0);


    /* update optimal parameter values */
    if ((sse < min_sse) || (min_sse == 0.0)) {
      min_sse = sse;
      for (i = 0; i < param.n; ++i) param.optimal_vals[i] = param.vals[i];
    }
  }
  /* update parameter values */
  for (i = 0; i < param.n; ++i) param.vals[i] = param.optimal_vals[i];
}


int main(void) {

  int tmp;
  float f = 0.0;

  init_dims();
  mem_alloc();
  fill_in();

  /* Simulate only */
  if (opt.tf_restarts + opt.fs_restarts == 0) {
    print_full_sim = 1;
    if (opt.init_state_fit == 1) param.n = opt.p;
    sim_model_CVODE(&(opt.n), &(param.n), param.vals, &(sim.success_flag), NULL, sim.buffer_2, NULL, NULL, NULL, NULL);
    /* print three asterisks to separate the full simulation output from the obs vs. sim output */
    printf("***\n");
    print_obs_vs_sim();
    if (sim.success_flag != 0) {
      printf("SSE, reMSE, r2, AIC, BIC: %f %f %f %f %f\n",
	     SSE(sim.buffer_2, 0), SSE(sim.buffer_2, 1), r2(sim.buffer_2),
	     AIC(sim.buffer_2, 0), BIC(sim.buffer_2, 0));
    }

    return(0);
  }

/* as implemented, teacher forcing does not apply when missing values are present/allowed */
#ifndef __ALLOW_MISSING_VALUES
  /* Teacher forcing is allowed, so try it first */
  if (sys_var.n == sys_var.nobserved) {
    run_fit(sim_model_TF, opt.tf_restarts);
    tmp = opt.init_state_fit;
    opt.init_state_fit = 0;
    sim_model_CVODE(&(opt.n), &(param.n), param.vals, &(sim.success_flag), NULL, sim.buffer_2, NULL, NULL, NULL, NULL);
    if (sim.success_flag != 0) {
      /* print_double_vector("PARAMS", param.vals, param.n); */
      /* printf("SSE, SSEn, r2: %f %f %f\n", SSE(sim.buffer_2, 0), SSE(sim.buffer_2, 1), r2(sim.buffer_2)); */
    }
    opt.init_state_fit = tmp;
  }
#endif
  /* Call the parameter optimization code over the full trajectory */
  if (opt.init_state_fit == 1) param.n = opt.p;
  if (opt.fs_restarts > 0) {
    run_fit(sim_model_CVODE, opt.fs_restarts);
    sim_model_CVODE(&(opt.n), &(param.n), param.vals, &(sim.success_flag), NULL, sim.buffer_2, NULL, NULL, NULL, NULL);
  }
  if (sim.success_flag != 0) {
    print_double_vector("PARAMS", param.vals, param.n);
    printf("SSE, reMSE, r2, AIC, BIC: %f %f %f %f %f\n", 
	   SSE(sim.buffer_2, 0), SSE(sim.buffer_2, 1), r2(sim.buffer_2),
	   AIC(sim.buffer_2, 0), BIC(sim.buffer_2, 0));
  }
  return(0);
}
