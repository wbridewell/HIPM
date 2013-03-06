				HIPM:
	       Software for inductive process modeling



Copyright (c) 2008 
Institute for the Study of Learning and Expertise
See the LICENSE file for details.



Contents
  * Introduction
  * System Requirements
  * Installation Instructions
  * Example Modeling Problems
  * HIPM Source Files
  * Web Pages
  * Relevant Citations
  * Acknowledgments



Introduction:

HIPM is an inductive process modeler that addresses the task defined
in Bridewell et al. (2008).  Unlike the earlier IPM and RPM programs,
HIPM supports entities that group associated variables and parameters
much like structs.  HIPM also organizes the generic process library
into a hierarchy that serves as a context-free grammar for generating
model structures.  Todorovski et al. (2005) introduced HIPM, and a 
more detailed technical report is in progress.

This version of HIPM does not support conditions on the processes does
not output the value of variables defined by algebraic equations.
These limitations are addressed in the SC-IPM program.



System Requirements:

We originally developed HIPM for use on the Linux operating system,
but we have successfully ported it to other operating systems and
hardware architectures.

In particular, the program is known to run, with some reconfiguration,
on 32- and 64-bit versions of (see instructions for TOMS-717 and
model-t.c for configuration details): 
  * Linux (Intel x86 processors)
  * OS X (Intel x86 and PowerPC processors)
  * SunOS (Sparc processors)

We expect that HIPM will run on any Unix-like system that has the
following support software:
  * C compiler (e.g., gcc, cc)
  * Fortran compiler (e.g., g77, gfortran)
  * Python 2.5 (other versions work, but may require minor code changes)
  * CVODE 2.2.1 (later versions use a slightly different API)
      CVODE is a numerical solver for systems of ordinary differential
      equations.  The source code, license, and build instructions are
      included with the HIPM software.
  * TOMS-717
      This Fortran code implements Algorithm-717 (1993), a parameter
      estimation routine that has a long development history.
      Specifically, it builds on Algorithm-573 (1981) and
      Algorithm-611 (1983), both published in the ACM Transactions on
      Mathematical Software.  The source code, license, and build
      instructions are included with the HIPM software.



Installation Instructions:

The HIPM archive contains a prebuilt version of the software that runs
on Linux with a recent C compiler (e.g., gcc 4.2) and gfortran.
Recent (2007/2008) versions of Fedora Core, CentOS, and Ubuntu are
known to meet these requirements.  In this case, the contents of the
"hipm" directory should work without further effort.  Otherwise, see
the build instructions below.

Note that other operating systems (specifically versions of OS X) may
require slight changes in this procedure.

Building TOMS-717:
  * unpack TOMS-717 with tar -zxf toms-717.tar.gz
  * cd TOMS-717
  * if you are using a non-x86 processor, be sure to edit dmdc.f0
    for your architecture and store the resulting file as dmdc.f.
    the included dmdc.f works for 32-bit x86 machines.
  * determine which of g77 or gfortran is installed
  * edit the following makefile line to reflect the Fortran environment
      F77=gfortran 
      F77=g77
  * make clean
  * make
  * copy dgletc.o dglfgb.o dmdc.o to the HIPM alg-717 directory

Building CVODE 2.2.1:
  * unpack CVODE with tar -zxf cvode-ser-2.2.1.tar.gz
  * cd sundials
  * ./configure --prefix=`pwd`
  * make; make install
  * copy include/* to the hipm/cvode directory
  * copy lib/* to the hipm/cvode directory  
  
HIPM Makefile
  * if using gfortran, ensure that the A717_LFLAGS includes -lgfortran
  * if using g77, ensure that the A717_FLAGS includes -lg77

HIPM model-t.c
  * if attempting to build for a 64-bit platform, determine whether
    "#define fortran_int long" or "#define fortran_int int" should be 
    used.


Example Modeling Problems:

To give a feel for using HIPM, we provide three example scenarios.
Each one consists of a generic process library, one or more data
files, an entity definition file, and an experiment run file. The file
LIBRARY-FORMAT.readme in the "examples" directory describes the format
of a generic process library.

Note that HIPM is primarily a Python program, and with the exception
of the data files, the input files are written in Python syntax.  Data
files are formatted with space-separated columns.  The first column is
assumed to be the "time" index, and the first row is assumed to
contain variable names.  The data sets are described by README files
in their respective directories.

Of the provided libraries, the one for population dynamics is the most
straightforward.  This library relates two species in a predatory
relationship.  Classical models include a growth process for the prey,
a death process for the predator, and a predation process that relates
the two.  Each of these process types may have alternative functional
forms (e.g., logistic and exponential versions of the growth process).
The provided library defines a space of 22 model structures and ranges
on their associated parameters.  In comparison, the aquatic ecosystems
and fjord dynamics libraries define over 1,000 structures depending
largely on the entities being modeled.  The file "protist_entities.py"
details the code for entity instantiation.  The file "run-protists.py"
shows how to search the space of models and simulate the best model.

Bridewell et al. (2008) gives a detailed description of the protist
data and the development of the population dynamics library as it
existed before the introduction of entities and hierarchical
constraints in HIPM.  That paper also describes modeling tasks
associated with the Ross Sea and Ringkobing Fjord, which are also
included as examples.

  Protist data 
    -- copy population_dynamics_lib.py, protist_entities.py, and
       run-protists.py from examples directory to hipm directory
    -- copy 1cs2.txt from data/protists directory to hipm directory
    -- cd to hipm directory
    -- ./run-protists.py

  Ross Sea data
    -- copy aquatic_ecosystems_lib.py, ross_entities.py, and 
       run-ross-sea.py from examples directory to hipm directory
    -- copy ross-sea-yr1.data from data/ross-sea to hipm directory
    -- cd to hipm directory
    -- ./run-ross-sea.py

  Ringkobing Fjord data
    -- copy fjord_lib.py, ringkobing_entities.py, and run-ringkobing.py
       from examples directory to hipm directory
    -- copy fjord.data from data/ringkobing to hipm directory
    -- cd to hipm directory
    -- ./run-ringkobing.py



Websites:

HIPM is part of a larger research program on the computational
induction of scientific process models at the Institute for the Study
of Learning and Expertise.  Information on this project, along with an
extensive list of related publications is at 

     http://www.isle.org/process.html

The Institute for the Study of Learning and Expertise (ISLE) is a
nonprofit organization dedicated to research and education in machine
learning, artificial intelligence, and cognitive science.  The ISLE
website is at

     http://www.isle.org/



HIPM Source Files:

The "hipm" directory contains the source files for HIPM.

  * data.py
      includes code that defines a data structure for storing time
      series and reading them from a file.

  * entities.py
      includes code that defines objects for generic entities and
      instantiated entities.

  * process.py 
      includes code that defines objects for generic processes and
      instantiated ones.  contains the code for fitting parameters to
      model structures, printing models in C format, and simulating
      models.

  * library.py
      defines the structure of a generic process library.

  * search.py
      contains code for carrying out exhaustive and beam search of the
      models defined by a generic process library.  describes the
      arguments to these functions.

  * simulate.py
      contains code that reads a model as printed by HIPM and stores
      it in an internal data structure.  this is useful for building
      models by hand and for simulating models generated in prior
      HIPM runs.
      
  * utilities.py
      contains alternate code for reading models from a text file and
      for generating a table of model constants.

  * model-t.c
      this is the C template for model simulation and parameter
      estimation.  to carry out these tasks, the Python code outputs
      the file ms.c.  this model specific code is appended to
      model-t.c, the resulting file model.c is compiled through a call
      to gcc and the resulting executable "model" is created according
      to the instructions in the Makefile.

  * Makefile
      the instructions for building a model for the purposes of
      simulation or parameter estimation.



Relevant Citations:

HIPM played a key role in the following articles.  The software was
introduced in the 2005 paper by Todorovski et al. and directly builds
on the IPM and RPM systems and Todorovski's Lagramge/Lagrange software
for equation discovery.

Bridewell, W., Langley, P., Todorovski, L., & Dzeroski, S. (2008).
Inductive process modeling. Machine Learning, 71, 1--32.

Bridewell, W., Borrett, S., & Todorovski, L. (2007). Extracting
constraints for process modeling. Proceedings of the Fourth
International Conference on Knowledge Capture (pp. 87--94). 
Whistler, BC.

Bridewell, W., & Todorovski, L. (2007). Learning declarative
bias. Proceedings of the Seventeenth International Conference on
Inductive Logic Programming. Corvallis, OR.

Bridewell, W., Langley P., Racunas, S., & Borrett, S. (2006). Learning
process models with missing data. Proceedings of the Seventeenth
European Conference on Machine Learning, 557--565.

Langley, P., Shiran, O., Shrager, J., Todorovski, L., & Pohorille, A.
(2006). Constructing explanatory process models from biological data
and knowledge. AI in Medicine, 37, 191--201.

Bridewell, W., Bani Asadi, N., Langley, P., & Todorovski, L. (2005).
Reducing overfitting in process model induction. Proceedings of the
Twenty-Second International Conference on Machine Learning, 81--88.

Todorovski, L., Bridewell, W., Shiran, O., & Langley, P. (2005).
Inducing hierarchical process models in dynamic domains. Proceedings
of the Twentieth National Conference on Artificial Intelligence,
892--897.



Acknowledgments:

HIPM is based upon work supported by the National Science Foundation
under Grant IIS-0326059.  Any opinions, findings, and conclusions or
recommendations expressed in this material do not necessarily reflect
the views of the National Science Foundation.  Contributors to HIPM
include Ljupco Todorovski, Will Bridewell, and Oren Shiran.  Stuart
Borrett, Ljupco Todorovski, and Will Bridewell compiled the generic
process libraries for the example problems.  Kevin Arrigo and Gert van
Dijken provided data and domain knowledge for the Ross Sea.
