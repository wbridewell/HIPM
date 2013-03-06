# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

CC = gcc
#CFLAGS = -ggdb -Wall -pedantic
CFLAGS = -O

A717_LFLAGS = ./alg-717/dgletc.o ./alg-717/dglfgb.o ./alg-717/dmdc.o -lgfortran -lm
CVODE_LFLAGS = -L./cvode -lsundials_cvode -lsundials_nvecserial -lsundials_shared
CVODE_IFLAGS = -I./cvode

all: clean model

clean:
	rm -f model model.o model.c

model: model.o
	$(CC) -o $@ model.o $(A717_LFLAGS) $(CVODE_LFLAGS)

model.o: model.c
	$(CC) $(CFLAGS) $(CVODE_IFLAGS) -c -o $@ model.c

model.c: model-t.c ms.c
	cat model-t.c ms.c > $@
