########## User-definable stuff ##########
#
###Compiler and compilation options
COMP_SER = cc
OPTIONS = -Wall -O3 -std=c99
#
### Behavioural flags
#Generate debug help. Only useful for development
DEFINEFLAGS += -D_DEBUG
#
###Path to libraries and headers
###If two or more of the dependencies reside in the same paths, only
###one instance is necessary.
#GSL
#GSL_INC = -I/add/path
#GSL_LIB = -L/add/path
GSL_INC = -I/usr/common/software/gsl/2.1/intel/include
GSL_LIB = -L/usr/common/software/gsl/2.1/intel/lib
#FFTW
FFTW_INC = -I/opt/cray/pe/fftw/3.3.6.2/haswell/include
FFTW_LIB = -L/opt/cray/pe/fftw/3.3.6.2/haswell/lib
#CCL
CCL_INC = -I/global/homes/d/damonge/include
CCL_LIB = -L/global/homes/d/damonge/lib
#
########## End of user-definable ##########

OPTIONS += $(DEFINEFLAGS)

INC_ALL = -I./src $(GSL_INC) $(FFTW_INC) $(CCL_INC)
LIB_ALL = $(CCL_LIB) $(GSL_LIB) $(FFTW_LIB) -lccl -lgsl -lgslcblas -lfftw3 -lm

FFTLOGO = src/fftlog.o
MAIN = src/main.c
OFILES = $(FFTLOGO)

EXEC = lnpred

default : $(EXEC)

%.o : %.c
	$(COMP_CC) $(OPTIONS) $(INC_ALL) -c $< -o $@

$(OFILES) : COMP_CC := $(COMP_SER)

$(EXEC) : $(OFILES)
	$(COMP_SER) $(OPTIONS) $(INC_ALL) $(OFILES) $(MAIN) -o $(EXEC) $(LIB_ALL)

clean :
	rm -f src/*.o

cleaner :
	rm -f *~ src/*.o src/*~  $(EXEC)
