#FC			=gfortran#ifort
FC			=ifort
FFLAGS	=-O3
INCLUDE	=-I$(mylib)/include -I/opt/local/include
LIBS		=-lheader -lInterpolation -L$(mylib) -lfftw3 -lm -L/opt/local/lib
MOD		=mod
#OBJS		=twoFAST.o special_functions.o cdgamma.o
OBJS_qc	=twoFAST_qc.o cqlgamma.o
d_OBJS_qc	=direct_twoFAST_qc.o cqlgamma.o
OBJS_c	=twoFAST_complex.o cdlgamma.o
OBJS		=twoFAST.o cdlgamma.o 
TARGET	=twoFAST.exe twoFAST_c.exe twoFAST_qc.exe direct_twoFAST_qc.exe

.SUFFIXES:	.f90
.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDE)

default : $(TARGET)

twoFAST.exe : $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LIBS) $(INCLUDE)

twoFAST_c.exe : $(OBJS_c)
	$(FC) $(FFLAGS) -o $@ $(OBJS_c) $(LIBS) $(INCLUDE)

twoFAST_qc.exe : $(OBJS_qc)
	$(FC) $(FFLAGS) -o $@ $(OBJS_qc) $(LIBS) $(INCLUDE)

direct_twoFAST_qc.exe : $(d_OBJS_qc)
	$(FC) $(FFLAGS) -o $@ $(d_OBJS_qc) $(LIBS) $(INCLUDE)

clean:
	rm -rf $(OBJS) $(TARGET)
	rm -f *.$(MOD)

