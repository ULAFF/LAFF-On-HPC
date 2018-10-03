#HOME      := /Users/rvdg
# Make sure you have BLIS installed in your home directory.
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O3 -I$(BLAS_INC) -m64 -msse3 -std=c99 -march=native -D_POSIX_C_SOURCE=200112L
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.  The best time is reported.
NFIRST     := 48#     smallest size to be timed.
NLAST_SMALL:= 500#    largest size to be timed for slow implementations.
NLAST      := 2000#   largest size to be timed for fast implementations.
NINC       := 48#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# ---------------------

OBJS_IJP := driver.o Gemm_IJP.o

driver_IJP.x: $(OBJS_IJP) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP.x $(LDFLAGS) 

IJP: driver_IJP.x
	echo "Name = 'IJP';" > data/output_IJP.m
	echo " " >> data/output_IJP.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJP.x >> data/output_IJP.m 
	more data/output_IJP.m 

# ---------------------

OBJS_IPJ := driver.o Gemm_IPJ.o

driver_IPJ.x: $(OBJS_IPJ) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IPJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IPJ.x $(LDFLAGS) 

IPJ: driver_IPJ.x
	echo "Name = 'IPJ';" > data/output_IPJ.m
	echo " " >> data/output_IPJ.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IPJ.x >>  data/output_IPJ.m
	more data/output_IPJ.m

# ---------------------

OBJS_IP_Axpy := driver.o Gemm_IP_Axpy.o Axpy.o

driver_IP_Axpy.x: $(OBJS_IP_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IP_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_Axpy.x $(LDFLAGS) 

IP_Axpy: driver_IP_Axpy.x
	echo "Name = 'IP\_Axpy';" > data/output_IP_Axpy.m
	echo " " >> data/output_IP_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_Axpy.x >>  data/output_IP_Axpy.m
	more data/output_IP_Axpy.m

# ---------------------

OBJS_IP_daxpy := driver.o Gemm_IP_daxpy.o

driver_IP_daxpy.x: $(OBJS_IP_daxpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IP_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_daxpy.x $(LDFLAGS) 

IP_daxpy: driver_IP_daxpy.x
	echo "Name = 'IP\_daxpy';" > data/output_IP_daxpy.m
	echo " " >> data/output_IP_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_daxpy.x >> data/output_IP_daxpy.m
	more data/output_IP_daxpy.m

# ---------------------

OBJS_IP_bli_daxpyv := driver.o Gemm_IP_bli_daxpyv.o

driver_IP_bli_daxpyv.x: $(OBJS_IP_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IP_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_bli_daxpyv.x $(LDFLAGS) 

IP_bli_daxpyv: driver_IP_bli_daxpyv.x
	echo "Name = 'IP\_bli\_daxpyv';" > data/output_IP_bli_daxpyv.m
	echo " " >> data/output_IP_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_bli_daxpyv.x >>  data/output_IP_bli_daxpyv.m
	more data/output_IP_bli_daxpyv.m

# ---------------------

OBJS_JIP := driver.o Gemm_JIP.o

driver_JIP.x: $(OBJS_JIP) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JIP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JIP.x $(LDFLAGS) 

JIP: driver_JIP.x
	echo "Name = 'JIP';" > data/output_JIP.m
	echo " " >> data/output_JIP.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JIP.x >> data/output_JIP.m 
	more data/output_JIP.m

# ---------------------

OBJS_JPI := driver.o Gemm_JPI.o

driver_JPI.x: $(OBJS_JPI) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JPI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JPI.x $(LDFLAGS) 

JPI: driver_JPI.x
	echo "Name = 'JPI';" > data/output_JPI.m
	echo " " >> data/output_JPI.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JPI.x >> data/output_JPI.m 
	more data/output_JPI.m

# ---------------------

OBJS_PIJ := driver.o Gemm_PIJ.o

driver_PIJ.x: $(OBJS_PIJ) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PIJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PIJ.x $(LDFLAGS) 

PIJ: driver_PIJ.x
	echo "Name = 'PIJ';" > data/output_PIJ.m
	echo " " >> data/output_PIJ.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PIJ.x >> data/output_PIJ.m 
	more data/output_PIJ.m

# ---------------------

OBJS_PI_Axpy := driver.o Gemm_PI_Axpy.o Axpy.o

driver_PI_Axpy.x: $(OBJS_PI_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PI_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_Axpy.x $(LDFLAGS) 

PI_Axpy: driver_PI_Axpy.x
	echo "Name = 'PI_Axpy';" > data/output_PI_Axpy.m
	echo " " >> data/output_PI_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_Axpy.x >>  data/output_PI_Axpy.m
	more data/output_PI_Axpy.m

# ---------------------

OBJS_PI_daxpy := driver.o Gemm_PI_daxpy.o

driver_PI_daxpy.x: $(OBJS_PI_daxpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PI_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_daxpy.x $(LDFLAGS) 

PI_daxpy: driver_PI_daxpy.x
	echo "Name = 'PI_daxpy';" > data/output_PI_daxpy.m
	echo " " >> data/output_PI_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_daxpy.x >> data/output_PI_daxpy.m 
	more data/output_PI_daxpy.m

# ---------------------

OBJS_PI_bli_daxpyv := driver.o Gemm_PI_bli_daxpyv.o

driver_PI_bli_daxpyv.x: $(OBJS_PI_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PI_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_bli_daxpyv.x $(LDFLAGS) 

PI_bli_daxpyv: driver_PI_bli_daxpyv.x
	echo "Name = 'PI_bli_daxpyv';" > data/output_PI_bli_daxpyv.m
	echo " " >> data/output_PI_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_bli_daxpyv.x >> data/output_PI_bli_daxpyv.m 
	more data/output_PI_bli_daxpyv.m

# ---------------------

OBJS_IJ_DOTS := driver.o Gemm_IJ_Dots.o Dots.o

driver_IJ_Dots.x: $(OBJS_IJ_DOTS) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IJ_DOTS) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_Dots.x $(LDFLAGS) 

IJ_Dots: driver_IJ_Dots.x
	echo "Name = 'IJ\_Dots';" > data/output_IJ_Dots.m
	echo " " >> data/output_IJ_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_Dots.x >> data/output_IJ_Dots.m 
	more data/output_IJ_Dots.m

# ---------------------

OBJS_IJ_ddot := driver.o Gemm_IJ_ddot.o

driver_IJ_ddot.x: $(OBJS_IJ_ddot) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IJ_ddot) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_ddot.x $(LDFLAGS) 

IJ_ddot: driver_IJ_ddot.x
	echo "Name = 'IJ\_ddot';" > data/output_IJ_ddot.m
	echo " " >> data/output_IJ_ddot.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_ddot.x >> data/output_IJ_ddot.m
	more data/output_IJ_ddot.m

# ---------------------

OBJS_IJ_bli_ddotxv := driver.o Gemm_IJ_bli_ddotxv.o

driver_IJ_bli_ddotxv.x: $(OBJS_IJ_bli_ddotxv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IJ_bli_ddotxv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_bli_ddotxv.x $(LDFLAGS) 

IJ_bli_ddotxv: driver_IJ_bli_ddotxv.x
	echo "Name = 'IJ\_bli\_ddotxv';" > data/output_IJ_bli_ddotxv.m
	echo " " >> data/output_IJ_bli_ddotxv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_bli_ddotxv.x >>  data/output_IJ_bli_ddotxv.m
	more data/output_IJ_bli_ddotxv.m

# ---------------------

OBJS_JI_DOTS := driver.o Gemm_JI_Dots.o Dots.o

driver_JI_Dots.x: $(OBJS_JI_DOTS) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JI_DOTS) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_Dots.x $(LDFLAGS) 

JI_Dots: driver_JI_Dots.x
	echo "Name = 'JI\_Dots';" > data/output_JI_Dots.m
	echo " " >> data/output_JI_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_Dots.x >> data/output_JI_Dots.m
	more data/output_JI_Dots.m

# ---------------------

OBJS_JI_ddot := driver.o Gemm_JI_ddot.o

driver_JI_ddot.x: $(OBJS_JI_ddot) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JI_ddot) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_ddot.x $(LDFLAGS) 

JI_ddot: driver_JI_ddot.x
	echo "Name = 'JI\_ddot';" > data/output_JI_ddot.m
	echo " " >> data/output_JI_ddot.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_ddot.x >> data/output_JI_ddot.m
	more data/output_JI_ddot.m

# ---------------------

OBJS_JI_bli_ddotxv := driver.o Gemm_JI_bli_ddotxv.o

driver_JI_bli_ddotxv.x: $(OBJS_JI_bli_ddotxv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JI_bli_ddotxv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_bli_ddotxv.x $(LDFLAGS) 

JI_bli_ddotxv: driver_JI_bli_ddotxv.x
	echo "Name = 'JI\_bli\_ddotxv';" > data/output_JI_bli_ddotxv.m
	echo " " >> data/output_JI_bli_ddotxv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_bli_ddotxv.x >>  data/output_JI_bli_ddotxv.m
	more data/output_JI_bli_ddotxv.m

# ---------------------

OBJS_JP_Axpy := driver.o Gemm_JP_Axpy.o Axpy.o

driver_JP_Axpy.x: $(OBJS_JP_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JP_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_Axpy.x $(LDFLAGS) 

JP_Axpy: driver_JP_Axpy.x
	echo "Name = 'JP\_Axpy';" > data/output_JP_Axpy.m
	echo " " >> data/output_JP_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_Axpy.x >> data/output_JP_Axpy.m
	more data/output_JP_Axpy.m

# ---------------------

OBJS_JP_daxpy := driver.o Gemm_JP_daxpy.o 

driver_JP_daxpy.x: $(OBJS_JP_daxpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JP_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_daxpy.x $(LDFLAGS) 

JP_daxpy: driver_JP_daxpy.x
	echo "Name = 'JP\_daxpy';" > data/output_JP_daxpy.m
	echo " " >> data/output_JP_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_daxpy.x >>  data/output_JP_daxpy.m
	more data/output_JP_daxpy.m

# ---------------------

OBJS_JP_bli_daxpyv := driver.o Gemm_JP_bli_daxpyv.o 

driver_JP_bli_daxpyv.x: $(OBJS_JP_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JP_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_bli_daxpyv.x $(LDFLAGS) 

JP_bli_daxpyv: driver_JP_bli_daxpyv.x
	echo "Name = 'JP\_axpyv';" > data/output_JP_bli_daxpyv.m
	echo " " >> data/output_JP_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_bli_daxpyv.x >>  data/output_JP_bli_daxpyv.m
	more data/output_JP_bli_daxpyv.m

# ---------------------

OBJS_PJI := driver.o Gemm_PJI.o

driver_PJI.x: $(OBJS_PJI) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PJI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJI.x $(LDFLAGS) 

PJI: driver_PJI.x
	echo "Name = 'PJI';" > data/output_PJI.m
	echo " " >> data/output_PJI.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJI.x>>  data/output_PJI.m
	more data/output_PJI.m

# ---------------------

OBJS_PJ_Axpy := driver.o Gemm_PJ_Axpy.o Axpy.o

driver_PJ_Axpy.x: $(OBJS_PJ_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PJ_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_Axpy.x $(LDFLAGS) 

PJ_Axpy: driver_PJ_Axpy.x
	echo "Name = 'PJ\_Axpy';" > data/output_PJ_Axpy.m
	echo " " >> data/output_PJ_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_Axpy.x >> data/output_PJ_Axpy.m
	more data/output_PJ_Axpy.m

# ---------------------

OBJS_PJ_daxpy := driver.o Gemm_PJ_daxpy.o

driver_PJ_daxpy.x: $(OBJS_PJ_daxpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PJ_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_daxpy.x $(LDFLAGS) 

PJ_daxpy: driver_PJ_daxpy.x
	echo "Name = 'PJ\_daxpy';" > data/output_PJ_daxpy.m
	echo " " >> data/output_PJ_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_daxpy.x >>  data/output_PJ_daxpy.m
	more data/output_PJ_daxpy.m

# ---------------------

OBJS_PJ_bli_daxpyv := driver.o Gemm_PJ_bli_daxpyv.o

driver_PJ_bli_daxpyv.x: $(OBJS_PJ_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PJ_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_bli_daxpyv.x $(LDFLAGS) 

PJ_bli_daxpyv: driver_PJ_bli_daxpyv.x
	echo "Name = 'PJ\_bli\_daxpyv';" > data/output_PJ_bli_daxpyv.m
	echo " " >> data/output_PJ_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_bli_daxpyv.x >>  data/output_PJ_bli_daxpyv.m
	more data/output_PJ_bli_daxpyv.m

# ---------------------

OBJS_I_dgemv := driver.o Gemm_I_dgemv.o

driver_I_dgemv.x: $(OBJS_I_dgemv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_I_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_I_dgemv.x $(LDFLAGS) 

I_dgemv: driver_I_dgemv.x
	echo "Name = 'I\_dgemv';" > data/output_I_dgemv.m
	echo " " >> data/output_I_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_I_dgemv.x >> data/output_I_dgemv.m data/output_I_dgemv.m
	more data/output_I_dgemv.m

# ---------------------

OBJS_I_bli_dgemv := driver.o Gemm_I_bli_dgemv.o

driver_I_bli_dgemv.x: $(OBJS_I_bli_dgemv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_I_bli_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_I_bli_dgemv.x $(LDFLAGS) 

I_bli_dgemv: driver_I_bli_dgemv.x
	echo "Name = 'I\_bli\_dgemv';" > data/output_I_bli_dgemv.m
	echo " " >> data/output_I_bli_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_I_bli_dgemv.x >> data/output_I_bli_dgemv.m data/output_I_bli_dgemv.m
	more data/output_I_bli_dgemv.m

# ---------------------

OBJS_P_GER_J_AXPY := driver.o Gemm_P_Ger.o Ger_J_Axpy.o Axpy.o

driver_P_Ger_J_Axpy.x: $(OBJS_P_GER_J_AXPY) $(UTIL_OBJS)
	$(LINKER) $(OBJS_P_GER_J_AXPY) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_Ger_J_Axpy.x $(LDFLAGS) 

P_Ger_J_Axpy: driver_P_Ger_J_Axpy.x
	echo "Name = 'P\_Ger\_J\_Axpy';" > data/output_P_Ger_J_Axpy.m
	echo " " >> data/output_P_Ger_J_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_Ger_J_Axpy.x >>  data/output_P_Ger_J_Axpy.m
	more data/output_P_Ger_J_Axpy.m

# ---------------------

OBJS_P_GER_I_AXPY := driver.o Gemm_P_Ger.o Ger_I_Axpy.o Axpy.o

driver_P_Ger_I_Axpy.x: $(OBJS_P_GER_I_AXPY) $(UTIL_OBJS)
	$(LINKER) $(OBJS_P_GER_I_AXPY) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_Ger_I_Axpy.x $(LDFLAGS) 

P_Ger_I_Axpy: driver_P_Ger_I_Axpy.x
	echo "Name = 'P\_Ger\_I\_Axpy';" > data/output_P_Ger_I_Axpy.m
	echo " " >> data/output_P_Ger_I_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_Ger_I_Axpy.x >>  data/output_P_Ger_I_Axpy.m
	more data/output_P_Ger_I_Axpy.m

# ---------------------

OBJS_P_DGER := driver.o Gemm_P_dger.o

driver_P_dger.x: $(OBJS_P_DGER) $(UTIL_OBJS)
	$(LINKER) $(OBJS_P_DGER) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_dger.x $(LDFLAGS) 

P_dger: driver_P_dger.x
	echo "Name = 'P\_dger';" > data/output_P_dger.m
	echo " " >> data/output_P_dger.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_dger.x >>  data/output_P_dger.m
	more data/output_P_dger.m

# ---------------------

OBJS_P_BLI_DGER := driver.o Gemm_P_bli_dger.o

driver_P_bli_dger.x: $(OBJS_P_BLI_DGER) $(UTIL_OBJS)
	$(LINKER) $(OBJS_P_BLI_DGER) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_bli_dger.x $(LDFLAGS) 

P_bli_dger: driver_P_bli_dger.x
	echo "Name = 'P\_bli\_dger';" > data/output_P_bli_dger.m
	echo " " >> data/output_P_bli_dger.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_bli_dger.x >>  data/output_P_bli_dger.m
	more data/output_P_bli_dger.m

# ---------------------

OBJS_J_GEMV_J_Axpy := driver.o Gemm_J_Gemv.o Gemv_J_Axpy.o Axpy.o

driver_J_Gemv_J_Axpy.x: $(OBJS_J_GEMV_J_Axpy) $(UTIL_OBJS)
	$(LINKER) $(OBJS_J_GEMV_J_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_Gemv_J_Axpy.x $(LDFLAGS) 

J_Gemv_J_Axpy: driver_J_Gemv_J_Axpy.x
	echo "Name = 'J\_Gemv\_J\_Axpy';" > data/output_J_Gemv_J_Axpy.m
	echo " " >> data/output_J_Gemv_J_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_Gemv_J_Axpy.x >> data/output_J_Gemv_J_Axpy.m data/output_J_Gemv_J_Axpy.m
	more data/output_J_Gemv_J_Axpy.m

# ---------------------

OBJS_J_dgemv := driver.o Gemm_J_dgemv.o

driver_J_dgemv.x: $(OBJS_J_dgemv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_J_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_dgemv.x $(LDFLAGS) 

J_dgemv: driver_J_dgemv.x
	echo "Name = 'J\_dgemv';" > data/output_J_dgemv.m
	echo " " >> data/output_J_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_dgemv.x >> data/output_J_dgemv.m data/output_J_dgemv.m
	more data/output_J_dgemv.m

# ---------------------

OBJS_J_bli_dgemv := driver.o Gemm_J_bli_dgemv.o

driver_J_bli_dgemv.x: $(OBJS_J_bli_dgemv) $(UTIL_OBJS)
	$(LINKER) $(OBJS_J_bli_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_bli_dgemv.x $(LDFLAGS) 

J_bli_dgemv: driver_J_bli_dgemv.x
	echo "Name = 'J\_dgemv';" > data/output_J_bli_dgemv.m
	echo " " >> data/output_J_bli_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_bli_dgemv.x >> data/output_J_bli_dgemv.m data/output_J_bli_dgemv.m
	more data/output_J_bli_dgemv.m

# ---------------------

Plot_All_Orderings:
	make IJP
	make IPJ
	make JIP
	make JPI
	make PIJ
	make PJI
# --------------------- 

Plot_Inner_I:
	make JPI
	make JP_Axpy
	make JP_daxpy
	make JP_bli_daxpyv
	make PJI
	make PJ_Axpy
	make PJ_daxpy
	make PJ_bli_daxpyv

# ---------------------                                                                    

Plot_Inner_J:
	make IPJ
	make IP_Axpy
	make IP_daxpy
	make IP_bli_daxpyv
	make PIJ
	make PI_Axpy
	make PI_daxpy
	make PI_bli_daxpyv

# ---------------------                                                               

Plot_Inner_P:
	make IJP
	make IJ_Dots
	make IJ_ddot
	make IJ_bli_ddotxv
	make JIP
	make JI_Dots
	make JI_ddot
	make JI_bli_ddotxv

# ---------------------                                                           
    
Plot_Outer_I:
	make IJP
	make IPJ
	make I_dgemv
	make I_bli_dgemv

# ---------------------                                                           
    
Plot_Outer_J:
	make JIP
	make JPI
	make J_dgemv
	make J_bli_dgemv

# ---------------------                                                               
Plot_Outer_P:
	make PIJ
	make PJI
	make P_Ger_J_Axpy
	make P_Ger_I_Axpy
	make P_dger
	make P_bli_dger

# ---------------------                                                               
NLAST_SMALL:= 1000 # Change the size of the largest problem tested

Plot_All_Outer:
	make I_dgemv
	make I_bli_dgemv
	make J_dgemv
	make J_bli_dgemv
	make P_dger
	make P_bli_dger

# ---------------------                                                               
clean:
	rm -f *.o *~ core *.x *.pdf
