# HOME      := /Users/rvdg
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
NLAST_SMALL:= 1000#    largest size to be timed for slow implementations.
NLAST      := 2000#   largest size to be timed for fast implementations.
NINC       := 48#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# ---------------------

TEST_OBJS_IJP := driver.o Gemm_IJP.o

driver_IJP.x: $(TEST_OBJS_IJP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP.x $(LDFLAGS) 

test_IJP: driver_IJP.x
	echo "Name = 'IJP';" > data/output_IJP.m
	echo " " >> data/output_IJP.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJP.x >> data/output_IJP.m

# ---------------------

TEST_OBJS_IPJ := driver.o Gemm_IPJ.o

driver_IPJ.x: $(TEST_OBJS_IPJ) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IPJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IPJ.x $(LDFLAGS) 

test_IPJ: driver_IPJ.x
	echo "Name = 'IPJ';" > data/output_IPJ.m
	echo " " >> data/output_IPJ.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IPJ.x >> data/output_IPJ.m

# ---------------------

TEST_OBJS_IP_Axpy := driver.o Gemm_IP_Axpy.o Axpy.o

driver_IP_Axpy.x: $(TEST_OBJS_IP_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IP_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_Axpy.x $(LDFLAGS) 

test_IP_Axpy: driver_IP_Axpy.x
	echo "Name = 'IP\_Axpy';" > data/output_IP_Axpy.m
	echo " " >> data/output_IP_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_Axpy.x >> data/output_IP_Axpy.m

# ---------------------

TEST_OBJS_IP_daxpy := driver.o Gemm_IP_daxpy.o

driver_IP_daxpy.x: $(TEST_OBJS_IP_daxpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IP_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_daxpy.x $(LDFLAGS) 

test_IP_daxpy: driver_IP_daxpy.x
	echo "Name = 'IP\_daxpy';" > data/output_IP_daxpy.m
	echo " " >> data/output_IP_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_daxpy.x >> data/output_IP_daxpy.m

# ---------------------

TEST_OBJS_IP_bli_daxpyv := driver.o Gemm_IP_bli_daxpyv.o

driver_IP_bli_daxpyv.x: $(TEST_OBJS_IP_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IP_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_bli_daxpyv.x $(LDFLAGS) 

test_IP_bli_daxpyv: driver_IP_bli_daxpyv.x
	echo "Name = 'IP\_bli\_daxpyv';" > data/output_IP_bli_daxpyv.m
	echo " " >> data/output_IP_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_bli_daxpyv.x >> data/output_IP_bli_daxpyv.m

# ---------------------

TEST_OBJS_JIP := driver.o Gemm_JIP.o

driver_JIP.x: $(TEST_OBJS_JIP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JIP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JIP.x $(LDFLAGS) 

test_JIP: driver_JIP.x
	echo "Name = 'JIP';" > data/output_JIP.m
	echo " " >> data/output_JIP.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JIP.x >> data/output_JIP.m

# ---------------------

TEST_OBJS_JPI := driver.o Gemm_JPI.o

driver_JPI.x: $(TEST_OBJS_JPI) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JPI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JPI.x $(LDFLAGS) 

test_JPI: driver_JPI.x
	echo "Name = 'JPI';" > data/output_JPI.m
	echo " " >> data/output_JPI.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JPI.x >> data/output_JPI.m

# ---------------------

TEST_OBJS_PIJ := driver.o Gemm_PIJ.o

driver_PIJ.x: $(TEST_OBJS_PIJ) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PIJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PIJ.x $(LDFLAGS) 

test_PIJ: driver_PIJ.x
	echo "Name = 'PIJ';" > data/output_PIJ.m
	echo " " >> data/output_PIJ.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PIJ.x >> data/output_PIJ.m

# ---------------------

TEST_OBJS_PI_Axpy := driver.o Gemm_PI_Axpy.o Axpy.o

driver_PI_Axpy.x: $(TEST_OBJS_PI_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PI_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_Axpy.x $(LDFLAGS) 

test_PI_Axpy: driver_PI_Axpy.x
	echo "Name = 'PI\_Axpy';" > data/output_PI_Axpy.m
	echo " " >> data/output_PI_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_Axpy.x >> data/output_PI_Axpy.m

# ---------------------

TEST_OBJS_PI_daxpy := driver.o Gemm_PI_daxpy.o

driver_PI_daxpy.x: $(TEST_OBJS_PI_daxpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PI_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_daxpy.x $(LDFLAGS) 

test_PI_daxpy: driver_PI_daxpy.x
	echo "Name = 'PI\_daxpy';" > data/output_PI_daxpy.m
	echo " " >> data/output_PI_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_daxpy.x >> data/output_PI_daxpy.m

# ---------------------

TEST_OBJS_PI_bli_daxpyv := driver.o Gemm_PI_bli_daxpyv.o

driver_PI_bli_daxpyv.x: $(TEST_OBJS_PI_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PI_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_bli_daxpyv.x $(LDFLAGS) 

test_PI_bli_daxpyv: driver_PI_bli_daxpyv.x
	echo "Name = 'PI\_bli\_daxpyv';" > data/output_PI_bli_daxpyv.m
	echo " " >> data/output_PI_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_bli_daxpyv.x >> data/output_PI_bli_daxpyv.m

# ---------------------

TEST_OBJS_IJ_DOTS := driver.o Gemm_IJ_Dots.o Dots.o

driver_IJ_Dots.x: $(TEST_OBJS_IJ_DOTS) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJ_DOTS) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_Dots.x $(LDFLAGS) 

test_IJ_Dots: driver_IJ_Dots.x
	echo "Name = 'IJ\_Dots';" > data/output_IJ_Dots.m
	echo " " >> data/output_IJ_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_Dots.x >> data/output_IJ_Dots.m

# ---------------------

TEST_OBJS_IJ_ddot := driver.o Gemm_IJ_ddot.o

driver_IJ_ddot.x: $(TEST_OBJS_IJ_ddot) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJ_ddot) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_ddot.x $(LDFLAGS) 

test_IJ_ddot: driver_IJ_ddot.x
	echo "Name = 'IJ\_ddot';" > data/output_IJ_ddot.m
	echo " " >> data/output_IJ_ddot.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_ddot.x >> data/output_IJ_ddot.m

# ---------------------

TEST_OBJS_IJ_bli_ddotxv := driver.o Gemm_IJ_bli_ddotxv.o

driver_IJ_bli_ddotxv.x: $(TEST_OBJS_IJ_bli_ddotxv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJ_bli_ddotxv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_bli_ddotxv.x $(LDFLAGS) 

test_IJ_bli_ddotxv: driver_IJ_bli_ddotxv.x
	echo "Name = 'IJ\_bli\_ddotxv';" > data/output_IJ_bli_ddotxv.m
	echo " " >> data/output_IJ_bli_ddotxv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_bli_ddotxv.x >> data/output_IJ_bli_ddotxv.m

# ---------------------

TEST_OBJS_JI_DOTS := driver.o Gemm_JI_Dots.o Dots.o

driver_JI_Dots.x: $(TEST_OBJS_JI_DOTS) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_DOTS) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_Dots.x $(LDFLAGS) 

test_JI_Dots: driver_JI_Dots.x
	echo "Name = 'JI\_Dots';" > data/output_JI_Dots.m
	echo " " >> data/output_JI_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_Dots.x >> data/output_JI_Dots.m

# ---------------------

TEST_OBJS_JI_ddot := driver.o Gemm_JI_ddot.o

driver_JI_ddot.x: $(TEST_OBJS_JI_ddot) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_ddot) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_ddot.x $(LDFLAGS) 

test_JI_ddot: driver_JI_ddot.x
	echo "Name = 'JI\_ddot';" > data/output_JI_ddot.m
	echo " " >> data/output_JI_ddot.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_ddot.x >> data/output_JI_ddot.m

# ---------------------

TEST_OBJS_JI_bli_ddotxv := driver.o Gemm_JI_bli_ddotxv.o

driver_JI_bli_ddotxv.x: $(TEST_OBJS_JI_bli_ddotxv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_bli_ddotxv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_bli_ddotxv.x $(LDFLAGS) 

test_JI_bli_ddotxv: driver_JI_bli_ddotxv.x
	echo "Name = 'JI\_bli\_ddotxv';" > data/output_JI_bli_ddotxv.m
	echo " " >> data/output_JI_bli_ddotxv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_bli_ddotxv.x >> data/output_JI_bli_ddotxv.m

# ---------------------

TEST_OBJS_JP_Axpy := driver.o Gemm_JP_Axpy.o Axpy.o

driver_JP_Axpy.x: $(TEST_OBJS_JP_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JP_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_Axpy.x $(LDFLAGS) 

test_JP_Axpy: driver_JP_Axpy.x
	echo "Name = 'JP\_Axpy';" > data/output_JP_Axpy.m
	echo " " >> data/output_JP_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_Axpy.x >> data/output_JP_Axpy.m

# ---------------------

TEST_OBJS_JP_daxpy := driver.o Gemm_JP_daxpy.o 

driver_JP_daxpy.x: $(TEST_OBJS_JP_daxpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JP_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_daxpy.x $(LDFLAGS) 

test_JP_daxpy: driver_JP_daxpy.x
	echo "Name = 'JP\_daxpy';" > data/output_JP_daxpy.m
	echo " " >> data/output_JP_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_daxpy.x >> data/output_JP_daxpy.m

# ---------------------

TEST_OBJS_JP_bli_daxpyv := driver.o Gemm_JP_bli_daxpyv.o 

driver_JP_bli_daxpyv.x: $(TEST_OBJS_JP_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JP_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_bli_daxpyv.x $(LDFLAGS) 

test_JP_bli_daxpyv: driver_JP_bli_daxpyv.x
	echo "Name = 'JP\_bli\_daxpyv';" > data/output_JP_bli_daxpyv.m
	echo " " >> data/output_JP_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_bli_daxpyv.x >> data/output_JP_bli_daxpyv.m

# ---------------------

TEST_OBJS_PJI := driver.o Gemm_PJI.o

driver_PJI.x: $(TEST_OBJS_PJI) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJI.x $(LDFLAGS) 

test_PJI: driver_PJI.x
	echo "Name = 'PJI';" > data/output_PJI.m
	echo " " >> data/output_PJI.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJI.x >> data/output_PJI.m

# ---------------------

TEST_OBJS_PJ_Axpy := driver.o Gemm_PJ_Axpy.o Axpy.o

driver_PJ_Axpy.x: $(TEST_OBJS_PJ_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJ_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_Axpy.x $(LDFLAGS) 

test_PJ_Axpy: driver_PJ_Axpy.x
	echo "Name = 'PJ\_Axpy';" > data/output_PJ_Axpy.m
	echo " " >> data/output_PJ_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_Axpy.x >> data/output_PJ_Axpy.m

# ---------------------

TEST_OBJS_PJ_daxpy := driver.o Gemm_PJ_daxpy.o

driver_PJ_daxpy.x: $(TEST_OBJS_PJ_daxpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJ_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_daxpy.x $(LDFLAGS) 

test_PJ_daxpy: driver_PJ_daxpy.x
	echo "Name = 'PJ\_daxpy';" > data/output_PJ_daxpy.m
	echo " " >> data/output_PJ_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_daxpy.x >> data/output_PJ_daxpy.m

# ---------------------

TEST_OBJS_PJ_bli_daxpyv := driver.o Gemm_PJ_bli_daxpyv.o

driver_PJ_bli_daxpyv.x: $(TEST_OBJS_PJ_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJ_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_bli_daxpyv.x $(LDFLAGS) 

test_PJ_bli_daxpyv: driver_PJ_bli_daxpyv.x
	echo "Name = 'PJ\_bli\_daxpyv';" > data/output_PJ_bli_daxpyv.m
	echo " " >> data/output_PJ_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_bli_daxpyv.x >> data/output_PJ_bli_daxpyv.m

# ---------------------

TEST_OBJS_J_Gemv_I_Dots := driver.o Gemm_J_Gemv_I_Dots.o Gemv_I_Dots.o Dots.o

driver_J_Gemv_I_Dots.x: $(TEST_OBJS_J_Gemv_I_Dots) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_J_Gemv_I_Dots) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_Gemv_I_Dots.x $(LDFLAGS) 

test_J_Gemv_I_Dots: driver_J_Gemv_I_Dots.x
	echo "Name = 'J\_Gemv\_I\_Dots';" > data/output_J_Gemv_I_Dots.m
	echo " " >> data/output_J_Gemv_I_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_Gemv_I_Dots.x >> data/output_J_Gemv_I_Dots.m

# ---------------------

TEST_OBJS_J_Gemv_J_Axpy := driver.o Gemm_J_Gemv_J_Axpy.o Gemv_J_Axpy.o Axpy.o

driver_J_Gemv_J_Axpy.x: $(TEST_OBJS_J_Gemv_J_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_J_Gemv_J_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_Gemv_J_Axpy.x $(LDFLAGS) 

test_J_Gemv_J_Axpy: driver_J_Gemv_J_Axpy.x
	echo "Name = 'J\_Gemv\_J\_Axpy';" > data/output_J_Gemv_J_Axpy.m
	echo " " >> data/output_J_Gemv_J_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_Gemv_J_Axpy.x >> data/output_J_Gemv_J_Axpy.m

# ---------------------

TEST_OBJS_J_dgemv := driver.o Gemm_J_dgemv.o

driver_J_dgemv.x: $(TEST_OBJS_J_dgemv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_J_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_dgemv.x $(LDFLAGS) 

test_J_dgemv: driver_J_dgemv.x
	echo "Name = 'J\_dgemv';" > data/output_J_dgemv.m
	echo " " >> data/output_J_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_dgemv.x >> data/output_J_dgemv.m

# ---------------------

TEST_OBJS_J_bli_dgemv := driver.o Gemm_J_bli_dgemv.o

driver_J_bli_dgemv.x: $(TEST_OBJS_J_bli_dgemv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_J_bli_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_bli_dgemv.x $(LDFLAGS) 

test_J_bli_dgemv: driver_J_bli_dgemv.x
	echo "Name = 'J\_bli\_dgemv';" > data/output_J_bli_dgemv.m
	echo " " >> data/output_J_bli_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_bli_dgemv.x >> data/output_J_bli_dgemv.m

# ---------------------

TEST_OBJS_P_Ger_J_Axpy := driver.o Gemm_P_Ger_J_Axpy.o Ger_J_Axpy.o Axpy.o

driver_P_Ger_J_Axpy.x: $(TEST_OBJS_P_Ger_J_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_P_Ger_J_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_Ger_J_Axpy.x $(LDFLAGS) 

test_P_Ger_J_Axpy: driver_P_Ger_J_Axpy.x
	echo "Name = 'P\_Ger\_J\_Axpy';" > data/output_P_Ger_J_Axpy.m
	echo " " >> data/output_P_Ger_J_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_Ger_J_Axpy.x >> data/output_P_Ger_J_Axpy.m

# ---------------------

TEST_OBJS_P_Ger_I_Axpy := driver.o Gemm_P_Ger_I_Axpy.o Ger_I_Axpy.o Axpy.o

driver_P_Ger_I_Axpy.x: $(TEST_OBJS_P_Ger_I_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_P_Ger_I_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_Ger_I_Axpy.x $(LDFLAGS) 

test_P_Ger_I_Axpy: driver_P_Ger_I_Axpy.x
	echo "Name = 'P\_Ger\_I\_Axpy';" > data/output_P_Ger_I_Axpy.m
	echo " " >> data/output_P_Ger_I_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_Ger_I_Axpy.x >> data/output_P_Ger_I_Axpy.m

# ---------------------

TEST_OBJS_P_dger := driver.o Gemm_P_dger.o 

driver_P_dger.x: $(TEST_OBJS_P_dger) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_P_dger) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_dger.x $(LDFLAGS) 

test_P_dger: driver_P_dger.x
	echo "Name = 'P\_dger';" > data/output_P_dger.m
	echo " " >> data/output_P_dger.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_dger.x >> data/output_P_dger.m

# ---------------------

TEST_OBJS_P_bli_dger := driver.o Gemm_P_bli_dger.o 

driver_P_bli_dger.x: $(TEST_OBJS_P_bli_dger) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_P_bli_dger) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_bli_dger.x $(LDFLAGS) 

test_P_bli_dger: driver_P_bli_dger.x
	echo "Name = 'P\_bli\_dger';" > data/output_P_bli_dger.m
	echo " " >> data/output_P_bli_dger.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_bli_dger.x >> data/output_P_bli_dger.m

# ---------------------

TEST_OBJS_I_dgemv := driver.o Gemm_I_dgemv.o

driver_I_dgemv.x: $(TEST_OBJS_I_dgemv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_I_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_I_dgemv.x $(LDFLAGS) 

test_I_dgemv: driver_I_dgemv.x
	echo "Name = 'I\_dgemv';" > data/output_I_dgemv.m
	echo " " >> data/output_I_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_I_dgemv.x >> data/output_I_dgemv.m

# ---------------------

TEST_OBJS_I_bli_dgemv := driver.o Gemm_I_bli_dgemv.o

driver_I_bli_dgemv.x: $(TEST_OBJS_I_bli_dgemv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_I_bli_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_I_bli_dgemv.x $(LDFLAGS) 

test_I_bli_dgemv: driver_I_bli_dgemv.x
	echo "Name = 'I\_bli\_dgemv';" > data/output_I_bli_dgemv.m
	echo " " >> data/output_I_bli_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_I_bli_dgemv.x >> data/output_I_bli_dgemv.m

# ---------------------

all:
	make test_IJP
	make test_JIP
	make test_IPJ
	make test_PIJ
	make test_JPI
	make test_PJI
	make test_IJ_Dots
	make test_JI_Dots
	make test_IP_Axpy
	make test_PI_Axpy
	make test_JP_Axpy
	make test_PJ_Axpy
	make test_IJ_ddot
	make test_JI_ddot
	make test_IP_daxpy
	make test_PI_daxpy
	make test_JP_daxpy
	make test_PJ_daxpy
	make test_IJ_bli_ddotxv
	make test_JI_bli_ddotxv
	make test_IP_bli_daxpyv
	make test_PI_bli_daxpyv
	make test_JP_bli_daxpyv
	make test_PJ_bli_daxpyv
	make test_J_Gemv_I_Dots
	make test_J_Gemv_J_Axpy
	make test_J_dgemv
	make test_J_bli_dgemv
	make test_P_Ger_J_Axpy
	make test_P_Ger_I_Axpy
	make test_P_dger
	make test_P_bli_dger
	make test_I_dgemv
	make test_I_bli_dgemv

clean:
	rm -f *.o *~ core *.x *.pdf
