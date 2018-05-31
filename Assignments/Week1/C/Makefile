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
NLAST_SMALL:= 1000#    largest size to be timed for slow implementations.
NLAST      := 2000#   largest size to be timed for fast implementations.
NINC       := 48#     increment between sizes.

LDFLAGS    := -fopenmp -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# ---------------------

TEST_OBJS_IJP := driver.o Gemm_IJP.o

driver_IJP.x: $(TEST_OBJS_IJP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP.x $(LDFLAGS) 

test_IJP: driver_IJP.x
	echo "your_version = 'IJP';" > data/output_IJP.m
	echo " " >> data/output_IJP.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJP.x >> data/output_IJP.m

# ---------------------

TEST_OBJS_IPJ := driver.o Gemm_IPJ.o

driver_IPJ.x: $(TEST_OBJS_IPJ) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IPJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IPJ.x $(LDFLAGS) 

test_IPJ: driver_IPJ.x
	echo "your_version = 'IPJ';" > data/output_IPJ.m
	echo " " >> data/output_IPJ.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IPJ.x >> data/output_IPJ.m

# ---------------------

TEST_OBJS_IP_Axpy := driver.o Gemm_IP_Axpy.o Axpy.o

driver_IP_Axpy.x: $(TEST_OBJS_IP_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IP_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_Axpy.x $(LDFLAGS) 

test_IP_Axpy: driver_IP_Axpy.x
	echo "your_version = 'IP\_Axpy';" > data/output_IP_Axpy.m
	echo " " >> data/output_IP_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_Axpy.x >> data/output_IP_Axpy.m

# ---------------------

TEST_OBJS_IP_daxpy := driver.o Gemm_IP_daxpy.o

driver_IP_daxpy.x: $(TEST_OBJS_IP_daxpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IP_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_daxpy.x $(LDFLAGS) 

test_IP_daxpy: driver_IP_daxpy.x
	echo "your_version = 'IP\_daxpy';" > data/output_IP_daxpy.m
	echo " " >> data/output_IP_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_daxpy.x >> data/output_IP_daxpy.m

# ---------------------

TEST_OBJS_IP_bli_daxpyv := driver.o Gemm_IP_bli_daxpyv.o

driver_IP_bli_daxpyv.x: $(TEST_OBJS_IP_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IP_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IP_bli_daxpyv.x $(LDFLAGS) 

test_IP_bli_daxpyv: driver_IP_bli_daxpyv.x
	echo "your_version = 'IP\_bli_daxpyv';" > data/output_IP_bli_daxpyv.m
	echo " " >> data/output_IP_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IP_bli_daxpyv.x >> data/output_IP_bli_daxpyv.m

# ---------------------

TEST_OBJS_JIP := driver.o Gemm_JIP.o

driver_JIP.x: $(TEST_OBJS_JIP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JIP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JIP.x $(LDFLAGS) 

test_JIP: driver_JIP.x
	echo "your_version = 'JIP';" > data/output_JIP.m
	echo " " >> data/output_JIP.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JIP.x >> data/output_JIP.m

# ---------------------

TEST_OBJS_JPI := driver.o Gemm_JPI.o

driver_JPI.x: $(TEST_OBJS_JPI) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JPI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JPI.x $(LDFLAGS) 

test_JPI: driver_JPI.x
	echo "your_version = 'JPI';" > data/output_JPI.m
	echo " " >> data/output_JPI.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JPI.x >> data/output_JPI.m

# ---------------------

TEST_OBJS_PIJ := driver.o Gemm_PIJ.o

driver_PIJ.x: $(TEST_OBJS_PIJ) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PIJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PIJ.x $(LDFLAGS) 

test_PIJ: driver_PIJ.x
	echo "your_version = 'PIJ';" > data/output_PIJ.m
	echo " " >> data/output_PIJ.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PIJ.x >> data/output_PIJ.m

# ---------------------

TEST_OBJS_PI_Axpy := driver.o Gemm_PI_Axpy.o Axpy.o

driver_PI_Axpy.x: $(TEST_OBJS_PI_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PI_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_Axpy.x $(LDFLAGS) 

test_PI_Axpy: driver_PI_Axpy.x
	echo "your_version = 'PI_Axpy';" > data/output_PI_Axpy.m
	echo " " >> data/output_PI_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_Axpy.x >> data/output_PI_Axpy.m

# ---------------------

TEST_OBJS_PI_daxpy := driver.o Gemm_PI_daxpy.o

driver_PI_daxpy.x: $(TEST_OBJS_PI_daxpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PI_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_daxpy.x $(LDFLAGS) 

test_PI_daxpy: driver_PI_daxpy.x
	echo "your_version = 'PI_daxpy';" > data/output_PI_daxpy.m
	echo " " >> data/output_PI_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_daxpy.x >> data/output_PI_daxpy.m

# ---------------------

TEST_OBJS_PI_bli_daxpyv := driver.o Gemm_PI_bli_daxpyv.o

driver_PI_bli_daxpyv.x: $(TEST_OBJS_PI_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PI_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_bli_daxpyv.x $(LDFLAGS) 

test_PI_bli_daxpyv: driver_PI_bli_daxpyv.x
	echo "your_version = 'PI_bli_daxpyv';" > data/output_PI_bli_daxpyv.m
	echo " " >> data/output_PI_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_bli_daxpyv.x >> data/output_PI_bli_daxpyv.m

# ---------------------

TEST_OBJS_IJ_DOTS := driver.o Gemm_IJ_Dots.o Dots.o

driver_IJ_Dots.x: $(TEST_OBJS_IJ_DOTS) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJ_DOTS) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_Dots.x $(LDFLAGS) 

test_IJ_Dots: driver_IJ_Dots.x
	echo "your_version = 'IJ\_Dots';" > data/output_IJ_Dots.m
	echo " " >> data/output_IJ_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_Dots.x >> data/output_IJ_Dots.m

# ---------------------

TEST_OBJS_IJ_ddot := driver.o Gemm_IJ_ddot.o

driver_IJ_ddot.x: $(TEST_OBJS_IJ_ddot) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJ_ddot) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_ddot.x $(LDFLAGS) 

test_IJ_ddot: driver_IJ_ddot.x
	echo "your_version = 'IJ\_ddot';" > data/output_IJ_ddot.m
	echo " " >> data/output_IJ_ddot.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_ddot.x >> data/output_IJ_ddot.m

# ---------------------

TEST_OBJS_IJ_bli_ddotxv := driver.o Gemm_IJ_bli_ddotxv.o

driver_IJ_bli_ddotxv.x: $(TEST_OBJS_IJ_bli_ddotxv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJ_bli_ddotxv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_bli_ddotxv.x $(LDFLAGS) 

test_IJ_bli_ddotxv: driver_IJ_bli_ddotxv.x
	echo "your_version = 'IJ\_bli\_ddotxv';" > data/output_IJ_bli_ddotxv.m
	echo " " >> data/output_IJ_bli_ddotxv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_bli_ddotxv.x >> data/output_IJ_bli_ddotxv.m

# ---------------------

TEST_OBJS_JI_DOTS := driver.o Gemm_JI_Dots.o Dots.o

driver_JI_Dots.x: $(TEST_OBJS_JI_DOTS) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_DOTS) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_Dots.x $(LDFLAGS) 

test_JI_Dots: driver_JI_Dots.x
	echo "your_version = 'JI\_Dots';" > data/output_JI_Dots.m
	echo " " >> data/output_JI_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_Dots.x >> data/output_JI_Dots.m

# ---------------------

TEST_OBJS_JI_ddot := driver.o Gemm_JI_ddot.o

driver_JI_ddot.x: $(TEST_OBJS_JI_ddot) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_ddot) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_ddot.x $(LDFLAGS) 

test_JI_ddot: driver_JI_ddot.x
	echo "your_version = 'JI\_ddot';" > data/output_JI_ddot.m
	echo " " >> data/output_JI_ddot.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_ddot.x >> data/output_JI_ddot.m

# ---------------------

TEST_OBJS_JI_bli_ddotxv := driver.o Gemm_JI_bli_ddotxv.o

driver_JI_bli_ddotxv.x: $(TEST_OBJS_JI_bli_ddotxv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_bli_ddotxv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_bli_ddotxv.x $(LDFLAGS) 

test_JI_bli_ddotxv: driver_JI_bli_ddotxv.x
	echo "your_version = 'JI\_bli\_ddotxv';" > data/output_JI_bli_ddotxv.m
	echo " " >> data/output_JI_bli_ddotxv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_bli_ddotxv.x >> data/output_JI_bli_ddotxv.m

# ---------------------

TEST_OBJS_JP_Axpy := driver.o Gemm_JP_Axpy.o Axpy.o

driver_JP_Axpy.x: $(TEST_OBJS_JP_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JP_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_Axpy.x $(LDFLAGS) 

test_JP_Axpy: driver_JP_Axpy.x
	echo "your_version = 'JP\_Axpy';" > data/output_JP_Axpy.m
	echo " " >> data/output_JP_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_Axpy.x >> data/output_JP_Axpy.m

# ---------------------

TEST_OBJS_JP_daxpy := driver.o Gemm_JP_daxpy.o 

driver_JP_daxpy.x: $(TEST_OBJS_JP_daxpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JP_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_daxpy.x $(LDFLAGS) 

test_JP_daxpy: driver_JP_daxpy.x
	echo "your_version = 'JP\_daxpy';" > data/output_JP_daxpy.m
	echo " " >> data/output_JP_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_daxpy.x >> data/output_JP_daxpy.m

# ---------------------

TEST_OBJS_JP_bli_daxpyv := driver.o Gemm_JP_bli_daxpyv.o 

driver_JP_bli_daxpyv.x: $(TEST_OBJS_JP_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JP_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JP_bli_daxpyv.x $(LDFLAGS) 

test_JP_bli_daxpyv: driver_JP_bli_daxpyv.x
	echo "your_version = 'JP\_bli_daxpyv';" > data/output_JP_bli_daxpyv.m
	echo " " >> data/output_JP_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JP_bli_daxpyv.x >> data/output_JP_bli_daxpyv.m

# ---------------------

TEST_OBJS_PJI := driver.o Gemm_PJI.o

driver_PJI.x: $(TEST_OBJS_PJI) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJI.x $(LDFLAGS) 

test_PJI: driver_PJI.x
	echo "your_version = 'PJI';" > data/output_PJI.m
	echo " " >> data/output_PJI.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJI.x >> data/output_PJI.m

# ---------------------

TEST_OBJS_PJ_Axpy := driver.o Gemm_PJ_Axpy.o Axpy.o

driver_PJ_Axpy.x: $(TEST_OBJS_PJ_Axpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJ_Axpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_Axpy.x $(LDFLAGS) 

test_PJ_Axpy: driver_PJ_Axpy.x
	echo "your_version = 'PJ\_Axpy';" > data/output_PJ_Axpy.m
	echo " " >> data/output_PJ_Axpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_Axpy.x >> data/output_PJ_Axpy.m

# ---------------------

TEST_OBJS_PJ_daxpy := driver.o Gemm_PJ_daxpy.o

driver_PJ_daxpy.x: $(TEST_OBJS_PJ_daxpy) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJ_daxpy) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_daxpy.x $(LDFLAGS) 

test_PJ_daxpy: driver_PJ_daxpy.x
	echo "your_version = 'PJ\_daxpy';" > data/output_PJ_daxpy.m
	echo " " >> data/output_PJ_daxpy.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_daxpy.x >> data/output_PJ_daxpy.m

# ---------------------

TEST_OBJS_PJ_bli_daxpyv := driver.o Gemm_PJ_bli_daxpyv.o

driver_PJ_bli_daxpyv.x: $(TEST_OBJS_PJ_bli_daxpyv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJ_bli_daxpyv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJ_bli_daxpyv.x $(LDFLAGS) 

test_PJ_bli_daxpyv: driver_PJ_bli_daxpyv.x
	echo "your_version = 'PJ\_bli_daxpyv';" > data/output_PJ_bli_daxpyv.m
	echo " " >> data/output_PJ_bli_daxpyv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJ_bli_daxpyv.x >> data/output_PJ_bli_daxpyv.m

# ---------------------

TEST_OBJS_P_GER := driver.o Gemm_P_Ger.o

driver_P_Ger.x: $(TEST_OBJS_P_GER) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_P_GER) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_Ger.x $(LDFLAGS) 

test_P_Ger: driver_P_Ger.x
	echo "your_version = 'P\_Ger';" > data/output_P_Ger.m
	echo " " >> data/output_P_Ger.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_Ger.x >> data/output_P_Ger.m

# ---------------------

TEST_OBJS_J_GEMV := driver.o Gemm_J_Gemv.o

driver_J_Gemv.x: $(TEST_OBJS_J_GEMV) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_J_GEMV) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_Gemv.x $(LDFLAGS) 

test_J_Gemv: driver_J_Gemv.x
	echo "your_version = 'J\_Gemv';" > data/output_J_Gemv.m
	echo " " >> data/output_J_Gemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_Gemv.x >> data/output_J_Gemv.m

# ---------------------

TEST_OBJS_J_dgemv := driver.o Gemm_J_dgemv.o

driver_J_dgemv.x: $(TEST_OBJS_J_dgemv) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_J_dgemv) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_dgemv.x $(LDFLAGS) 

test_J_dgemv: driver_J_dgemv.x
	echo "your_version = 'J\_dgemv';" > data/output_J_dgemv.m
	echo " " >> data/output_J_dgemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_dgemv.x >> data/output_J_dgemv.m

# ---------------------

clean:
	rm -f *.o *~ core *.x *.pdf
