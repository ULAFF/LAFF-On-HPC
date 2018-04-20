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
NLAST      := 1500#   largest size to be timed for fast implementations.
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

TEST_OBJS_IJP_MRxNRxKR := driver.o Gemm_IJP_MRxNRxKR.o

driver_IJP_MRxNRxKR.x: $(TEST_OBJS_IJP_MRxNRxKR) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP_MRxNRxKR) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP_MRxNRxKR.x $(LDFLAGS) 

test_IJP_MRxNRxKR: driver_IJP_MRxNRxKR.x
	echo "Name = 'IJP\_MRxNRxKR';" > data/output_IJP_MRxNRxKR.m
	echo " " >> data/output_IJP_MRxNRxKR.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJP_MRxNRxKR.x >> data/output_IJP_MRxNRxKR.m

# ---------------------

TEST_OBJS_IJ_4x4Kernel := driver.o Gemm_IJ_4x4Kernel.o

driver_IJ_4x4Kernel.x: $(TEST_OBJS_IJ_4x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJ_4x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_4x4Kernel.x $(LDFLAGS) 

test_IJ_4x4Kernel: driver_IJ_4x4Kernel.x
	echo "Name = 'IJ\_4x4Kernel';" > data/output_IJ_4x4Kernel.m
	echo " " >> data/output_IJ_4x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_4x4Kernel.x >> data/output_IJ_4x4Kernel.m

# ---------------------

TEST_OBJS_JI_4x4Kernel := driver.o Gemm_JI_4x4Kernel.o

driver_JI_4x4Kernel.x: $(TEST_OBJS_JI_4x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_4x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_4x4Kernel.x $(LDFLAGS) 

test_JI_4x4Kernel: driver_JI_4x4Kernel.x
	echo "Name = 'JI\_4x4Kernel';" > data/output_JI_4x4Kernel.m
	echo " " >> data/output_JI_4x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_4x4Kernel.x >> data/output_JI_4x4Kernel.m

# ---------------------

TEST_OBJS_JI_8x4Kernel := driver.o Gemm_JI_8x4Kernel.o

driver_JI_8x4Kernel.x: $(TEST_OBJS_JI_8x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_8x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_8x4Kernel.x $(LDFLAGS) 

test_JI_8x4Kernel: driver_JI_8x4Kernel.x
	echo "Name = 'JI\_8x4Kernel';" > data/output_JI_8x4Kernel.m
	echo " " >> data/output_JI_8x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_8x4Kernel.x >> data/output_JI_8x4Kernel.m

# ---------------------

TEST_OBJS_JI_12x4Kernel := driver.o Gemm_JI_12x4Kernel.o

driver_JI_12x4Kernel.x: $(TEST_OBJS_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_12x4Kernel.x $(LDFLAGS) 

test_JI_12x4Kernel: driver_JI_12x4Kernel.x
	echo "Name = 'JI\_12x4Kernel';" > data/output_JI_12x4Kernel.m
	echo " " >> data/output_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_12x4Kernel.x >> data/output_JI_12x4Kernel.m

# ---------------------

TEST_OBJS_JI_8x6Kernel := driver.o Gemm_JI_8x6Kernel.o

driver_JI_8x6Kernel.x: $(TEST_OBJS_JI_8x6Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_8x6Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_8x6Kernel.x $(LDFLAGS) 

test_JI_8x6Kernel: driver_JI_8x6Kernel.x
	echo "Name = 'JI\_8x6Kernel';" > data/output_JI_8x6Kernel.m
	echo " " >> data/output_JI_8x6Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_8x6Kernel.x >> data/output_JI_8x6Kernel.m

# ---------------------

TEST_OBJS_JI_4x8Kernel := driver.o Gemm_JI_4x8Kernel.o

driver_JI_4x8Kernel.x: $(TEST_OBJS_JI_4x8Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_4x8Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_4x8Kernel.x $(LDFLAGS) 

test_JI_4x8Kernel: driver_JI_4x8Kernel.x
	echo "Name = 'JI\_4x8Kernel';" > data/output_JI_4x8Kernel.m
	echo " " >> data/output_JI_4x8Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_4x8Kernel.x >> data/output_JI_4x8Kernel.m

# ---------------------

TEST_OBJS_PI_JI_MCxKC := driver_MCxKC.o Gemm_PI_JI_12x4Kernel.o

driver_PI_JI_MCxKC.x: $(TEST_OBJS_PI_JI_MCxKC) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PI_JI_MCxKC) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_JI_MCxKC.x $(LDFLAGS) 

test_PI_JI_MCxKC: driver_PI_JI_MCxKC.x
	echo "Name = 'PI_JI\_MCxKC';" > data/output_PI_JI_MCxKC.m
	echo " " >> data/output_PI_JI_MCxKC.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_JI_MCxKC.x >> data/output_PI_JI_MCxKC.m

# ---------------------

TEST_OBJS_IJP_JI_4x4Kernel := driver.o Gemm_IJP_JI_4x4Kernel.o

driver_IJP_JI_4x4Kernel.x: $(TEST_OBJS_IJP_JI_4x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP_JI_4x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP_JI_4x4Kernel.x $(LDFLAGS) 

test_IJP_JI_4x4Kernel: driver_IJP_JI_4x4Kernel.x
	echo "Name = 'IJP\_JI\_4x4Kernel';" > data/output_IJP_JI_4x4Kernel.m
	echo " " >> data/output_IJP_JI_4x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJP_JI_4x4Kernel.x >> data/output_IJP_JI_4x4Kernel.m

# ---------------------

TEST_OBJS_IJP_JI_12x4Kernel := driver.o Gemm_IJP_JI_12x4Kernel.o

driver_IJP_JI_12x4Kernel.x: $(TEST_OBJS_IJP_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP_JI_12x4Kernel.x $(LDFLAGS) 

test_IJP_JI_12x4Kernel: driver_IJP_JI_12x4Kernel.x
	echo "Name = 'IJP\_JI\_12x4Kernel';" > data/output_IJP_JI_12x4Kernel.m
	echo " " >> data/output_IJP_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJP_JI_12x4Kernel.x >> data/output_IJP_JI_12x4Kernel.m

# ---------------------

TEST_OBJS_JPI_JI_12x4Kernel := driver.o Gemm_JPI_JI_12x4Kernel.o

driver_JPI_JI_12x4Kernel.x: $(TEST_OBJS_JPI_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JPI_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JPI_JI_12x4Kernel.x $(LDFLAGS) 

test_JPI_JI_12x4Kernel: driver_JPI_JI_12x4Kernel.x
	echo "Name = 'JPI\_JI\_12x4Kernel';" > data/output_JPI_JI_12x4Kernel.m
	echo " " >> data/output_JPI_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JPI_JI_12x4Kernel.x >> data/output_JPI_JI_12x4Kernel.m

# ---------------------

TEST_OBJS_PIJ_JI_12x4Kernel := driver.o Gemm_PIJ_JI_12x4Kernel.o

driver_PIJ_JI_12x4Kernel.x: $(TEST_OBJS_PIJ_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PIJ_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PIJ_JI_12x4Kernel.x $(LDFLAGS) 

test_PIJ_JI_12x4Kernel: driver_PIJ_JI_12x4Kernel.x
	echo "Name = 'PIJ\_JI\_12x4Kernel';" > data/output_PIJ_JI_12x4Kernel.m
	echo " " >> data/output_PIJ_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PIJ_JI_12x4Kernel.x >> data/output_PIJ_JI_12x4Kernel.m

# ---------------------

TEST_OBJS_PI_JI_12x4Kernel := driver.o Gemm_PI_JI_12x4Kernel.o

driver_PI_JI_12x4Kernel.x: $(TEST_OBJS_PI_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PI_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_JI_12x4Kernel.x $(LDFLAGS) 

test_PI_JI_12x4Kernel: driver_PI_JI_12x4Kernel.x
	echo "Name = 'PI\_JI\_12x4Kernel';" > data/output_PI_JI_12x4Kernel.m
	echo " " >> data/output_PI_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_JI_12x4Kernel.x >> data/output_PI_JI_12x4Kernel.m

# ---------------------

TEST_OBJS_Five_Loops_12x4Kernel := driver.o Gemm_Five_Loops_12x4Kernel.o

driver_Five_Loops_12x4Kernel.x: $(TEST_OBJS_Five_Loops_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_Five_Loops_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_12x4Kernel.x $(LDFLAGS) 

test_Five_Loops_12x4Kernel: driver_Five_Loops_12x4Kernel.x
	echo "Name = 'Five\_Loops\_12x4Kernel';" > data/output_Five_Loops_12x4Kernel.m
	echo " " >> data/output_Five_Loops_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_Five_Loops_12x4Kernel.x >> data/output_Five_Loops_12x4Kernel.m

# ---------------------

TEST_OBJS_Five_Loops_Pack_4x4Kernel := driver.o Gemm_Five_Loops_Pack_4x4Kernel.o

driver_Five_Loops_Pack_4x4Kernel.x: $(TEST_OBJS_Five_Loops_Pack_4x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_Five_Loops_Pack_4x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_Pack_4x4Kernel.x $(LDFLAGS) 

test_Five_Loops_Pack_4x4Kernel: driver_Five_Loops_Pack_4x4Kernel.x
	echo "Name = 'Five\_Loops\_Pack\_4x4Kernel';" > data/output_Five_Loops_Pack_4x4Kernel.m
	echo " " >> data/output_Five_Loops_Pack_4x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_Five_Loops_Pack_4x4Kernel.x >> data/output_Five_Loops_Pack_4x4Kernel.m

# ---------------------

TEST_OBJS_Five_Loops_Pack_12x4Kernel := driver.o Gemm_Five_Loops_Pack_12x4Kernel.o

driver_Five_Loops_Pack_12x4Kernel.x: $(TEST_OBJS_Five_Loops_Pack_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_Five_Loops_Pack_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_Pack_12x4Kernel.x $(LDFLAGS) 

test_Five_Loops_Pack_12x4Kernel: driver_Five_Loops_Pack_12x4Kernel.x
	echo "Name = 'Five\_Loops\_Pack\_12x4Kernel';" > data/output_Five_Loops_Pack_12x4Kernel.m
	echo " " >> data/output_Five_Loops_Pack_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_Five_Loops_Pack_12x4Kernel.x >> data/output_Five_Loops_Pack_12x4Kernel.m

# ---------------------

clean:
	rm -f *.o *~ core *.x *.pdf
