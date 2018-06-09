# HOME      := /Users/rvdg
# Make sure you have BLIS installed in your home directory.
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O3 -I$(BLAS_INC) -m64 -msse3 -std=c99 -march=native -D_POSIX_C_SOURCE=200112L -fopenmp
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.  The best time is reported.
NFIRST     := 48#     smallest size to be timed.
NLAST_SMALL:= 1500#    largest size to be timed for slow implementations.
NLAST      := 1500#   largest size to be timed for fast implementations.
NINC       := 48#     increment between sizes.

LDFLAGS    := -fopenmp -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# ---------------------

HelloWorld:    HelloWorld.o
	$(LINKER) $(CFLAGS) $(LDFLAGS) HelloWorld.o -o HelloWorld.x


HelloWorld2:    HelloWorld2.o
	$(LINKER) $(CFLAGS) $(LDFLAGS) HelloWorld2.o -o HelloWorld2.x

# ---------------------

TEST_OBJS_Parallel_Loop1_12x4Kernel := driver.o Gemm_Parallel_Loop1_12x4Kernel.o

driver_Parallel_Loop1_12x4Kernel.x: $(TEST_OBJS_Parallel_Loop1_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_Parallel_Loop1_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Parallel_Loop1_12x4Kernel.x $(LDFLAGS) 

test_Parallel_Loop1_12x4Kernel: driver_Parallel_Loop1_12x4Kernel.x
	echo "Name = 'Parallel\_Loop1\_12x4Kernel';" > data/output_Parallel_Loop1_12x4Kernel.m
	echo " " >> data/output_Parallel_Loop1_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Parallel_Loop1_12x4Kernel.x >> data/output_Parallel_Loop1_12x4Kernel.m

# ---------------------

TEST_OBJS_Parallel_Loop2_12x4Kernel := driver.o Gemm_Parallel_Loop2_12x4Kernel.o

driver_Parallel_Loop2_12x4Kernel.x: $(TEST_OBJS_Parallel_Loop2_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_Parallel_Loop2_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Parallel_Loop2_12x4Kernel.x $(LDFLAGS) 

test_Parallel_Loop2_12x4Kernel: driver_Parallel_Loop2_12x4Kernel.x
	echo "Name = 'Parallel\_Loop2\_12x4Kernel';" > data/output_Parallel_Loop2_12x4Kernel.m
	echo " " >> data/output_Parallel_Loop2_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Parallel_Loop2_12x4Kernel.x >> data/output_Parallel_Loop2_12x4Kernel.m

# ---------------------

TEST_OBJS_Parallel_Loop3_12x4Kernel := driver.o Gemm_Parallel_Loop3_12x4Kernel.o

driver_Parallel_Loop3_12x4Kernel.x: $(TEST_OBJS_Parallel_Loop3_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_Parallel_Loop3_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Parallel_Loop3_12x4Kernel.x $(LDFLAGS) 

test_Parallel_Loop3_12x4Kernel: driver_Parallel_Loop3_12x4Kernel.x
	echo "Name = 'Parallel\_Loop3\_12x4Kernel';" > data/output_Parallel_Loop3_12x4Kernel.m
	echo " " >> data/output_Parallel_Loop3_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Parallel_Loop3_12x4Kernel.x >> data/output_Parallel_Loop3_12x4Kernel.m

# ---------------------

clean:
	rm -f *.o *~ core *.x *.pdf
