# HOME      := /Users/rvdg
# Make sure you have BLIS installed in your home directory.
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O3 -I$(BLAS_INC) -m64 -msse3 -std=c99 -march=native -D_POSIX_C_SOURCE=200112L # -fopenmp
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.  The best time is reported.
NFIRST     := 48#     smallest size to be timed.
NLAST_SMALL:= 1000#    largest size to be timed for slow implementations.
NLAST      := 1000#   largest size to be timed for fast implementations.
NINC       := 48#     increment between sizes.

LDFLAGS    := -fopenmp -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# ---------------------

OBJS_Five_Loops_4x4Kernel := driver.o Gemm_Five_Loops_4x4Kernel.o

driver_Five_Loops_4x4Kernel.x: $(OBJS_Five_Loops_4x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Five_Loops_4x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_4x4Kernel.x $(LDFLAGS) 

Five_Loops_4x4Kernel: driver_Five_Loops_4x4Kernel.x
	echo "Name = 'Five\_Loops\_4x4Kernel';" > data/output_Five_Loops_4x4Kernel.m
	echo " " >> data/output_Five_Loops_4x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Five_Loops_4x4Kernel.x >> data/output_Five_Loops_4x4Kernel.m
	more data/output_Five_Loops_4x4Kernel.m

# ---------------------

OBJS_Five_Loops_12x4Kernel := driver.o Gemm_Five_Loops_12x4Kernel.o

driver_Five_Loops_12x4Kernel.x: $(OBJS_Five_Loops_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Five_Loops_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_12x4Kernel.x $(LDFLAGS) 

Five_Loops_12x4Kernel: driver_Five_Loops_12x4Kernel.x
	echo "Name = 'Five\_Loops\_12x4Kernel';" > data/output_Five_Loops_12x4Kernel.m
	echo " " >> data/output_Five_Loops_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Five_Loops_12x4Kernel.x >> data/output_Five_Loops_12x4Kernel.m
	more data/output_Five_Loops_12x4Kernel.m

# ---------------------

OBJS_Five_Loops_Pack_4x4Kernel := driver.o Gemm_Five_Loops_Pack_4x4Kernel.o

driver_Five_Loops_Pack_4x4Kernel.x: $(OBJS_Five_Loops_Pack_4x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Five_Loops_Pack_4x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_Pack_4x4Kernel.x $(LDFLAGS) 

Five_Loops_Pack_4x4Kernel: driver_Five_Loops_Pack_4x4Kernel.x
	echo "Name = 'Five\_Loops\_Pack\_4x4Kernel';" > data/output_Five_Loops_Pack_4x4Kernel.m
	echo " " >> data/output_Five_Loops_Pack_4x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Five_Loops_Pack_4x4Kernel.x >> data/output_Five_Loops_Pack_4x4Kernel.m
	more data/output_Five_Loops_Pack_4x4Kernel.m

# ---------------------

OBJS_Five_Loops_Pack_12x4Kernel := driver.o Gemm_Five_Loops_Pack_12x4Kernel.o

driver_Five_Loops_Pack_12x4Kernel.x: $(OBJS_Five_Loops_Pack_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Five_Loops_Pack_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_Pack_12x4Kernel.x $(LDFLAGS) 

Five_Loops_Pack_12x4Kernel: driver_Five_Loops_Pack_12x4Kernel.x
	echo "Name = 'Five\_Loops\_Pack\_12x4Kernel';" > data/output_Five_Loops_Pack_12x4Kernel.m
	echo " " >> data/output_Five_Loops_Pack_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Five_Loops_Pack_12x4Kernel.x >> data/output_Five_Loops_Pack_12x4Kernel.m
	more data/output_Five_Loops_Pack_12x4Kernel.m

# ---------------------

OBJS_Five_Loops_Pack_12x4Kernel_OPT := driver.o Gemm_Five_Loops_Pack_12x4Kernel_OPT.o

driver_Five_Loops_Pack_12x4Kernel_OPT.x: $(OBJS_Five_Loops_Pack_12x4Kernel_OPT) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Five_Loops_Pack_12x4Kernel_OPT) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_Pack_12x4Kernel_OPT.x $(LDFLAGS) 

Five_Loops_Pack_12x4Kernel_OPT: driver_Five_Loops_Pack_12x4Kernel_OPT.x
	echo "Name = 'Five\_Loops\_Pack\_12x4Kernel\_OPT';" > data/output_Five_Loops_Pack_12x4Kernel_OPT.m
	echo " " >> data/output_Five_Loops_Pack_12x4Kernel_OPT.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Five_Loops_Pack_12x4Kernel_OPT.x >> data/output_Five_Loops_Pack_12x4Kernel_OPT.m
	more data/output_Five_Loops_Pack_12x4Kernel_OPT.m

# ---------------------

clean:
	rm -f *.o *~ core *.x *.pdf
