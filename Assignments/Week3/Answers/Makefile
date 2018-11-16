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
NLAST_SMALL:= 512#    largest size to be timed for slow implementations.
NLAST      := 1000#   largest size to be timed for fast implementations.
NINC       := 48#     increment between sizes.

LDFLAGS    := -fopenmp -lpthread -m64 -lm  

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# ---------------------

OBJS_IJP_JI_4x4Kernel := driver.o Gemm_IJP_JI_4x4Kernel.o

driver_IJP_JI_4x4Kernel.x: $(OBJS_IJP_JI_4x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IJP_JI_4x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP_JI_4x4Kernel.x $(LDFLAGS) 

IJP_JI_4x4Kernel: driver_IJP_JI_4x4Kernel.x
	echo "Name = 'IJP\_JI\_4x4Kernel';" > data/output_IJP_JI_4x4Kernel.m
	echo " " >> data/output_IJP_JI_4x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_IJP_JI_4x4Kernel.x >> data/output_IJP_JI_4x4Kernel.m
	more data/output_IJP_JI_4x4Kernel.m

# ---------------------

OBJS_IJP_JI_12x4Kernel := driver.o Gemm_IJP_JI_12x4Kernel.o

driver_IJP_JI_12x4Kernel.x: $(OBJS_IJP_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IJP_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP_JI_12x4Kernel.x $(LDFLAGS) 

IJP_JI_12x4Kernel: driver_IJP_JI_12x4Kernel.x
	echo "Name = 'IJP\_JI\_12x4Kernel';" > data/output_IJP_JI_12x4Kernel.m
	echo " " >> data/output_IJP_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_IJP_JI_12x4Kernel.x >> data/output_IJP_JI_12x4Kernel.m
	more data/output_IJP_JI_12x4Kernel.m

# ---------------------

OBJS_IPJ_JI_12x4Kernel := driver.o Gemm_IPJ_JI_12x4Kernel.o

driver_IPJ_JI_12x4Kernel.x: $(OBJS_IPJ_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_IPJ_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IPJ_JI_12x4Kernel.x $(LDFLAGS) 

IPJ_JI_12x4Kernel: driver_IPJ_JI_12x4Kernel.x
	echo "Name = 'IPJ\_JI\_12x4Kernel';" > data/output_IPJ_JI_12x4Kernel.m
	echo " " >> data/output_IPJ_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_IPJ_JI_12x4Kernel.x >> data/output_IPJ_JI_12x4Kernel.m
	more data/output_IPJ_JI_12x4Kernel.m

# ---------------------

OBJS_JIP_JI_12x4Kernel := driver.o Gemm_JIP_JI_12x4Kernel.o

driver_JIP_JI_12x4Kernel.x: $(OBJS_JIP_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JIP_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JIP_JI_12x4Kernel.x $(LDFLAGS) 

JIP_JI_12x4Kernel: driver_JIP_JI_12x4Kernel.x
	echo "Name = 'JIP\_JI\_12x4Kernel';" > data/output_JIP_JI_12x4Kernel.m
	echo " " >> data/output_JIP_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_JIP_JI_12x4Kernel.x >> data/output_JIP_JI_12x4Kernel.m
	more data/output_JIP_JI_12x4Kernel.m

# ---------------------

OBJS_JPI_JI_12x4Kernel := driver.o Gemm_JPI_JI_12x4Kernel.o

driver_JPI_JI_12x4Kernel.x: $(OBJS_JPI_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_JPI_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JPI_JI_12x4Kernel.x $(LDFLAGS) 

JPI_JI_12x4Kernel: driver_JPI_JI_12x4Kernel.x
	echo "Name = 'JPI\_JI\_12x4Kernel';" > data/output_JPI_JI_12x4Kernel.m
	echo " " >> data/output_JPI_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_JPI_JI_12x4Kernel.x >> data/output_JPI_JI_12x4Kernel.m
	more data/output_JPI_JI_12x4Kernel.m

# ---------------------

OBJS_PIJ_JI_12x4Kernel := driver.o Gemm_PIJ_JI_12x4Kernel.o

driver_PIJ_JI_12x4Kernel.x: $(OBJS_PIJ_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PIJ_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PIJ_JI_12x4Kernel.x $(LDFLAGS) 

PIJ_JI_12x4Kernel: driver_PIJ_JI_12x4Kernel.x
	echo "Name = 'PIJ\_JI\_12x4Kernel';" > data/output_PIJ_JI_12x4Kernel.m
	echo " " >> data/output_PIJ_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_PIJ_JI_12x4Kernel.x >> data/output_PIJ_JI_12x4Kernel.m
	more data/output_PIJ_JI_12x4Kernel.m

# ---------------------

OBJS_PJI_JI_12x4Kernel := driver.o Gemm_PJI_JI_12x4Kernel.o

driver_PJI_JI_12x4Kernel.x: $(OBJS_PJI_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PJI_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJI_JI_12x4Kernel.x $(LDFLAGS) 

PJI_JI_12x4Kernel: driver_PJI_JI_12x4Kernel.x
	echo "Name = 'PJI\_JI\_12x4Kernel';" > data/output_PJI_JI_12x4Kernel.m
	echo " " >> data/output_PJI_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_PJI_JI_12x4Kernel.x >> data/output_PJI_JI_12x4Kernel.m
	more data/output_PJI_JI_12x4Kernel.m

# ---------------------

OBJS_PI_JI_12x4Kernel := driver.o Gemm_PI_JI_12x4Kernel.o

driver_PI_JI_12x4Kernel.x: $(OBJS_PI_JI_12x4Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PI_JI_12x4Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_JI_12x4Kernel.x $(LDFLAGS) 

PI_JI_12x4Kernel: driver_PI_JI_12x4Kernel.x
	echo "Name = 'PI\_JI\_12x4Kernel';" > data/output_PI_JI_12x4Kernel.m
	echo " " >> data/output_PI_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_PI_JI_12x4Kernel.x >> data/output_PI_JI_12x4Kernel.m
	more data/output_PI_JI_12x4Kernel.m

# ---------------------

OBJS_PI_JI_12x4Kernel_MCxKC := driver_MCxKC.o Gemm_PI_JI_12x4Kernel_MCxKC.o

driver_PI_JI_12x4Kernel_MCxKC.x: $(OBJS_PI_JI_12x4Kernel_MCxKC) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PI_JI_12x4Kernel_MCxKC) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_JI_12x4Kernel_MCxKC.x $(LDFLAGS) 

PI_JI_12x4Kernel_MCxKC: driver_PI_JI_12x4Kernel_MCxKC.x
	echo "Name = 'PI\_JI\_12x4Kernel\_MCxKC';" > data/output_PI_JI_12x4Kernel_MCxKC.m
	echo " " >> data/output_PI_JI_12x4Kernel_MCxKC.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PI_JI_12x4Kernel_MCxKC.x >> data/output_PI_JI_12x4Kernel_MCxKC.m
	more data/output_PI_JI_12x4Kernel_MCxKC.m

# ---------------------

OBJS_PI_JI_12x4Kernel_ldim := driver_ldim.o Gemm_PI_JI_12x4Kernel.o

driver_PI_JI_12x4Kernel_ldim.x: $(OBJS_PI_JI_12x4Kernel_ldim) $(UTIL_OBJS)
	$(LINKER) $(OBJS_PI_JI_12x4Kernel_ldim) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PI_JI_12x4Kernel_ldim.x $(LDFLAGS) 

PI_JI_12x4Kernel_ldim: driver_PI_JI_12x4Kernel_ldim.x
	echo "Name = 'PI\_JI\_12x4Kernel_ldim';" > data/output_PI_JI_12x4Kernel_ldim.m
	echo " " >> data/output_PI_JI_12x4Kernel_ldim.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_PI_JI_12x4Kernel_ldim.x >> data/output_PI_JI_12x4Kernel_ldim.m
	more data/output_PI_JI_12x4Kernel_ldim.m

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

OBJS_Five_Loops_12x4Kernel_ldim := driver_ldim.o Gemm_Five_Loops_12x4Kernel.o

driver_Five_Loops_12x4Kernel_ldim.x: $(OBJS_Five_Loops_12x4Kernel_ldim) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Five_Loops_12x4Kernel_ldim) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_12x4Kernel_ldim.x $(LDFLAGS) 

Five_Loops_12x4Kernel_ldim: driver_Five_Loops_12x4Kernel_ldim.x
	echo "Name = 'Five\_Loops\_12x4Kernel\_ldim';" > data/output_Five_Loops_12x4Kernel_ldim.m
	echo " " >> data/output_Five_Loops_12x4Kernel_ldim.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Five_Loops_12x4Kernel_ldim.x >> data/output_Five_Loops_12x4Kernel_ldim.m
	more data/output_Five_Loops_12x4Kernel_ldim.m

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

OBJS_Five_Loops_Pack_6x8Kernel := driver.o Gemm_Five_Loops_Pack_6x8Kernel.o

driver_Five_Loops_Pack_6x8Kernel.x: $(OBJS_Five_Loops_Pack_6x8Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Five_Loops_Pack_6x8Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_Pack_6x8Kernel.x $(LDFLAGS) 

Five_Loops_Pack_6x8Kernel: driver_Five_Loops_Pack_6x8Kernel.x
	echo "Name = 'Five\_Loops\_Pack\_6x8Kernel';" > data/output_Five_Loops_Pack_6x8Kernel.m
	echo " " >> data/output_Five_Loops_Pack_6x8Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Five_Loops_Pack_6x8Kernel.x >> data/output_Five_Loops_Pack_6x8Kernel.m
	more data/output_Five_Loops_Pack_6x8Kernel.m

# ---------------------

OBJS_Five_Loops_Pack_8x6Kernel := driver.o Gemm_Five_Loops_Pack_8x6Kernel.o

driver_Five_Loops_Pack_8x6Kernel.x: $(OBJS_Five_Loops_Pack_8x6Kernel) $(UTIL_OBJS)
	$(LINKER) $(OBJS_Five_Loops_Pack_8x6Kernel) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Five_Loops_Pack_8x6Kernel.x $(LDFLAGS) 

Five_Loops_Pack_8x6Kernel: driver_Five_Loops_Pack_8x6Kernel.x
	echo "Name = 'Five\_Loops\_Pack\_8x6Kernel';" > data/output_Five_Loops_Pack_8x6Kernel.m
	echo " " >> data/output_Five_Loops_Pack_8x6Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Five_Loops_Pack_8x6Kernel.x >> data/output_Five_Loops_Pack_8x6Kernel.m
	more data/output_Five_Loops_Pack_8x6Kernel.m

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
