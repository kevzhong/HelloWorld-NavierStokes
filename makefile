FC = gfortran

SRCDIR = src
OBJDIR = obj
OUTDIR = outputdir

TARGET = a.out

FFLAGS = -O3 -cpp -fdefault-real-8 -fdefault-double-8 -fcheck=bounds -fallow-argument-mismatch -funroll-loops -std=f2008 -Wall -Wextra -Wpedantic -fopenmp -I$(OBJDIR) -J$(OBJDIR)

#----FFTW library------------------- (adjust prefix if needed)--------------
# For personal Linux machine
FFTW_PREFIX := $(shell pkg-config --variable=prefix fftw3)

# For personal MAC machine
#FFTW_PREFIX := $(shell brew --prefix fftw)

#--- End prefix modification -----------------------------------------------

FFLAGS += -I$(FFTW_PREFIX)/include 
LIBS = -L$(FFTW_PREFIX)/lib  -lfftw3_omp -lfftw3 -lm

#----- END FFTW library-------------------

# Source files
SRC =	$(SRCDIR)/main.f90 \
		$(SRCDIR)/memAlloc.f90 \
		$(SRCDIR)/grid.f90 \
		$(SRCDIR)/dt.f90 \
		$(SRCDIR)/init_FFT.f90 \
		$(SRCDIR)/initialCondition.f90 \
		$(SRCDIR)/updateFields.f90 \
		$(SRCDIR)/fullyExplicit_update.f90 \
		$(SRCDIR)/rhsVelocity.f90 \
		$(SRCDIR)/rhsScalar.f90 \
		$(SRCDIR)/ADI_implicitUpdate.f90 \
		$(SRCDIR)/pressurePoisson.f90 \
		$(SRCDIR)/helmholtz_implicitUpdate.f90 \
		$(SRCDIR)/solve_helmholtz.f90 \
		$(SRCDIR)/tridiag.f90 \
		$(SRCDIR)/stats.f90 \
		$(SRCDIR)/fileIO.f90 \

MOD_FILES = $(SRCDIR)/my_fftw.f90 \
			$(SRCDIR)/types.f90 \
			$(SRCDIR)/params.f90 \
			$(SRCDIR)/rk3.f90 \
			$(SRCDIR)/ghost.f90 \

# Object files
OBJ = $(addprefix $(OBJDIR)/, $(notdir $(SRC:.f90=.o)))
MOBJ = $(addprefix $(OBJDIR)/, $(notdir $(MOD_FILES:.f90=.o)))

# Default target
all: $(OUTDIR) $(OBJDIR) $(TARGET)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# Pattern rule to compile .f90 files into .o files in SRCDIR
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

## In src/adi/
#$(OBJDIR)/%.o: $(SRCDIR)/adi/%.f90
#	$(FC) $(FFLAGS) -c $< -o $@

# Linking the final executable
$(TARGET): $(MOBJ) $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(MOBJ) $(OBJ) $(LIBS)

# Clean up generated files
clean:
	@rm -f $(OBJDIR)/*.o $(TARGET)
	@rm -f $(OBJDIR)/*.mod

veryclean: clean
	@rm -rf $(OUTDIR)
	@rm -rf $(OBJDIR)

.PHONY : clean veryclean