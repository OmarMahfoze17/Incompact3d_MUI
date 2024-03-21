# Fortran compiler
FC = mpif90

# Compiler flags
FFLAGS = -O3 -cpp -std=f2003 -fcheck=all -Wall

# Library path
MUI_PATH = /home/omar/WORK/codes/MUI/MUI_2.0/wrappers/Fortran

# Output executable
OUTPUT = dummyWRF

# Source files
SOURCES = $(MUI_PATH)/mui_f_wrapper_3d.f90 dummyWRF.f90

# Libraries
LIBRARIES = -lMUI_Fortran_wrapper -lstdc++ -lmpi_cxx -I/usr/include -L/usr/lib -lnetcdff

# Compilation rule
$(OUTPUT): $(SOURCES)
	$(FC) $(FFLAGS) -o $@ $(SOURCES) -L$(MUI_PATH) -I$(MUI_PATH) $(LIBRARIES) -Wl,-rpath,$(MUI_PATH)

# Phony target to clean compiled files
clean:
	rm -f $(OUTPUT)
