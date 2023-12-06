# Fortran compiler
FC = mpifort

# Compiler flags
FFLAGS = -O3 -cpp -std=f2003 -fcheck=all -Wall

# Library path
LIB_PATH = /home/omar/WORK/codes/MUI/MUI_2.0/wrappers/Fortran

# Output executable
OUTPUT = dummySendReceive

# Source files
SOURCES = $(LIB_PATH)/mui_f_wrapper_3d.f90 dummySendReceive.f90

# Libraries
LIBRARIES = -lMUI_Fortran_wrapper -lstdc++ -lmpi_cxx

# Compilation rule
$(OUTPUT): $(SOURCES)
	$(FC) $(FFLAGS) -o $@ $(SOURCES) -L$(LIB_PATH) -I$(LIB_PATH) $(LIBRARIES) -Wl,-rpath,$(LIB_PATH)

# Phony target to clean compiled files
clean:
	rm -f $(OUTPUT)