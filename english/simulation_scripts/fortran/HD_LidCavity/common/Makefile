PROJECT_DIR := $(shell cd .. && pwd)
include $(PROJECT_DIR)/Makefile_inc

TARGET = $(PROJECT_DIR)/libcommon.a
OBJS = mpi_setup.o operator.o getNewdt.o flux_calc.o boundary.o dac_header.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	ar -r $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

# Dependencies
boundary.o : mpi_setup.o
# getNewdt.o : mpi_setup.o
flux_calc.o : operator.o

clean:
	rm -f $(OBJS) *.mod *.i90 $(TARGET)
