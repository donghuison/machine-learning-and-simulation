PROJECT_DIR := $(shell cd .. && pwd)
include $(PROJECT_DIR)/Makefile_inc

TARGET = a.out
OBJS =  openfile.o const.o model.o init.o \
		integrate.o bnd.o main.o rhs.o

COMMON_LIB = $(PROJECT_DIR)/libcommon.a
COMMON_INC = $(PROJECT_DIR)/common
LDFLAGS = -L$(PROJECT_DIR) -lcommon

.PHONY: all clean clean_all clean_data

all: $(TARGET)

$(TARGET): $(OBJS) $(COMMON_LIB)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.f90
	$(FC) $(FFLAGS) -I$(COMMON_INC) -c $<

$(COMMON_LIB):
	$(MAKE) -C $(COMMON_INC)

#Dependencies
init.o: $(COMMON_LIB) const.o model.o bnd.o
model.o: $(COMMON_LIB) const.o
bnd.o: $(COMMON_LIB)
integrate.o: $(COMMON_LIB) bnd.o rhs.o
main.o: $(COMMON_LIB) init.o const.o openfile.o integrate.o
openfile.o: $(COMMON_LIB) const.o
rhs.o: $(COMMON_LIB)

clean:
	rm -f *.mod *.out *.i90 $(TARGET) $(OBJS)

# clean_data:
# 	find . -name "*.dac" -delete
clean_data_300:
	find ./data_300x300 -name "*.dac" -delete
	find ./data_300x300 -name "*.txt" -delete
	@echo "Cleaned .dac and .txt files from data_300x300 directory"

clean_data_600:
	find ./data_600x600 -name "*.dac" -delete
	find ./data_600x600 -name "*.txt" -delete
	@echo "Cleaned .dac and .txt files from data_600x600 directory"

clean_data_1200:
	find ./data_1200x1200 -name "*.dac" -delete
	find ./data_1200x1200 -name "*.txt" -delete
	@echo "Cleaned .dac and .txt files from data_1200x1200 directory"

clean_data_2400:
	find ./data_2400x2400 -name "*.dac" -delete
	find ./data_2400x2400 -name "*.txt" -delete
	@echo "Cleaned .dac and .txt files from data_2400x2400 directory"