SHELL = /bin/sh

CC = g++
FLAGS = -std=c++11 -g -pedantic \
	-Werror=int-conversion -Werror=implicit -L/usr/lib/openmpi/  \
	-Werror=return-type -Werror=uninitialized -Wall

LNK = -L./ -llammpstools
INC = -I./ -I../../../c_lib/

EXT = cpp
SRC = $(wildcard *.$(EXT))

COMP = $(CC) $(FLAGS) $(INC)
LINK = $(CC) $(FLAGS) $(INC) $(LNK)
AR   = ar rcs

OBJ_DIR = Obj
OBJ = $(SRC:%.$(EXT)=$(OBJ_DIR)$(S)%.o)
OBJ_DIRS = $(dir $(OBJ))
DEPS = $(OBJ:%.o=%.d)

.PHONY: dirs all help clean print-flags

OBJ_GUARD = $(OBJ_DIR)$(S).guard
# For windows:
#MAKE_DIR = $(if exist $(1),,mkdir $(1))
#S=\\
# Linux and Unix-like:
MAKE_DIR = mkdir -p $(1)
S=/

EXE = test_gsd

all : dirs $(EXE)

dirs : $(OBJ_DIR)

$(OBJ_DIR) :
	$(call $(MAKE_DIR),$@)

help :
	@echo "SRC    is   $(SRC)"
	@echo "OBJ    is   $(OBJ)"
	@echo "DEPS   is   $(DEPS)"
	@echo "INC    is   $(INC)"
	@echo "LIB    is   $(LIB)"


$(EXE) : $(OBJ)
	$(LINK) $(OBJ) -o $@

$(OBJ_DIR)$(S)%.o : %.$(EXT)
	$(call MAKE_DIR,$(dir $@))
	$(COMP) -c $< -o $@
	$(COMP) -M -MT '$@' $< -MF $(@:%.o=%.d)

clean:
	rm -r $(OBJ_DIR)
	rm -f $(EXE)


install: $(EXE)
	$(call MAKE_DIR,$(LIB_DIR))
	cp $(EXE) $(LIB_DIR)


uninstall:
	rm -rf $(LIB_DIR)/$(subst ../,,$(EXE))

print-flags:
	@echo "Compilation flags: $(COMP)"

-include $(DEPS)
