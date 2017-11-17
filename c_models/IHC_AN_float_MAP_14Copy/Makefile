ifndef NEURAL_MODELLING_DIRS
    $(error NEURAL_MODELLING_DIRS is not set.  Please define NEURAL_MODELLING_DIRS (possibly by running "source setup" in the neural_modelling folder within the sPyNNaker source folder))
endif

APP := IHC_AN_float_MAP
BUILD_DIR = build/
SOURCES = src/src/spinn_IHC_AN_softfloat_MAP_pipe_MC.c#spinn_IHC_AN_softfloat_MAP_single.c#src/spinn_IHC_AN_softfloat_MAP.c
MAKEFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CURRENT_DIR := $(dir $(MAKEFILE_PATH))
EXTRA_SRC_DIR := $(abspath $(CURRENT_DIR)/src) 
SOURCE_DIRS += $(EXTRA_SRC_DIR)
SOURCE_DIRS += $(SPINN_DIRS)/spinn_common/src
APP_OUTPUT_DIR := $(abspath $(CURRENT_DIR)../../../OME_SpiNN/model_binaries/)/#$(abspath $(CURRENT_DIR)../../model_binaries/)/#
CFLAGS +=-lm -I$(NEURAL_MODELLING_DIRS)/src -Isrc/ -Iinclude/ -I$(SPINN_DIRS)/spinn_common/include -ffast-math

include $(NEURAL_MODELLING_DIRS)/src/Makefile.common
