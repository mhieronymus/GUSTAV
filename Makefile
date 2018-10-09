CXX = nvcc 
CXXFLAGS = -std=c++11 -arch=sm_50 -D_FORCE_INLINES

BUILD = build
OBJ_DIR = $(BUILD)/objects 
APP_DIR = $(BUILD)/apps
TARGET = GUSTAV 
INCLUDE = -L/usr/include/ -lcfitsio
SRC = \
	$(wildcard src/*.cpp) \
	$(wildcard src/*.cu) \
	$(wildcard src/helper/photospline/*.c) \

OBJECTS = $(SRC:%.*=$(OBJ_DIR)/%.o)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.*
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ 

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $(APP_DIR)/$(TARGET) $(OBJECTS) 


.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -g
debug: all

release: CXXFLAGS += -O2 
release: all


