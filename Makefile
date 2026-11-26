CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
EIGEN_INC ?= /usr/include/eigen3

TARGET := NucMutInt
SRCS := src/HelixFromTmpl.cpp third_party/tinyxml2/tinyxml2.cpp
INCLUDES := -Isrc -Ithird_party/tinyxml2 -I$(EIGEN_INC)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRCS) src/PNI_Calculate.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET) $(TARGET).exe
