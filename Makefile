G++ = g++
G++_FLAGS = -std=c++17 -c -Wall -Ideps/eigen-3.3.8 -g

LD_FLAGS = -lgtest -lpthread

SOURCE_FILE = $(wildcard *.cpp) $(wildcard ./*/*.cpp)
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCE_FILE))
TARGET = solver

all: $(TARGET)

$(TARGET): $(OBJECTS)
	g++ -o $(TARGET) $(OBJECTS) $(LD_FLAGS)

%.o : %.cpp
	$(G++) -o $@ $(G++_FLAGS) $<

clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: all clean
