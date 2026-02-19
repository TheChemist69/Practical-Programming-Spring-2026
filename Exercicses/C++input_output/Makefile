CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
TARGET = main
OBJS = main.o input_output.o

all: $(TARGET)
Out.txt : main
	./main > Out.txt
	@cat Out.txt

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) -lm

main.o: main.cc input_output.h
	$(CXX) $(CXXFLAGS) -c main.cc
	

input_output.o: input_output.cc input_output.h
	$(CXX) $(CXXFLAGS) -c input_output.cc

clean:
	rm -f $(TARGET) $(OBJS)
