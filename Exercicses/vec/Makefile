CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++14
TARGET = main
OBJS = main.o


all: run
run: $(TARGET)
	./$(TARGET) > Out.txt
	@cat Out.txt
Out.txt : main
	./main > Out.txt
	@cat Out.txt

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

main.o: main.cc vec.h
	$(CXX) $(CXXFLAGS) -c main.cc

clean:
	rm -f $(TARGET) $(OBJS)

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
