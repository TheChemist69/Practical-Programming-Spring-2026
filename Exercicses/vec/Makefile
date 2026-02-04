 CXX = g++
CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm

Out.txt : main
	./main > Out.txt
	@cat Out.txt

main : main.o vec.o
	$(CXX) $^ -o $@ $(LDLIBS)

main.o: main.cc vec.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vec.o: vec.cc vec.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) main main.o Out.txt *.o

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
