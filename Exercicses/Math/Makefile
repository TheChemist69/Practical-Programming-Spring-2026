CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm


Out.txt : main
	./main > Out.txt
	@cat Out.txt

main : main.o task_1.o task_2.o task_3.o
	$(CXX) main.o task_1.o task_2.o task_3.o -o main $(LDLIBS)

main.o: task_1.h task_2.h task_3.h main.cc
	$(CXX) $(CXXFLAGS) -c main.cc #-o main.o

task_1.o: task_1.h task_1.cc
	$(CXX) $(CXXFLAGS) -c task_1.cc #-o task_1.o

task_2.o: task_2.h task_2.cc
	$(CXX) $(CXXFLAGS) -c task_2.cc #-o task_2.o

task_3.o: task_3.h task_3.cc 
	$(CXX) $(CXXFLAGS) -c task_3.cc #-o task_3.o




clean:
	$(RM)  main main.o Out.txt *.o

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
