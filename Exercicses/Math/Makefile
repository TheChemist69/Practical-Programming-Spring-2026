CXXFLAGS = -Wall -Wextra -Werror -std=c++23 -O
LDLIBS = -lstdc++ -lm


Out.txt : main
	./main > Out.txt
	@cat Out.txt

main : main.o task_1.o task_2.o task_3.o
	$(CXX) $^ -o $@ $(LDLIBS)

main.o: main.cc task_1.h task_2.h task_3.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

task_1.o: task_1.cc task_1.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

task_2.o: task_2.cc task_2.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

task_3.o: task_3.cc task_3.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) main main.o Out.txt *.o

test:
	@echo CXX
	@echo $(CXX)
	@echo $(RM)
