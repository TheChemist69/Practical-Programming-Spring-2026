CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm

MAIN_OBJS = main.o mc.o quasi.o problems.o reporting.o vector.o

all: Out.txt plot_plain_error.svg plot_quasi_vs_plain.svg plot_stratified_vs_plain.svg plot_difficult.svg

Out.txt: main
	./main > Out.txt
	@cat Out.txt

two_dim_examples.data plain_scaling.data quasi_scaling.data stratified_scaling.data difficult_integral.data: Out.txt

plot_plain_error.svg: plain_scaling.data plot_plain_error.gnuplot
	gnuplot plot_plain_error.gnuplot

plot_quasi_vs_plain.svg: quasi_scaling.data plot_quasi_vs_plain.gnuplot
	gnuplot plot_quasi_vs_plain.gnuplot

plot_stratified_vs_plain.svg: stratified_scaling.data plot_stratified_vs_plain.gnuplot
	gnuplot plot_stratified_vs_plain.gnuplot

plot_difficult.svg: difficult_integral.data plot_difficult.gnuplot
	gnuplot plot_difficult.gnuplot

main: $(MAIN_OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

main.o: main.cc mc.h quasi.h problems.h reporting.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

mc.o: mc.cc mc.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

quasi.o: quasi.cc quasi.h mc.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

problems.o: problems.cc problems.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

reporting.o: reporting.cc reporting.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vector.o: vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) main *.o Out.txt *.data *.svg fit.log

.PHONY: all clean
