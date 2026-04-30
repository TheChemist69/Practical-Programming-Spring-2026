CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm

MAIN_OBJS = main.o problems.o rkstep.o vector.o

all: Out.txt plot_oscillator.svg plot_lotka_volterra.svg plot_relativistic.svg plot_three_body.svg

Out.txt: main
	./main > Out.txt
	@cat Out.txt

oscillator.data lotka_volterra.data orbit_newtonian_circle.data orbit_newtonian_ellipse.data orbit_relativistic.data figure8.data: Out.txt

plot_oscillator.svg: oscillator.data plot_oscillator.gnuplot
	gnuplot plot_oscillator.gnuplot

plot_lotka_volterra.svg: lotka_volterra.data plot_lotka_volterra.gnuplot
	gnuplot plot_lotka_volterra.gnuplot

plot_relativistic.svg: orbit_newtonian_circle.data orbit_newtonian_ellipse.data orbit_relativistic.data plot_relativistic.gnuplot
	gnuplot plot_relativistic.gnuplot

plot_three_body.svg: figure8.data plot_three_body.gnuplot
	gnuplot plot_three_body.gnuplot

main: $(MAIN_OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

main.o: main.cc problems.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

problems.o: problems.cc problems.h rkstep.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

rkstep.o: rkstep.cc rkstep.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vector.o: ../vector.cc ../vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) main *.o Out.txt *.data *.svg

.PHONY: clean all
