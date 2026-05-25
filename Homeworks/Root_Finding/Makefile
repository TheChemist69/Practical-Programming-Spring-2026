# Compiler configuration for this homework project.
CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm

MAIN_OBJS = main.o vector.o matrix.o qr.o newton.o test_problems.o ode.o hydrogen.o reporting.o

all: Out.txt plot_wavefunction.svg plot_conv_rmax.svg plot_conv_rmin.svg plot_conv_acc.svg plot_conv_eps.svg

Out.txt: main
	./main > Out.txt
	@cat Out.txt

hydrogen_wavefunction.data hydrogen_convergence_rmax.data hydrogen_convergence_rmin.data hydrogen_convergence_acc.data hydrogen_convergence_eps.data line_search_comparison.data: Out.txt

plot_wavefunction.svg: hydrogen_wavefunction.data plot_wavefunction.gnuplot
	gnuplot plot_wavefunction.gnuplot

plot_conv_rmax.svg: hydrogen_convergence_rmax.data plot_conv_rmax.gnuplot
	gnuplot plot_conv_rmax.gnuplot

plot_conv_rmin.svg: hydrogen_convergence_rmin.data plot_conv_rmin.gnuplot
	gnuplot plot_conv_rmin.gnuplot

plot_conv_acc.svg: hydrogen_convergence_acc.data plot_conv_acc.gnuplot
	gnuplot plot_conv_acc.gnuplot

plot_conv_eps.svg: hydrogen_convergence_eps.data plot_conv_eps.gnuplot
	gnuplot plot_conv_eps.gnuplot

main: $(MAIN_OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

main.o: main.cc newton.h test_problems.h hydrogen.h reporting.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vector.o: vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

matrix.o: matrix.cc matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

qr.o: qr.cc qr.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

newton.o: newton.cc newton.h matrix.h qr.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

test_problems.o: test_problems.cc test_problems.h newton.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

ode.o: ode.cc ode.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

hydrogen.o: hydrogen.cc hydrogen.h ode.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

reporting.o: reporting.cc reporting.h hydrogen.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) main *.o Out.txt *.data *.svg

.PHONY: all clean
