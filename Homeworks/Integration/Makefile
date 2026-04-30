CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm

MAIN_OBJS = main.o integrator.o benchmarks.o erf_integral.o reporting.o

all: Out.txt plot_erf.svg plot_erf_acc.svg scipy_compare

Out.txt: main
	./main > Out.txt
	@cat Out.txt

benchmark_results.data erf_curve.data erf_tabulated.data erf_acc_sweep.data: Out.txt

plot_erf.svg: erf_curve.data erf_tabulated.data plot_erf.gnuplot
	gnuplot plot_erf.gnuplot

plot_erf_acc.svg: erf_acc_sweep.data plot_erf_acc.gnuplot
	gnuplot plot_erf_acc.gnuplot

scipy_compare: compare_scipy.py benchmark_results.data
	@if command -v python3 >/dev/null 2>&1; then \
		python3 compare_scipy.py; \
	else \
		echo "python3 not found, skipping SciPy comparison"; \
	fi

		$(CXX) $(CXXFLAGS) -c $< -o $@
	$(CXX) $^ -o $@ $(LDLIBS)
	matrix.o: ../matrix.cc ../matrix.h ../vector.h
main.o: main.cc benchmarks.h erf_integral.h reporting.h integrator.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

integrator.o: integrator.cc integrator.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

benchmarks.o: benchmarks.cc benchmarks.h integrator.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

erf_integral.o: erf_integral.cc erf_integral.h integrator.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

reporting.o: reporting.cc reporting.h benchmarks.h erf_integral.h integrator.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) main *.o Out.txt *.data *.svg

.PHONY: all clean scipy_compare
