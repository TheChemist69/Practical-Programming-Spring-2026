CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm

LIB_OBJS = vector.o matrix.o evd.o
MAIN_OBJS = main.o $(LIB_OBJS)

all : Out.txt plot_convergence_dr.svg plot_convergence_rmax.svg plot_wavefunctions.svg

Out.txt : main
	./main -rmax 25 -dr 0.1 > Out.txt
	@cat Out.txt

main : $(MAIN_OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

main.o : main.cc vector.h matrix.h evd.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vector.o : vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

matrix.o : matrix.cc matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

evd.o : evd.cc evd.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ── Convergence studies ──────────────────────────────────────────
convergence_dr.data : main
	./main -rmax 25 -mode convergence_dr

convergence_rmax.data : main
	./main -dr 0.1 -mode convergence_rmax

# ── Wavefunctions (produced alongside Out.txt) ──────────────────
wavefunctions.data : Out.txt

# ── Plots ────────────────────────────────────────────────────────
plot_convergence_dr.svg : convergence_dr.data plot_convergence_dr.gnuplot
	gnuplot plot_convergence_dr.gnuplot

plot_convergence_rmax.svg : convergence_rmax.data plot_convergence_rmax.gnuplot
	gnuplot plot_convergence_rmax.gnuplot

plot_wavefunctions.svg : wavefunctions.data plot_wavefunctions.gnuplot
	gnuplot plot_wavefunctions.gnuplot

# ── Timing / scaling (Task C) ───────────────────────────────────
Ns := $(shell seq 20 5 80)
out.times.data : main
	>$@
	for N in $(Ns); do \
		./main -size $$N >> $@ ; \
	done
	sort -n $@ -o $@

plot_times.svg : out.times.data plot_times.gnuplot
	gnuplot plot_times.gnuplot

# ── Cleanup ──────────────────────────────────────────────────────
clean :
	$(RM) main *.o Out.txt *.data *.svg

.PHONY : clean all
