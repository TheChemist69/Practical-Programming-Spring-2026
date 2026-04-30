
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++23 -O2

LIB_OBJS =
MAIN_OBJS = main.o minimization.o

all : Out.txt plot_fit.svg

# ── Run program (produces Out.txt, data.out, fit.out) ────────────
Out.txt : main higgs.data.txt
	./main < higgs.data.txt > Out.txt
	@cat Out.txt

# ── Link ─────────────────────────────────────────────────────────
main : $(MAIN_OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# ── Compile ──────────────────────────────────────────────────────
main.o : main.cc minimization.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

minimization.o : minimization.cc minimization.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ── Data files (side effects of Out.txt) ─────────────────────────
data.out fit.out : Out.txt

# ── Plot ─────────────────────────────────────────────────────────
plot_fit.svg : data.out fit.out plot.gnuplot
	gnuplot plot.gnuplot


# ── Cleanup ──────────────────────────────────────────────────────
clean :
	rm -f *.o main Out.txt data.out fit.out *.svg
