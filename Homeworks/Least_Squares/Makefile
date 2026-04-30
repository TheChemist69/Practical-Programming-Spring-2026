CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm

LIB_OBJS = vector.o matrix.o qr.o lsfit.o
MAIN_OBJS = main.o $(LIB_OBJS)

all : Out.txt plot_fit.svg plot_logfit.svg

Out.txt : main
	./main > Out.txt
	@cat Out.txt

main : $(MAIN_OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

main.o : main.cc vector.h matrix.h lsfit.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vector.o : ../vector.cc ../vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

matrix.o : ../matrix.cc ../matrix.h ../vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

qr.o : qr.cc qr.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

lsfit.o : lsfit.cc lsfit.h qr.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ── Data files (produced by main alongside Out.txt) ──────────────
data.txt fit.txt logdata.txt logfit.txt : Out.txt

# ── Plots ────────────────────────────────────────────────────────
plot_fit.svg : data.txt fit.txt plot_fit.gnuplot
	gnuplot plot_fit.gnuplot

plot_logfit.svg : logdata.txt logfit.txt plot_logfit.gnuplot
	gnuplot plot_logfit.gnuplot

# ── Cleanup ──────────────────────────────────────────────────────
clean :
	rm -f *.o main Out.txt data.txt fit.txt logdata.txt logfit.txt *.svg
