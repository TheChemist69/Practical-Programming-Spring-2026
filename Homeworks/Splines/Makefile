CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm

LIB_OBJS = vector.o spline.o
MAIN_OBJS = main.o $(LIB_OBJS)

all : Out.txt plot_linear.svg plot_quadratic.svg plot_cubic.svg

Out.txt : main
	./main > Out.txt
	@cat Out.txt

main : $(MAIN_OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

main.o : main.cc vector.h spline.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vector.o : vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

spline.o : spline.cc spline.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# -- Data files (produced by main alongside Out.txt) ----------------
cos_points.txt linear_interp.txt sin_points.txt quad_interp.txt cubic_interp.txt : Out.txt

# -- Plots -----------------------------------------------------------
plot_linear.svg : cos_points.txt linear_interp.txt plot_linear.gnuplot
	gnuplot plot_linear.gnuplot

plot_quadratic.svg : sin_points.txt quad_interp.txt plot_quadratic.gnuplot
	gnuplot plot_quadratic.gnuplot

plot_cubic.svg : cos_points.txt cubic_interp.txt plot_cubic.gnuplot
	gnuplot plot_cubic.gnuplot

# -- Cleanup ---------------------------------------------------------
clean :
	rm -f *.o main Out.txt cos_points.txt sin_points.txt linear_interp.txt quad_interp.txt cubic_interp.txt *.svg
