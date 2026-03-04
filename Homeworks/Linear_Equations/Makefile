CXXFLAGS = -Wall -Wextra -std=c++23 -O2
LDLIBS = -lstdc++ -lm
TIME = /usr/bin/time --format "$$N %e" --append --output $@

OBJS = main.o vector.o matrix.o qr.o lu.o cholesky.o householder.o givens.o

Out.txt : main
	./main > Out.txt
	@cat Out.txt

main : $(OBJS)
	$(CXX) $^ -o $@ $(LDLIBS)

main.o : main.cc vector.h matrix.h qr.h lu.h cholesky.h householder.h givens.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vector.o : vector.cc vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

matrix.o : matrix.cc matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

qr.o : qr.cc qr.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

lu.o : lu.cc lu.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

cholesky.o : cholesky.cc cholesky.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

householder.o : householder.cc householder.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

givens.o : givens.cc givens.h matrix.h vector.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

out.times.data : main
	>$@
	for N in $$(seq 100 50 1000); do \
		$(TIME) ./main -size $$N 1>out 2>err ;\
	done
	@cat $@

plot.svg : out.times.data plot_times.gnuplot
	gnuplot plot_times.gnuplot

clean :
	$(RM) main *.o Out.txt out.times.data out err plot.svg

.PHONY : clean
