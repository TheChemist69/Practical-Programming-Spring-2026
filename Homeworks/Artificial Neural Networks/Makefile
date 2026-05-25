CXX := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic -I..

TARGET := main
OBJ := main.o ann.o activation.o optimizer.o output_writer.o vector.o

all: $(TARGET) run

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

main.o: main.cc ann.h activation.h optimizer.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

ann.o: ann.cc ann.h activation.h optimizer.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

activation.o: activation.cc activation.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

output_writer.o: output_writer.cc output_writer.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

optimizer.o: optimizer.cc optimizer.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

vector.o: ../vector.cc ../vector.h
	$(CXX) $(CXXFLAGS) -c ../vector.cc -o $@

clean:
	rm -f $(OBJ) $(TARGET) Out.txtx ann_fit.data ann_derivatives.data ann_ode.data ann_error.data ann_ode_residual.data ann_train_history.data ann_ode_history.data \
		plot_fit.svg plot_derivatives.svg plot_ode.svg plot_error.svg plot_ode_residual.svg plot_training.svg

run: $(TARGET)
	./$(TARGET)
	gnuplot plot_fit.gnuplot
	gnuplot plot_derivatives.gnuplot
	gnuplot plot_ode.gnuplot
	gnuplot plot_error.gnuplot
	gnuplot plot_ode_residual.gnuplot
	gnuplot plot_training.gnuplot

.PHONY: all clean run
