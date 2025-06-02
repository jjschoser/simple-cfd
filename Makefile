CXX=g++
CXXFLAGS=-O3 -Wall -Wextra

simple-cfd : main.o EquationOfState.o Euler.o FluxSolver.o Mesh.o Solver.o
	$(CXX) $^ -o $@

main.o : main.cpp EquationOfState.H Euler.H Macros.H Mesh.H Solver.H
	$(CXX) -c $< -o $@ $(CXXFLAGS)

EquationOfState.o : EquationOfState.cpp EquationOfState.H Macros.H
	$(CXX) -c $< -o $@ $(CXXFLAGS)

Euler.o : Euler.cpp Euler.H EquationOfState.H Macros.H
	$(CXX) -c $< -o $@ $(CXXFLAGS)

FluxSolver.o : FluxSolver.cpp FluxSolver.H Euler.H Macros.H
	$(CXX) -c $< -o $@ $(CXXFLAGS)

Mesh.o : Mesh.cpp Mesh.H Macros.H
	$(CXX) -c $< -o $@ $(CXXFLAGS)

Solver.o : Solver.cpp Solver.H Euler.H FluxSolver.H Macros.H Mesh.H
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f *.o simple-cfd