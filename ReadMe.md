Implementations of MISOCP formulation of `COMPUTING EXACT D-OPTIMAL DESIGNS BY MIXED INTEGER SECOND-ORDER CONE PROGRAMMING`
Note: there is an epsilon paramerer to add the information matrix M with the term epsion * I, which can be modelled as a slack variable w_0, A_0 = epsion * I

1. data.py generates a set of instances with data matrix
2. Two implementations: One is implemented julia (with solver options: CPLEX.jl, GUROBI.jl, SCIP.jl. SCIP.jl does not support geometric means, but it should be easy). The other is implemented in SCIP. Now, SCIP seems to have numerical problems, and CPLEX is the most stable solver. Using Gurobi's log display, you can find the numerical condition of the problems.
3. Install SCIP solver to read instance:  "cd build / cmake .. -DSCIP_DIR=$SCIP_DR", where $SCIP_DR is the location of SCIP installation containing the directory "build" and "src".
