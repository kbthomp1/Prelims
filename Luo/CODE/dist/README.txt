This code explicitly solves the Laplace equation for a channel with a spherical
bump.  The namelist "fe_input.nml" controls the user-definable parameters for the
solver.  The namelist parameters used are the same as those used to achieve the
results presented in the report.  Given that a large number of timesteps are
required to converge the solution to steady state, a robust restart capability
was added.  If at any point the the user wishes to stop the solver gracefully,
simply run

  touch stop.dat

This creates a file called "stop.dat" which the solver looks for at each
iteration.  If found, the iteration loop is halted and all restart/output files
are written.  When restarting, it is important to run

  rm stop.dat

and also add/change the namelist option

  read_restart = .true.

which will look for the restart file "flow.restart", which contains the solution
nodal coefficients at each degree of freedom.
