# Waves
Wave simulator for the nonlinear wave equation, u_tt=c^2*(u_xx+u_yy)-V'(u), using Lobatto IIIA-IIIB discretisation in space (with 2, 3, or 4 stages for now) to get explicit ODEs in time (see my PhD thesis or Section 3 of https://doi.org/10.1137/140958050 for details), which are then integrated using leapfrog to get a high-order fully explicit multisymplectic integrator. The observed numerical order in space for an r-stage discretisation is 2 for r=2 or r+1 for r>2, while the order in time is 2 for leapfrog (other time-stepping methods may give higher order in time).

A Qt frontend to the integrator can be found in the Qt folder.
As of version 3.0, the integrator has been implemented in OpenCL using boost.compute. This allows for a significant performance gain sufficient to run a 1 million element simulation using 3 stages in x and y.
