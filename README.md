# Waves
Wave simulator for the nonlinear wave equation, u_tt=c^2*(u_xx+u_yy)-V'(u), using Lobatto IIIA-IIIB discretisation in space (order 2, 3, or 4 for now) to get explicit ODEs in time (see my PhD thesis for details), which are then integrated using leapfrog to get a fully explicit multisymplectic integrator.
