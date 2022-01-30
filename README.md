# Waves
Wave simulator for the nonlinear wave equation, u_tt=c^2*(u_xx+u_yy)-V'(u), using Lobatto IIIA-IIIB discretisation in space (with 2, 3, or 4 stages for now) to get explicit ODEs in time (see my PhD thesis or Section 3 of https://doi.org/10.1137/140958050 for details), which are then integrated using leapfrog to get a high-order fully explicit multisymplectic integrator.
The observed global numerical order in space for an r-stage discretisation is 2 for r=2 or r+1 for r>2, while the order in time is 2 for leapfrog (other time-stepping methods may give higher order in time).
Additionally, the explicit ODEs are local (only depending on neighbouring cells), which allows various boundary conditions to be handled simply.

A Qt frontend to the integrator can be found in the Qt folder.
As of version 3.0, the integrator has been implemented in OpenCL using boost.compute. This allows for a significant performance gain sufficient to run a 1 million element simulation using 3 stages in x and y.

## Build Instructions
Run `make` in `Qt/QtWaves`.
Hopefully everything required is there.

## Running
Run `./QtWaves` in the `Qt/QtWaves` folder and choose the desired configuration in the GUI.

Some configurations have a much higher number of cells than others, see the configurations below.

The `SPF` number in the `Time` box indicates the number of steps being taken per frame, so adjust the `Timestep` to get something approaching real-time.
Initial configuration was built for an older laptop with lower specs, so real-time simulations will generally require much lower time-steps.

## Configurations:
### Boundary Conditions
- Periodic square:
	- 200x200 cells square domain, centered at the origin with periodic boundary conditions.
	- dx = dy = 0.1
	- wave speed = 1
- Dirichlet, square:
	- 200x200 cells square domain with Dirichlet boundary conditions.
	- dx = dy = 0.1
	- wave speed = 1
- Dirichlet, circle:
	- Circle of radius 100 cells, centered at the origin with Dirichlet boundary conditions.
	- dx = dy = 0.1
	- wave speed = 1
- Dirichlet, circle with a cusp:
	- Circle of radius 1303 cells (totalling ~1 million cells) with Dirichlet boundar conditions.
	- dx = dy = 0.1
	- wave speed = 1
- Dirichlet, intersecting circles:
	- Two slightly intersecting circles within a 200x200 cell square with Dirichlet boundary conditions.
	- dx = dy = 0.1
	- wave speed = 10
- Double slit (mixed BCs)
	- A double slit in a 1000x1000 cell square with Dirichlet boundary conditions at the slit and two edges, periodic boundary conditions on the side edges.
	- dx = dy = 0.025
	- wave speed = 0.3

### Initial Conditions
- Single frequency:
	- A traveling cosine wave.
- Continuous spectrum:
	- A wide Gaussian hump with a continuous spectrum of frequencies centred at the origin.
- Localised hump:
	- A narrow Gaussian hump limited to an area around the origin.
	- `u=exp(f(x,y))-1`, where `f(x,y)=a*(1-x^2/bx^2-y^2/by^2)`
- Two localised humps:
	- As above, but two of them at the centre of two intersecting circles boundary condition.
- Off-centre localised hump:
	- Just one of the two localised humps.
- A localised wave:
	- A Gaussian hump in one dimension only with the other dimension being constant to make a wave.
	- Good starting condition for the double slit boundary condition.
