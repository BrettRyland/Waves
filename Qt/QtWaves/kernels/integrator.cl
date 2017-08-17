// The precision used by the integrator. This must be set in all .cl files to mirror that in integratorCL.h.
typedef float integrator_precision;

// This should actually by unsigned to match the definition in integratorCL.h, but C99 doesn't allow this. Fortunately, the number of types is small enough to not be an issue.
typedef enum { Normal, Periodic_left, Periodic_right, Dirichlet_left, Dirichlet_right, Neumann_left, Neumann_right } Cell_Type;

// Perform a step of size dt in the U variable
__kernel void step_U(__global integrator_precision *U, __global const unsigned int *masks, __global const integrator_precision *V, __global const unsigned int *cell_type_x, __global const unsigned int *cell_type_y, unsigned int cell_type_count, unsigned int nodes_per_cell, integrator_precision dt) {
	size_t cell = get_global_id(0);
	for (unsigned int node = 0; node != nodes_per_cell; ++node)
		if (masks[(cell_type_y[cell] * cell_type_count + cell_type_x[cell])*nodes_per_cell + node] == 1) {
			U[cell*nodes_per_cell + node] += dt * V[cell*nodes_per_cell + node]; // No dissipation
			//U[id] *= 1.0 - dt * m_artificial_dissipation; // Dampens everything evenly.
		}
}

// V'(u)=2*u+4*u^3
inline integrator_precision derivative_of_the_potential(integrator_precision u) { return 2.0 * u + 4.0 * u * u * u; }

// Perform a step of size dt in the V variable
__kernel void step_V(__global integrator_precision *V, __global const unsigned int *masks, __global const integrator_precision *U, __global const integrator_precision *U_xx, __global const integrator_precision *U_yy, __global const unsigned int *cell_type_x, __global const unsigned int *cell_type_y, unsigned int cell_type_count, unsigned int nodes_per_cell, integrator_precision dt, integrator_precision wavespeed) {
	size_t cell = get_global_id(0);
	for (unsigned int node = 0; node != nodes_per_cell; ++node)
		if (masks[(cell_type_y[cell] * cell_type_count + cell_type_x[cell])*nodes_per_cell + node] == 1)
			V[cell*nodes_per_cell + node] += dt * (wavespeed * (U_xx[cell*nodes_per_cell + node] + U_yy[cell*nodes_per_cell + node]) - derivative_of_the_potential(U[cell*nodes_per_cell + node]));
}

// A custom "inner product" on first and second that operates over count values with a custom start and step size for each vector and a given initial value.
// Note: no bounds checking.
inline integrator_precision inner_product_n(__global const integrator_precision *first, size_t start1, int step1, __global const integrator_precision* second, size_t start2, int step2, int count, integrator_precision initial_value) {
	for (int n = 0; n != count; ++n, start1 += step1, start2 += step2)
		initial_value += first[start1] * second[start2];
	return initial_value;
}

// Apply the main_node stencil to the given main node
inline integrator_precision main_node_xx(__global const integrator_precision *U, __global const integrator_precision *coefficients, size_t cell, size_t left_cell, size_t right_cell, unsigned int stages, integrator_precision dx) {
	return inner_product_n(U, left_cell, 1, coefficients, 0, 1, stages,
		inner_product_n(U, cell, 1, coefficients, stages, 1, stages,
			U[right_cell] * coefficients[2 * stages])) / dx / dx;
}

// Apply the appropriate internal_node stencil to the given internal node
inline integrator_precision internal_node_xx(__global const integrator_precision *U, __global const integrator_precision *coefficients, size_t cell, size_t right_cell, unsigned int node, unsigned int stages, integrator_precision dx) {
	return inner_product_n(U, cell, 1, coefficients, 2 * stages + 1 + (node - 1) * (stages + 1), 1, stages,
		U[right_cell] * coefficients[2 * stages + node * (stages + 1)]) / dx / dx;
}

// Compute the second derivative in x at each node in the cell
__kernel void compute_U_xx(__global integrator_precision *U_xx, __global const integrator_precision *U, __global const unsigned int *cell_type, __global const integrator_precision *coefficients, __global const unsigned int *adjacency_information, unsigned int stages_x, unsigned int stages_y, integrator_precision dx) {
	size_t cell = get_global_id(0);
	unsigned int nodes_per_cell = stages_x * stages_y;
	// FIXME the update mask should be implemented in here
	switch (cell_type[cell]) {
		case Normal:
		case Periodic_left:
		case Periodic_right: {
			for (unsigned int main_node = 0; main_node != nodes_per_cell; main_node += stages_x) {
				U_xx[cell * nodes_per_cell + main_node] = main_node_xx(U, coefficients, cell * nodes_per_cell + main_node, adjacency_information[cell * 4] * nodes_per_cell + main_node, adjacency_information[cell * 4 + 1] * nodes_per_cell + main_node, stages_x, dx);
				for (unsigned int internal_node = 1; internal_node != stages_x; ++internal_node)
					U_xx[cell * nodes_per_cell + main_node + internal_node] = internal_node_xx(U, coefficients, cell * nodes_per_cell + main_node, adjacency_information[cell * 4 + 1] * nodes_per_cell + main_node, internal_node, stages_x, dx);
			}
			break;
		}
		case Neumann_left: { // FIXME: not implemented yet
			break;
		}
		case Neumann_right: { // FIXME: not implemented yet
			break;
		}
		case Dirichlet_left: { // We only need the internal nodes
			for (unsigned int main_node = 0; main_node != nodes_per_cell; main_node += stages_x)
				for (unsigned int internal_node = 1; internal_node != stages_x; ++internal_node)
					U_xx[cell * nodes_per_cell + main_node + internal_node] = internal_node_xx(U, coefficients, cell * nodes_per_cell + main_node, adjacency_information[cell * 4 + 1] * nodes_per_cell + main_node, internal_node, stages_x, dx);
			break;
		}
		case Dirichlet_right: { // Do nothing. We don't need the derivatives for any nodes in this case.
			break;
		}
	}
}

// Apply the main_node stencil to the given main node
inline integrator_precision main_node_yy(__global const integrator_precision *U, __global const integrator_precision *coefficients, size_t cell, size_t left_cell, size_t right_cell, unsigned int stages_x, unsigned int stages_y, integrator_precision dy) {
	return inner_product_n(U, left_cell, stages_x, coefficients, 0, 1, stages_y,
		inner_product_n(U, cell, stages_x, coefficients, stages_y, 1, stages_y,
			U[right_cell] * coefficients[2 * stages_y])) / dy / dy;
}

// Apply the appropriate internal_node stencil to the given internal node
inline integrator_precision internal_node_yy(__global const integrator_precision *U, __global const integrator_precision *coefficients, size_t cell, size_t right_cell, unsigned int node, unsigned int stages_x, unsigned int stages_y, integrator_precision dy) {
	return inner_product_n(U, cell, stages_x, coefficients, 2 * stages_y + 1 + (node - 1) * (stages_y + 1), 1, stages_y,
		U[right_cell] * coefficients[2 * stages_y + node * (stages_y + 1)]) / dy / dy;
}

// Compute the second derivative in y at each node in the cell
__kernel void compute_U_yy(__global integrator_precision *U_yy, __global const integrator_precision *U, __global const unsigned int *cell_type, __global const integrator_precision *coefficients, __global const unsigned int *adjacency_information, unsigned int stages_x, unsigned int stages_y, integrator_precision dy) {
	size_t cell = get_global_id(0);
	unsigned int nodes_per_cell = stages_x * stages_y;
	// FIXME the update mask should be implemented in here
	switch (cell_type[cell]) {
		case Normal:
		case Periodic_left:
		case Periodic_right: {
			for (unsigned int main_node = 0; main_node != stages_x; ++main_node) {
				U_yy[cell * nodes_per_cell + main_node] = main_node_yy(U, coefficients, cell * nodes_per_cell + main_node, adjacency_information[cell * 4 + 2] * nodes_per_cell + main_node, adjacency_information[cell * 4 + 3] * nodes_per_cell + main_node, stages_x, stages_y, dy);
				for (unsigned int internal_node = 1; internal_node != stages_y; ++internal_node)
					U_yy[cell * nodes_per_cell + main_node + internal_node * stages_x] = internal_node_yy(U, coefficients, cell * nodes_per_cell + main_node, adjacency_information[cell * 4 + 3] * nodes_per_cell + main_node, internal_node, stages_x, stages_y, dy);
			}
			break;
		}
		case Neumann_left: { // We need all nodes, but we have to use phantom points in the left cell for the main node. FIXME: not implemented yet
			break;
		}
		case Neumann_right: { // We only need the main nodes, but we have to use phantom points in the current cell. FIXME: not implemented yet.
			break;
		}
		case Dirichlet_left: { // We only need the internal nodes.
			for (unsigned int main_node = 0; main_node != stages_x; ++main_node)
				for (unsigned int internal_node = 1; internal_node != stages_y; ++internal_node)
					U_yy[cell * nodes_per_cell + main_node + internal_node * stages_x] = internal_node_yy(U, coefficients, cell * nodes_per_cell + main_node, adjacency_information[cell * 4 + 3] * nodes_per_cell + main_node, internal_node, stages_x, stages_y, dy);
			break;
		}
		case Dirichlet_right: { // Do nothing. We don't need the derivatives for any nodes in this case.
			break;
		}
	}
}
