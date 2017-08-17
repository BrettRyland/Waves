///@file

#include <cmath> // for std::sin(), std::cos()
#include <algorithm> // for std::max
#include <numeric> // for std::inner_product
#include <limits> // for std::numeric_limits::max()
#include <cassert> // for assert()
#include "integratorCL.h"
#include "br/common/mathconstants.h"

namespace Waves {
	// Initialise the Integrator.
	void Integrator::initialise(int rx, int ry, integrator_precision dt, int ic, int bc, boost::compute::context &shared_context, boost::compute::command_queue &queue) {
		stages_x = rx - 1; // We use rx-1 here since the last node of each cell is the same as the first node of the next cell.
		stages_y = ry - 1; // Similarly for ry.
		nodes_per_cell = stages_x*stages_y;
		step_size_time = dt;
		initial_conditions = ic;
		boundary_conditions = bc;
		coords_x = setup_coords(stages_x + 1);
		coords_y = setup_coords(stages_y + 1);

		// Set up device side variables
		m_context = shared_context; // Obtain the shared context and command queue from the OpenGL widget
		m_queue = queue;
		{ // Coefficients
			[&]() {
				auto coeffs = get_coefficients(rx);
				std::vector<integrator_precision> tmp;
				for (auto c : coeffs)
					tmp.insert(tmp.end(), c.begin(), c.end());
				coefficients_x = boost::compute::vector<integrator_precision>(tmp.size(), m_context);
				boost::compute::copy(tmp.begin(), tmp.end(), coefficients_x.begin(), m_queue);
			}();
			[&]() {
				auto coeffs = get_coefficients(ry);
				std::vector<integrator_precision> tmp;
				for (auto c : coeffs)
					tmp.insert(tmp.end(), c.begin(), c.end());
				coefficients_y = boost::compute::vector<integrator_precision>(tmp.size(), m_context);
				boost::compute::copy(tmp.begin(), tmp.end(), coefficients_y.begin(), m_queue);
			}();
		}
		U = boost::compute::vector<integrator_precision>(m_context);
		V = boost::compute::vector<integrator_precision>(m_context);
		U_xx = boost::compute::vector<integrator_precision>(m_context);
		U_yy = boost::compute::vector<integrator_precision>(m_context);
		adjacency_information = boost::compute::vector<cl_uint>(m_context);
		//position_information = boost::compute::vector<integrator_precision>(m_context);
		cell_type_x = boost::compute::vector<cl_uint>(m_context);
		cell_type_y = boost::compute::vector<cl_uint>(m_context);
		generate_kernels();
		generate_update_masks();
		change_boundary_conditions(boundary_conditions); // Sets the initial conditions and copies everything to the compute device too.
	}

	// Define and compile the kernels required for the simulation
	void Integrator::generate_kernels() {
		auto step_program = boost::compute::program::create_with_source_file("kernels/integrator.cl", m_context);
		step_program.build();
		kernel_step_U = step_program.create_kernel("step_U");
		kernel_step_V = step_program.create_kernel("step_V");
		kernel_compute_U_xx = step_program.create_kernel("compute_U_xx");
		kernel_compute_U_yy = step_program.create_kernel("compute_U_yy");
	}

	/** Set the buffers and other kernel args.
	This needs to be called whenever the configuration of cells is changed.
	*/
	void Integrator::configure_kernels_static_vars(std::vector<cl_uint> const& host_cell_type_x, std::vector<cl_uint> const& host_cell_type_y, std::vector<cl_uint> const& host_adjacency_information, std::vector<integrator_precision> const& host_position_information) {
		// Allocate and copy to the device the static (i.e. don't change while the integrator is running) variables.
		adjacency_information.resize(host_adjacency_information.size(), m_queue);
		boost::compute::copy(host_adjacency_information.begin(), host_adjacency_information.end(), adjacency_information.begin(), m_queue);
		//position_information.resize(host_position_information.size(), m_queue);
		//boost::compute::copy(host_position_information.begin(), host_position_information.end(), position_information.begin(), m_queue);
		cell_type_x.resize(host_cell_type_x.size(), m_queue);
		boost::compute::copy(host_cell_type_x.begin(), host_cell_type_x.end(), cell_type_x.begin(), m_queue);
		cell_type_y.resize(host_cell_type_y.size(), m_queue);
		boost::compute::copy(host_cell_type_y.begin(), host_cell_type_y.end(), cell_type_y.begin(), m_queue);
	}

	/** Copy new initial conditions to the simulation.
	This needs to be called whenever variables related to the PDE are updated on the host side.
	*/
	void Integrator::configure_kernels_dynamic_vars(std::vector<integrator_precision> const& host_U, std::vector<integrator_precision> const& host_V) {
		// Allocate and copy to the device the dynamic (i.e. change every step of the integrator) variables.
		U.resize(host_U.size(), m_queue);
		boost::compute::copy(host_U.begin(), host_U.end(), U.begin(), m_queue);
		V.resize(host_V.size(), m_queue);
		boost::compute::copy(host_V.begin(), host_V.end(), V.begin(), m_queue);
		U_xx.resize(host_U.size(), m_queue);
		U_yy.resize(host_U.size(), m_queue);
		m_queue.finish(); // Wait for the queue to finish before changing kernel args
		kernel_step_U.set_args(U.get_buffer(), masks.get_buffer(), V.get_buffer(), cell_type_x.get_buffer(), cell_type_y.get_buffer(), cell_type_count, nodes_per_cell, step_size_time);
		kernel_step_V.set_args(V.get_buffer(), masks.get_buffer(), U.get_buffer(), U_xx.get_buffer(), U_yy.get_buffer(), cell_type_x.get_buffer(), cell_type_y.get_buffer(), cell_type_count, nodes_per_cell, step_size_time, wave_speed);
		kernel_compute_U_xx.set_args(U_xx.get_buffer(), U.get_buffer(), cell_type_x.get_buffer(), coefficients_x.get_buffer(), adjacency_information.get_buffer(), stages_x, stages_y, step_size_x);
		kernel_compute_U_yy.set_args(U_yy.get_buffer(), U.get_buffer(), cell_type_y.get_buffer(), coefficients_y.get_buffer(), adjacency_information.get_buffer(), stages_x, stages_y, step_size_y);
		update_U_dependent_variables(); // U_xx and U_yy depend on U
		update_V_dependent_variables();
	}

	void Integrator::generate_update_masks() {
		std::vector<cl_uint> host_masks(cell_type_count*cell_type_count*nodes_per_cell, 1); // Initialise each mask to true (1) and then set various nodes to false (0)
		for (cl_uint type_y = 0; type_y != cell_type_count; ++type_y)
			for (cl_uint type_x = 0; type_x != cell_type_count; ++type_x) {
				auto cell_type = (type_y*cell_type_count + type_x)*nodes_per_cell;
				switch (type_x) { // x-direction.
					case Cell_Type::Normal:
					case Cell_Type::Periodic_left:
					case Cell_Type::Periodic_right:
					case Cell_Type::Neumann_left:
						break;
					case Cell_Type::Neumann_right:
						for (auto ndx = 0, i = 0; ndx != nodes_per_cell; ++ndx, ++i)
							if (i % stages_x != 0)
								host_masks[cell_type + ndx] = 0;
						break;
					case Cell_Type::Dirichlet_left:
						for (auto ndx = 0, i = 0; ndx != nodes_per_cell; ++ndx, ++i)
							if (i % stages_x == 0)
								host_masks[cell_type + ndx] = 0;
						break;
					case Cell_Type::Dirichlet_right:
						for (auto ndx = 0; ndx != nodes_per_cell; ++ndx)
							host_masks[cell_type + ndx] = 0;
						break;
					default:
						throw std::runtime_error("Invalid Cell_Type");
				}
				switch (type_y) { // y-direction.
					case Cell_Type::Normal:
					case Cell_Type::Periodic_left:
					case Cell_Type::Periodic_right:
					case Cell_Type::Neumann_left:
						break;
					case Cell_Type::Neumann_right:
						for (auto ndx = stages_x; ndx != nodes_per_cell; ++ndx)
							host_masks[cell_type + ndx] = 0;
						break;
					case Cell_Type::Dirichlet_left:
						for (auto ndx = 0; ndx != stages_x; ++ndx)
							host_masks[cell_type + ndx] = 0;
						break;
					case Cell_Type::Dirichlet_right:
						for (auto ndx = 0; ndx != nodes_per_cell; ++ndx)
							host_masks[cell_type + ndx] = 0;
						break;
					default:
						throw std::runtime_error("Invalid Cell_Type");
				}
			}
		masks = boost::compute::vector<cl_uint>(host_masks.size(), m_context);
		boost::compute::copy(host_masks.begin(), host_masks.end(), masks.begin(), m_queue);
	}

	/** Step forwards in time using 2-stage Lobatto IIIA-IIIB in time (i.e., leapfrog) and rx-stage Lobatto IIIA-IIIB in x and y directions.
		It is explicit, local, multisymplectic and handles boundary conditions easily. See my PhD thesis for the details.
		Note: we perform an initial half-step during the initialisation in Change_Initial_Conditions() (which is called by Change_Boundary_Conditions()) so that we can do full steps in V.
	*/
	void Integrator::step()
	{
		m_queue.enqueue_1d_range_kernel(kernel_step_U, 0, num_cells, 0);
		update_U_dependent_variables();
		m_queue.enqueue_1d_range_kernel(kernel_step_V, 0, num_cells, 0);
		update_V_dependent_variables();
		time += step_size_time;
		m_queue.finish(); // wait for the queue to finish or else we queue up steps faster than we can process them
	}

	// Apply a half-step of the integrator (i.e. a step of dt/2 in V)
	void Integrator::half_step() {
		m_queue.finish(); // Wait for the queue to finish before changing kernel args
		kernel_step_V.set_arg(9, static_cast<integrator_precision>(step_size_time / 2.0));
		m_queue.enqueue_1d_range_kernel(kernel_step_V, 0, num_cells, 0);
		update_V_dependent_variables();
		m_queue.finish(); // Wait for the queue to finish before changing kernel args back
		kernel_step_V.set_arg(9, step_size_time);
	}

	// Update variables that depend on U (i.e. U_xx and U_yy).
	void Integrator::update_U_dependent_variables()
	{
		m_queue.enqueue_1d_range_kernel(kernel_compute_U_xx, 0, num_cells, 0);
		m_queue.enqueue_1d_range_kernel(kernel_compute_U_yy, 0, num_cells, 0);
	}

	// Update variables that depend on V (i.e. nothing).
	void Integrator::update_V_dependent_variables()
	{ // Nothing to do. Should be optimised away, but useful placeholder for future simulations that may use it.
	}

	/* Switch to a new set of initial conditions.
	This calls configure_kernels_dynamic_vars with the new values.
	TODO: allow reading of initial conditions from a file
	*/
	void Integrator::change_initial_conditions(int ic) { // Change the initial conditions
		if (ic == -1)
			++initial_conditions;
		else
			initial_conditions = ic;

		const integrator_precision zero = 0.0, half = 0.5, one = 1.0, two = 2.0; // Convenience values for when we change integrator_precision
		std::vector<integrator_precision> host_U(num_cells*nodes_per_cell),
			host_V(num_cells*nodes_per_cell);
		assert(host_position_information.size() == num_cells * 2);
		switch (initial_conditions) {
			case 0: { // Single frequency
				for (cl_uint c = 0; c < num_cells; ++c)
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							integrator_precision x = (host_position_information[c * 2 + 0] + coords_x[i] * step_size_x) * two / domain_scaling_factor[0]; // [-1,1)
							integrator_precision y = (host_position_information[c * 2 + 1] + coords_y[j] * step_size_y) * two / domain_scaling_factor[1]; // [-1,1)
							// u=cos(br::tau<integrator_precision>*(t+c1*x+c2*y)), which has a wave speed of sqrt(c1^2+c2^2)
							host_U[c*nodes_per_cell + j*stages_x + i] = cos(br::tau<integrator_precision>*(two*x + one*y));
							// set V to be the derivative of U to get a travelling wave train
							host_V[c*nodes_per_cell + j*stages_x + i] = -half*br::pi<integrator_precision>*br::pi<integrator_precision> * std::sin(br::tau<integrator_precision>*(two*x + one*y));
						}
				break;
			}
			case 1: { // Continuous spectrum of frequencies
				for (cl_uint c = 0; c < num_cells; ++c)
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							integrator_precision x = (host_position_information[c * 2 + 0] + coords_x[i] * step_size_x) * two / domain_scaling_factor[0];
							integrator_precision y = (host_position_information[c * 2 + 1] + coords_y[j] * step_size_y) * two / domain_scaling_factor[1];
							// this choice should be zero along the edge of the domain but non-zero everywhere else
							host_U[c*nodes_per_cell + j*stages_x + i] = static_cast<integrator_precision>(0.0625)*(std::exp(std::sin((one*x + half)*br::pi<integrator_precision>) + one) - one)*(std::exp(std::sin((one*y + half)*br::pi<integrator_precision>) + one) - one);
						}
				// start with stationary surface
				std::fill(host_V.begin(), host_V.end(), zero);
				break;
			}
			case 2: { // Continuous spectrum of frequencies, but limited to localised hump in the middle
				// u=exp(f(x,y))-1, where f(x,y)=a*(1-x^2/bx^2-y^2/by^2)
				// note: if c/bx^2 or c/by^2 is too large then the integrator will break due to a CFL condition!
				integrator_precision bx = static_cast<integrator_precision>(0.3), by = static_cast<integrator_precision>(0.3), a = static_cast<integrator_precision>(1.5), x0 = static_cast<integrator_precision>(0.0), y0 = static_cast<integrator_precision>(0.0);
				for (cl_uint c = 0; c < num_cells; ++c)
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							integrator_precision x = (host_position_information[c * 2 + 0] + coords_x[i] * step_size_x) * two / domain_scaling_factor[0];
							integrator_precision y = (host_position_information[c * 2 + 1] + coords_y[j] * step_size_y) * two / domain_scaling_factor[1];
							host_U[c*nodes_per_cell + j*stages_x + i] = std::max(zero, std::exp(a - a / bx / bx*(x - x0)*(x - x0) - a / by / by*(y - y0)*(y - y0)) - one);
						}
				// start with stationary surface
				std::fill(host_V.begin(), host_V.end(), zero);
				break;
			}
			case 3: { // Continuous spectrum of frequencies, but limited to 2 localised humps.
				// u=exp(f(x,y))-1, where f(x,y)=c*(1-x^2/bx^2-y^2/by^2)
				// note: if c/bx^2 or c/by^2 is too large then the integrator will break due to a CFL condition!
				integrator_precision bx = static_cast<integrator_precision>(0.3), by = static_cast<integrator_precision>(0.3), a = static_cast<integrator_precision>(1.5), x0 = sqrt(half)*half, y0 = x0;
				for (cl_uint c = 0; c < num_cells; ++c)
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							integrator_precision x = (host_position_information[c * 2 + 0] + coords_x[i] * step_size_x) * two / domain_scaling_factor[0];
							integrator_precision y = (host_position_information[c * 2 + 1] + coords_y[j] * step_size_y) * two / domain_scaling_factor[1];
							if (x > 0)
								host_U[c*nodes_per_cell + j*stages_x + i] = std::max(zero, std::exp(a - a / bx / bx*(x - x0)*(x - x0) - a / by / by*(y - y0)*(y - y0)) - one);
							else
								host_U[c*nodes_per_cell + j*stages_x + i] = std::max(zero, std::exp(a - a / bx / bx*(x + x0)*(x + x0) - a / by / by*(y + y0)*(y + y0)) - one);
						}
				// start with stationary surface
				std::fill(host_V.begin(), host_V.end(), zero);
				break;
			}
			case 4: { // Continuous spectrum of frequencies, but limited to an off-centre localised hump.
				// u=exp(f(x,y))-1, where f(x,y)=c*(1-x^2/bx^2-y^2/by^2)
				// note: if c/bx^2 or c/by^2 is too large then the integrator will break due to a CFL condition!
				integrator_precision bx = static_cast<integrator_precision>(0.3), by = static_cast<integrator_precision>(0.3), a = static_cast<integrator_precision>(1.5), x0 = sqrt(half)*half, y0 = x0;
				for (cl_uint c = 0; c < num_cells; ++c)
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							integrator_precision x = (host_position_information[c * 2 + 0] + coords_x[i] * step_size_x) * two / domain_scaling_factor[0];
							integrator_precision y = (host_position_information[c * 2 + 1] + coords_y[j] * step_size_y) * two / domain_scaling_factor[1];
							if (x > 0)
								host_U[c*nodes_per_cell + j*stages_x + i] = std::max(zero, std::exp(a - a / bx / bx*(x - x0)*(x - x0) - a / by / by*(y - y0)*(y - y0)) - one);
							else
								host_U[c*nodes_per_cell + j*stages_x + i] = zero;// std::max(zero, std::exp(a - a / bx / bx*(x + x0)*(x + x0) - a / by / by*(y + y0)*(y + y0)) - one);
						}
				// start with stationary surface
				std::fill(host_V.begin(), host_V.end(), zero);
				break;
			}
			case 5: { // Continuous spectrum of frequencies in the y-direction
				// u=exp(f(x,y))-1, where f(x,y)=a*(1-y^2/by^2)
				// note: if a/by^2 is too large then the integrator will break due to a CFL condition!
				integrator_precision b = static_cast<integrator_precision>(0.2), a = static_cast<integrator_precision>(1.5), y0 = static_cast<integrator_precision>(-0.8);
				for (cl_uint c = 0; c < num_cells; ++c)
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							integrator_precision x = (host_position_information[c * 2 + 0] + coords_x[i] * step_size_x) * two / domain_scaling_factor[0];
							integrator_precision y = (host_position_information[c * 2 + 1] + coords_y[j] * step_size_y) * two / domain_scaling_factor[1];
							host_U[c*nodes_per_cell + j*stages_x + i] = std::max(zero, std::exp(a - a / b / b * (y - y0) * (y - y0)) - one);
						}
				// start with stationary surface
				std::fill(host_V.begin(), host_V.end(), zero);
				break;
			}
			default: // We've reached the end or been given a wrong value, reset to the first initial condition.
				return Integrator::change_initial_conditions(0);
		}
		configure_kernels_dynamic_vars(host_U, host_V);
		time = zero;
		half_step(); // Perform a half-step to initialise the leapfrog algorithm
	}

	// Helper function to create adjacency and cell type information for non-periodic domains with Dirichlet boundaries.
	// Note: this could be rolled into a __kernel too.
	std::tuple<std::vector<cl_uint>, std::vector<cl_uint>, std::vector<cl_uint>> Integrator::find_neighbours(std::vector<integrator_precision> const& position_information)
	{
		const integrator_precision half = 0.5, one_half = 1.5; // Convenience values for when we change integrator_precision
		num_cells = static_cast<cl_uint>(position_information.size() / 2);
		std::vector<cl_uint> adjacency_information(num_cells * 4, missing_index),
			cell_type_x(num_cells),
			cell_type_y(num_cells);
		for (cl_uint c = 0; c < num_cells; ++c) {
			// Find up to 4 neighbouring cells, then assign them to top, left, right, bottom and assign Cell_Type information
			int i = c - 1;
			for (; i != -1; --i) { // look for left neighbour
				if (fabs(position_information[i * 2 + 1] - position_information[c * 2 + 1]) > half*step_size_y) // abort if we're no longer on the same row
					break;
				if ((position_information[c * 2 + 0] - position_information[i * 2 + 0] > half*step_size_x) && (position_information[c * 2 + 0] - position_information[i * 2 + 0] < one_half*step_size_x)) {
					adjacency_information[c * 4 + 0] = i;
					break;
				}
			}
			for (; i != -1; --i) { // continue looking for top neighbour
				if (position_information[c * 2 + 1] - position_information[i * 2 + 1] > one_half*step_size_y) // abort if we're more than 1 row back
					break;
				if ((fabs(position_information[i * 2 + 0] - position_information[c * 2 + 0]) < half*step_size_x) && (position_information[c * 2 + 1] - position_information[i * 2 + 1] > half*step_size_y) && (position_information[c * 2 + 1] - position_information[i * 2 + 1] < one_half*step_size_y)) {
					adjacency_information[c * 4 + 2] = i;
					break;
				}
			}
			i = c + 1;
			for (; i != num_cells; ++i) { // look for right neighbour
				if (fabs(position_information[i * 2 + 1] - position_information[c * 2 + 1]) > half*step_size_y) // abort if we're no longer on the same row
					break;
				if ((position_information[i * 2 + 0] - position_information[c * 2 + 0] > half*step_size_x) && (position_information[i * 2 + 0] - position_information[c * 2 + 0] < one_half*step_size_x)) {
					adjacency_information[c * 4 + 1] = i;
					break;
				}
			}
			for (; i != num_cells; ++i) { // continue looking for bottom neighbour
				if (position_information[i * 2 + 1] - position_information[c * 2 + 1] > one_half*step_size_y) // abort if we're more than 1 row ahead
					break;
				if ((fabs(position_information[i * 2 + 0] - position_information[c * 2 + 0]) < half*step_size_x) && (position_information[i * 2 + 1] - position_information[c * 2 + 1] > half*step_size_y) && (position_information[c * 2 + 1] - position_information[i * 2 + 1] < one_half*step_size_y)) {
					adjacency_information[c * 4 + 3] = i;
					break;
				}
			}

			// FIXME: If a cell is both Dirchlet_right and Dirichlet_left, then we should remove it and change the neighbour (in the other direction) appropriately.
			if (adjacency_information[c * 4 + 1] == missing_index)
				cell_type_x[c] = Cell_Type::Dirichlet_right;
			else if (adjacency_information[c * 4 + 0] == missing_index)
				cell_type_x[c] = Cell_Type::Dirichlet_left;
			else
				cell_type_x[c] = Cell_Type::Normal;

			if (adjacency_information[c * 4 + 3] == missing_index)
				cell_type_y[c] = Cell_Type::Dirichlet_right;
			else if (adjacency_information[c * 4 + 2] == missing_index)
				cell_type_y[c] = Cell_Type::Dirichlet_left;
			else
				cell_type_y[c] = Cell_Type::Normal;
		}
		return std::make_tuple(adjacency_information, cell_type_x, cell_type_y);
	}

	/* Switch to a new set of boundary conditions.
	This calls configure_kernels_static_vars with the new configuration.
	TODO: allow reading of boundary conditions from a file
	*/
	void Integrator::change_boundary_conditions(int bc)
	{ // Change the boundary conditons
		if (bc == -1)
			++boundary_conditions;
		else
			boundary_conditions = bc;

		const integrator_precision half = 0.5, one = 1.0, two = 2.0; // Convenience values for when we change integrator_precision
		//std::vector<cl_uint> host_cell_type_x, host_cell_type_y;
		host_position_information.clear();
		switch (boundary_conditions) {
			case 0: { // Periodic square boundary.
				/*
				200 * 200 square domain centered at the origin with periodic boundary conditions
				dx = dy = 0.1
				wave_speed = 1^2
				FIXME: for consistency, we should label some of these Periodic_left, though this has no effect on anything except in constructing the adjacency information
				*/
				host_cell_type_x.clear();
				host_cell_type_y.clear();
				host_adjacency_information.clear();
				cl_uint n = 200;
				step_size_x = static_cast<integrator_precision>(0.1);
				step_size_y = static_cast<integrator_precision>(0.1);
				wave_speed = static_cast<integrator_precision>(10.0);
				domain_scaling_factor = { n*step_size_x,n*step_size_y };
				for (cl_uint j = 0; j != n - 1; ++j) {
					for (cl_uint i = 0; i != n - 1; ++i) {
						host_cell_type_x.emplace_back(Cell_Type::Normal);
						host_cell_type_y.emplace_back(Cell_Type::Normal);
						host_adjacency_information.insert(host_adjacency_information.end(), { n * j + (i + n - 1) % n, n * j + i + 1, n * ((j + n - 1) % n) + i, n * (j + 1) + i });
						host_position_information.insert(host_position_information.end(), { step_size_x * (static_cast<integrator_precision>(i) - static_cast<integrator_precision>(n / 2)), step_size_y * (static_cast<integrator_precision>(j) - static_cast<integrator_precision>(n / 2)) });
					}
					// last column
					host_cell_type_x.emplace_back(Cell_Type::Periodic_right);
					host_cell_type_y.emplace_back(Cell_Type::Normal);
					host_adjacency_information.insert(host_adjacency_information.end(), { n * j + n - 2, n * j, n * ((j + n - 1) % n) + n - 1, n * (j + 1) + n - 1 });
					host_position_information.insert(host_position_information.end(), { step_size_x * (static_cast<integrator_precision>(n - 1) - static_cast<integrator_precision>(n / 2)), step_size_y * (static_cast<integrator_precision>(j) - static_cast<integrator_precision>(n / 2)) });
				}
				// last row
				for (cl_uint i = 0; i != n - 1; ++i) {
					host_cell_type_x.emplace_back(Cell_Type::Normal);
					host_cell_type_y.emplace_back(Cell_Type::Periodic_right);
					host_adjacency_information.insert(host_adjacency_information.end(), { n * (n - 1) + (i + n - 1) % n, n * (n - 1) + i + 1, n * (n - 2) + i, i });
					host_position_information.insert(host_position_information.end(), { step_size_x * (static_cast<integrator_precision>(i) - static_cast<integrator_precision>(n / 2)), step_size_y * (static_cast<integrator_precision>(n - 1) - static_cast<integrator_precision>(n / 2)) });
				}
				// last corner
				host_cell_type_x.emplace_back(Cell_Type::Periodic_right);
				host_cell_type_y.emplace_back(Cell_Type::Periodic_right);
				host_adjacency_information.insert(host_adjacency_information.end(), { n * n - 2, n * n - n, n * n - 1 - n, n - 1 });
				host_position_information.insert(host_position_information.end(), { step_size_x * (static_cast<integrator_precision>(n - 1) - static_cast<integrator_precision>(n / 2)), step_size_y * (static_cast<integrator_precision>(n - 1) - static_cast<integrator_precision>(n / 2)) });
				break;
			}
			case 1: { // Square with Dirichlet boundaries.
				/*
				200 * 200 square domain with Dirichlet boundary conditions
				dx = dy = 0.1
				wave_speed = 1^2
				*/
				int n = 200;
				step_size_x = static_cast<integrator_precision>(0.1);
				step_size_y = static_cast<integrator_precision>(0.1);
				wave_speed = static_cast<integrator_precision>(10.0);
				domain_scaling_factor = { n*step_size_x, n*step_size_y };
				for (int j = 0; j != n; ++j)
					for (int i = 0; i != n; ++i) {
						auto x = step_size_x * (i - n / 2);
						auto y = step_size_y * (j - n / 2);
						host_position_information.insert(host_position_information.end(), { x, y });
					}
				std::tie(host_adjacency_information, host_cell_type_x, host_cell_type_y) = find_neighbours(host_position_information);
				break;
			}
			case 2: { // Circle with Dirichlet boundaries.
				/*
				Circle of radius n/2 centred at the origin.
				dx = dy = 0.1
				wave_speed = 1^2
				*/

				int n = 200;
				step_size_x = static_cast<integrator_precision>(0.1);
				step_size_y = static_cast<integrator_precision>(0.1);
				wave_speed = static_cast<integrator_precision>(10.0);
				domain_scaling_factor = { n*step_size_x,n*step_size_y };
				integrator_precision r = static_cast<integrator_precision>(0.49)*domain_scaling_factor[0];
				for (int j = 0; j != n; ++j)
					for (int i = 0; i != n; ++i) {
						auto x = step_size_x * (i - n / 2);
						auto y = step_size_y * (j - n / 2);
						if (x*x + y*y < r*r)
							host_position_information.insert(host_position_information.end(), { x, y });
					}
				std::tie(host_adjacency_information, host_cell_type_x, host_cell_type_y) = find_neighbours(host_position_information);
				break;
			}
			case 3: { // Circle with a cusp, Dirichlet boundaries
				/*
				As with case 2, but with a cusp.
				dx = dy = 0.1
				wave_speed = 1^2
				*/

				int n = 1303; // This gives a million element simulation.
				step_size_x = static_cast<integrator_precision>(0.023);
				step_size_y = static_cast<integrator_precision>(0.023);
				wave_speed = static_cast<integrator_precision>(2.0);
				domain_scaling_factor = { n*step_size_x, n*step_size_y };
				integrator_precision r = half*domain_scaling_factor[0], r2 = half*r, s = std::sqrt(half)*half*r;
				for (int j = 0; j < n; ++j)
					for (int i = 0; i < n; ++i) {
						auto x = step_size_x * (i - n / 2);
						auto y = step_size_y * (j - n / 2);
						if ((x*x + y*y < r*r) && ((y > x) || (((y - s)*(y - s) + (x - s)*(x - s) < r2*r2) || ((y + s)*(y + s) + (x + s)*(x + s) < r2*r2))))
							host_position_information.insert(host_position_information.end(), { x, y });
					}
				std::tie(host_adjacency_information, host_cell_type_x, host_cell_type_y) = find_neighbours(host_position_information);
				break;
			}
			case 4: { // 2 intersecting circles, Dirichlet boundaries
				int n = 300;
				step_size_x = static_cast<integrator_precision>(0.1);
				step_size_y = static_cast<integrator_precision>(0.1);
				wave_speed = static_cast<integrator_precision>(10.0);
				domain_scaling_factor = { n*step_size_x,n*step_size_y };
				integrator_precision r = static_cast<integrator_precision>(0.26)*domain_scaling_factor[0], s = std::sqrt(half)*half*half*domain_scaling_factor[0];
				for (int j = 0; j < n; ++j)
					for (int i = 0; i < n; ++i) {
						auto x = step_size_x * (i - n / 2);
						auto y = step_size_y * (j - n / 2);
						if (((y - s)*(y - s) + (x - s)*(x - s) < r*r) || ((y + s)*(y + s) + (x + s)*(x + s) < r*r))
							host_position_information.insert(host_position_information.end(), { x, y });
					}
				std::tie(host_adjacency_information, host_cell_type_x, host_cell_type_y) = find_neighbours(host_position_information);
				break;
			}
			case 5: { // Double slit in a square with some periodic and some Dirichlet boundaries
				int n{ 1000 }, slit_position_y{ n / 5 }, slit_position_x{ n / 10 }, slit_width_x{ n / 100 }, slit_width_y{ 1 };
				step_size_x = static_cast<integrator_precision>(0.025);
				step_size_y = static_cast<integrator_precision>(0.025);
				wave_speed = static_cast<integrator_precision>(0.3);
				domain_scaling_factor = { n*step_size_x, n*step_size_y };
				// Source chamber
				for (int j = 0; j != slit_position_y; ++j)
					for (int i = 0; i != n; ++i) {
						auto x = step_size_x * (i - n / 2);
						auto y = step_size_y * (j - n / 2);
						host_position_information.insert(host_position_information.end(), { x, y });
					}
				// Slits
				for (int j = slit_position_y; j != slit_position_y + slit_width_y; ++j) {
					for (int i = n / 2 - slit_position_x - slit_width_x; i != n / 2 - slit_position_x + slit_width_x; ++i) {
						auto x = step_size_x * (i - n / 2);
						auto y = step_size_y * (j - n / 2);
						host_position_information.insert(host_position_information.end(), { x, y });
					}
					for (int i = n / 2 + slit_position_x - slit_width_x; i != n / 2 + slit_position_x + slit_width_x; ++i) {
						auto x = step_size_x * (i - n / 2);
						auto y = step_size_y * (j - n / 2);
						host_position_information.insert(host_position_information.end(), { x, y });
					}
				}
				// Free chamber
				for (int j = slit_position_y + slit_width_y; j != n; ++j)
					for (int i = 0; i != n; ++i) {
						auto x = step_size_x * (i - n / 2);
						auto y = step_size_y * (j - n / 2);
						host_position_information.insert(host_position_information.end(), { x, y });
					}
				std::tie(host_adjacency_information, host_cell_type_x, host_cell_type_y) = find_neighbours(host_position_information);

				// Correct left and right boundaries to be periodic
				for (cl_uint cell = 0; cell != host_cell_type_x.size(); ++cell) {
					if ((host_cell_type_x[cell] == Cell_Type::Dirichlet_left) && (host_position_information[cell * 2] < step_size_x*(2.0 - n / 2))) {
						host_cell_type_x[cell] = Cell_Type::Periodic_left;
						host_adjacency_information[cell * 4 + 0] = cell + n - 1;
					}
					else if ((host_cell_type_x[cell] == Cell_Type::Dirichlet_right) && (host_position_information[cell * 2] > step_size_x*(n / 2 - 2.0))) {
						host_cell_type_x[cell] = Cell_Type::Periodic_right;
						host_adjacency_information[cell * 4 + 1] = cell + 1 - n;
					}
				}
				break;
			}
			default: // We've reached the end or been given a wrong value, reset to the first boundary condition.
				return Integrator::change_boundary_conditions(0);
		}
		// Record the number of cells in the simulation
		num_cells = static_cast<cl_uint>(host_cell_type_x.size());

		// Transfer the new configuration to the compute device.
		configure_kernels_static_vars(host_cell_type_x, host_cell_type_y, host_adjacency_information, host_position_information);

		// Report on new setup
		std::cout << *this << std::endl;
		return change_initial_conditions(initial_conditions); // Reset the initial conditions for the new boundary conditions
	}

	// Set the coefficients of the Lobatto IIIA-IIIB stencils.
	// For r stages, the first 2*r-1 values are the coefficients of the main node and the remaining (r-2)*r values are the coefficients of the internal nodes.
	std::vector<std::vector<integrator_precision>> Integrator::get_coefficients(int stages)
	{
		std::vector<std::vector<integrator_precision>> coefficients(stages - 1); // The last node in a cell is the first node of the next cell, so we drop it.
		switch (stages) {
			case 2:
				// Main node
				coefficients[0].emplace_back(static_cast<integrator_precision>(1.0));
				coefficients[0].emplace_back(static_cast<integrator_precision>(-2.0));
				coefficients[0].emplace_back(static_cast<integrator_precision>(1.0));
				break;

			case 3:
				// Main node
				coefficients[0].emplace_back(static_cast<integrator_precision>(-1.0));
				coefficients[0].emplace_back(static_cast<integrator_precision>(8.0));
				coefficients[0].emplace_back(static_cast<integrator_precision>(-14.0));
				coefficients[0].emplace_back(static_cast<integrator_precision>(8.0));
				coefficients[0].emplace_back(static_cast<integrator_precision>(-1.0));

				// Internal node
				coefficients[1].emplace_back(static_cast<integrator_precision>(4.0));
				coefficients[1].emplace_back(static_cast<integrator_precision>(-8.0));
				coefficients[1].emplace_back(static_cast<integrator_precision>(4.0));
				break;

			case 4:
				// Main node
				coefficients[0].emplace_back(static_cast<integrator_precision>(1.0));
				coefficients[0].emplace_back(static_cast<integrator_precision>(0.5*(25.0 - 15.0*sqrt(5.0))));
				coefficients[0].emplace_back(static_cast<integrator_precision>(0.5*(25.0 + 15.0*sqrt(5.0))));
				coefficients[0].emplace_back(static_cast<integrator_precision>(-52.0));
				coefficients[0].emplace_back(static_cast<integrator_precision>(0.5*(25.0 + 15.0*sqrt(5.0))));
				coefficients[0].emplace_back(static_cast<integrator_precision>(0.5*(25.0 - 15.0*sqrt(5.0))));
				coefficients[0].emplace_back(static_cast<integrator_precision>(1.0));

				// First internal node
				coefficients[1].emplace_back(static_cast<integrator_precision>(5.0 + 3.0*sqrt(5.0)));
				coefficients[1].emplace_back(static_cast<integrator_precision>(-20.0));
				coefficients[1].emplace_back(static_cast<integrator_precision>(10.0));
				coefficients[1].emplace_back(static_cast<integrator_precision>(5.0 - 3.0*sqrt(5.0)));

				// Second internal node
				coefficients[2].emplace_back(static_cast<integrator_precision>(5.0 - 3.0*sqrt(5.0)));
				coefficients[2].emplace_back(static_cast<integrator_precision>(10.0));
				coefficients[2].emplace_back(static_cast<integrator_precision>(-20.0));
				coefficients[2].emplace_back(static_cast<integrator_precision>(5.0 + 3.0*sqrt(5.0)));
				break;

			default: // FIXME: Calculate the higher order coefficients or set up a function to generate them.
				throw std::runtime_error(stages + " stages not implemented yet.");
		}
		return coefficients;
	}

	// Get the coordinates of the nodes within a unit cell.
	std::vector<integrator_precision> Integrator::setup_coords(int stages)
	{
		std::vector<integrator_precision> coords;
		coords.reserve(stages);
		switch (stages)
		{
			case 2:
				coords.emplace_back(static_cast<integrator_precision>(0.0));
				coords.emplace_back(static_cast<integrator_precision>(1.0));
				break;
			case 3:
				coords.emplace_back(static_cast<integrator_precision>(0.0));
				coords.emplace_back(static_cast<integrator_precision>(0.5));
				coords.emplace_back(static_cast<integrator_precision>(1.0));
				break;
			case 4:
				coords.emplace_back(static_cast<integrator_precision>(0.0));
				coords.emplace_back(static_cast<integrator_precision>(0.5 - std::sqrt(5.0) / 10.0));
				coords.emplace_back(static_cast<integrator_precision>(0.5 + std::sqrt(5.0) / 10.0));
				coords.emplace_back(static_cast<integrator_precision>(1.0));
				break;
			default: // FIXME: Calculate the higher order coordinates or set up a function to generate them.
				throw std::runtime_error(stages + " stages not implemented yet.");
		}
		return coords;
	}

	//auto Integrator::get_position_information()->std::vector<integrator_precision> {
	//	std::vector<integrator_precision> p(2 * num_cells);
	//	boost::compute::copy(position_information.begin(), position_information.end(), p.begin(), m_queue);
	//	return p;
	//}
	auto Integrator::get_adjacency_information()->std::vector<cl_uint> {
		std::vector<cl_uint> a(4 * num_cells);
		boost::compute::copy(adjacency_information.begin(), adjacency_information.end(), a.begin(), m_queue);
		return a;
	}
	auto Integrator::get_cell_type_x()->std::vector<cl_uint> {
		std::vector<cl_uint> ctx(num_cells);
		boost::compute::copy(cell_type_x.begin(), cell_type_x.end(), ctx.begin(), m_queue);
		return ctx;
	}
	auto Integrator::get_cell_type_y()->std::vector<cl_uint> {
		std::vector<cl_uint> cty(num_cells);
		boost::compute::copy(cell_type_y.begin(), cell_type_y.end(), cty.begin(), m_queue);
		return cty;
	}

	/// Friendly ostream operator<<
	std::ostream& operator<<(std::ostream& os, const Integrator& integrator)
	{
		os << "A " << integrator.num_cells << " element simulation using Lobatto IIIA-IIIB discretisation in space with " << integrator.stages_x + 1 << " stages in x, " << integrator.stages_y + 1 << " stages in y (giving " << integrator.stages_x*integrator.stages_y << " nodes per element), and stepsizes " << "dx=" << integrator.step_size_x << ", dy=" << integrator.step_size_y << " and dt=" << integrator.step_size_time << ".";
		os << "\nRunning on the device: " << integrator.m_context.get_device().name() << "\n";
		return os;
	}
}
