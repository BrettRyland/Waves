///@file

#include <cmath> // for std::sin(), std::cos()
#include <algorithm> // for std::max
#include <numeric> // for std::inner_product
#include <limits> // for std::numeric_limits::max()
#include <cassert> // for assert()
#include "integrator.h"

namespace Waves {
	/// Our global instance of the Integrator. This needs to be global for OpenGL to access it.
	Integrator g_waves;

	// Initialise the Integrator.
	// We return a hint for the z-scaling factor for use later in the renderer.
	float Integrator::Initialise(int rx, int ry, double dt, int ic, int bc) {
		stages_x = rx - 1; // We use rx-1 here since the last node of each cell is the same as the first node of the next cell.
		stages_y = ry - 1; // Similarly for ry.
		step_size_time = dt;
		initial_conditions = ic;
		boundary_conditions = bc;
		coefficients_x = Setup_Coefficients(stages_x + 1);
		coefficients_y = Setup_Coefficients(stages_y + 1);
		coords_x = Setup_Coords(stages_x + 1);
		coords_y = Setup_Coords(stages_y + 1);
		return Change_Boundary_Conditions(boundary_conditions); // Sets the initial conditions too.
	}

	/** Generate a true/false mask for a given cell for determining which nodes in the cell need updating.
	\param[out] mask is the mask
	\param[in] cell is the cell to examine
	\param[in] stages_x is the number of stages in the x direction
	\param[in] stages_y is the number of stages in the y direction
	*/
	void Generate_Update_Mask(std::vector<bool>&  mask, const Cell&  cell, unsigned int stages_x, unsigned int stages_y)
	{
		auto it = mask.begin();
		int i;
		switch (cell.cell_type.first) { // x-direction.
			case Cell_Type::Normal:
			case Cell_Type::Periodic:
			case Cell_Type::Neumann_left:
				mask.assign(stages_x*stages_y, true);
				break;
			case Cell_Type::Neumann_right:
				for (it = mask.begin(), i = 0; it != mask.end(); ++it, ++i)
					if (i % stages_x != 0)
						*it = false;
					else
						*it = true;
				break;
			case Cell_Type::Dirichlet_left:
				for (it = mask.begin(), i = 0; it != mask.end(); ++it, ++i)
					if (i % stages_x != 0)
						*it = true;
					else
						*it = false;
				break;
			case Cell_Type::Dirichlet_right:
				mask.assign(stages_x*stages_y, false);
				break;
			default:
				throw std::runtime_error("Invalid cell type.");
		}
		switch (cell.cell_type.second) { // y-direction. Here we only set extra false values.
			case Cell_Type::Normal:
			case Cell_Type::Periodic:
			case Cell_Type::Neumann_left:
				break;
			case Cell_Type::Neumann_right:
				for (it = mask.begin() + stages_x; it != mask.end(); ++it)
					*it = false;
				break;
			case Cell_Type::Dirichlet_left:
				for (it = mask.begin(); it != mask.begin() + stages_x; ++it)
					*it = false;
				break;
			case Cell_Type::Dirichlet_right:
				mask.assign(stages_x*stages_y, false);
				break;
			default:
				throw std::runtime_error("Invalid cell type.");
		}
	}

	// Part of the integrator. Perform a step of size step_size in the V variables.
	void Integrator::Step_V(double step_size)
	{
		auto derivative_of_the_potential = [](double i) {return 2 * i + 4 * i*i*i; }; // V'(u)=2*u+4*u^3
#pragma omp parallel for
		for (auto c = 0; c < cells.size(); ++c) {
			std::vector<bool> mask(stages_x*stages_y, true);
			//for (auto& cell : cells) {
			Generate_Update_Mask(mask, cells[c], stages_x, stages_y);
			auto it_mask = mask.begin();
			auto it_V = cells[c].V.begin();
			auto it_U = cells[c].U.begin();
			auto it_U_xx = cells[c].U_xx.begin();
			auto it_U_yy = cells[c].U_yy.begin();
			for (; it_V != cells[c].V.end(); ++it_V, ++it_U_xx, ++it_U_yy, ++it_U, ++it_mask)
				if (*it_mask)
					*it_V += step_size * (wave_speed * (*it_U_xx + *it_U_yy) - derivative_of_the_potential(*it_U));
		}
		update_V_dependent_variables();
	}

	// Part of the integrator. Perform a step of size step_size in the U variables.
	void Integrator::Step_U(double step_size)
	{
#pragma omp parallel for
		for (auto c = 0; c < cells.size(); ++c) {
			std::vector<bool> mask(stages_x*stages_y, true);
			//for (auto& cell : cells) {
			Generate_Update_Mask(mask, cells[c], stages_x, stages_y);
			auto it_mask = mask.begin();
			auto it_U = cells[c].U.begin();
			auto it_V = cells[c].V.begin();
			for (; it_U != cells[c].U.end(); ++it_U, ++it_V)
				if (*it_mask)
					*it_U += step_size * (*it_V);
		}
		update_U_dependent_variables();
	}

	// Apply one step of the integrator.
	void Integrator::Step()
	{ // The integrator.
		// Step forwards in time using 2-stage Lobatto IIIA-IIIB in time (i.e., leapfrog) and rx-stage Lobatto IIIA-IIIB in x and y directions.
		// It is explicit, local, multisymplectic and handles boundary conditions easily. See my PhD thesis for the details.
		// Note: we perform an initial half-step during the initialisation in Change_Initial_Conditions() (which is called by Change_Boundary_Conditions()) so that we can do full steps in V.

		//Step_V(0.5*step_size_time);
		Step_U(step_size_time);
		//Step_V(0.5*step_size_time);
		Step_V(step_size_time);
		Time += step_size_time;
	}

	// Compute the second order derivative of U in the x-direction using the Lobatto IIIA-IIIB stencils.
	void Integrator::Compute_U_xx()
	{
		// Function to apply the first coefficient stencil to the local cell centred on the main node (includes all values in cell to the left and first value in cell to right).
		auto main_node_xx = [this](std::vector<double>::const_iterator it_U, std::vector<double>::const_iterator it_U_left, std::vector<double>::const_iterator it_U_right, std::vector<double>::const_iterator it_coefs) {
			return std::inner_product(it_U_left, it_U_left + stages_x, it_coefs,
				std::inner_product(it_U, it_U + stages_x, it_coefs + stages_x, // *it_coefs has size 2*stages_x+1, so this selects the middle value
					*it_U_right * *(it_coefs + 2 * stages_x))) / step_size_x / step_size_x;
		};
		// Funtion to apply the remaining coefficient stencils to the local cell (includes first value in cell to the right).
		auto internal_node_xx = [this](std::vector<double>::const_iterator it_U, std::vector<double>::const_iterator it_U_right, std::vector<double>::const_iterator it_coefs) {
			return std::inner_product(it_U, it_U + stages_x, it_coefs, *it_U_right * *(it_coefs + stages_x)) / step_size_x / step_size_x; // *it_coefs has size stages_x+1, so this selects the last value
		};

		//#pragma omp parallel for
		for (int c = 0; c < cells.size(); ++c) {
			// Apply the appropriate stencils based on the cell type.
			switch (cells[c].cell_type.first) {
				case Cell_Type::Normal:
				case Cell_Type::Periodic:
					assert((adjacency_information[c][0] != missing_index) && (adjacency_information[c][1] != missing_index));
					for (auto offset = 0; offset != stages_x*stages_y; offset += stages_x) {
						auto it_coefs = coefficients_x.begin();
						cells[c].U_xx[offset] = main_node_xx(cells[c].U.begin() + offset, cells[adjacency_information[c][0]].U.begin() + offset, cells[adjacency_information[c][1]].U.begin() + offset, (*it_coefs).begin());
						for (++it_coefs; it_coefs != coefficients_x.end(); ++it_coefs) {
							auto offset2 = offset + std::distance(coefficients_x.begin(), it_coefs);
							cells[c].U_xx[offset2] = internal_node_xx(cells[c].U.begin() + offset, cells[adjacency_information[c][1]].U.begin() + offset, (*it_coefs).begin());
						}
					}
					break;
				case Cell_Type::Neumann_left:
					// We need all nodes, but we have to use phantom points in the left cell for the main node.
					throw "Function not implemented yet.";
					assert(adjacency_information[c][1] != missing_index);
					for (auto offset = 0; offset != stages_x*stages_y; offset += stages_x) {
						auto it_coefs = coefficients_x.begin();
						//main_node_xx((*it).U_xx.begin() + offset, (*it).U.begin() + offset, (*x_left).U.begin() + offset, (*x_right).U.begin() + offset, (*it_coefs).begin());
						for (++it_coefs; it_coefs != coefficients_x.end(); ++it_coefs) {
							auto offset2 = offset + std::distance(coefficients_x.begin(), it_coefs);
							cells[c].U_xx[offset2] = internal_node_xx(cells[c].U.begin() + offset, cells[adjacency_information[c][1]].U.begin() + offset, (*it_coefs).begin());
						}
					}
					break;
				case Cell_Type::Neumann_right:
					// We only need the main nodes, but we have to use phantom points in the current cell.
					throw "Function not implemented yet.";
					assert(adjacency_information[c][0] != missing_index);
					//for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_y*stages_x; offset += stages_x) {
					//	auto it_coefs = coefficients_x.begin();
					//	//main_node_xx((*it).U_xx.begin() + offset, (*it).U.begin() + offset, (*x_left).U.begin() + offset, (*x_right).U.begin() + offset, (*it_coefs).begin());
					//}
					break;
				case Cell_Type::Dirichlet_left:
					// We only need the internal nodes.
					assert(adjacency_information[c][1] != missing_index);
					for (auto offset = 0; offset != stages_x*stages_y; offset += stages_x) {
						auto it_coefs = coefficients_x.begin() + 1;
						for (; it_coefs != coefficients_x.end(); ++it_coefs) {
							auto offset2 = offset + std::distance(coefficients_x.begin(), it_coefs);
							cells[c].U_xx[offset2] = internal_node_xx(cells[c].U.begin() + offset, cells[adjacency_information[c][1]].U.begin() + offset, (*it_coefs).begin());
						}
					}
					break;
				case Cell_Type::Dirichlet_right:
					// Do nothing. We don't need the derivatives for any nodes in this case.
					break;
				default:
					throw std::runtime_error("Invalid cell type.");
			}
		}
	}

	// Compute the second order derivative of U in the y-direction using the Lobatto IIIA-IIIB stencils.
	void Integrator::Compute_U_yy()
	{
		// Unfortunately, we can't use std::inner_product here as the relevant elements of U, etc., are not adjacent. So, here's a version that has non-unity step size in the first argument, but with no bounds checking!
		auto inner_product_n = [](std::vector<double>::const_iterator input1_start, std::vector<double>::const_iterator input2_start, size_t n, std::iterator_traits<std::vector<double>::const_iterator>::difference_type adv, double value) {
			auto ret = value;
			for (size_t i = 0; i != n - 1; ++i, std::advance(input1_start, adv), std::advance(input2_start, 1))
				ret += *input1_start * *input2_start;
			return ret + *input1_start * *input2_start;
		};

		// Function to apply the first coefficient stencil to the local cell centred on the main node (includes all values in cell to the left and first value in cell to right).
		auto main_node_yy = [this,& inner_product_n](std::vector<double>::const_iterator it_U, std::vector<double>::const_iterator it_U_left, std::vector<double>::const_iterator it_U_right, std::vector<double>::const_iterator it_coefs) {
			return inner_product_n(it_U_left, it_coefs, stages_y, stages_x,
				inner_product_n(it_U, it_coefs + stages_y, stages_y, stages_x, // *it_coefs has size 2*stages_x+1, so this selects the middle value
					*it_U_right * *(it_coefs + 2 * stages_y))) / step_size_y / step_size_y;
		};

		// Function to apply the remaining coefficient stencils to the local cell (includes first value in cell to the right).
		auto internal_node_yy = [this,& inner_product_n](std::vector<double>::const_iterator it_U, std::vector<double>::const_iterator it_U_right, std::vector<double>::const_iterator it_coefs) {
			return inner_product_n(it_U, it_coefs, stages_y, stages_x, *it_U_right * *(it_coefs + stages_y)) / step_size_y / step_size_y; // *it_coefs has size stages_x+1, so this selects the last value
		};

		//#pragma omp parallel for
		for (int c = 0; c < cells.size(); ++c) {
			// Apply the appropriate stencils based on the cell type.
			switch (cells[c].cell_type.second) {
				case Cell_Type::Normal:
				case Cell_Type::Periodic:
					assert((adjacency_information[c][2] != missing_index) && (adjacency_information[c][3] != missing_index));
					for (auto offset = 0; offset != stages_x; ++offset) {
						auto it_coefs = coefficients_y.begin();
						cells[c].U_yy[offset] = main_node_yy(cells[c].U.begin() + offset, cells[adjacency_information[c][2]].U.begin() + offset, cells[adjacency_information[c][3]].U.begin() + offset, (*it_coefs).begin());
						for (++it_coefs; it_coefs != coefficients_y.end(); ++it_coefs) {
							auto offset2 = offset + std::distance(coefficients_y.begin(), it_coefs)*stages_x;
							cells[c].U_yy[offset2] = internal_node_yy(cells[c].U.begin() + offset, cells[adjacency_information[c][3]].U.begin() + offset, (*it_coefs).begin());
						}
					}
					break;
				case Cell_Type::Neumann_left:
					// We need all nodes, but we have to use phantom points in the left cell for the main node.
					throw "Function not implemented yet.";
					assert(adjacency_information[c][3] != missing_index);
					for (auto offset = 0; offset != stages_x; ++offset) {
						auto it_coefs = coefficients_y.begin();
						//main_node_xx((*it).U_yy.begin() + offset, (*it).U.begin() + offset, (*y_left).U.begin() + offset, (*y_right).U.begin() + offset, (*it_coefs).begin());
						for (++it_coefs; it_coefs != coefficients_y.end(); ++it_coefs) {
							auto offset2 = offset + std::distance(coefficients_y.begin(), it_coefs)*stages_x;
							cells[c].U_yy[offset2] = internal_node_yy(cells[c].U.begin() + offset, cells[adjacency_information[c][3]].U.begin() + offset, (*it_coefs).begin());
						}
					}
					break;
				case Cell_Type::Neumann_right:
					// We only need the main nodes, but we have to use phantom points in the current cell.
					throw "Function not implemented yet.";
					for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_x; ++offset) {
						auto it_coefs = coefficients_y.begin();
						//main_node_yy((*it).Uyy.begin() + offset, (*it).U.begin() + offset, (*x_left).U.begin() + offset, (*x_right).U.begin() + offset, (*it_coefs).begin());
					}
					break;
				case Cell_Type::Dirichlet_left:
					// We only need the internal nodes.
					assert(adjacency_information[c][3] != missing_index);
					for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_x; ++offset) {
						auto it_coefs = coefficients_y.begin() + 1;
						for (; it_coefs != coefficients_y.end(); ++it_coefs) {
							auto offset2 = offset + std::distance(coefficients_y.begin(), it_coefs)*stages_x;
							cells[c].U_yy[offset2] = internal_node_yy(cells[c].U.begin() + offset, cells[adjacency_information[c][3]].U.begin() + offset, (*it_coefs).begin());
						}
					}
					break;
				case Cell_Type::Dirichlet_right:
					// Do nothing. We don't need the derivatives for any nodes in this case.
					break;
				default:
					throw std::runtime_error("Invalid cell type.");
			}
		}
	}

	// Update variables that depend on U (i.e. U_xx and U_yy).
	void Integrator::update_U_dependent_variables()
	{
		Compute_U_xx();
		Compute_U_yy();
	}

	// Update variables that depend on V (i.e. nothing).
	void Integrator::update_V_dependent_variables()
	{ // Nothing to do. Should be optimised away, but useful placeholder for future simulations that may use it.
	}

	// Switch to a new set of initial conditions.
	// TODO: allow reading of initial conditions from a file
	float Integrator::Change_Initial_Conditions(int ic)
	{ // Change the initial conditions
		float hint = 1.0f; // Hint for the z-scaling factor in the renderer.
		if (ic == -1)
			++initial_conditions;
		else
			initial_conditions = ic;

		// Re-initialise the cells.
		for (auto& cell : cells) {
			cell.U.resize(stages_x*stages_y);
			cell.V.resize(stages_x*stages_y);
			cell.U_xx.resize(stages_x*stages_y);
			cell.U_yy.resize(stages_x*stages_y);
			// minimise memory usage
			cell.U.shrink_to_fit();
			cell.V.shrink_to_fit();
			cell.U_xx.shrink_to_fit();
			cell.U_yy.shrink_to_fit();
		}

		switch (initial_conditions) {
			case 0: // Single frequency
#pragma omp parallel for
				for (int c = 0; c < position_information.size(); ++c) {
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							auto x = (position_information[c][0] + coords_x[i] * step_size_x) * 2.0 / domain_scaling_factor[0]; // [-1,1)
							auto y = (position_information[c][1] + coords_y[j] * step_size_y) * 2.0 / domain_scaling_factor[1]; // [-1,1)
							// u=cos(2*PI*(t+c1*x+c2*y)), which has a wave speed of sqrt(c1^2+c2^2)
							cells[c].U[j*stages_x + i] = cos(2.0*PI*(2.0*x + 1.0*y));
							// set V to be the derivative of U to get a travelling wave train
							cells[c].V[j*stages_x + i] = -0.5*PI*PI * std::sin(2.0*PI*(2.0*x + 1.0*y));
						}
				}
				hint = 0.25f;
				break;
			case 1: // Continuous spectrum of frequencies
#pragma omp parallel for
				for (int c = 0; c < position_information.size(); ++c) {
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							auto x = (position_information[c][0] + coords_x[i] * step_size_x) * 2.0 / domain_scaling_factor[0];
							auto y = (position_information[c][1] + coords_y[j] * step_size_y) * 2.0 / domain_scaling_factor[1];
							// this choice should be zero along the edge of the domain but non-zero everywhere else
							cells[c].U[j*stages_x + i] = 0.0625*(std::exp(std::sin((1.0*x + 0.5)*PI) + 1.0) - 1.0)*(std::exp(std::sin((1.0*y + 0.5)*PI) + 1.0) - 1.0);
						}
					// start with stationary surface
					std::fill(cells[c].V.begin(), cells[c].V.end(), 0.0);
				}
				hint = 0.7f;
				break;
			case 2: // Continuous spectrum of frequencies, but limited to localised hump in the middle
			{
				// u=exp(f(x,y))-1, where f(x,y)=a*(1-x^2/bx^2-y^2/by^2)
				// note: if c/bx^2 or c/by^2 is too large then the integrator will break due to a CFL condition!
				double bx = 0.3, by = 0.3, a = 1.5, x0 = 0.0, y0 = 0.0;
#pragma omp parallel for
				for (int c = 0; c < position_information.size(); ++c) {
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							auto x = (position_information[c][0] + coords_x[i] * step_size_x) * 2.0 / domain_scaling_factor[0];
							auto y = (position_information[c][1] + coords_y[j] * step_size_y) * 2.0 / domain_scaling_factor[1];
							cells[c].U[j*stages_x + i] = std::max(0.0, std::exp(a - a / bx / bx*(x - x0)*(x - x0) - a / by / by*(y - y0)*(y - y0)) - 1.0);
						}
					// start with stationary surface
					std::fill(cells[c].V.begin(), cells[c].V.end(), 0.0);
				}
			}
			break;
			case 3: // Continuous spectrum of frequencies, but limited to 2 localised humps.
			{
				// u=exp(f(x,y))-1, where f(x,y)=c*(1-x^2/bx^2-y^2/by^2)
				// note: if c/bx^2 or c/by^2 is too large then the integrator will break due to a CFL condition!
				double bx = 0.3, by = 0.3, a = 1.5, x0 = sqrt(0.5)*0.5, y0 = x0;
#pragma omp parallel for
				for (int c = 0; c < position_information.size(); ++c) {
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							auto x = (position_information[c][0] + coords_x[i] * step_size_x) * 2.0 / domain_scaling_factor[0];
							auto y = (position_information[c][1] + coords_y[j] * step_size_y) * 2.0 / domain_scaling_factor[1];
							if (x > 0)
								cells[c].U[j*stages_x + i] = std::max(0.0, std::exp(a - a / bx / bx*(x - x0)*(x - x0) - a / by / by*(y - y0)*(y - y0)) - 1.0);
							else
								cells[c].U[j*stages_x + i] = std::max(0.0, std::exp(a - a / bx / bx*(x + x0)*(x + x0) - a / by / by*(y + y0)*(y + y0)) - 1.0);
						}
					// start with stationary surface
					std::fill(cells[c].V.begin(), cells[c].V.end(), 0.0);
				}
			}
			break;
			case 4: // Continuous spectrum of frequencies, but limited to an off-centre localised hump.
			{
				// u=exp(f(x,y))-1, where f(x,y)=c*(1-x^2/bx^2-y^2/by^2)
				// note: if c/bx^2 or c/by^2 is too large then the integrator will break due to a CFL condition!
				double bx = 0.3, by = 0.3, a = 1.5, x0 = sqrt(0.5)*0.5, y0 = x0;
#pragma omp parallel for
				for (int c = 0; c < position_information.size(); ++c) {
					for (unsigned int j = 0; j < stages_y; ++j)
						for (unsigned int i = 0; i < stages_x; ++i) {
							auto x = (position_information[c][0] + coords_x[i] * step_size_x) * 2.0 / domain_scaling_factor[0];
							auto y = (position_information[c][1] + coords_y[j] * step_size_y) * 2.0 / domain_scaling_factor[1];
							if (x > 0)
								cells[c].U[j*stages_x + i] = std::max(0.0, std::exp(a - a / bx / bx*(x - x0)*(x - x0) - a / by / by*(y - y0)*(y - y0)) - 1.0);
							else
								cells[c].U[j*stages_x + i] = 0.0;// std::max(0.0, std::exp(a - a / bx / bx*(x + x0)*(x + x0) - a / by / by*(y + y0)*(y + y0)) - 1.0);
						}
					// start with stationary surface
					std::fill(cells[c].V.begin(), cells[c].V.end(), 0.0);
				}
			}
			break;
			default: // We've reached the end or been given a wrong value, reset to the first initial condition.
				return Integrator::Change_Initial_Conditions(0);
		}
		update_U_dependent_variables();
		update_V_dependent_variables();
		Half_Step();

		return hint;
	}

	// Helper function to create adjacency and cell type information for non-periodic domains.
	void Integrator::Find_Neighbours()
	{
		adjacency_information.resize(position_information.size());
		cells.resize(position_information.size());
#pragma omp parallel for
		for (int c = 0; c < position_information.size(); ++c) {
			// Find up to 4 neighbouring cells, then assign them to top, left, right, bottom and assign Cell_Type information
			std::array<unsigned int, 4> neighbours{ missing_index,missing_index,missing_index,missing_index };
			int i = c - 1;
			for (; i != -1; --i) { // look for left neighbour
				if (fabs(position_information[i][1] - position_information[c][1]) > 0.5*step_size_y) // abort if we're no longer on the same row
					break;
				if ((position_information[c][0] - position_information[i][0] > 0.5*step_size_x) && (position_information[c][0] - position_information[i][0] < 1.5*step_size_x)) {
					neighbours[0] = i;
					break;
				}
			}
			for (; i != -1; --i) { // continue looking for top neighbour
				if (position_information[c][1] - position_information[i][1] > 1.5*step_size_y) // abort if we're more than 1 row back
					break;
				if ((fabs(position_information[i][0] - position_information[c][0]) < 0.5*step_size_x) && (position_information[c][1] - position_information[i][1] > 0.5*step_size_y) && (position_information[c][1] - position_information[i][1] < 1.5*step_size_y)) {
					neighbours[2] = i;
					break;
				}
			}
			i = c + 1;
			for (; i != position_information.size(); ++i) { // look for right neighbour
				if (fabs(position_information[i][1] - position_information[c][1]) > 0.5*step_size_y) // abort if we're no longer on the same row
					break;
				if ((position_information[i][0] - position_information[c][0] > 0.5*step_size_x) && (position_information[i][0] - position_information[c][0] < 1.5*step_size_x)) {
					neighbours[1] = i;
					break;
				}
			}
			for (; i != position_information.size(); ++i) { // continue looking for bottom neighbour
				if (position_information[i][1] - position_information[c][1] > 1.5*step_size_y) // abort if we're more than 1 row ahead
					break;
				if ((fabs(position_information[i][0] - position_information[c][0]) < 0.5*step_size_x) && (position_information[i][1] - position_information[c][1] > 0.5*step_size_y) && (position_information[c][1] - position_information[i][1] < 1.5*step_size_y)) {
					neighbours[3] = i;
					break;
				}
			}
			adjacency_information[c] = neighbours;

			Cell_Type type_x, type_y;
			if (neighbours[1] == missing_index)
				type_x = Cell_Type::Dirichlet_right; // right => do nothing cell, so it should be ok to mistake a cell for Dirichlet_right when it could also be Dirichlet_left
			else if (neighbours[0] == missing_index)
				type_x = Cell_Type::Dirichlet_left;
			else
				type_x = Cell_Type::Normal;
			if (neighbours[3] == missing_index)
				type_y = Cell_Type::Dirichlet_right;
			else if (neighbours[2] == missing_index)
				type_y = Cell_Type::Dirichlet_left;
			else
				type_y = Cell_Type::Normal;
			cells[c] = Cell(std::make_pair(type_x, type_y));
		}
	}

	// Switch to a new set of boundary conditions.
	// TODO: allow reading of boundary conditions from a file
	float Integrator::Change_Boundary_Conditions(int bc)
	{ // Change the boundary conditons
		if (bc == -1)
			++boundary_conditions;
		else
			boundary_conditions = bc;

		switch (boundary_conditions) {
			case 0: // Periodic square boundary.
			{
				/*
				200 * 200 square domain centered at the origin with periodic boundary conditions
				dx = dy = 0.1
				wave_speed = 1^2
				*/
				unsigned int n = 200;
				cells.clear();
				adjacency_information.clear();
				position_information.clear();
				step_size_x = 0.1;
				step_size_y = 0.1;
				wave_speed = 10.0;
				domain_scaling_factor = { n*step_size_x,n*step_size_y };
				for (unsigned int j = 0; j < n - 1; ++j) {
					for (unsigned int i = 0; i < n - 1; ++i) {
						cells.emplace_back(std::make_pair(Cell_Type::Normal, Cell_Type::Normal));
						adjacency_information.emplace_back(std::array<unsigned int, 4>{ n * j + (i + n - 1) % n, n * j + i + 1, n * ((j + n - 1) % n) + i, n * (j + 1) + i });
						position_information.emplace_back(std::array<double, 2>{ step_size_x * static_cast<int>(i - n / 2), step_size_y * static_cast<int>(j - n / 2)});
					}
					// last column
					cells.emplace_back(std::make_pair(Cell_Type::Periodic, Cell_Type::Normal));
					adjacency_information.emplace_back(std::array<unsigned int, 4>{ n * j + n - 2, n * j, n * ((j + n - 1) % n) + n - 1, n * (j + 1) + n - 1 });
					position_information.emplace_back(std::array<double, 2>{ step_size_x * static_cast<int>(n - 1 - n / 2), step_size_y * static_cast<int>(j - n / 2)});
				}
				// last row
				for (unsigned int i = 0; i < n - 1; ++i) {
					cells.emplace_back(std::make_pair(Cell_Type::Normal, Cell_Type::Periodic));
					adjacency_information.emplace_back(std::array<unsigned int, 4>{ n * (n - 1) + (i + n - 1) % n, n * (n - 1) + i + 1, n * (n - 2) + i, i });
					position_information.emplace_back(std::array<double, 2>{ step_size_x * static_cast<int>(i - n / 2), step_size_y * static_cast<int>(n - 1 - n / 2)});
				}
				// last corner
				cells.emplace_back(std::make_pair(Cell_Type::Periodic, Cell_Type::Periodic));
				adjacency_information.emplace_back(std::array<unsigned int, 4>{ n * n - 2, n * n - n, n * n - 1 - n, n - 1 });
				position_information.emplace_back(std::array<double, 2>{ step_size_x * static_cast<int>(n - 1 - n / 2), step_size_y * static_cast<int>(n - 1 - n / 2)});
				break;
			}
			case 1: // Square with Dirichlet boundaries.
			{
				/*
				200 * 200 square domain with Dirichlet boundary conditions
				dx = dy = 0.1
				wave_speed = 1^2
				*/
				unsigned int n = 200;
				cells.clear();
				adjacency_information.clear();
				position_information.clear();
				step_size_x = 0.1;
				step_size_y = 0.1;
				wave_speed = 10.0;
				domain_scaling_factor = { n*step_size_x,n*step_size_y };
				for (unsigned int j = 0; j < n; ++j)
					for (unsigned int i = 0; i < n; ++i)
						position_information.emplace_back(std::array<double, 2>{step_size_x * static_cast<int>(i - n / 2), step_size_y * static_cast<int>(j - n / 2)});
				Find_Neighbours();
				break;
			}
			case 2: // Circle with Dirichlet boudaries.
			{
				/*
				Circle of radius n/2 centred at the origin.
				dx = dy = 0.1
				wave_speed = 1^2
				*/
				cells.clear();
				adjacency_information.clear();
				position_information.clear();
				unsigned int n = 200;
				step_size_x = 0.1;
				step_size_y = 0.1;
				wave_speed = 10.0;
				domain_scaling_factor = { n*step_size_x,n*step_size_y };
				double r = 0.49*domain_scaling_factor[0];
				for (unsigned int j = 0; j < n; ++j) {
					for (unsigned int i = 0; i < n; ++i) {
						auto x = step_size_x * static_cast<int>(i - n / 2);
						auto y = step_size_y * static_cast<int>(j - n / 2);
						if (x*x + y*y < r*r) {
							position_information.emplace_back(std::array<double, 2>{x, y});
						}
					}
				}
				Find_Neighbours();
			}
			break;
			case 3: // Circle with a cusp, Dirichlet boundaries
			{
				/*
				As with case 2, but with a cusp.
				dx = dy = 0.1
				wave_speed = 1^2
				*/
				cells.clear();
				adjacency_information.clear();
				position_information.clear();
				unsigned int n = 300;
				step_size_x = 0.1;
				step_size_y = 0.1;
				wave_speed = 10.0;
				domain_scaling_factor = { n*step_size_x,n*step_size_y };
				double r = 0.5*domain_scaling_factor[0], r2 = 0.5*r, s = std::sqrt(0.5)*0.5*r;
				for (unsigned int j = 0; j < n; ++j) {
					for (unsigned int i = 0; i < n; ++i) {
						auto x = step_size_x * static_cast<int>(i - n / 2);
						auto y = step_size_y * static_cast<int>(j - n / 2);
						if ((x*x + y*y < r*r) && ((y > x) || (((y - s)*(y - s) + (x - s)*(x - s) < r2*r2) || ((y + s)*(y + s) + (x + s)*(x + s) < r2*r2))))
							position_information.emplace_back(std::array<double, 2>{x, y});
					}
				}
				Find_Neighbours();
				break;
			}
			case 4: // 2 intersecting circles, Dirichlet boundaries
			{
				cells.clear();
				adjacency_information.clear();
				position_information.clear();
				unsigned int n = 300;
				step_size_x = 0.1;
				step_size_y = 0.1;
				wave_speed = 10.0;
				domain_scaling_factor = { n*step_size_x,n*step_size_y };
				double r = 0.26*domain_scaling_factor[0], s = std::sqrt(0.5)*0.25*domain_scaling_factor[0];
				for (unsigned int j = 0; j < n; ++j) {
					for (unsigned int i = 0; i < n; ++i) {
						auto x = step_size_x * static_cast<int>(i - n / 2);
						auto y = step_size_y * static_cast<int>(j - n / 2);
						if (((y - s)*(y - s) + (x - s)*(x - s) < r*r) || ((y + s)*(y + s) + (x + s)*(x + s) < r*r))
							position_information.emplace_back(std::array<double, 2>{x, y});
					}
				}
				Find_Neighbours();
				break;
			}
			default: // We've reached the end or been given a wrong value, reset to the first boundary condition.
				return Integrator::Change_Boundary_Conditions(0);
		}
		// minimise meory usage
		cells.shrink_to_fit();

		// Report on new setup
		std::cout << g_waves << std::endl;
		return Change_Initial_Conditions(initial_conditions);
	}

	// Set the coefficients of the Lobatto IIIA-IIIB stencils.
	std::vector<std::vector<double>> Integrator::Setup_Coefficients(int stages)
	{
		std::vector<std::vector<double>> coefficients(stages - 1); // The last node in a cell is the first node of the next cell, so we drop it.
		switch (stages) {
			case 2:
				// Main node
				coefficients[0].emplace_back(1.0);
				coefficients[0].emplace_back(-2.0);
				coefficients[0].emplace_back(1.0);
				break;

			case 3:
				// Main node
				coefficients[0].emplace_back(-1.0);
				coefficients[0].emplace_back(8.0);
				coefficients[0].emplace_back(-14.0);
				coefficients[0].emplace_back(8.0);
				coefficients[0].emplace_back(-1.0);

				// Internal node
				coefficients[1].emplace_back(4.0);
				coefficients[1].emplace_back(-8.0);
				coefficients[1].emplace_back(4.0);
				break;

			case 4:
				// Main node
				coefficients[0].emplace_back(1.0);
				coefficients[0].emplace_back(0.5*(25.0 - 15.0*sqrt(5.0)));
				coefficients[0].emplace_back(0.5*(25.0 + 15.0*sqrt(5.0)));
				coefficients[0].emplace_back(-52.0);
				coefficients[0].emplace_back(0.5*(25.0 + 15.0*sqrt(5.0)));
				coefficients[0].emplace_back(0.5*(25.0 - 15.0*sqrt(5.0)));
				coefficients[0].emplace_back(1.0);

				// First internal node
				coefficients[1].emplace_back(5.0 + 3.0*sqrt(5.0));
				coefficients[1].emplace_back(-20.0);
				coefficients[1].emplace_back(10.0);
				coefficients[1].emplace_back(5.0 - 3.0*sqrt(5.0));

				// Second internal node
				coefficients[2].emplace_back(5.0 - 3.0*sqrt(5.0));
				coefficients[2].emplace_back(10.0);
				coefficients[2].emplace_back(-20.0);
				coefficients[2].emplace_back(5.0 + 3.0*sqrt(5.0));
				break;

			default: // FIXME: Calculate the higher order coefficients or set up a function to generate them.
				throw std::runtime_error(stages + " stages not implemented yet.");
		}
		return coefficients;
	}

	// Get the coordinates of the nodes within a unit cell.
	std::vector<double> Integrator::Setup_Coords(int stages)
	{
		std::vector<double> coords;
		coords.reserve(stages);
		switch (stages)
		{
			case 2:
				coords.emplace_back(0);
				coords.emplace_back(1);
				break;
			case 3:
				coords.emplace_back(0);
				coords.emplace_back(0.5);
				coords.emplace_back(1);
				break;
			case 4:
				coords.emplace_back(0);
				coords.emplace_back(0.5 - std::sqrt(5.0) / 10.0);
				coords.emplace_back(0.5 + std::sqrt(5.0) / 10.0);
				coords.emplace_back(1);
				break;
			default: // FIXME: Calculate the higher order coordinates or set up a function to generate them.
				throw std::runtime_error(stages + " stages not implemented yet.");
		}
		return coords;
	}

	/// Our friendly ostream operator<<
	std::ostream&  operator<<(std::ostream&  os, const Integrator&  integrator)
	{
		os << "A " << integrator.cells.size() << " element simulation using Lobatto IIIA-IIIB discretisation in space with " << integrator.stages_x + 1 << " stages in x, " << integrator.stages_y + 1 << " stages in y, and stepsizes " << "dx=" << integrator.step_size_x << ", dy=" << integrator.step_size_y << " and dt=" << integrator.step_size_time << ".";
		return os;
	}

}
