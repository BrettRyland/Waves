#include <cmath> // for std::sin(), std::cos()
#include <algorithm> // for std::max
#include <numeric> // for std::inner_product
#include "integrator.h"
#include "renderer.h"

namespace Waves {
	// Our global instance of the Integrator. This needs to be global for OpenGL to access it.
	Integrator g_waves;
}

using namespace Waves;

const double PI = 3.141592653589793;

// Initialise the Integrator.
// We return a hint for the z-scaling factor for use later in the renderer.
float Waves::Integrator::Initialise(int rx, int ry, double dt, int ic, int bc) {
	stages_x = rx - 1; // We use rx-1 here since the last node of each cell is the same as the first node of the next cell.
	stages_y = ry - 1; // Similarly for ry.
	step_size_time = dt;
	initial_conditions = ic;
	boundary_conditions = bc;
	coefficients_x = Setup_Coefficients(stages_x + 1);
	coefficients_y = Setup_Coefficients(stages_y + 1);
	coords_x = Setup_Coords(stages_x + 1);
	coords_y = Setup_Coords(stages_y + 1);
	Change_Boundary_Conditions(boundary_conditions);
	return Change_Initial_Conditions(initial_conditions);
}

// Generate a true/false mask for a given cell for determining which nodes in the cell need updating.
void Generate_Update_Mask(std::vector<bool> & mask, const Waves::Cell & cell, unsigned int stages_x, unsigned int stages_y)
{
	auto it = mask.begin();
	int i;
	switch (cell.cell_type.first) { // x-direction.
		case Waves::Cell_Type::Normal:
		case Waves::Cell_Type::Neumann_left:
			mask.assign(stages_x*stages_y, true);
			break;
		case Waves::Cell_Type::Neumann_right:
			for (it = mask.begin(), i = 0; it != mask.end(); ++it, ++i)
				if (i % stages_x != 0)
					*it = false;
				else
					*it = true;
			break;
		case Waves::Cell_Type::Dirichlet_left:
			for (it = mask.begin(), i = 0; it != mask.end(); ++it, ++i)
				if (i % stages_x != 0)
					*it = true;
				else
					*it = false;
			break;
		case Waves::Cell_Type::Dirichlet_right:
			mask.assign(stages_x*stages_y, false);
			break;
		default:
			throw std::runtime_error("Invalid cell type.");
	}
	switch (cell.cell_type.second) { // y-direction. Here we only set extra false values.
		case Waves::Cell_Type::Normal:
		case Waves::Cell_Type::Neumann_left:
			break;
		case Waves::Cell_Type::Neumann_right:
			for (it = mask.begin() + stages_x; it != mask.end(); ++it)
				*it = false;
			break;
		case Waves::Cell_Type::Dirichlet_left:
			for (it = mask.begin(); it != mask.begin() + stages_x; ++it)
				*it = false;
			break;
		case Waves::Cell_Type::Dirichlet_right:
			mask.assign(stages_x*stages_y, false);
			break;
		default:
			throw std::runtime_error("Invalid cell type.");
	}
}

// Part of the integrator. Perform a step of size step_size in the V variables.
void Waves::Integrator::Step_V(double step_size)
{
	auto derivative_of_the_potential = [](double i) {return 2 * i + 4 * i*i*i; }; // V'(u)=2*u+4*u^3
	// Initialised the mask outside the for loop to save on reallocations as we're going to use it lots.
	std::vector<bool> mask(stages_x*stages_y, true);
	for (auto& cell : cells)
		if (cell.cell_type.first != Waves::Cell_Type::Inactive) {
			Generate_Update_Mask(mask, cell, stages_x, stages_y);
			auto it_mask = mask.begin();
			auto it_V = cell.V.begin();
			auto it_U = cell.U.begin();
			auto it_U_xx = cell.U_xx.begin();
			auto it_U_yy = cell.U_yy.begin();
			for (; it_V != cell.V.end(); ++it_V, ++it_U_xx, ++it_U_yy, ++it_U, ++it_mask)
				if (*it_mask)
					*it_V += step_size * (wave_speed * (*it_U_xx + *it_U_yy) - derivative_of_the_potential(*it_U));
		}
	this->update_V_dependent_variables();
}

// Part of the integrator. Perform a step of size step_size in the U variables.
void Waves::Integrator::Step_U(double step_size)
{
	// Initialised the mask outside the for loop to save on reallocations as we're going to use it lots.
	std::vector<bool> mask(stages_x*stages_y, true);
	for (auto& cell : cells)
		if (cell.cell_type.first != Waves::Cell_Type::Inactive) {
			Generate_Update_Mask(mask, cell, stages_x, stages_y);
			auto it_mask = mask.begin();
			auto it_U = cell.U.begin();
			auto it_V = cell.V.begin();
			for (; it_U != cell.U.end(); ++it_U, ++it_V)
				if (*it_mask)
					*it_U += step_size * (*it_V);
		}
	this->update_U_dependent_variables();
}

// Apply one step of the integrator.
void Waves::Integrator::Step()
{ // The integrator.
	// Step forwards in time using 2-stage Lobatto IIIA-IIIB in time (i.e., leapfrog) and rx-stage Lobatto IIIA-IIIB in x and y directions.
	// It is explicit, local, multisymplectic and handles boundary conditions easily. See my PhD thesis for the details.
	// FIXME: We don't need the values of V between steps, so we should perform a half-step in V during initialisation (or after any modification to the set-up) and then could just use whole-steps in U and V here.

	Step_V(0.5*step_size_time);
	Step_U(step_size_time);
	Step_V(0.5*step_size_time);
	Time += step_size_time;
}

// Compute the second order derivative of U in the x-direction using the Lobatto IIIA-IIIB stencils.
void Waves::Integrator::Compute_U_xx(Waves::Cell& cell, const Waves::Cell& x_left, const Waves::Cell& x_right)
{
	// Function to apply the first coefficient stencil to the local cell centred on the main node (includes all values in cell to the left and first value in cell to right).
	auto main_node_xx = [this](std::vector<double>::iterator it_xx, std::vector<double>::const_iterator it_U, std::vector<double>::const_iterator it_U_left, std::vector<double>::const_iterator it_U_right, std::vector<double>::const_iterator it_coefs) {
		*it_xx = std::inner_product(it_U_left, it_U_left + stages_x, it_coefs,
			std::inner_product(it_U, it_U + stages_x, it_coefs + stages_x, // *it_coefs has size 2*stages_x+1, so this selects the middle value
				*it_U_right * *(it_coefs + 2 * stages_x))) / step_size_x / step_size_x;
	};
	// Funtion to apply the remaining coefficient stencils to the local cell (includes first value in cell to the right).
	auto internal_node_xx = [this](std::vector<double>::iterator it_xx, std::vector<double>::const_iterator it_U, std::vector<double>::const_iterator it_U_right, std::vector<double>::const_iterator it_coefs) {
		*it_xx = std::inner_product(it_U, it_U + stages_x, it_coefs, *it_U_right * *(it_coefs + stages_x)) / step_size_x / step_size_x; // *it_coefs has size stages_x+1, so this selects the last value
	};

	// Apply the appropriate stencils based on the cell type.
	switch (cell.cell_type.first) {
		case Waves::Cell_Type::Normal:
			for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_y*stages_x; offset += stages_x) {
				auto it_coefs = coefficients_x.begin();
				main_node_xx(cell.U_xx.begin() + offset, cell.U.begin() + offset, x_left.U.begin() + offset, x_right.U.begin() + offset, (*it_coefs).begin());
				for (++it_coefs; it_coefs != coefficients_x.end(); ++it_coefs) {
					auto offset2 = offset + std::distance(coefficients_x.begin(), it_coefs);
					internal_node_xx(cell.U_xx.begin() + offset2, cell.U.begin() + offset2, x_right.U.begin() + offset2, (*it_coefs).begin());
				}
			}
			break;
		case Waves::Cell_Type::Inactive:
			// Do nothing. Actually, we shouldn't get to here as this function should never be called for the Inactive case.
			break;
		case Waves::Cell_Type::Neumann_left:
			// We need all nodes, but we have to use phantom points in the left cell for the main node.
			throw "Function not implemented yet.";
			for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_y*stages_x; offset += stages_x) {
				auto it_coefs = coefficients_x.begin();
				//main_node_xx((*it).U_xx.begin() + offset, (*it).U.begin() + offset, (*x_left).U.begin() + offset, (*x_right).U.begin() + offset, (*it_coefs).begin());
				for (++it_coefs; it_coefs != coefficients_x.end(); ++it_coefs) {
					auto offset2 = offset + std::distance(coefficients_x.begin(), it_coefs);
					internal_node_xx(cell.U_xx.begin() + offset2, cell.U.begin() + offset2, x_right.U.begin() + offset2, (*it_coefs).begin());
				}
			}
			break;
		case Waves::Cell_Type::Neumann_right:
			// We only need the main nodes, but we have to use phantom points in the current cell.
			throw "Function not implemented yet.";
			for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_y*stages_x; offset += stages_x) {
				auto it_coefs = coefficients_x.begin();
				//main_node_xx((*it).U_xx.begin() + offset, (*it).U.begin() + offset, (*x_left).U.begin() + offset, (*x_right).U.begin() + offset, (*it_coefs).begin());
			}
			break;
		case Waves::Cell_Type::Dirichlet_left:
			// We only need the internal nodes.
			for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_y*stages_x; offset += stages_x) {
				auto it_coefs = coefficients_x.begin() + 1;
				for (; it_coefs != coefficients_x.end(); ++it_coefs) {
					auto offset2 = offset + std::distance(coefficients_x.begin(), it_coefs);
					internal_node_xx(cell.U_xx.begin() + offset2, cell.U.begin() + offset2, x_right.U.begin() + offset2, (*it_coefs).begin());
				}
			}
			break;
		case Waves::Cell_Type::Dirichlet_right:
			// Do nothing. We don't need the derivatives for any nodes in this case.
			break;
		default:
			throw std::runtime_error("Invalid cell type.");
	}
}

// Compute the second order derivative of U in the y-direction using the Lobatto IIIA-IIIB stencils.
void Waves::Integrator::Compute_U_yy(Waves::Cell& cell, const Waves::Cell& y_left, const Waves::Cell& y_right)
{
	// Unfortunately, we can't use std::inner_product here as the relevant elements of U, etc., are not adjacent. So, here's a version that has non-unity step size, but with no bounds checking!
	auto inner_product_n = [](std::vector<double>::const_iterator input1_start, std::vector<double>::const_iterator input2_start, size_t n, std::iterator_traits<std::vector<double>::const_iterator>::difference_type adv, double value) {
		auto ret = value;
		for (size_t i = 0; i != n - 1; ++n, std::advance(input1_start, adv), std::advance(input2_start, adv))
			ret += *input1_start * *input2_start;
		return ret + *input1_start * *input2_start;
	};

	// Function to apply the first coefficient stencil to the local cell centred on the main node (includes all values in cell to the left and first value in cell to right).
	auto main_node_yy = [this, &inner_product_n](std::vector<double>::iterator it_xx, std::vector<double>::const_iterator it_U, std::vector<double>::const_iterator it_U_left, std::vector<double>::const_iterator it_U_right, std::vector<double>::const_iterator it_coefs) {
		*it_xx = inner_product_n(it_U_left, it_coefs, stages_y, stages_x,
			inner_product_n(it_U, it_coefs + stages_y, stages_y, stages_x, // *it_coefs has size 2*stages_x+1, so this selects the middle value
				*it_U_right * *(it_coefs + 2 * stages_y))) / step_size_y / step_size_y;
	};

	// Function to apply the remaining coefficient stencils to the local cell (includes first value in cell to the right).
	auto internal_node_yy = [this, &inner_product_n](std::vector<double>::iterator it_xx, std::vector<double>::const_iterator it_U, std::vector<double>::const_iterator it_U_right, std::vector<double>::const_iterator it_coefs) {
		*it_xx = inner_product_n(it_U, it_coefs, stages_y, stages_x, *it_U_right * *(it_coefs + stages_y)) / step_size_y / step_size_y; // *it_coefs has size stages_x+1, so this selects the last value
	};

	// Apply the appropriate stencils based on the cell type.
	switch (cell.cell_type.second) {
		case Waves::Cell_Type::Normal:
			for (std::iterator_traits<std::vector<double>::const_iterator>::difference_type offset = 0; offset != stages_x; ++offset) {
				auto it_coefs = coefficients_y.begin();
				main_node_yy(cell.U_yy.begin() + offset, cell.U.begin() + offset, y_left.U.begin() + offset, y_right.U.begin() + offset, (*it_coefs).begin());
				for (++it_coefs; it_coefs != coefficients_y.end(); ++it_coefs) {
					auto offset2 = offset + std::distance(coefficients_y.begin(), it_coefs)*stages_x;
					internal_node_yy(cell.U_yy.begin() + offset2, cell.U.begin() + offset2, y_right.U.begin() + offset2, (*it_coefs).begin());
				}
			}
			break;
		case Waves::Cell_Type::Inactive:
			// Do nothing. Actually, we shouldn't get to here as this function should never be called for the Inactive case.
			break;
		case Waves::Cell_Type::Neumann_left:
			// We need all nodes, but we have to use phantom points in the left cell for the main node.
			throw "Function not implemented yet.";
			for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_x; ++offset) {
				auto it_coefs = coefficients_y.begin();
				//main_node_xx((*it).U_yy.begin() + offset, (*it).U.begin() + offset, (*y_left).U.begin() + offset, (*y_right).U.begin() + offset, (*it_coefs).begin());
				for (++it_coefs; it_coefs != coefficients_y.end(); ++it_coefs) {
					auto offset2 = offset + std::distance(coefficients_y.begin(), it_coefs)*stages_x;
					internal_node_yy(cell.U_yy.begin() + offset2, cell.U.begin() + offset2, y_right.U.begin() + offset2, (*it_coefs).begin());
				}
			}
			break;
		case Waves::Cell_Type::Neumann_right:
			// We only need the main nodes, but we have to use phantom points in the current cell.
			throw "Function not implemented yet.";
			for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_x; ++offset) {
				auto it_coefs = coefficients_y.begin();
				//main_node_yy((*it).Uyy.begin() + offset, (*it).U.begin() + offset, (*x_left).U.begin() + offset, (*x_right).U.begin() + offset, (*it_coefs).begin());
			}
			break;
		case Waves::Cell_Type::Dirichlet_left:
			// We only need the internal nodes.
			for (std::iterator_traits<std::vector<double>::iterator>::difference_type offset = 0; offset != stages_x; ++offset) {
				auto it_coefs = coefficients_y.begin() + 1;
				for (; it_coefs != coefficients_y.end(); ++it_coefs) {
					auto offset2 = offset + std::distance(coefficients_y.begin(), it_coefs)*stages_x;
					internal_node_yy(cell.U_yy.begin() + offset2, cell.U.begin() + offset2, y_right.U.begin() + offset2, (*it_coefs).begin());
				}
			}
			break;
		case Waves::Cell_Type::Dirichlet_right:
			// Do nothing. We don't need the derivatives for any nodes in this case.
			break;
		default:
			throw std::runtime_error("Invalid cell type.");
	}
}

// Update variables that depend on U (i.e. U_xx and U_yy).
void Waves::Integrator::update_U_dependent_variables()
{
	for (long c = 0; c < adjacency_information.size(); ++c) {
		Compute_U_xx(cells[c], cells[adjacency_information[c][0]], cells[adjacency_information[c][1]]);
		Compute_U_yy(cells[c], cells[adjacency_information[c][2]], cells[adjacency_information[c][3]]);
	}

	//// Note: A slight waste of processing power here as we only exclude Inactive cells. We could use the update mask from Generate_Update_Mask(), but this would probably be more expensive.

	//// We use an external iterator instead of a range-for loop as we need to be aware of the target cell's location and whom its neighbours are.
	//auto it = cells.begin();
	//// Some local iterators for referencing the neighbouring cells (starting in the left,left corner).
	//auto x_left{ it + (domain_size_x - 1) };
	//auto x_right{ it + 1 };
	//auto y_left{ it + (domain_size_y - 1)*domain_size_x };
	//auto y_right{ it + domain_size_x };

	//// left,left corner
	//if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//	Compute_U_xx(it, x_left, x_right);
	//	Compute_U_yy(it, y_left, y_right);
	//}
	//// centre,left
	//for (x_left = it, it = x_right++, ++y_left, ++y_right; it != cells.begin() + domain_size_x - 1; x_left = it, it = x_right++, ++y_left, ++y_right)
	//{
	//	if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//		Compute_U_xx(it, x_left, x_right);
	//		Compute_U_yy(it, y_left, y_right);
	//	}
	//}
	//// right,left corner
	//x_right = it - (domain_size_x - 1);
	//if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//	Compute_U_xx(it, x_left, x_right);
	//	Compute_U_yy(it, y_left, y_right);
	//}
	//for (auto row = 2; row != domain_size_y; ++row) {
	//	// left, centre
	//	++it;
	//	x_left = it + domain_size_x - 1;
	//	x_right = it + 1;
	//	y_left = it - domain_size_x;
	//	y_right = it + domain_size_x;
	//	if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//		Compute_U_xx(it, x_left, x_right);
	//		Compute_U_yy(it, y_left, y_right);
	//	}
	//	// centre, centre
	//	for (x_left = it, it = x_right++, ++y_left, ++y_right; it != cells.begin() + row*domain_size_x - 1; x_left = it, it = x_right++, ++y_left, ++y_right) {
	//		if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//			Compute_U_xx(it, x_left, x_right);
	//			Compute_U_yy(it, y_left, y_right);
	//		}
	//	}
	//	// right, centre
	//	x_right = it - (domain_size_x - 1);
	//	if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//		Compute_U_xx(it, x_left, x_right);
	//		Compute_U_yy(it, y_left, y_right);
	//	}
	//}
	//// left,right corner
	//++it;
	//x_left = it + (domain_size_x - 1);
	//x_right = it + 1;
	//y_left = it - domain_size_x;
	//y_right = it - domain_size_x*(domain_size_y - 1);
	//if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//	Compute_U_xx(it, x_left, x_right);
	//	Compute_U_yy(it, y_left, y_right);
	//}
	//// centre,right
	//for (x_left = it, it = x_right++, ++y_left, ++y_right; it != cells.end() - 1; x_left = it, it = x_right++, ++y_left, ++y_right)
	//{
	//	if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//		Compute_U_xx(it, x_left, x_right);
	//		Compute_U_yy(it, y_left, y_right);
	//	}
	//}
	//// right,right
	//x_right = it - domain_size_x + 1;
	//if ((*it).cell_type.first != Waves::Cell_Type::Inactive) {
	//	Compute_U_xx(it, x_left, x_right);
	//	Compute_U_yy(it, y_left, y_right);
	//}
}

// Update variables that depend on V (i.e. nothing).
void Waves::Integrator::update_V_dependent_variables()
{ // Nothing to do. Should be optimised away, but useful placeholder for future simulations that may use it.
}

// Switch to a new set of initial conditions.
float Waves::Integrator::Change_Initial_Conditions(int ic)
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
	}

	switch (initial_conditions) {
		case 0: // Single frequency
			for (long c = 0; c < position_information.size(); ++c) {
				for (int j = 0; j < stages_y; ++j)
					for (int i = 0; i < stages_x; ++i) {
						auto x = (position_information[c][0] + coords_x[i] * step_size_x) / 100; // [-1,1)
						auto y = (position_information[c][1] + coords_y[j] * step_size_y) / 100; // [-1,1)
						// u=cos(2*PI*(t+c1*x+c2*y)), which has a wave speed of sqrt(c1^2+c2^2)
						cells[c].U[j*stages_x + i] = cos(4.0*PI*(x + y));
						// set V to be the derivative of U to get a travelling wave train (actually -ve direction)
						cells[c].V[j*stages_x + i] = -4.0*PI / g_waves.wave_speed * std::sqrt(static_cast<double>(200 * step_size_x) + static_cast<double>(200 * step_size_y))*std::sin(4.0*PI*(x + y));
					}
			}
			//for (auto it = cells.begin(); it != cells.end(); ++it) {
			//	auto it_U = (*it).U.begin();
			//	auto it_V = (*it).V.begin();
			//	for (; it_U != (*it).U.end(); ++it_U, ++it_V) {
			//		auto x = (static_cast<double>(std::distance(cells.begin(), it) % domain_size_x) + coords_x[std::distance((*it).U.begin(), it_U) % stages_x]) / domain_size_x; // [0,1)
			//		auto y = (static_cast<double>(std::distance(cells.begin(), it) / domain_size_x) + coords_x[std::distance((*it).U.begin(), it_U) / stages_x]) / domain_size_y; // [0,1)
			//		// u=cos(2*PI*(t+c1*x+c2*y)), which has a wave speed of sqrt(c1^2+c2^2)
			//		*it_U = cos(4.0*PI*(x + y));
			//		// set V to be the derivative of U to get a travelling wave train (actually -ve direction)
			//		*it_V = -4.0*PI / g_waves.wave_speed * std::sqrt(static_cast<double>(domain_size_x*step_size_x) + static_cast<double>(domain_size_y*step_size_y))*std::sin(4.0*PI*(x + y));
			//	}
			//}
			hint = 0.5f;
			break;
		case 1: // Continuous spectrum of frequencies
			for (long c = 0; c < position_information.size(); ++c) {
				for (int j = 0; j < stages_y; ++j)
					for (int i = 0; i < stages_x; ++i) {
						auto x = (position_information[c][0] + coords_x[i] * step_size_x) / 200; // [-1,1)
						auto y = (position_information[c][1] + coords_y[j] * step_size_y) / 200; // [-1,1)
						// this choice should be zero along the edge of the domain but non-zero everywhere else
						cells[c].U[j*stages_x + i] = 0.0625*(std::exp(std::sin((4.0*x - 0.5)*PI) + 1.0) - 1.0)*(std::exp(std::sin((2.0*y - 0.5)*PI) + 1.0) - 1.0);
						// start with stationary surface
						cells[c].V[j*stages_x + i] = 0.0;
					}
			}
			//for (auto it = cells.begin(); it != cells.end(); ++it) {
			//	auto it_U = (*it).U.begin();
			//	auto it_V = (*it).V.begin();
			//	for (; it_U != (*it).U.end(); ++it_U, ++it_V) {
			//		auto x = (static_cast<double>(std::distance(cells.begin(), it) % domain_size_x) + coords_x[std::distance((*it).U.begin(), it_U) % stages_x]) / domain_size_x; // [0,1)
			//		auto y = (static_cast<double>(std::distance(cells.begin(), it) / domain_size_x) + coords_y[std::distance((*it).U.begin(), it_U) / stages_x]) / domain_size_y; // [0,1)
			//		// this choice should be zero along the edge of the domain but non-zero everywhere else
			//		*it_U = 0.0625*(std::exp(std::sin((4.0*x - 0.5)*PI) + 1.0) - 1.0)*(std::exp(std::sin((2.0*y - 0.5)*PI) + 1.0) - 1.0);
			//		// start with stationary surface
			//		*it_V = 0.0;
			//	}
			//}
			break;
		case 2: // Continuous spectrum of frequencies, but limited to localised hump in the middle
		{
			// u=exp(f(x,y))-1, where f(x,y)=c*(1-x^2/bx^2-y^2/by^2)
			// note: if c/bx^2 or c/by^2 is too large then the integrator will break due to a CFL condition!
			double bx = 0.3*step_size_x, by = 0.3*step_size_y, c = 1.5, x0 = 0.0, y0 = 0.0;
			for (long c = 0; c < position_information.size(); ++c) {
				for (int j = 0; j < stages_y; ++j)
					for (int i = 0; i < stages_x; ++i) {
						auto x = (position_information[c][0] + coords_x[i] * step_size_x) / 200 * step_size_x; // [-dx,dx)
						auto y = (position_information[c][1] + coords_y[j] * step_size_y) / 200 * step_size_y; // [-dy,dy)
						// this choice should be zero along the edge of the domain but non-zero everywhere else
						cells[c].U[j*stages_x + i] = std::max(0.0, std::exp(c - c / bx / bx*(x - x0)*(x - x0) - c / by / by*(y - y0)*(y - y0)) - 1.0);
						// start with stationary surface
						cells[c].V[j*stages_x + i] = 0.0;
					}
			}
			//for (auto it = cells.begin(); it != cells.end(); ++it) {
			//	auto it_U = (*it).U.begin();
			//	auto it_V = (*it).V.begin();
			//	for (; it_U != (*it).U.end(); ++it_U, ++it_V) {
			//		auto x = (static_cast<double>(std::distance(cells.begin(), it) % domain_size_x) + coords_x[std::distance((*it).U.begin(), it_U) % stages_x]) / domain_size_x * step_size_x * 2.0 - step_size_x; // [-dx,dx)
			//		auto y = (static_cast<double>(std::distance(cells.begin(), it) / domain_size_x) + coords_y[std::distance((*it).U.begin(), it_U) / stages_x]) / domain_size_y * step_size_y * 2.0 - step_size_y; // [-dy,dy)
			//		// this choice should be zero along the edge of the domain but non-zero everywhere else
			//		*it_U = std::max(0.0, std::exp(c - c / bx / bx*(x - x0)*(x - x0) - c / by / by*(y - y0)*(y - y0)) - 1.0);
			//		// start with stationary surface
			//		*it_V = 0.0;
			//	}
			//}
		}
		break;
		case 3: // Continuous spectrum of frequencies, but limited to 2 localised humps. (Best used with BC=4 once it has been implemented.)
		{
			// u=exp(f(x,y))-1, where f(x,y)=c*(1-x^2/bx^2-y^2/by^2)
			// note: if c/bx^2 or c/by^2 is too large then the integrator will break due to a CFL condition!
			double bx = 0.2*step_size_x, by = 0.2*step_size_y, c = 1.5, x0 = 0.49, y0 = 0.0;
			for (long c = 0; c < position_information.size(); ++c) {
				for (int j = 0; j < stages_y; ++j)
					for (int i = 0; i < stages_x; ++i) {
						auto x = (position_information[c][0] + coords_x[i] * step_size_x) / 200 * step_size_x; // [-dx,dx)
						auto y = (position_information[c][1] + coords_y[j] * step_size_y) / 200 * step_size_y; // [-dy,dy)
						// this choice should be zero along the edge of the domain but non-zero everywhere else
						if(x>0)
							cells[c].U[j*stages_x + i] = std::max(0.0, std::exp(c - c / bx / bx*(x - x0)*(x - x0) - c / by / by*(y - y0)*(y - y0)) - 1.0);
						else
							cells[c].U[j*stages_x + i] = std::max(0.0, std::exp(c - c / bx / bx*(x + x0)*(x + x0) - c / by / by*(y - y0)*(y - y0)) - 1.0);
						// start with stationary surface
						cells[c].V[j*stages_x + i] = 0.0;
					}
			}
			//for (auto it = cells.begin(); it != cells.end(); ++it) {
			//	auto it_U = (*it).U.begin();
			//	auto it_V = (*it).V.begin();
			//	for (; it_U != (*it).U.end(); ++it_U, ++it_V) {
			//		auto x = (static_cast<double>(std::distance(cells.begin(), it) % domain_size_x) + coords_x[std::distance((*it).U.begin(), it_U) % stages_x]) / domain_size_x * step_size_x * 2.0 - step_size_x; // [-dx,dx)
			//		auto y = (static_cast<double>(std::distance(cells.begin(), it) / domain_size_x) + coords_y[std::distance((*it).U.begin(), it_U) / stages_x]) / domain_size_y * step_size_y * 2.0 - step_size_y; // [-dy,dy)
			//		// this choice should be zero along the edge of the domain but non-zero everywhere else
			//		if (x > 0)
			//			*it_U = std::max(0.0, std::exp(c - c / bx / bx*(x - x0)*(x - x0) - c / by / by*(y - y0)*(y - y0)) - 1.0);
			//		else
			//			*it_U = std::max(0.0, std::exp(c - c / bx / bx*(x + x0)*(x + x0) - c / by / by*(y - y0)*(y - y0)) - 1.0);
			//		// start with stationary surface
			//		*it_V = 0.0;
			//	}
			//}
		}
		break;
		default: // We've reached the end or been given a wrong value, reset to the first initial condition.
			return Waves::Integrator::Change_Initial_Conditions(0);
	}
	this->update_U_dependent_variables();
	this->update_V_dependent_variables();

	return hint;
}

// Switch to a new set of boundary conditions. (FIXME: Only the first two have been implemented so far.)
void Waves::Integrator::Change_Boundary_Conditions(int bc)
{ // Change the boundary conditons
	if (bc == -1)
		++boundary_conditions;
	else
		boundary_conditions = bc;

	switch (boundary_conditions) {
		case 0: // Periodic square boundary.
			/*
			200 * 200 square domain with periodic boundary conditions
			dx = dy = 1
			wave_speed = 10^2
			*/
			long n = 200;
			cells.clear();
			adjacency_information.clear();
			position_information.clear();
			step_size_x = 1.0;
			step_size_y = 1.0;
			wave_speed = 10.0;
			for (long j = 0; j < n; ++j) {
				for (long i = 0; i < n; ++i) {
					cells.emplace_back(std::make_pair(Waves::Cell_Type::Normal, Waves::Cell_Type::Normal));
					adjacency_information.emplace_back(n * j + (i + n - 1) % n, n * j + (i + 1) % n, n * ((j + n - 1) % n) + i, n * ((j + 1) % n) + i);
					position_information.emplace_back(i*step_size_x - n / 2 * step_size_x, j*step_size_y - n / 2 * step_size_y);
				}
			}
			break;
		case 1: // Square with Dirichlet boundaries.
			/*
			200 * 200 square domain with Dirichlet boundary conditions
			dx = dy = 1
			wave_speed = 10^2
			*/
			long n = 200;
			cells.clear();
			adjacency_information.clear();
			position_information.clear();
			step_size_x = 1.0;
			step_size_y = 1.0;
			wave_speed = 10.0;

			// left,left corner
			cells.emplace_back(std::make_pair(Waves::Cell_Type::Dirichlet_left, Waves::Cell_Type::Dirichlet_left));
			adjacency_information.emplace_back(std::array<long, 4>{-1, 1, -1, n});
			position_information.emplace_back(std::array<double, 2>{-n / 2 * step_size_x, -n / 2 * step_size_y});
			// center,left
			for (long i = 1; i < n - 1; ++i) {
				cells.emplace_back(std::make_pair(Waves::Cell_Type::Normal, Waves::Cell_Type::Dirichlet_left));
				adjacency_information.emplace_back(std::array<long, 4>{i - 1, i + 1, -1, n + i});
				position_information.emplace_back(std::array<double, 2>{i*step_size_x - n / 2 * step_size_x, -n / 2 * step_size_y});
			}
			// right,left corner
			cells.emplace_back(std::make_pair(Waves::Cell_Type::Dirichlet_right, Waves::Cell_Type::Dirichlet_left));
			adjacency_information.emplace_back(std::array<long, 4>{n - 2, -1, -1, 2 * n - 1});
			position_information.emplace_back(std::array<double, 2>{(n - 1)*step_size_x - n / 2 * step_size_x, -n / 2 * step_size_y});
			for (long j = 1; j < n - 1; ++j) {
				// left,center
				cells.emplace_back(std::make_pair(Waves::Cell_Type::Dirichlet_left, Waves::Cell_Type::Normal));
				adjacency_information.emplace_back(std::array<long, 4>{-1, n * j + 1, n * (j - 1), n * (j + 1)});
				position_information.emplace_back(std::array<double, 2>{-n / 2 * step_size_x, j*step_size_y - n / 2 * step_size_y});
				// center,center
				for (long i = 1; i < n - 1; ++i) {
					cells.emplace_back(std::make_pair(Waves::Cell_Type::Normal, Waves::Cell_Type::Normal));
					adjacency_information.emplace_back(std::array<long, 4>{n * j + (i + n - 1) % n, n * j + (i + 1) % n, n * ((j + n - 1) % n) + i, n * ((j + 1) % n) + i});
					position_information.emplace_back(std::array<double, 2>{i*step_size_x - n / 2 * step_size_x, j*step_size_y - n / 2 * step_size_y});
				}
				// right,center
				cells.emplace_back(std::make_pair(Waves::Cell_Type::Dirichlet_right, Waves::Cell_Type::Normal));
				adjacency_information.emplace_back(std::array<long, 4>{n * j + n - 2, -1, n * (j - 1), n * (j + 1)});
				position_information.emplace_back(std::array<double, 2>{(n - 1)*step_size_x - n / 2 * step_size_x, j*step_size_y - n / 2 * step_size_y});
			}
			// left,right corner
			cells.emplace_back(std::make_pair(Waves::Cell_Type::Dirichlet_left, Waves::Cell_Type::Dirichlet_right));
			adjacency_information.emplace_back(std::array<long, 4>{-1, n * (n - 1) + 1, n * (n - 2), -1});
			position_information.emplace_back(std::array<double, 2>{-n / 2 * step_size_x, (n - 1)*step_size_y - n / 2 * step_size_y});
			// center,right
			for (long i = 1; i < n - 1; ++i) {
				cells.emplace_back(std::make_pair(Waves::Cell_Type::Normal, Waves::Cell_Type::Dirichlet_right));
				adjacency_information.emplace_back(std::array<long, 4>{i - 1, i + 1, n * (n - 2) + i, -1});
				position_information.emplace_back(std::array<double, 2>{i*step_size_x - n / 2 * step_size_x, (n - 1)*step_size_y - n / 2 * step_size_y});
			}
			// right,right corner
			cells.emplace_back(std::make_pair(Waves::Cell_Type::Dirichlet_right, Waves::Cell_Type::Dirichlet_right));
			adjacency_information.emplace_back(std::array<long, 4>{n - 2, -1, n * (n - 2), -1});
			position_information.emplace_back(std::array<double, 2>{(n - 1)*step_size_x - n / 2 * step_size_x, (n - 1)*step_size_y - n / 2 * step_size_y});
			break;
			//case 3: // Circle with a cusp, Dirichlet boundaries
			//{
			//	double r = 1, r2 = 0.51;
			//	// if( (2.0*j/(double)((stages_y-1)*domain_size_y)-1.0 < 0.0) && ((2.0*j/(double)((stages_y-1)*domain_size_y)-1.0)*(2.0*j/(double)((stages_y-1)*domain_size_y)-1.0)+(2.0*i/(double)((stages_x-1)*domain_size_x)-(2.0-radius2))*(2.0*i/(double)((stages_x-1)*domain_size_x)-(2.0-radius2)) >= radius2*radius2) && ((2.0*j/(double)((stages_y-1)*domain_size_y)-1.0)*(2.0*j/(double)((stages_y-1)*domain_size_y)-1.0)+(2.0*i/(double)((stages_x-1)*domain_size_x)-radius2)*(2.0*i/(double)((stages_x-1)*domain_size_x)-radius2) >= radius2*radius2) )
			//	// 	// domain[i*(stages_y-1)*domain_size_y+j] = false;
			//	break;
			//}
			//case 4: // 2 intersecting circles, Dirichlet boundaries
			//{
			//	double r = 1, r2 = 0.51;
			//	// if( ((2.0*j/(double)((stages_y-1)*domain_size_y)-1.0)*(2.0*j/(double)((stages_y-1)*domain_size_y)-1.0)+(2.0*i/(double)((stages_x-1)*domain_size_x)-(2.0-radius2))*(2.0*i/(double)((stages_x-1)*domain_size_x)-(2.0-radius2)) >= radius2*radius2) && ((2.0*j/(double)((stages_y-1)*domain_size_y)-1.0)*(2.0*j/(double)((stages_y-1)*domain_size_y)-1.0)+(2.0*i/(double)((stages_x-1)*domain_size_x)-radius2)*(2.0*i/(double)((stages_x-1)*domain_size_x)-radius2) >= radius2*radius2) )
			//	// 	// domain[i*(stages_y-1)*domain_size_y+j] = false;
			//	break;
			//}
		default: // We've reached the end or been given a wrong value, reset to the first boundary condition.
			return Waves::Integrator::Change_Boundary_Conditions(0);
	}
}

// Set the coefficients of the Lobatto IIIA-IIIB stencils.
std::vector<std::vector<double>> Waves::Integrator::Setup_Coefficients(int stages)
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
std::vector<double> Waves::Integrator::Setup_Coords(int stages)
{
	std::vector<double> coords(stages);
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

std::ostream& Waves::operator<< (std::ostream& os, const Waves::Integrator& integrator)
{
	os << "dx=" << integrator.step_size_x << ", dy=" << integrator.step_size_y << ", " << integrator.stages_x << " stages in x, " << integrator.stages_y << " stages in y, and time step " << integrator.step_size_time << ".";
	return os;
}
