#pragma once

#include <vector>
#include <array>
#include <iostream>

namespace Waves {
	/* A differentiation of cell types for use in describing boundary conditions.
		Normal          = A cell with no special boundary conditions.
		Periodic		= The same as normal, but for a top or right edge cell that wraps around periodically (needed for rendering purposes).
		Dirichlet_left  = Main node is constant, other nodes are normal.
		Dirichlet_right = Main node is constant, other nodes are outside the domain.
		Neumann_left    = Main node uses phantom points on the left in computations, other nodes are normal.
		Neumann_right   = Main node uses phantom points on the right in computations, other nodes are outside the domain.
	*/
	enum class Cell_Type { Normal, Periodic, Dirichlet_left, Dirichlet_right, Neumann_left, Neumann_right };
	// The basic cell within the domain resulting from applying Lobatto IIIA-IIIB discretisation in space (x,y).
	// Each cell has stages_x * stages_y nodes.
	class Cell
	{
	public:
		std::pair<Cell_Type, Cell_Type> cell_type; // boundary condition types in x (first) and y (second) directions
		std::vector<double> U, V, U_xx, U_yy;
		Cell() = default;
		Cell(std::pair<Cell_Type, Cell_Type> type) : cell_type(type) {};
	};

	// The integrator class.
	class Integrator
	{ // Functions and variables related to the PDE and integrator.
	private:
		unsigned int stages_x, stages_y; // Number of stages per cell.
		double wave_speed;
		int initial_conditions;
		int boundary_conditions;

		// Coefficient stencils used in the approximation of U_{xx} and U_{yy}.
		std::vector<std::vector<double>> coefficients_x, coefficients_y;
		std::vector<std::vector<double>> Setup_Coefficients(int r);

		// Position within each cell (from 0 to 1) of each node {main node, inner node1, inner node2, ...}
		std::vector<double> coords_x, coords_y;
		std::vector<double> Setup_Coords(int r);

		// Compute second-order spatial derivatives of U.
		void Compute_U_xx();
		void Compute_U_yy();

		// Components of the integrator.
		void Step_U(double step_size);
		void Step_V(double step_size);
		void update_U_dependent_variables();
		void update_V_dependent_variables();

		// Helper function for constructing boundary conditions
		std::array<double,2> domain_scaling_factor;
		void Find_Neighbours();

	public:
		std::vector<Waves::Cell> cells; // Variables in the PDE (including derivatives and other useful values) for each cell.
		// We keep the following information arrays separate from the cells to minimise the memory footprint of the cells.
		std::vector<std::array<unsigned int, 4>> adjacency_information; // std::numeric_limits<unsigned int>::max() indicates no adjacent cell (i.e. for boundary cells)
		std::vector<std::array<double, 2>> position_information;
		double Time = 0.0; // Time the system has evolved for.
		double step_size_time;
		double step_size_x;
		double step_size_y;
		void Step(); // Step forwards in time, i.e., the integrator.
		void Half_Step() { Step_V(0.5*step_size_time); }; // Perform an initial half-step in V to optimise the leapfrog algorithm.
		float Change_Initial_Conditions(int ic = -1); // Defaults to switching to the next initial condition.
		float Change_Boundary_Conditions(int bc = -1); // Defaults to switching to the next initial condition. Also resets the initial conditions and returns the initial conditions hint.

		// Initialise the Integrator.
		float Initialise(int rx, int ry, double dt, int ic, int bc);

		// Friendly ostream operator<< so it can access private variables.
		friend std::ostream& operator<< (std::ostream& os, const Integrator& integrator);
	};

	// We need to use a single global instance of this class for OpenGL to have access to it. We declare it extern here and construct it in integrator.cpp.
	extern Integrator g_waves;

	//// ostream operator for describing the Integator.
	//std::ostream& operator<< (std::ostream& os, const Integrator& integrator);
}
