#pragma once

#include <vector>
#include <array>
#include <iostream>

namespace Waves {
	/* A differentiation of cell types for use in describing boundary conditions.
		Normal          = A cell with no special boundary conditions.
		Inactive        = A cell outside the domain (so we can ignore it).
		Dirichlet_left  = Main node is constant, other nodes are normal.
		Dirichlet_right = Main node is constant, other nodes are outside the domain.
		Neumann_left    = Main node uses phantom points on the left in computations, other nodes are normal.
		Neumann_right   = Main node uses phantom points on the right in computations, other nodes are outside the domain.
	*/
	enum class Cell_Type { Normal, Inactive, Dirichlet_left, Dirichlet_right, Neumann_left, Neumann_right };
	// The basic cell within the domain resulting from applying Lobatto IIIA-IIIB discretisation in space (x,y).
	// Each cell has stages_x * stages_y nodes.
	class Cell
	{
	public:
		std::pair<Cell_Type, Cell_Type> cell_type; // boundary condition types in x (first) and y (second) directions
		std::vector<double> U, V, U_xx, U_yy;
		Cell(std::pair<Cell_Type, Cell_Type> type) : cell_type(type) {};
	};

	// The integrator class.
	class Integrator
	{ // Functions and variables related to the PDE and integrator.
	private:
		unsigned int stages_x, stages_y; // Number of stages per cell.
		double step_size_x;
		double step_size_y;
		double step_size_time;
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
		inline void Compute_U_xx(Waves::Cell& cell, const Waves::Cell& x_left, const Waves::Cell& x_right);
		inline void Compute_U_yy(Waves::Cell& cell, const Waves::Cell& y_left, const Waves::Cell& y_right);

		// Components of the integrator.
		void Step_U(double step_size);
		void Step_V(double step_size);
		void update_U_dependent_variables();
		void update_V_dependent_variables();

	public:
		std::vector<Waves::Cell> cells; // Variables in the PDE (including derivatives and other useful values) for each cell.
		// We keep the following information arrays separate from the cells to minimise the memory footprint of the cells.
		std::vector<std::array<long, 4>> adjacency_information; // -1 indicates no adjacent cell (i.e. for boundary cells)
		std::vector<std::array<double, 2>> position_information;
		double Time = 0.0; // Time the system has evolved for.
		void Step(); // Step forwards in time, i.e., the integrator.
		float Change_Initial_Conditions(int ic = -1); // Defaults to switching to the next initial condition.
		void Change_Boundary_Conditions(int bc = -1); // Defaults to switching to the next initial condition.

		// Initialise the Integrator.
		float Initialise(int rx, int ry, double dt, int ic, int bc);

		// Friendly ostream operator<< so it can access private variables.
		friend std::ostream& operator<< (std::ostream& os, const Integrator& integrator);
	};

	// We need to use a single global instance of this class for OpenGL to have access to it. We declare it extern here and construct it in integrator.cpp.
	extern Integrator g_waves;

	// ostream operator for describing the Integator.
	std::ostream& operator<< (std::ostream& os, const Integrator& integrator);

}
