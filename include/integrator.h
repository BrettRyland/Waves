#pragma once
///@file

#include <vector>
#include <array>
#include <iostream>

/** Our base namespace for the wave simulation
*/
namespace Waves {
	/** A differentiation of cell types for use in describing boundary conditions.
		Normal          = A cell with no special boundary conditions.
		Periodic		= The same as normal, but for a top or right edge cell that wraps around periodically (needed for rendering purposes).
		Dirichlet_left  = Main node is constant, other nodes are normal.
		Dirichlet_right = Main node is constant, other nodes are outside the domain.
		Neumann_left    = Main node uses phantom points on the left in computations, other nodes are normal.
		Neumann_right   = Main node uses phantom points on the right in computations, other nodes are outside the domain.
	*/
	enum class Cell_Type { Normal, Periodic, Dirichlet_left, Dirichlet_right, Neumann_left, Neumann_right };
	
	const double PI = 3.141592653589793; ///< Numeric constant
	const unsigned int missing_index = std::numeric_limits<unsigned int>::max(); ///< Value to use for missing indices
	
	/** The basic cell within the domain resulting from applying Lobatto IIIA-IIIB discretisation in space (x,y).
	Each cell has stages_x * stages_y nodes.
	*/
	class Cell
	{
	public:
		std::pair<Cell_Type, Cell_Type> cell_type; ///< boundary condition types in x (first) and y (second) directions
		std::vector<double> U, ///< Wave surface height
			V, ///< dU/dt
			U_xx, ///< d^2U/dx^2
			U_yy; ///< d^2U/dy^2
		/// Constructor
		Cell(std::pair<Cell_Type, Cell_Type> type) : cell_type(type) {};
		Cell() = default; ///< @overload
	};

	/// The integrator class.
	class Integrator
	{ // Functions and variables related to the PDE and integrator.
	private:
		unsigned int stages_x, ///< Number of stages per cell in the x direction.
			stages_y; ///< Number of stages per cell in the y direction.
		double wave_speed; ///< Wave speed
		int initial_conditions; ///< Initial conditions indicator
		int boundary_conditions; ///< Boudary conditions indicator

		/** @name Stencils
		@{*/
		std::vector<std::vector<double>> coefficients_x, ///< Coefficient stencil used in the approximation of U_{xx}
			coefficients_y; ///< Coefficient stencil used in the approximation of U_{yy}
		/** Construct the coefficients of the stencil
		\param r is the number of stages in the stencil
		*/
		std::vector<std::vector<double>> Setup_Coefficients(int r);
		///@}

		/** @name Coordinates
		Position within each cell (from 0 to 1) of each node {main node, inner node1, inner node2, ...}
		@{ */
		std::vector<double> coords_x, ///< x-coordinates
			coords_y; ///< y-coordinates
		/** Setup the coordinates for a given number of stage values
		\param r is the number of stage values
		*/
		std::vector<double> Setup_Coords(int r);
		///@}

		void Step_U(double step_size); ///< Take a step in U
		void Step_V(double step_size); ///< Take a step in V
		void update_U_dependent_variables(); ///< Update the U-dependent variables
		void update_V_dependent_variables(); ///< Update the V-dependent variables
		void Compute_U_xx(); ///< Compute second-order spatial derivative of U in the x direction.
		void Compute_U_yy(); ///< Compute second-order spatial derivative of U in the y direction.

		std::array<double, 2> domain_scaling_factor; ///< A scaling factor to get the graphics in the right place
		void Find_Neighbours(); ///< Helper function for constructing boundary conditions

	public:
		std::vector<Waves::Cell> cells; ///< Variables in the PDE (including derivatives and other useful values) for each cell.
		// We keep the following information arrays separate from the cells to minimise the memory footprint of the cells.
		std::vector<std::array<unsigned int, 4>> adjacency_information; ///< Adjacency information of each cell (left, right, down, up). missing_index indicates no adjacent cell (i.e. for boundary cells)
		std::vector<std::array<double, 2>> position_information; ///< Absolute position of the lower left corner of each cell
		double Time = 0.0; ///< Time the system has evolved for.
		double step_size_time; ///< Step size in time
		double step_size_x; ///< Step size in the x direction
		double step_size_y; ///< Step size in the y direction

		/// Step forwards in time, i.e., the integrator.
		void Step();
		/// Perform an initial half-step in V to optimise the leapfrog algorithm.
		void Half_Step() { Step_V(0.5*step_size_time); };

		/** Initialise the Integrator.
		\param rx is the number of stages in the x direction
		\param ry is the number of stages in the y direction
		\param dt is the step size in time
		\param ic is the initial condition
		\param bc is the boundary condition
		*/
		float Initialise(int rx, int ry, double dt, int ic, int bc);

		float Change_Initial_Conditions(int ic = -1); ///< Defaults to switching to the next initial condition.
		float Change_Boundary_Conditions(int bc = -1); ///< Defaults to switching to the next initial condition. Also resets the initial conditions and returns the initial conditions hint.

		/// Friendly ostream operator<<
		friend std::ostream& operator<< (std::ostream& os, const Integrator& integrator);
	};

	/// We need to use a single global instance of this class for OpenGL to have access to it. We declare it extern here and construct it in integrator.cpp.
	extern Integrator g_waves;

}
