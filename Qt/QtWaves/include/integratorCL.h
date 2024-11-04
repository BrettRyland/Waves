#pragma once
///@file

#include <vector>
#include <array>
#include <iostream>
#include "boost_compute.h"

/** Our base namespace for the wave simulation
@note: we use cl_uint for indexing throughout as that is the type we're passing to glDrawElements in OGLWidget::paintGL()
*/
namespace Waves
{
	/** A differentiation of cell types for use in describing boundary conditions.
		Normal          0 = A cell with no special boundary conditions.
		Periodic_left	1 = The same as normal, but for a bottom or left edge cell that wraps around periodically
		Periodic_right	2 = The same as normal, but for a top or right edge cell that wraps around periodically.
		Dirichlet_left  3 = Main node is constant, other nodes are normal.
		Dirichlet_right 4 = Main node is constant, other nodes are outside the domain.
		Neumann_left    5 = Main node uses phantom points on the left in computations, other nodes are normal.
		Neumann_right   6 = Main node uses phantom points on the right in computations, other nodes are outside the domain.

		Note: we can't make this an enum class due to interoperability issues with boost.compute. Also, changes made here must be mirrored in the C99 version in integrator.cl.
	*/
	enum Cell_Type : cl_uint
	{
		Normal,
		Periodic_left,
		Periodic_right,
		Dirichlet_left,
		Dirichlet_right,
		Neumann_left,
		Neumann_right
	};
	const cl_uint cell_type_count = 7; ///< The number of different cell types in the enum

	const cl_uint missing_index = std::numeric_limits<cl_uint>::max(); ///< Value to use for missing indices

	typedef cl_float integrator_precision; ///< The precision to use in the integrator. @note This must be set here and a non-cl version in the OpenCL kernels in all .cl files!

	/// The integrator class.
	class Integrator
	{ // Functions and variables related to the PDE and integrator.
	private:
		/** @name Compute device context and queue
		 */
		boost::compute::context m_context;		 ///< The context of the compute device
		boost::compute::command_queue m_queue; ///< The queue associated with the context
		///@}

		/** @name Core integrator components
		These are the main variables used in the simulation (PDE).
		Size of U, V, U_xx and U_yy are stages_x*stages_y*num_cells
		@{*/
		boost::compute::vector<integrator_precision> U, ///< Surface height
				V,																					///< dU/dt
				U_xx,																				///< d^2U/dx^2
				U_yy;																				///< d^2U/dy^2
		void update_U_dependent_variables();						///< Update the U-dependent variables
		void update_V_dependent_variables();						///< Update the V-dependent variables
		///@}

		/** @name Integrator control variables
		@{*/
		cl_uint num_cells;																	///< The number of cells in the simulation
		cl_uint stages_x,																		///< Number of stages per cell in the x direction
				stages_y,																				///< Number of stages per cell in the y direction
				nodes_per_cell;																	///< The number of nodes per cell (rx-1)*(ry-1)
		integrator_precision time,													///< Time the system has evolved for.
				step_size_time,																	///< Step size in time
				step_size_x,																		///< Step size in the x direction
				step_size_y,																		///< Step size in the y direction
				wave_speed,																			///< Wave speed
				wave_speed_scale,																///< Wave speed scale
				height_scale;																		///< Height scale for the initial conditions
		int initial_conditions;															///< Initial conditions indicator
		int boundary_conditions;														///< Boudary conditions indicator
		integrator_precision m_artificial_dissipation{0.0}; ///< Artificial dissipation amount (default: 0.0 = no dissipation). Anything other than 0.0 breaks multisymplecticity and most likely stability of time integrator. Note: currently not enabled.
		///@}

		/** @name Integrator structural components
		Adjacency information, position information, cell types, update masks.
		@{*/
		boost::compute::vector<cl_uint> adjacency_information; ///< Adjacency information of each cell (left, right, down, up quadruplets). missing_index indicates no adjacent cell (i.e. for boundary cells)
		std::vector<cl_uint> host_adjacency_information;			 ///< Host copy of adjacency_information
		boost::compute::vector<cl_uint> masks;								 ///< An update mask dependent on the type of cell. Note: we can't use bool as the underlying type due to opencl issues with bools.
		void generate_update_masks();													 ///< Generate the update masks (different cell types only update certain nodes in the cell)
		boost::compute::vector<cl_uint> cell_type_x,					 ///< Cell type in the x direction
				cell_type_y;																			 ///< Cell type in the y direction
		std::vector<cl_uint> host_cell_type_x,								 ///< Host copy of cell_type_x
				host_cell_type_y;																	 ///< Host copy of cell_type_y
		// boost::compute::vector<integrator_precision> position_information; ///< Absolute position (x,y pairs) of the lower left corner of each cell
		std::vector<integrator_precision> host_position_information;																																																 ///< Host copy of position_information
		std::tuple<std::vector<cl_uint>, std::vector<cl_uint>, std::vector<cl_uint>> find_neighbours(std::vector<integrator_precision> const &position_information); ///< Helper function for constructing boundary conditions
		std::array<integrator_precision, 2> domain_scaling_factor;																																																	 ///< A scaling factor to get the graphics in the right place
		///@}

		/** @name Stencils
		@{*/
		boost::compute::vector<integrator_precision> coefficients_x, ///< Coefficient stencil used in the approximation of U_{xx}
				coefficients_y;																					 ///< Coefficient stencil used in the approximation of U_{yy}
		/** Construct the coefficients of the stencil
		\param r is the number of stages in the stencil
		*/
		std::vector<std::vector<integrator_precision>> get_coefficients(int r);
		///@}

		/** @name Coordinates
		Position within each cell (from 0 to 1) of each node {main node, inner node1, inner node2, ...}
		@{ */
		std::vector<integrator_precision> coords_x, ///< x-coordinates
				coords_y;																///< y-coordinates
		/** Setup the coordinates for a given number of stage values
		\param r is the number of stage values
		*/
		std::vector<integrator_precision> setup_coords(int r);
		///@}

		/** @name OpenCL kernels
		@{*/
		boost::compute::kernel kernel_step_U, ///< Take a step in the U variable
				kernel_step_V,										///< Take a step in the V variable
				kernel_compute_U_xx,							///< Compute the second spatial derivative of U in the x-direction
				kernel_compute_U_yy;							///< Compute the second spatial derivative of U in the y-direction

		void generate_kernels(); ///< Generate the kernels required for the simulation

		void configure_kernels_static_vars(std::vector<cl_uint> const &host_cell_type_x, std::vector<cl_uint> const &host_cell_type_y, std::vector<cl_uint> const &host_adjacency_information, std::vector<integrator_precision> const &host_position_information); ///< Update the kernel configuration (arguments) whenever the structure of the simulation is changed. E.g. a change in boundary conditions.

		void configure_kernels_dynamic_vars(std::vector<integrator_precision> const &host_U, std::vector<integrator_precision> const &host_V); ///< Update the kernel dynamic variables (arguments). E.g. a change in initial conditions.

		///@}

	public:
		/** @name Core integrator methods
		@{*/
		void step();			///< Step forwards in time, i.e., the integrator
		void half_step(); ///< Perform an initial half-step in V to optimise the leapfrog algorithm
		///@}

		/** @name Integrator configuration
		@{*/
		/** Initialise the Integrator.
		@param rx is the number of stages in the x direction
		@param ry is the number of stages in the y direction
		@param dt is the step size in time
		@param ic is the initial condition
		@param bc is the boundary condition
		@param shared_context is the context of the OpenCL/OpenGL compute device
		@param queue is the queue associated with the shared_context
		*/
		void initialise(int rx, int ry, integrator_precision dt, int ic, int bc, boost::compute::context &shared_context, boost::compute::command_queue &queue);
		void change_boundary_conditions(int bc = -1);																				///< Defaults to switching to the next initial condition. Also resets the initial conditions and returns the initial conditions hint.
		void change_initial_conditions(int ic = -1);																				///< Defaults to switching to the next initial condition
		void reset_initial_conditiones() { change_initial_conditions(initial_conditions); } ///< Reset the simulation to the current initial conditions
		void change_time_step(integrator_precision dt)
		{
			step_size_time = dt;
			kernel_step_U.set_arg(7, step_size_time);
			kernel_step_V.set_arg(9, step_size_time);
		} ///< Change the time step of the integrator
		void set_artificial_dissipation(integrator_precision value)
		{
			m_artificial_dissipation = value /*std::clamp(value, 0.0, 1.0) std::clamp is not part of gcc-6 in linux*/;
			kernel_step_U.set_arg(8, m_artificial_dissipation);
		} ///< Set the amount of artificial dissipation in the integrator
		void change_wave_speed(integrator_precision scale)
		{
			wave_speed_scale = scale;
			kernel_step_V.set_arg(10, wave_speed * wave_speed_scale);
		} ///< Change the wave speed in the integrator
		void change_height_scale(integrator_precision scale)
		{
			height_scale = scale;
			if (time == static_cast<integrator_precision>(0.0))
				change_initial_conditions(initial_conditions);
		} ///< Set the height scale for initial conditions
		///@}

		/** @name Accessors
		@{*/
		auto get_IC() const { return initial_conditions; }		 ///< Get the current initial conditions
		auto get_BC() const { return boundary_conditions; }		 ///< Get the current boundary conditions
		auto get_time_step() const { return step_size_time; }	 ///< Get the current time step of the integrator
		auto get_time() const { return time; }								 ///< Get the current time reported by the integrator
		auto get_number_of_cells() const { return num_cells; } ///< Get the number of cells in the simulation
		auto get_wave_speed() const { return wave_speed; }		 ///< Get the wave speed of the current boundary conditions
		///@}

		/** @name Device to host accessors
		@note The cell types returned are vectors of cl_uint instead of Cell_Type due to copying from the device not allowing direct copying to type Cell_Type
		@{*/
		auto get_cell_type_x() -> std::vector<cl_uint>; ///< Get the cell types in the x-direction
		auto get_cell_type_y() -> std::vector<cl_uint>; ///< Get the cell types in the y-direction
		// auto get_position_information()->std::vector<integrator_precision>; ///< Get the position information
		auto get_adjacency_information() -> std::vector<cl_uint>; ///< Get the adjacency information
		///@}

		/// Friendly ostream operator<<
		friend std::ostream &operator<<(std::ostream &os, const Integrator &integrator);

		/// Friend class Mesh so that it can access the raw data
		friend class Mesh;
	};
}
