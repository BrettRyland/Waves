#pragma once
///@file

#include <QThread>
#include "integratorCL.h"
#include "mesh.h"

namespace Waves
{
	/** A wrapper around the Integrator class that allows the integrator to run in its own thread and provides controls for pausing, modifying and quiting.
	@note The slots in this class are executed in the parent thread, not the same thread as run().
	*/
	class Integrator_Wrapper : public QThread
	{
		Q_OBJECT
		/// Our modified run() function. Called by QThread::thread.start().
		void run() Q_DECL_OVERRIDE;

	public:
		/** Initialise the integrator and mesh
		@param dt is the time step
		@param rx is the number of stages in the integrator in the x-direction
		@param ry is the number of stages in the integrator in the y-direction
		@param ic is the initial conditions
		@param bc is the boundary conditions
		@param shared_context is the shared context between OpenGL and OpenCL
		@param queue is the queue associated with the shared context
		*/
		void initialise(double dt, int rx, int ry, int ic, int bc, boost::compute::context &shared_context, boost::compute::command_queue &queue);
		// integrator
		void change_boundary_conditions(int bc);   ///< Change the boundary conditions and generate the static components of a new surface mesh. The vertex buffer and mesh will need to be re-connected and updated after calling this.
		void change_initial_conditions(int ic);	   ///< Change the initial conditions and update the surface mesh (dynamic components)
		void set_artificial_dissipation(double d); ///< Set the artificial dissipation in the integrator
		void change_time_step(double dt);		   ///< Change the time step of the integrator
		void reverse_time();					   ///< Take 2 half-steps to set up the integrator for going backwards
		// mesh
		void connect_mesh_with_vertex_buffer(QOpenGLBuffer &vertex_buffer); ///< Connect the mesh with the OpenGL vertex buffer
		void update_surface_mesh();											///< Update the surface mesh ("dynamic" components

		/** Accessors
		@{*/
		inline bool const is_paused() const { return m_paused; } ///< Check our paused state
		auto const get_steps_taken()
		{
			auto steps_taken = m_step_counter;
			m_step_counter = 0;
			return steps_taken;
		} ///< Returns the number of steps taken since the last time this function was called.
		// integrator pass-through accessors
		auto const get_IC() const { return m_integrator.get_IC(); }		///< Get the current initial conditions
		auto const get_BC() const { return m_integrator.get_BC(); }		///< Get the current boundary conditions
		auto const get_time() const { return m_integrator.get_time(); } ///< Get the time reported by the integrator
		// mesh pass-through accessors
		auto const get_mesh_vertex_data() const { return m_mesh.get_vertex_data(); }		   ///< Get the mesh vertex data
		auto const get_mesh_vertex_data_size() const { return m_mesh.get_vertex_data_size(); } ///< Get the size of the mesh vertex data
		auto const get_mesh_index_data() const { return m_mesh.get_index_data(); }			   ///< Get the mesh index data
		auto const get_mesh_index_data_size() const { return m_mesh.get_index_data_size(); }   ///< Get the size of the mesh index data
		///@}

	public slots:
		inline void pause() { m_paused = true; }			  ///< Pause the integrator
		inline void unpause() { m_paused = false; }			  ///< Unpause the integrator
		inline void toggle_paused() { m_paused = !m_paused; } ///< Toggle paused state
		inline void quit() { m_quit = true; }				  ///< Quit
		void modify_integrator();							  ///< Signal that we are going to modify the integrator. modify_integrator() then waits until the integrator is in a safe state for modifying before returning.

	signals:
		void result_ready(); ///< Signal that an iteration of the integrator has been performed and is ready to be displayed.

	private:
		Integrator m_integrator;  ///< Our integrator
		Mesh m_mesh;			  ///< Our mesh
		bool m_paused{true};	  ///< Flag to indicate that we should pause on the next loop of run()
		bool m_quit{false};		  ///< Flag to indicate that we should quit on the next loop of run()
		bool m_modifiable{true};  ///< Flag to indicate when it is safe to modify the integrator
		size_t m_step_counter{0}; ///< A step counter for measuring the number of steps taken per frame rendered
	};
}
