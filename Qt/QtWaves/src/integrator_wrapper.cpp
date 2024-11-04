///@file

#include "integrator_wrapper.h"

namespace Waves
{
	void Integrator_Wrapper::initialise(double dt, int rx, int ry, int ic, int bc, boost::compute::context &shared_context, boost::compute::command_queue &queue)
	{
		m_integrator.initialise(rx, ry, static_cast<integrator_precision>(dt), ic, bc, shared_context, queue);
		m_mesh.generate_surface_mesh(m_integrator);
		m_mesh.generate_mesh_update_kernel(m_integrator);
	}

	// Our modified run() function.
	void Integrator_Wrapper::run()
	{
		while (!m_quit)
		{
			if (!m_paused)
			{
				// Integrator is running, so it's not safe to modify it.
				m_modifiable = false;
				// Actually perform one step of the integrator.
				m_integrator.step();
				++m_step_counter;
				// Signal that we're done.
				emit result_ready();
			}
			else
			{
				// Integrator is paused, so it's safe to modify it now.
				m_modifiable = true;
				// Also, we can sleep for a bit to conserve CPU resources. (The run() function runs in this thread.)
				QThread::msleep(100);
			}
		}
	}

	// Signal that we're going to modify the integrator.
	void Integrator_Wrapper::modify_integrator()
	{
		// Pause the integrator.
		m_paused = true;
		// Wait until it's safe to modify the integrator by telling the calling thread to sleep for a bit (this is a slot in a single run QThread class, so it runs in the parent thread).
		while (!m_modifiable) // FIXME: this should probably now be done with m_queue.finish()
			QThread::msleep(100);
	}

	void Integrator_Wrapper::change_initial_conditions(int ic)
	{
		m_integrator.change_initial_conditions(ic);
		m_mesh.update_surface_mesh_CL(m_integrator);
	}

	void Integrator_Wrapper::change_boundary_conditions(int bc)
	{
		m_integrator.change_boundary_conditions(bc);
		m_mesh.generate_surface_mesh(m_integrator);
	}

	void Integrator_Wrapper::set_artificial_dissipation(double d)
	{
		m_integrator.set_artificial_dissipation(static_cast<integrator_precision>(d));
	}

	void Integrator_Wrapper::change_time_step(double dt)
	{
		m_integrator.change_time_step(static_cast<integrator_precision>(dt));
	}

	void Integrator_Wrapper::change_wave_speed(double speed)
	{
		m_integrator.change_wave_speed(static_cast<integrator_precision>(speed));
	}

	void Integrator_Wrapper::change_height_scale(double scale)
	{
		m_integrator.change_height_scale(static_cast<integrator_precision>(scale));
	}

	void Integrator_Wrapper::reverse_time()
	{
		m_integrator.change_time_step(-m_integrator.get_time_step());
		m_integrator.half_step();
		m_integrator.half_step();
	}

	void Integrator_Wrapper::connect_mesh_with_vertex_buffer(QOpenGLBuffer &vertex_buffer)
	{
		m_mesh.configure_mesh_update_kernel(m_integrator, vertex_buffer);
		m_mesh.update_surface_mesh_CL(m_integrator);
	}

	void Integrator_Wrapper::update_surface_mesh()
	{
		m_mesh.update_surface_mesh_CL(m_integrator);
	}
}
