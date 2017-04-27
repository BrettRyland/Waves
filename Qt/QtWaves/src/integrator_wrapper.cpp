///@file

#include "integrator_wrapper.h"

namespace Waves
{
	Integrator_Wrapper::Integrator_Wrapper(Mesh& mesh, QObject *parent) : m_mesh(&mesh)
	{
		// A default initialisation with some nice values for the integrator.
		// TODO: These values should be settable through the mainwindow interface and shouldn't need setting here.
		double dt = 0.01667; // 1/60 (aiming for 60fps once timing is implemented)
		int rx = 2; // Number of stages in x.
		int ry = 2; // Number of stages in y.
		int IC = 4; // Continuous spectrum of frequencies in a central localised hump.
		int BC = 3; // Square Dirichlet (i.e. fixed) boundary.
		g_waves.Initialise(rx, ry, dt, IC, BC);
	}

	Integrator_Wrapper::~Integrator_Wrapper()
	{
	}

	// Our modified run() function.
	void Integrator_Wrapper::run()
	{
		while (!m_quit) {
			if (!m_paused) {
				// Integrator is running, so it's not safe to modify it.
				m_modifiable = false;
				// Actually perform one step of the integrator.
				g_waves.Step();
				// Update the surface mesh heights and normals.
				m_mesh->update_surface_mesh();
				// Signal that we're done.
				emit resultReady();
			}
			else {
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
		while (!m_modifiable)
			QThread::msleep(100);
	}

}
