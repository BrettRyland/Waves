#pragma once
#include <QThread>
#include "integrator.h"
#include "mesh.h"

namespace Waves {
	/** Integrator_Wrapper class.
	This class allows the integrator to run in its own thread and provides controls for pausing, modifying and quiting.
	*/
	class Integrator_Wrapper : public QThread
	{
		Q_OBJECT
		// Our modified run() function. Called by QThread::thread.start().
		void run() Q_DECL_OVERRIDE;

	public:
		// A modifiable Mesh object is required in the constructor so that we can update it.
		Integrator_Wrapper(Mesh &mesh, QObject *parent = Q_NULLPTR);
		~Integrator_Wrapper();
		inline bool is_paused() { return m_paused; }

	public slots: // Note: these slots are executed in the parent thread, not the same thread as run().
		inline void pause() { m_paused = true; }
		inline void unpause() { m_paused = false; }
		inline void toggle_paused() { m_paused = !m_paused; }
		inline void quit() { m_quit = true; }
		// Signal that we are going to modify the integrator. modify_integrator() then waits until the integrator is in a safe state for modifying before returning.
		void modify_integrator();

	signals:
		// Signal that an iteration of the integrator has been performed and is ready to be displayed.
		void resultReady();

	private:
		bool m_paused{ true };
		bool m_quit{ false };
		bool m_modifiable{ true };
		Mesh *m_mesh; // Our reference to the mesh.
	};
}
