#pragma once
///@file

#include <QThread>
#include "integrator.h"
#include "mesh.h"

namespace Waves {
	/** A wrapper around the Integrator class that allows the integrator to run in its own thread and provides controls for pausing, modifying and quiting.
	@note The slots in this class are executed in the parent thread, not the same thread as run().
	*/
	class Integrator_Wrapper : public QThread
	{
		Q_OBJECT
			/// Our modified run() function. Called by QThread::thread.start().
			void run() Q_DECL_OVERRIDE;

	public:
		/** Constructor
		A modifiable Mesh object is required in the constructor so that we can update it.
		*/
		Integrator_Wrapper(Mesh &mesh, QObject *parent = Q_NULLPTR);
		~Integrator_Wrapper(); ///< Destructor
		inline bool is_paused() { return m_paused; } ///< Check our paused state

		public slots:
		inline void pause() { m_paused = true; } ///< Pause the integrator
		inline void unpause() { m_paused = false; } ///< Unpause the integrator
		inline void toggle_paused() { m_paused = !m_paused; } ///< Toggle paused state
		inline void quit() { m_quit = true; } ///< Quit
		/// Signal that we are going to modify the integrator. modify_integrator() then waits until the integrator is in a safe state for modifying before returning.
		void modify_integrator();

	signals:
		/// Signal that an iteration of the integrator has been performed and is ready to be displayed.
		void resultReady();

	private:
		bool m_paused{ true }; ///< Flag to indicate that we should pause on the next loop of run()
		bool m_quit{ false }; ///< Flag to indicate that we should quit on the next loop of run()
		bool m_modifiable{ true }; ///< Flag to indicate when it is safe to modify the integrator
		Mesh *m_mesh; ///< Our reference to the mesh.
	};
}
