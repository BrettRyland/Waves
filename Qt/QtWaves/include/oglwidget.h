#pragma once
///@file

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QMatrix4x4>

#include <chrono>

#include "ui_oglwidget.h"
#include "integrator_wrapper.h"

// Fix #define bug in xlib.h header.
#undef Bool

namespace Waves
{

	/** Modelview helper class
	Translation, rotation and scaling are stored in separate matrices and only combined with the base matrix when update() is called.
	translate(), rotate() and scale() each update their own matrix.
	reset() sets the translation, rotation and scaling matrices to the identity and the modelview matrix to the base matrix.
	*/
	class Modelview
	{
	public:
		/** @name 4x4 transformation matrices
		These default to the identity
		@{*/
		QMatrix4x4 base, ///< The intial viewing position
				translation, ///< Translation
				rotation,		 ///< Rotation
				scaling,		 ///< Scaling
				modelview;	 ///< The current viewing position given by base*translation*rotation*scaling
		///@}
		/// Combine the various components to get the final modelview matrix
		inline void update() { modelview = base * translation * rotation * scaling; }
		/// Translate by the amount (x,y,z)
		inline void translate(float x, float y, float z) { translation.translate(x, y, z); }
		/// Rotate about the axis (x,y,z) by the amount angle
		inline void rotate(float angle, float x, float y, float z) { rotation.rotate(angle, x, y, z); }
		/// Scale by the amount (x,y,z)
		inline void scale(float x, float y, float z) { scaling.scale(x, y, z); }
		/// Reset to the initial view
		inline void reset()
		{
			translation.setToIdentity();
			rotation.setToIdentity();
			scaling.setToIdentity();
			update();
		}
	};

	/** OpenGLWidget class
	 */
	class OGLWidget : public QOpenGLWidget, protected QOpenGLFunctions
	{
		Q_OBJECT

	public:
		/// Constructor
		OGLWidget(QWidget *parent = Q_NULLPTR);
		/// Destructor
		~OGLWidget();

	public slots:
		void update(); ///< Update the display
		inline void reset_view()
		{
			m_modelview.reset();
			QOpenGLWidget::update();
		} ///< Revert to the default view
		void change_initial_conditions(int ic = -1);																					///< Switch to the next set of initial conditions
		void change_boundary_conditions(int bc = -1);																					///< Switch to the next set of boundary conditions
		inline bool is_paused() { return m_integrator_wrapper.is_paused(); }									///< Return the paused state
		void quit();																																					///< Quit
		void change_dissipation(double value);																								///< Change the value of the dissipation
		void change_timestep(double timestep);																								///< Change the value of the timestep
		void change_wave_speed(double speed);																									///< Change the wave speed
		void change_height_scale(double scale);																								///< Change the height scale
		void reset_integrator() { change_initial_conditions(m_integrator_wrapper.get_IC()); } ///< Reset the integrator to the initial conditions
		double get_time() const { return m_integrator_wrapper.get_time(); }										///< Get the time reported by the simulation (i.e. the time step * number of steps)
		double get_spf() { return m_integrator_wrapper.get_steps_taken(); }										///< Get the number of steps taken since last time this function was called
		double get_fps() { return m_fps; }																										///< Get the smoothed FPS.

	signals:
		void pause_integrator();					///< Pause the integrator
		void unpause_integrator();				///< Unpause the integrator
		void toggle_paused_integrator();	///< Toggle the pause state
		void notify_paused_state();				///< Notify the main window of the paused state
		void notify_IC_changed(int ic);		///< Notify the main window of the IC change
		void notify_BC_changed(int bc);		///< Notify the main window of the BC change
		void modify_integrator();					///< Used to signal to the integrator thread that we are going to change something important (i.e. avoids concurrency problems)
		void toggle_fullscreen();					///< Toggle fullscreen mode @todo{This is not implemented yet.}
		void update_time();								///< Signal to the main ui to update the displayed time
		void notify_time_step(double dt); ///< Notify the main window of the time step

	private:
		/**@name Boost compute information
		@{ */
		boost::compute::context m_shared_context; ///< The compute context
		boost::compute::command_queue m_queue;		///< The compute queue
		///@}

		/**@name Integrator information
		@{ */
		Integrator_Wrapper m_integrator_wrapper; ///< A wrapper around the integrator to allow it to run in its own thread
		void run_integrator();									 ///< Actually set the integrator running in its thread
		///@}

		/**@name OpenGL State Information
		@{ */
		QOpenGLBuffer m_vertex_buffer{QOpenGLBuffer::VertexBuffer}, ///< The vertex buffer
				m_index_buffer{QOpenGLBuffer::IndexBuffer};							///< The index buffer
		QOpenGLVertexArrayObject m_vertex_array_object;							///< The VAO
		QOpenGLShaderProgram *m_shader_program;											///< The shader program
		///@}

		/**@name Shader information
		@{ */
		int u_projection;				 ///< Projection matrix shader handle
		int u_modelview;				 ///< Modelview shader handle
		QMatrix4x4 m_projection; ///< Projection matrix
		Modelview m_modelview;	 ///< Modelview instance
		///@}

		/**@name OpenGL Helpers
		@{ */
		void initializeGL();									///< Set global OpenGL flags and initialise variables required for OpenGL
		void resizeGL(int width, int height); ///< Update the projection matrix when the widget is resized
		void paintGL();												///< Actually render the surface using our shader program
		void teardownGL();										///< Destroy our OpenGL information
		void make_shader_program();						///< Build and link the shader programs from the glsl files
		void update_vertex_buffer();					///< Update the vertex buffer based on the mesh data
		void rebuild_vertex_array_object();		///< Rebuild the entire vertex array object after the structure of the mesh has been modified (not just updated values)
		bool m_render_frame{true};						/// Ignore frames that are generated too fast
		///@}

		/**@name Mouse and key event helpers and handlers
		@{ */
		QPoint m_mouse_press_location;									///< Location of a mouse click
		Qt::MouseButton m_mouse_button_pressed;					///< Which mouse button was pressed
		void keyPressEvent(QKeyEvent *event);						///< Handle key press events
		void mousePressEvent(QMouseEvent *event);				///< Handle mouse press events
		void mouseReleaseEvent(QMouseEvent *event);			///< Hanlde mouse release events
		void mouseMoveEvent(QMouseEvent *event);				///< Handle mouse drag events
		void wheelEvent(QWheelEvent *event);						///< Handle mouse wheel events
		void mouseDoubleClickEvent(QMouseEvent *event); ///< Handle mouse double-click events
		///@}

		double m_max_fps{60.0};																												 ///< The maximum framerate at which to update the OpenGL widget
		double m_fps{m_max_fps};																											 ///< Smoothed framerate.
		std::chrono::steady_clock::time_point m_tic{std::chrono::steady_clock::now()}; ///< Time tracker.
	};
}