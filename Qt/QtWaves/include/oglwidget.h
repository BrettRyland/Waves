#pragma once
///@file

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
#include <QMatrix4x4>

#include "ui_oglwidget.h"
#include "mesh.h"
#include "integrator_wrapper.h"

namespace Waves {

	/** Modelview helper class
	Translation, rotation and scaling are stored in separate matrices and only combined with the base matrix when update() is called.
	translate(), rotate() and scale() each update their own matrix.
	reset() sets the translation, rotation and scaling matrices to the identity and the modelview matrix to the base matrix.
	*/
	class Modelview
	{
	public:
		/**@name 4x4 transformation matrices
		These default to the identity
		@{*/
		QMatrix4x4 base, ///< The intial viewing position
			translation, ///< Translation
			rotation,    ///< Rotation
			scaling,     ///< Scaling
			modelview;   ///< The current viewing position given by base*translation*rotation*scaling
		///@}
		/// Combine the various components to get the final modelview matrix
		inline void update() { modelview = base*translation*rotation*scaling; }
		/// Translate by the amount (x,y,z)
		inline void translate(float x, float y, float z) { translation.translate(x, y, z); }
		/// Rotate about the axis (x,y,z) by the amount angle
		inline void rotate(float angle, float x, float y, float z) { rotation.rotate(angle, x, y, z); }
		/// Scale by the amount (x,y,z)
		inline void scale(float x, float y, float z) { scaling.scale(x, y, z); }
		/// Reset to the initial view
		inline void reset() { translation.setToIdentity(); rotation.setToIdentity(); scaling.setToIdentity(); update(); }
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
		inline void reset_view() { m_modelview.reset(); QOpenGLWidget::update(); } ///< Revert to the default view
		void change_initial_conditions(); ///< Switch to the next set of initial conditions
		void change_boundary_conditions(); ///< Switch to the next set of boundary conditions
		inline bool is_paused() { return integrator_wrapper.is_paused(); } ///< Return the paused state
		void quit(); ///< Quit

	signals:
		void pause_integrator(); ///< Pause the integrator
		void unpause_integrator(); ///< Unpause the integrator
		void toggle_paused_integrator(); ///< Toggle the pause state
		void modify_integrator(); ///< Used to signal to the integrator thread that we are going to change something important (i.e. avoids concurrency problems)
		void toggle_fullscreen(); ///< Toggle fullscreen mode @todo{This is not implemented yet.}

	private:
		/**@name Integrator information
		@{ */
		Waves::Mesh mesh; ///< The mesh for rendering
		Integrator_Wrapper integrator_wrapper{ mesh, this }; ///< A wrapper around the integrator to allow it to run in its own thread
		void run_integrator(); ///< Actually set the integrator running in its thread
		///@}

		/**@name OpenGL State Information
		@{ */
		QOpenGLBuffer m_vertex_buffer{ QOpenGLBuffer::VertexBuffer }, ///< The vertex buffer
			m_index_buffer{ QOpenGLBuffer::IndexBuffer }; ///< The index buffer
		QOpenGLVertexArrayObject m_vertex_array_object; ///< The VAO
		QOpenGLShaderProgram *m_shader_program; ///< The shader program
		///@}

		/**@name Shader information
		@{ */
		int u_projection; ///< Projection matrix shader handle
		int u_modelview; ///< Modelview shader handle
		QMatrix4x4 m_projection; ///< Projection matrix
		Modelview m_modelview; ///< Modelview instance
		///@}

		/**@name OpenGL Helpers
		@{ */
		void initializeGL(); ///< Set global OpenGL flags and initialise variables required for OpenGL
		void resizeGL(int width, int height); ///< Update the projection matrix when the widget is resized
		void paintGL(); ///< Actually render the surface using our shader program
		void teardownGL(); ///< Destroy our OpenGL information
		void make_shader_program(); ///< Build and link the shader programs from the glsl files
		void update_vertex_buffer(); ///< Update the vertex buffer based on the mesh data
		void rebuild_vertex_array_object(); ///< Rebuild the entire vertex array object after the structure of the mesh has been modified (not just updated values)
		///@}

		/**@name Mouse and key event helpers and handlers
		@{ */
		QPoint mouse_press_location; ///< Location of a mouse click
		Qt::MouseButton mouse_button_pressed; ///< Which mouse button was pressed
		void keyPressEvent(QKeyEvent *event); ///< Handle key press events
		void mousePressEvent(QMouseEvent *event); ///< Handle mouse press events
		void mouseReleaseEvent(QMouseEvent *event); ///< Hanlde mouse release events
		void mouseMoveEvent(QMouseEvent *event); ///< Handle mouse drag events
		void wheelEvent(QWheelEvent *event); ///< Handle mouse wheel events
		void mouseDoubleClickEvent(QMouseEvent *event); ///< Handle mouse double-click events
		///@}
	};
}