#pragma once

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
		QMatrix4x4 base, translation, rotation, scaling, modelview; // defaults to identity
		// Combine the various components to get the final modelview matrix
		inline void update() { modelview = base*translation*rotation*scaling; }
		// translate, rotate and scale
		inline void translate(float x, float y, float z) { translation.translate(x, y, z); }
		inline void rotate(float angle, float x, float y, float z) { rotation.rotate(angle, x, y, z); }
		inline void scale(float x, float y, float z) { scaling.scale(x, y, z); }
		// Reset to the initial view
		inline void reset() { translation.setToIdentity(); rotation.setToIdentity(); scaling.setToIdentity(); update(); }
	};

	/** OpenGLWidget class
	*/
	class OGLWidget : public QOpenGLWidget, protected QOpenGLFunctions
	{
		Q_OBJECT

	public:
		OGLWidget(QWidget *parent = Q_NULLPTR);
		~OGLWidget();

	public slots:
		void update(); // Update the display
		inline void reset_view() { m_modelview.reset(); QOpenGLWidget::update(); } // Revert to the default view
		void change_initial_conditions(); // Switch to the next set of initial conditions
		void change_boundary_conditions(); // Switch to the next set of boundary conditions
		inline bool is_paused() { return integrator_wrapper.is_paused(); } // Return the paused state
		void quit();

	signals:
		void pause_integrator();
		void unpause_integrator();
		void toggle_paused_integrator();
		void modify_integrator(); // Used to signal to the integrator thread that we are going to change something important (i.e. avoids concurrency problems)
		void toggle_fullscreen();

	private:
		// *** Integrator information *** //
		Waves::Mesh mesh; // The mesh for rendering
        Integrator_Wrapper integrator_wrapper{ mesh }; // A wrapper around the integrator to allow it to run in its own thread
		void run_integrator(); // Actually set the integrator running in its thread

		// *** OpenGL State Information *** //
		QOpenGLBuffer m_vertex_buffer{ QOpenGLBuffer::VertexBuffer }, m_index_buffer{ QOpenGLBuffer::IndexBuffer };
		QOpenGLVertexArrayObject m_vertex_array_object;
		QOpenGLShaderProgram *m_shader_program;

		// *** Shader information *** //
		int u_projection;
		int u_modelview;
		QMatrix4x4 m_projection;
		Modelview m_modelview;

		// *** OpenGL Helpers *** //
		void initializeGL();
		void resizeGL(int width, int height);
		void paintGL();
		void teardownGL();
		void make_shader_program();
		void update_vertex_buffer();
		void rebuild_vertex_array_object();

		// *** Mouse and key event helpers and handlers *** //
		QPoint mouse_press_location;
		Qt::MouseButton mouse_button_pressed;
		void keyPressEvent(QKeyEvent *event);
		void mousePressEvent(QMouseEvent *event);
		void mouseReleaseEvent(QMouseEvent *event);
		void mouseMoveEvent(QMouseEvent *event);
		void wheelEvent(QWheelEvent *event);
		void mouseDoubleClickEvent(QMouseEvent *event);
	};
}
