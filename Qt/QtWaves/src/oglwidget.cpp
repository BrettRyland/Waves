///@file

#include <QKeyEvent>
#include <QThread>
#include "oglwidget.h"
#include "mesh.h"
#include "integrator_wrapper.h"

namespace Waves {
	OGLWidget::OGLWidget(QWidget *parent)
		: QOpenGLWidget(parent)
	{
		// Set the default modelview.
		m_modelview.base.setToIdentity();
		m_modelview.base.translate(0.0f, 0.0f, -30.0f);
		m_modelview.base.rotate(-60.0f, 1.0f, 0.0f, 0.0f);
		m_modelview.update();

		// Start the integrator (starts paused).
		run_integrator();
	}

	OGLWidget::~OGLWidget()
	{
		makeCurrent();
		teardownGL();
	}

	// *** Integrator information *** //

	// Run the simulation in a separate thread
	void OGLWidget::run_integrator()
	{
		QObject::connect(this, &OGLWidget::pause_integrator, &integrator_wrapper, &Integrator_Wrapper::pause);
		QObject::connect(this, &OGLWidget::unpause_integrator, &integrator_wrapper, &Integrator_Wrapper::unpause);
		QObject::connect(this, &OGLWidget::toggle_paused_integrator, &integrator_wrapper, &Integrator_Wrapper::toggle_paused);
		QObject::connect(this, &OGLWidget::modify_integrator, &integrator_wrapper, &Integrator_Wrapper::modify_integrator);
		QObject::connect(&integrator_wrapper, &Integrator_Wrapper::resultReady, this, &OGLWidget::update);
		QObject::connect(&integrator_wrapper, &Integrator_Wrapper::finished, &integrator_wrapper, &QObject::deleteLater);
		integrator_wrapper.start();
	}

	// *** OpenGL Helpers *** //

	void OGLWidget::initializeGL()
	{
		// We need to generate the surface mesh before we can render it
		mesh.init_surface_mesh();

		// Initialize OpenGL Backend
		initializeOpenGLFunctions();

		// Set global information
		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);
		
		// If we enable transparency, then we should draw something below the surface and sort the surface triangles from back to front before drawing.
		//glDisable(GL_DEPTH_TEST);
		//glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

		// Make the shader program
		make_shader_program();

		// Set some Qt behaviours
		setFocusPolicy(Qt::WheelFocus);
		setUpdateBehavior(QOpenGLWidget::NoPartialUpdate);
	}

	void OGLWidget::resizeGL(int width, int height)
	{
		// Update the perspective matrix to reflect the new window size.
		m_projection.setToIdentity();
		m_projection.perspective(45.0f, width / float(height), 0.1f, 1000.0f);
	}

	void OGLWidget::paintGL()
	{
		// Clear
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Render using our shader
		m_shader_program->bind();
		m_shader_program->setUniformValue(u_projection, m_projection);
		{
			m_vertex_array_object.bind();
			m_shader_program->setUniformValue(u_modelview, m_modelview.modelview);
			glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(mesh.indices.size()), GL_UNSIGNED_INT, mesh.indices.data());
			m_vertex_array_object.release();
		}
		m_shader_program->release();
	}

	void OGLWidget::teardownGL()
	{
		// Actually destroy our OpenGL information
		m_vertex_array_object.destroy();
		m_vertex_buffer.destroy();
		delete m_shader_program;
	}

	void OGLWidget::make_shader_program()
	{
		// Create Shader (Do not release until VAO is created)
		m_shader_program = new QOpenGLShaderProgram();
		m_shader_program->addShaderFromSourceFile(QOpenGLShader::Vertex, "shaders/surface.vert.glsl");
		m_shader_program->addShaderFromSourceFile(QOpenGLShader::Fragment, "shaders/surface.frag.glsl");
		m_shader_program->link();
		m_shader_program->bind();

		// Cache Uniform Locations
		u_projection = m_shader_program->uniformLocation("projection");
		u_modelview = m_shader_program->uniformLocation("modelview");

		// Index buffer
		m_index_buffer.create();
		m_index_buffer.bind();
		m_index_buffer.setUsagePattern(QOpenGLBuffer::StaticDraw);
		m_index_buffer.allocate(mesh.indices.data(), static_cast<int>(mesh.indices.size() * sizeof(unsigned int)));

		// Vertex Buffer
		m_vertex_buffer.create();
		m_vertex_buffer.bind();
		m_vertex_buffer.setUsagePattern(QOpenGLBuffer::DynamicDraw);
		m_vertex_buffer.allocate(mesh.vertices.data(), static_cast<int>(mesh.vertices.size() * sizeof(Vertex)));

		// Create Vertex Array Object
		m_vertex_array_object.create();
		m_vertex_array_object.bind();
		m_shader_program->enableAttributeArray(0);
		m_shader_program->enableAttributeArray(1);
		m_shader_program->enableAttributeArray(2);
		m_shader_program->enableAttributeArray(3);
		m_shader_program->setAttributeBuffer(0, GL_FLOAT, Vertex::positionOffset(), Vertex::PositionTupleSize, Vertex::stride());
		m_shader_program->setAttributeBuffer(1, GL_FLOAT, Vertex::normalOffset(), Vertex::NormalTupleSize, Vertex::stride());
		m_shader_program->setAttributeBuffer(2, GL_FLOAT, Vertex::shininessOffset(), Vertex::ShininessTupleSize, Vertex::stride());
		m_shader_program->setAttributeBuffer(3, GL_FLOAT, Vertex::specularOffset(), Vertex::SpecularTupleSize, Vertex::stride());

		// Release (unbind) all
		m_vertex_array_object.release();
		m_vertex_buffer.release();
		m_index_buffer.release();
		m_shader_program->release();
	}

	// Pass the modified mesh positions to OpenGL. Must be called whenever the mesh values are updated.
	void OGLWidget::update_vertex_buffer()
	{
		// Update the vertex buffer based on the mesh data.
		m_vertex_buffer.bind();
		m_vertex_buffer.write(0, mesh.vertices.data(), static_cast<int>(mesh.vertices.size() * sizeof(Vertex)));
		m_vertex_buffer.release();
	}

	// Rebuild the entire vertex array object after the structure of the mesh has been modified (not just updated values).
	void OGLWidget::rebuild_vertex_array_object()
	{
		m_vertex_array_object.bind();
		m_vertex_buffer.bind();
		m_index_buffer.bind();
		m_vertex_buffer.allocate(mesh.vertices.data(), static_cast<int>(mesh.vertices.size() * sizeof(Vertex)));
		m_index_buffer.allocate(mesh.indices.data(), static_cast<int>(mesh.indices.size() * sizeof(unsigned int)));
		m_vertex_buffer.release();
		m_index_buffer.release();
		m_vertex_array_object.release();
	}

	// *** Mouse and key event helpers and handlers *** //

	void OGLWidget::keyPressEvent(QKeyEvent *event)
	{
		switch (event->key()) {
			case Qt::Key_Up:
				m_modelview.translate(0.0f, 1.0f, 0.0f);
				m_modelview.update();
				QOpenGLWidget::update();
				break;
			case Qt::Key_Down:
				m_modelview.translate(0.0f, -1.0f, 0.0f);
				m_modelview.update();
				QOpenGLWidget::update();
				break;
			case Qt::Key_Left:
				m_modelview.translate(-1.0f, 0.0f, 0.0f);
				m_modelview.update();
				QOpenGLWidget::update();
				break;
			case Qt::Key_Right:
				m_modelview.translate(1.0f, 0.0f, 0.0f);
				m_modelview.update();
				QOpenGLWidget::update();
				break;
			case Qt::Key_Escape:
			case Qt::Key_Q: // Quit
				quit();
				break;
			case Qt::Key_I: // Change initial conditions
			{
				change_initial_conditions();
				break;
			}
			case Qt::Key_R: // Reverse time
				g_waves.step_size_time = -g_waves.step_size_time;
				g_waves.Half_Step(); // Take 2 half-steps to set up the integrator for gonig backwards
				g_waves.Half_Step();
				break;
			case Qt::Key_B: // Change boundary conditions
			{
				change_boundary_conditions();
				break;
			}
			case Qt::Key_Space: // Space, toggle pause
				emit toggle_paused_integrator();
				break;
			case Qt::Key_Z: // Zero the velocity to get a time-mirror effect
				for (auto& cell : g_waves.cells)
					std::fill(cell.V.begin(), cell.V.end(), 0.0);
				g_waves.Half_Step();
				break;
			default:
				event->ignore();
				return;
		}
		event->accept();
	}

	void OGLWidget::mousePressEvent(QMouseEvent *event)
	{
		mouse_button_pressed = event->button();
		mouse_press_location = event->pos();
		switch (event->button())
		{
			case Qt::LeftButton:
				break;
			case Qt::RightButton:
				break;
			default:
				event->ignore();
				return;
		}
		event->accept();
	}

	void OGLWidget::mouseReleaseEvent(QMouseEvent *event)
	{
		mouse_button_pressed = Qt::NoButton;
		event->accept();
	}

	void OGLWidget::mouseMoveEvent(QMouseEvent *event)
	{
		auto location = event->pos();
		auto rotation = m_modelview.rotation.transposed();
		switch (mouse_button_pressed) {
			case Qt::LeftButton:
			{
				rotation.rotate(-0.2f*(location.y() - mouse_press_location.y()), 1.0f, 0.0f, 0.0f);
				float distance_to_origin = std::sqrt(std::pow(location.x() - this->width() / 2, 2) + std::pow(location.y() - this->height() / 2, 2));
				if (distance_to_origin != 0) {
					rotation.rotate(0.2*(location.x() - mouse_press_location.x()) * (this->height() / 2 - location.y()) / distance_to_origin, 0.0f, 0.0f, 1.0f);
				}
				m_modelview.rotation = rotation.transposed();
			}
			break;
			case Qt::RightButton:
				rotation.rotate(-0.2f*(location.x() - mouse_press_location.x()), 0.0f, 1.0f, 0.0f);
				rotation.rotate(-0.2f*(location.y() - mouse_press_location.y()), 1.0f, 0.0f, 0.0f);
				m_modelview.rotation = rotation.transposed();
				break;
			case Qt::MiddleButton:
				m_modelview.translate((location.x() - mouse_press_location.x()) / 10.0f, -(location.y() - mouse_press_location.y()) / 10.0f, 0.0);
				break;
			default:
				event->ignore();
				return;
		}
		m_modelview.update();
		mouse_press_location = location;
		event->accept();
		QOpenGLWidget::update();
	}

	void OGLWidget::wheelEvent(QWheelEvent *event)
	{
		auto angle = event->angleDelta() / 120;
		m_modelview.translate(0.0f, -angle.y(), angle.y() * 0.5f);
		m_modelview.update();
		event->accept();
		QOpenGLWidget::update();
	}

	void OGLWidget::mouseDoubleClickEvent(QMouseEvent *event)
	{
		emit toggle_fullscreen();
		event->accept();
	}

	// *** public slots *** //

	// Refresh the OpenGLWidget
	void OGLWidget::update()
	{
		if (!is_paused())
			update_vertex_buffer(); // A possible race condition here if integrator_wrapper updates the surface mesh concurrently. The result would just be visual tearing though, so not bad.
		QOpenGLWidget::update();
	}

	// Switch the integrator to the next set of initial conditions.
	void OGLWidget::change_initial_conditions()
	{ // Note: we must signal to the integrator_wrapper thread before modifying the integrator
		emit modify_integrator();
		g_waves.Time = 0.0;
		float hint = g_waves.Change_Initial_Conditions();
		// Update the surface mesh and display it.
		mesh.update_surface_mesh();
		update_vertex_buffer();
		QOpenGLWidget::update();
	}

	// Switch the integrator to the next set of boudnary conditions.
	void OGLWidget::change_boundary_conditions()
	{ // Note: we must signal to the integrator_wrapper thread before modifying the integrator
		emit modify_integrator();
		g_waves.Time = 0.0;
		float hint = g_waves.Change_Boundary_Conditions();
		// Rebuild the surface mesh and display it.
		mesh.init_surface_mesh();
		rebuild_vertex_array_object();
		QOpenGLWidget::update();
	}

	// Quit
	void OGLWidget::quit()
	{
		// Tell the integrator_wrapper thread to quit and wait for it before actually quiting.
		integrator_wrapper.quit();
		integrator_wrapper.wait();
		QApplication::quit();
	}

}