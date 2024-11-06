///@file

#include <QKeyEvent>
#include <QThread>
#include <QTimer>
#include <cmath>
#include "oglwidget.h"
#include "integrator_wrapper.h"
#include <boost/compute/interop/opengl.hpp>

namespace Waves
{
	OGLWidget::OGLWidget(QWidget *parent)
			: QOpenGLWidget(parent)
	{
		// Set the default modelview.
		m_modelview.base.setToIdentity();
		m_modelview.base.translate(0.0f, 2.0f, -30.0f);
		m_modelview.base.rotate(-35.0f, 1.0f, 0.0f, 0.0f);
		m_modelview.update();

		// Limit the framerate
		auto *frame_rate_limiter = new QTimer(this);
		QObject::connect(frame_rate_limiter, &QTimer::timeout, [&]()
										 { m_render_frame = true; });
		frame_rate_limiter->start(static_cast<int>(1000.0 / m_max_fps));

		// Start the integrator (starts paused).
		run_integrator();
	}

	OGLWidget::~OGLWidget()
	{
		m_integrator_wrapper.quit(); // Tell the integrator thread to quit first, then wait for it.
		m_integrator_wrapper.wait();
		makeCurrent();
		teardownGL();
	}

	// *** Integrator information *** //

	// Run the simulation in a separate thread
	void OGLWidget::run_integrator()
	{
		// Connections to integrator_wrapper
		QObject::connect(this, &OGLWidget::pause_integrator, &m_integrator_wrapper, &Integrator_Wrapper::pause);
		QObject::connect(this, &OGLWidget::unpause_integrator, &m_integrator_wrapper, &Integrator_Wrapper::unpause);
		QObject::connect(this, &OGLWidget::toggle_paused_integrator, &m_integrator_wrapper, &Integrator_Wrapper::toggle_paused);
		QObject::connect(this, &OGLWidget::modify_integrator, &m_integrator_wrapper, &Integrator_Wrapper::modify_integrator);
		QObject::connect(&m_integrator_wrapper, &Integrator_Wrapper::result_ready, this, &OGLWidget::update);
		QObject::connect(&m_integrator_wrapper, &Integrator_Wrapper::finished, &m_integrator_wrapper, &QObject::deleteLater);
		m_integrator_wrapper.start();

		// Connections to ui
		QObject::connect(this, &OGLWidget::toggle_paused_integrator, this, &OGLWidget::notify_paused_state);
		QObject::connect(this, &OGLWidget::pause_integrator, this, &OGLWidget::notify_paused_state);
		QObject::connect(this, &OGLWidget::unpause_integrator, this, &OGLWidget::notify_paused_state);
	}

	// *** OpenGL Helpers *** //

	void OGLWidget::initializeGL()
	{
		// Initialise the boost.compute context and queue
		m_shared_context = boost::compute::opengl_create_shared_context();
		m_queue = boost::compute::command_queue(m_shared_context, m_shared_context.get_device());

		// Initialise the integrator wrapper
		double dt = 0.001;
		m_integrator_wrapper.initialise(dt, 2, 2, 5, 5, m_shared_context, m_queue);
		emit notify_time_step(dt);
		emit notify_IC_changed(m_integrator_wrapper.get_IC());
		emit notify_BC_changed(m_integrator_wrapper.get_BC());
		emit update_time();

		// Initialize OpenGL Backend
		initializeOpenGLFunctions();

		// Set global information
		glEnable(GL_DEPTH_TEST);
		glDisable(GL_CULL_FACE);

		// If we enable transparency, then we should draw something below the surface and sort the surface triangles from back to front before drawing.
		// glDisable(GL_DEPTH_TEST);
		// glEnable(GL_BLEND);
		// glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

		// Make the shader program
		make_shader_program(); // <- m_vertex_buffer is created here! m_mesh needs to be generated before this point, but cannot be updated (using OpenCL) until after.

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
			glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(m_integrator_wrapper.get_mesh_index_data_size()), GL_UNSIGNED_INT, m_integrator_wrapper.get_mesh_index_data());
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
		m_index_buffer.allocate(m_integrator_wrapper.get_mesh_index_data(), static_cast<int>(m_integrator_wrapper.get_mesh_index_data_size() * sizeof(unsigned int)));

		// Vertex Buffer
		m_vertex_buffer.create();
		m_vertex_buffer.bind();
		m_vertex_buffer.setUsagePattern(QOpenGLBuffer::DynamicDraw);
		m_vertex_buffer.allocate(m_integrator_wrapper.get_mesh_vertex_data(), static_cast<int>(m_integrator_wrapper.get_mesh_vertex_data_size() * sizeof(Vertex)));

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

		// Re-connect the mesh and vertex_buffer
		m_integrator_wrapper.connect_mesh_with_vertex_buffer(m_vertex_buffer);
	}

	// Pass the modified mesh positions to OpenGL. Must be called whenever the mesh values are updated.
	// With OpenCL directly copying values to the vertex buffer, this function is no longer necessary.
	void OGLWidget::update_vertex_buffer()
	{
		// Update the vertex buffer based on the mesh data.
		m_vertex_buffer.bind();
		m_vertex_buffer.write(0, m_integrator_wrapper.get_mesh_vertex_data(), static_cast<int>(m_integrator_wrapper.get_mesh_vertex_data_size() * sizeof(Vertex)));
		m_vertex_buffer.release();
	}

	// Rebuild the entire vertex array object after the structure of the mesh has been modified (not just updated values).
	void OGLWidget::rebuild_vertex_array_object()
	{
		m_vertex_array_object.bind();
		m_vertex_buffer.bind();
		m_index_buffer.bind();
		m_vertex_buffer.allocate(m_integrator_wrapper.get_mesh_vertex_data(), static_cast<int>(m_integrator_wrapper.get_mesh_vertex_data_size() * sizeof(Vertex)));
		m_index_buffer.allocate(m_integrator_wrapper.get_mesh_index_data(), static_cast<int>(m_integrator_wrapper.get_mesh_index_data_size() * sizeof(unsigned int)));
		m_vertex_buffer.release();
		m_index_buffer.release();
		m_vertex_array_object.release();

		// Connect the mesh to the vertex buffer
		m_integrator_wrapper.connect_mesh_with_vertex_buffer(m_vertex_buffer);
	}

	// *** Mouse and key event helpers and handlers *** //

	void OGLWidget::keyPressEvent(QKeyEvent *event)
	{
		switch (event->key())
		{
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
			m_integrator_wrapper.reverse_time();
			break;
		case Qt::Key_B: // Change boundary conditions
		{
			change_boundary_conditions();
			break;
		}
		case Qt::Key_Space: // Space, toggle pause
			emit toggle_paused_integrator();
			break;
			// case Qt::Key_Z: // Zero the velocity to get a time-mirror effect
			//	for (auto & cell : g_waves.cells)
			//		for (auto & v : cell.V)
			//			v = 0.0;
			//	g_waves.Half_Step();
			//	break;
		default:
			event->ignore();
			return;
		}
		event->accept();
	}

	void OGLWidget::mousePressEvent(QMouseEvent *event)
	{
		m_mouse_button_pressed = event->button();
		m_mouse_press_location = event->pos();
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
		m_mouse_button_pressed = Qt::NoButton;
		event->accept();
	}

	void OGLWidget::mouseMoveEvent(QMouseEvent *event)
	{
		auto location = event->pos();
		auto rotation = m_modelview.rotation.transposed();
		switch (m_mouse_button_pressed)
		{
		case Qt::LeftButton:
		{
			rotation.rotate(-0.2f * (location.y() - m_mouse_press_location.y()), 1.0f, 0.0f, 0.0f);
			float distance_to_origin = std::sqrt(std::pow(location.x() - this->width() / 2, 2) + std::pow(location.y() - this->height() / 2, 2));
			if (distance_to_origin != 0)
			{
				rotation.rotate(0.2 * (location.x() - m_mouse_press_location.x()) * (this->height() / 2 - location.y()) / distance_to_origin, 0.0f, 0.0f, 1.0f);
			}
			m_modelview.rotation = rotation.transposed();
		}
		break;
		case Qt::RightButton:
			rotation.rotate(-0.2f * (location.x() - m_mouse_press_location.x()), 0.0f, 1.0f, 0.0f);
			rotation.rotate(-0.2f * (location.y() - m_mouse_press_location.y()), 1.0f, 0.0f, 0.0f);
			m_modelview.rotation = rotation.transposed();
			break;
		case Qt::MiddleButton:
			m_modelview.translate((location.x() - m_mouse_press_location.x()) / 10.0f, -(location.y() - m_mouse_press_location.y()) / 10.0f, 0.0);
			break;
		default:
			event->ignore();
			return;
		}
		m_modelview.update();
		m_mouse_press_location = location;
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
		if (m_render_frame)
		{ // Ignore updates that have happened too soon since the previous one
			m_render_frame = false;
			auto tic = std::chrono::steady_clock::now();
			if (!is_paused())
			{
				m_integrator_wrapper.update_surface_mesh();
				m_fps = 0.95 * m_fps + 0.05 * 1000.0 / std::chrono::duration_cast<std::chrono::milliseconds>(tic - m_tic).count();
			}
			m_tic = tic;
			QOpenGLWidget::update();
			emit update_time(); // Notify the main window that it should update it's displayed time value
		}
	}

	// Switch the integrator to the next set of initial conditions.
	void OGLWidget::change_initial_conditions(int ic)
	{ // Note: we must signal to the integrator_wrapper thread before modifying the integrator
		emit modify_integrator();
		m_integrator_wrapper.change_initial_conditions(ic);
		QOpenGLWidget::update();
		emit notify_paused_state();
		emit notify_IC_changed(m_integrator_wrapper.get_IC());
		emit update_time();
	}

	// Switch the integrator to the next set of boudnary conditions.
	void OGLWidget::change_boundary_conditions(int bc)
	{ // Note: we must signal to the integrator_wrapper thread before modifying the integrator
		emit modify_integrator();
		m_integrator_wrapper.change_boundary_conditions(bc);
		rebuild_vertex_array_object();
		QOpenGLWidget::update();
		emit notify_paused_state();
		emit notify_BC_changed(m_integrator_wrapper.get_BC());
		emit update_time();
	}

	// Quit
	void OGLWidget::quit()
	{
		// Tell the integrator_wrapper thread to quit and wait for it before actually quiting.
		m_integrator_wrapper.quit();
		m_integrator_wrapper.wait();
		QApplication::quit();
	}

	void OGLWidget::change_dissipation(double value)
	{
		auto was_paused = is_paused();
		emit modify_integrator();
		m_integrator_wrapper.set_artificial_dissipation(value);
		if (!was_paused)
			emit unpause_integrator(); // integrator gets paused by modify_integrator
	}

	void OGLWidget::change_timestep(double timestep)
	{
		auto was_paused = is_paused();
		emit modify_integrator();
		m_integrator_wrapper.change_time_step(timestep);
		if (!was_paused)
			emit unpause_integrator(); // integrator gets paused by modify_integrator
	}

	void OGLWidget::change_wave_speed(double speed)
	{
		auto was_paused = is_paused();
		emit modify_integrator();
		m_integrator_wrapper.change_wave_speed(speed);
		if (!was_paused)
			emit unpause_integrator(); // integrator gets paused by modify_integrator
	}

	void OGLWidget::change_height_scale(double scale)
	{
		// auto was_paused = is_paused();
		// emit modify_integrator();
		m_integrator_wrapper.change_height_scale(scale);
		if (m_integrator_wrapper.get_time() == 0)
			change_initial_conditions(m_integrator_wrapper.get_IC());
		// if (!was_paused)
		// 	emit unpause_integrator(); // integrator gets paused by modify_integrator
	}
}
