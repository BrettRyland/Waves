#include <stdlib.h>
#include <GL/glew.h>
#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include "renderer.h"
#include "integrator.h"

namespace Waves {
	namespace Renderer {
		// Our global instance of the OpenGL_resources.
		OpenGL_Resources g_resources;
	}
}

using namespace Waves::Renderer;

// Read in the fragment or vertex shader files. (Unmodified from original C code.)
void* file_contents(const char *filename, GLint *length)
{
	FILE *f = fopen(filename, "r");
	void *buffer;

	if (!f) {
		fprintf(stderr, "Unable to open %s for reading\n", filename);
		return NULL;
	}

	fseek(f, 0, SEEK_END);
	*length = ftell(f);
	fseek(f, 0, SEEK_SET);

	buffer = malloc(*length + 1);
	*length = fread(buffer, 1, *length, f);
	fclose(f);
	((char*)buffer)[*length] = '\0';

	return buffer;
}

// Show why the shader program failed to compile or link. (Unmodified from original C code.)
void show_info_log(GLuint object, PFNGLGETSHADERIVPROC glGet__iv, PFNGLGETSHADERINFOLOGPROC glGet__InfoLog)
{
	GLint log_length;
	char *log;

	glGet__iv(object, GL_INFO_LOG_LENGTH, &log_length);
	log = (char *)malloc(log_length);
	glGet__InfoLog(object, log_length, NULL, log);
	fprintf(stderr, "%s", log);
	free(log);
}

// Create the shader components. (Unmodified from original C code.)
GLuint make_shader(GLenum type, const char *filename)
{
	GLint length;
	GLchar *source = (GLchar*)file_contents(filename, &length);
	GLuint shader;
	GLint shader_ok;

	if (!source)
		return 0;

	shader = glCreateShader(type);
	glShaderSource(shader, 1, (const GLchar**)&source, &length);
	free(source);
	glCompileShader(shader);

	glGetShaderiv(shader, GL_COMPILE_STATUS, &shader_ok);
	if (!shader_ok) {
		fprintf(stderr, "Failed to compile %s:\n", filename);
		show_info_log(shader, glGetShaderiv, glGetShaderInfoLog);
		glDeleteShader(shader);
		return 0;
	}
	return shader;
}

// Link the vertex and fragment shader components together to create the shader program. (Unmodified from original C code.)
GLuint make_program(GLuint vertex_shader, GLuint fragment_shader)
{
	GLint program_ok;

	GLuint program = glCreateProgram();

	glAttachShader(program, vertex_shader);
	glAttachShader(program, fragment_shader);
	glLinkProgram(program);

	glGetProgramiv(program, GL_LINK_STATUS, &program_ok);
	if (!program_ok) {
		fprintf(stderr, "Failed to link shader program:\n");
		show_info_log(program, glGetProgramiv, glGetProgramInfoLog);
		glDeleteProgram(program);
		return 0;
	}
	return program;
}

// Try to get a good estimate of the z-scaling for nice viewing.
glm::mat4 Waves::Renderer::rescale_z(float hint)
{
	double max_Z = 0.0;
	for (auto it = g_resources.surface_vertex_array.begin(); it != g_resources.surface_vertex_array.end(); ++it)
		if (fabs((*it).position[2]) > max_Z)
			max_Z = fabs((*it).position[2]);
	return glm::scale(glm::mat4(1.0f), glm::vec3{ 1.0f, 1.0f, 2.0f * hint / max_Z });
}

// Pass the surface mesh data to the shader program.
void Waves::Renderer::render_mesh(const Waves::Mesh::surface_mesh& mesh)
{
	glBindBuffer(GL_ARRAY_BUFFER, mesh.vertex_buffer);
	glVertexAttribPointer(
		g_resources.surface_program.attributes.position,
		3, GL_FLOAT, GL_FALSE, sizeof(Waves::Mesh::surface_vertex),
		(void*)offsetof(Waves::Mesh::surface_vertex, position)
	);
	glVertexAttribPointer(
		g_resources.surface_program.attributes.normal,
		3, GL_FLOAT, GL_FALSE, sizeof(Waves::Mesh::surface_vertex),
		(void*)offsetof(Waves::Mesh::surface_vertex, normal)
	);
	glVertexAttribPointer(
		g_resources.surface_program.attributes.shininess,
		1, GL_FLOAT, GL_FALSE, sizeof(Waves::Mesh::surface_vertex),
		(void*)offsetof(Waves::Mesh::surface_vertex, shininess)
	);
	glVertexAttribPointer(
		g_resources.surface_program.attributes.specular,
		4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(Waves::Mesh::surface_vertex),
		(void*)offsetof(Waves::Mesh::surface_vertex, specular)
	);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.element_buffer);
	glDrawElements(
		GL_TRIANGLES,
		mesh.element_count,
		GL_UNSIGNED_INT,
		(void*)0
	);
}

// Compile the shaders and make the shader program.
void Waves::Renderer::make_surface_program(GLuint& vertex_shader, GLuint& fragment_shader, GLuint& program)
{
	vertex_shader = make_shader(GL_VERTEX_SHADER, "shaders\\surface.v.glsl");
	if (vertex_shader == 0)
		throw std::runtime_error("Failed to compile vertex shader.");

	fragment_shader = make_shader(GL_FRAGMENT_SHADER, "shaders\\surface.f.glsl");
	if (fragment_shader == 0)
		throw std::runtime_error("Failed to compile fragment shader.");

	program = make_program(vertex_shader, fragment_shader);
	if (program == 0)
		throw std::runtime_error("Failed to create shader program.");
}

// Pass certain g_resource variables to the shader program.
void Waves::Renderer::enact_surface_program(GLuint vertex_shader, GLuint fragment_shader, GLuint program)
{
	g_resources.surface_program.vertex_shader = vertex_shader;
	g_resources.surface_program.fragment_shader = fragment_shader;
	g_resources.surface_program.program = program;

	g_resources.surface_program.uniforms.p_matrix = glGetUniformLocation(program, "p_matrix");
	g_resources.surface_program.uniforms.mv_matrix = glGetUniformLocation(program, "mv_matrix");

	g_resources.surface_program.attributes.position = glGetAttribLocation(program, "position");
	g_resources.surface_program.attributes.normal = glGetAttribLocation(program, "normal");
	g_resources.surface_program.attributes.shininess = glGetAttribLocation(program, "shininess");
	g_resources.surface_program.attributes.specular = glGetAttribLocation(program, "specular");
}

// Delete the shader program.
void Waves::Renderer::delete_surface_program(void)
{
	glDetachShader(
		g_resources.surface_program.program,
		g_resources.surface_program.vertex_shader
	);
	glDetachShader(
		g_resources.surface_program.program,
		g_resources.surface_program.fragment_shader
	);
	glDeleteProgram(g_resources.surface_program.program);
	glDeleteShader(g_resources.surface_program.vertex_shader);
	glDeleteShader(g_resources.surface_program.fragment_shader);
}

// Initialise g_resources and create the shader program.
void Waves::Renderer::make_resources(float hint)
{
	GLuint vertex_shader, fragment_shader, program;

	g_resources.surface_vertex_array = init_surface_mesh(g_resources.surface);

	make_surface_program(vertex_shader, fragment_shader, program);
	enact_surface_program(vertex_shader, fragment_shader, program);

	g_resources.p_matrix = glm::perspective(45.0f, static_cast<GLfloat>(g_resources.window_size[0]) / static_cast<GLfloat>(g_resources.window_size[1]), 0.0625f, 256.0f);
	g_resources.mv_base_matrix = glm::lookAt(g_resources.eye, g_resources.centre, g_resources.up);
	g_resources.scale_matrix = rescale_z(hint);
	g_resources.mv_matrix = g_resources.mv_base_matrix * g_resources.scale_matrix * g_resources.rotation_matrix;
}

// Set some OpenGL condtions.
void Waves::Renderer::init_gl_state(void)
{
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
}

// glutIdleFunc callback.
void Waves::Renderer::update(void)
{
	if (!g_resources.pause) {
		g_waves.Step();
		update_surface_mesh(g_resources.surface, g_resources.surface_vertex_array);
	}
	glutPostRedisplay();
}

// glutMotionFunc callback.
void Waves::Renderer::drag(int x, int y)
{
	if (g_resources.mouse_grab) {
		switch (g_resources.mouse_click[2]) {
			case 0: // left click
				g_resources.rotation_matrix = glm::rotate(g_resources.rotation_matrix, 0.005f*(y - g_resources.mouse_click[1]), glm::vec3{ 1.0f, 0.0f, 0.0f }); // rotate about x-axis
				g_resources.rotation_matrix = glm::rotate(g_resources.rotation_matrix, 0.005f*(x - g_resources.mouse_click[0]), glm::vec3{ 0.0f, 0.0f, 1.0f }); // rotate about z-axis
				break;
			case 2: // right click
				g_resources.rotation_matrix = glm::rotate(g_resources.rotation_matrix, 0.005f*(y - g_resources.mouse_click[1]), glm::vec3{ 1.0f, 0.0f, 0.0f }); // rotate about x-axis
				g_resources.rotation_matrix = glm::rotate(g_resources.rotation_matrix, 0.005f*(x - g_resources.mouse_click[0]), glm::vec3{ 0.0f, 1.0f, 0.0f }); // rotate about y-axis
				break;
		}
		g_resources.mv_matrix = g_resources.mv_base_matrix * g_resources.scale_matrix * g_resources.rotation_matrix;
	}
	g_resources.mouse_click[0] = x;
	g_resources.mouse_click[1] = y;
}

// glutMouseFunc callback.
void Waves::Renderer::mouse(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN) { // a mouse button has been pressed
		g_resources.mouse_click[2] = button; // register which button
		g_resources.mouse_click[0] = x; // and where
		g_resources.mouse_click[1] = y;
		switch (button) { // the button corresponds to the scroll wheel
			case 0: // left
			case 1: // middle
			case 2: // right
				g_resources.mouse_grab = true;
				break;
			case 4: // scroll down => zoom out
				g_resources.eye *= 1.1;
				break;
			case 3: // scroll up => zoom in
				g_resources.eye *= 0.9;
				break;
		}
		g_resources.mv_base_matrix = glm::lookAt(g_resources.eye, g_resources.centre, g_resources.up);
		g_resources.mv_matrix = g_resources.mv_base_matrix * g_resources.scale_matrix * g_resources.rotation_matrix;
	}
	else
		g_resources.mouse_grab = false;
}

// glutKeyboardFunc callback.
void Waves::Renderer::keyboard(unsigned char key, int x, int y)
{
	switch (key) {
		case 27: // Escape key
		case 'q': // Quit
		case 'Q':
			if (glutGameModeGet(GLUT_GAME_MODE_ACTIVE))
				glutLeaveGameMode();
			else
				glutDestroyWindow(g_resources.Window_ID);
			glutLeaveMainLoop();
			delete_surface_program();
			break;
		case 'i': // Change initial conditions
		case 'I':
		{
			float hint = g_waves.Change_Initial_Conditions();
			g_waves.Time = 0.0;
			g_resources.pause = true;
			update_surface_mesh(g_resources.surface, g_resources.surface_vertex_array);
			g_resources.scale_matrix = rescale_z(hint);
			g_resources.mv_matrix = g_resources.mv_base_matrix * g_resources.scale_matrix * g_resources.rotation_matrix;
			break;
		}
		case 'b': // Change boundary conditions
		case 'B':
			g_waves.Change_Boundary_Conditions();
			g_waves.Time = 0.0;
			g_resources.pause = true;
			update_surface_mesh(g_resources.surface, g_resources.surface_vertex_array);
			break;
		case ' ': // Space, toggle pause
			g_resources.pause = !g_resources.pause;
			break;
		default:
			std::cout << "Unknown key\n";
	}
}

// glutReshapeFunc callback.
void Waves::Renderer::reshape(int w, int h)
{
	g_resources.window_size[0] = w;
	g_resources.window_size[1] = h;
	g_resources.p_matrix = glm::perspective(90.0f, static_cast<GLfloat>(w) / static_cast<GLfloat>(h), 0.0625f, 256.0f);
	glViewport(0, 0, w, h);
}

// glutDisplayFunc callback.
void Waves::Renderer::render(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(g_resources.surface_program.program);

	//glActiveTexture(GL_TEXTURE0);
	//glUniform1i(g_resources.surface_program.uniforms.texture, 0);

	glUniformMatrix4fv(
		g_resources.surface_program.uniforms.p_matrix,
		1, GL_FALSE,
		glm::value_ptr(g_resources.p_matrix)
	);

	glUniformMatrix4fv(
		g_resources.surface_program.uniforms.mv_matrix,
		1, GL_FALSE,
		glm::value_ptr(g_resources.mv_matrix)
	);

	glEnableVertexAttribArray(g_resources.surface_program.attributes.position);
	glEnableVertexAttribArray(g_resources.surface_program.attributes.normal);
	//glEnableVertexAttribArray(g_resources.surface_program.attributes.texcoord);
	glEnableVertexAttribArray(g_resources.surface_program.attributes.shininess);
	glEnableVertexAttribArray(g_resources.surface_program.attributes.specular);

	render_mesh(g_resources.surface);

	glDisableVertexAttribArray(g_resources.surface_program.attributes.position);
	glDisableVertexAttribArray(g_resources.surface_program.attributes.normal);
	//glDisableVertexAttribArray(g_resources.surface_program.attributes.texcoord);
	glDisableVertexAttribArray(g_resources.surface_program.attributes.shininess);
	glDisableVertexAttribArray(g_resources.surface_program.attributes.specular);
	glutSwapBuffers();
}