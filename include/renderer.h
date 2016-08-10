#pragma once

#include <gl/glew.h>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <vector>
#include "mesh.h"

// Surface rendering engine, based on http://duriansoftware.com/joe/An-intro-to-modern-OpenGL.-Chapter-4:-Rendering-a-Dynamic-3D-Scene-with-Phong-Shading.html
// Modified to use glm and modern C++ instead of raw C pointers.
// Also, the texturing component is removed as we don't need it.

namespace Waves::Renderer {

	// Global resources for OpenGL.
	struct OpenGL_Resources {
		Mesh::surface_mesh surface;
		std::vector<Mesh::surface_vertex> surface_vertex_array;

		// Shader resources.
		struct {
			GLuint vertex_shader, fragment_shader, program;

			struct {
				GLint p_matrix, mv_matrix;
			} uniforms;

			struct {
				GLint position, normal, shininess, specular;
			} attributes;
		} surface_program;

		// Projection and modelview matrices.
		glm::mat4 p_matrix, mv_matrix, mv_base_matrix, rotation_matrix{ 1.0f }, scale_matrix{ 1.0f };
		glm::vec3 eye{ 0.0f, -4.0f, 4.0f }, centre{ 0.0f, 0.0f, 0.0f }, up = { 0.0f, 0.0f, 1.0f };

		int Window_ID;
		GLsizei window_size[2]{ 1200, 800 };
		bool mouse_grab{ false };
		bool pause{ true };
		int mouse_click[3]{ 0, 0, 0 };
	};

	// We need to use a single global instance of this class for OpenGL to have access to it. We declare it extern here and construct it in renderer.cpp.
	extern OpenGL_Resources g_resources;

	// Initialise g_resources.
	void make_resources(float hint);

	// Set some OpenGL condtions.
	void init_gl_state();

	// Shader helper functions.
	void make_surface_program(GLuint& vertex_shader, GLuint& fragment_shader, GLuint& program);
	void enact_surface_program(GLuint vertex_shader, GLuint fragment_shader, GLuint program);
	void delete_surface_program();

	// Functions for glut callbacks.
	void update();
	void drag(int x, int y);
	void mouse(int button, int state, int x, int y);
	void keyboard(unsigned char key, int x, int y);
	void reshape(int w, int h);
	void render();
	// render callback helper function.
	void render_mesh(const Mesh::surface_mesh& mesh);

	// Helper function for scaling the axes to make for better viewing.
	glm::mat4 rescale(float hint = 1.0f);
}