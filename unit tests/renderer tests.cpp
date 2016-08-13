#include "catch.hpp"
#include "integrator.h"
#include "renderer.h"

using namespace Waves::Renderer;

// Some basic sample unit testing using Catch https://github.com/philsquared/Catch

TEST_CASE("GL setup", "Renderer") {
	// Some initial setup to get the test case up and running
	bool Fullscreen = false;
	float hint = Waves::g_waves.Initialise(2, 2, 1e-2, 1, 1);
	int argc = 0;
	char** argv = NULL;
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(g_resources.window_size[0], g_resources.window_size[1]);
	if (Fullscreen == true) {
		char buf[20];
		g_resources.window_size[0] = glutGet(GLUT_SCREEN_WIDTH);
		g_resources.window_size[1] = glutGet(GLUT_SCREEN_HEIGHT);
		sprintf_s(buf, "%dx%d:32@60", g_resources.window_size[0], g_resources.window_size[1]);
		glutGameModeString(buf);
		glutEnterGameMode();
	}
	else
		g_resources.Window_ID = glutCreateWindow("Waves : u_tt=c^2*(u_xx+u_yy)-V'(u)");
	glutIdleFunc(&update);
	glutDisplayFunc(&render);
	glutReshapeFunc(&reshape);
	glutMotionFunc(&drag);
	glutMouseFunc(&mouse);
	glutKeyboardFunc(&keyboard);

	glewInit();
	if (!GLEW_VERSION_2_0) {
		throw std::runtime_error("OpenGL 2.0 not available.");
	}

	init_gl_state();
	make_resources(hint);

	// shader program resources
	REQUIRE(g_resources.surface_program.uniforms.mv_matrix != -1);
	REQUIRE(g_resources.surface_program.uniforms.p_matrix != -1);
	REQUIRE(g_resources.surface_program.attributes.normal != -1);
	REQUIRE(g_resources.surface_program.attributes.position != -1);
	REQUIRE(g_resources.surface_program.attributes.shininess != -1);
	REQUIRE(g_resources.surface_program.attributes.specular != -1);

	// projection matrix
	for (auto i = 0; i < 4; ++i)
		for (auto j = 0; j < 4; ++j)
			REQUIRE(std::isinf(g_resources.p_matrix[i][j]) == false);
	// modelview matrix
	for (auto i = 0; i < 4; ++i)
		for (auto j = 0; j < 4; ++j)
			REQUIRE(std::isinf(g_resources.mv_matrix[i][j]) == false);
	// modelview base matrix
	for (auto i = 0; i < 4; ++i)
		for (auto j = 0; j < 4; ++j)
			REQUIRE(std::isinf(g_resources.mv_base_matrix[i][j]) == false);
	// rotation matrix
	for (auto i = 0; i < 4; ++i)
		for (auto j = 0; j < 4; ++j)
			REQUIRE(std::isinf(g_resources.rotation_matrix[i][j]) == false);
	// rescale matrix
	for (auto i = 0; i < 4; ++i)
		for (auto j = 0; j < 4; ++j)
			REQUIRE(std::isinf(g_resources.scale_matrix[i][j]) == false);
}
