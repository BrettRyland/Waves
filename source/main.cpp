//#define PROGRAM_TITLE "u_tt=c^2*(u_xx+u_yy)-V'(u) -- Brett Ryland"

/*
 * Integrate the non-linear wave equation with 2 spatial and 1 temporal dimensions.
 * This is done using rx-stage Lobatto IIIA-IIIB in space and 2-stage Lobatto IIIA-IIIB in time.
 */

#include <iostream>
#include <string>
#include "renderer.h"
#include "integrator.h"

using namespace Waves::Renderer;

int main(int argc, char *argv[])
{
	// Initialise our simulation with some nice defaults.
	bool Fullscreen = false; // Default to not being fullscreen.
	int m = 200; // Note that any values of m and n (>1) are possible (within the processing power available).
	int n = 200; // The classic box scheme needs m and n to be odd... woot! yay for my thesis results.
	double dx = 1.0; // Some nice values for dx, dy and dt to satisfy CFL condition are: dx=1, dy=1, dt=0.01667, ws=10.
	double dy = 1.0;
	double dt = 0.01667; // 1/60 (aiming for 60fps once timing is implemented)
	double ws = 10.0;
	int rx = 2; // Number of stages in x.
	int ry = 2; // Number of stages in y.
	int IC = 2; // Continuous spectrum of frequencies in a central localised hump.
	int BC = 1; // Square Dirichlet (i.e. fixed) boundary.

	std::vector<std::string> args(argv + 1, argv + argc); // Convert input arguments to strings (ignoring program name in argv[0]).
	try { // Parse the input strings.
		for (auto arg : args) {
			auto pos = arg.find_first_of("=");
			auto var = arg.substr(0, pos);
			if (pos != std::string::npos) {
				if (var.compare("m") == 0)
					m = std::stoul(arg.substr(pos + 1));
				else if (var.compare("n") == 0)
					n = std::stoul(arg.substr(pos + 1));
				else if (var.compare("rx") == 0)
					rx = std::stoul(arg.substr(pos + 1));
				else if (var.compare("ry") == 0)
					ry = std::stoul(arg.substr(pos + 1));
				else if (var.compare("dx") == 0)
					dx = std::stod(arg.substr(pos + 1));
				else if (var.compare("dy") == 0)
					dy = std::stod(arg.substr(pos + 1));
				else if (var.compare("dt") == 0)
					dt = std::stod(arg.substr(pos + 1));
				else if (var.compare("ws") == 0)
					ws = std::stod(arg.substr(pos + 1));
				else if (var.compare("IC") == 0)
					IC = std::stoul(arg.substr(pos + 1));
				else if (var.compare("BC") == 0)
					BC = std::stoul(arg.substr(pos + 1));
				else
					throw std::runtime_error("Unknown argument " + arg);
			}
			else {
				if (var.compare("-fs") == 0) {
					Fullscreen = true;
				}
				else if (var.compare("-help") == 0) { // Give some useful help.
					std::cout << "Usage example: Wave m=100 n=100 rx=2 ry=2 dx=1e-2 dy=1e-2 dt=1e-3 ws=1 IC=1 BC=1 -fs\n";
					std::cout << "m = number of gridpoints in the x direction (default: 100)\n";
					std::cout << "n = number of gridpoints in the y direction (default: 100)\n";
					std::cout << "rx = number of stages in x (default: 2)\n";
					std::cout << "ry = number of stages in y (default: 2)\n";
					std::cout << "dx = x stepsize (default: 0.01)\n";
					std::cout << "dy = y stepsize (default: 0.01)\n";
					std::cout << "dt = time stepsize (default: 0.001)\n";
					std::cout << "ws = the speed of the wave (default: 1). Choose dx, dy, dt and ws wisely to satisfy CFL condition\n";
					std::cout << "IC = initial conditions [0 = single frequency, 1 = continuous spectrum, 2 = continuous spectrum with one hump, 3 = continuous spectrum with two humps] (default: 1)\n";
					std::cout << "BC = boundary conditions [0 = Periodic square, 1 = Dirichlet on a square, 2 = Dirichlet on a circle, 3 = Dirichlet on a circle with a cusp, 4 = Dirichlet on a double cusp (i.e., intersecting circles)] (default: 1)\n";
					std::cout << "-fs = fullscreen flag (not included by default)";
					std::cout << std::endl;
				}
				else
					throw std::runtime_error("Unknown flag " + var);
			}
		}
	}
	catch (...)
	{
		throw std::runtime_error("Unable to parse input arguments.");
	}

	float hint = Waves::g_waves.Initialise(m, n, rx, ry, dx, dy, dt, ws, IC, BC);
	std::cout << Waves::g_waves << '\n';
	std::cout << "Use the mouse to look around.\n";
	std::cout << "Press 'i' to change initial conditions.\n";
	std::cout << "Press 'b' to change boundary conditions.\n";
	std::cout << "Press <space> to toggle pause state.\n";
	std::cout << "Press <Esc> or 'q' to quit.\n";

	// Standard glut initialisation and binding of callback functions.
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

	glutMainLoop();

	return 0;
}
