#pragma once

#include <vector>
#include <array>

#include "vertex.h"

namespace Waves {

	class Mesh {
	public:
		std::vector<Vertex> vertices;
		std::vector<unsigned int> indices;
		std::vector<std::array<unsigned int, 3>> adjacency_information;

		// Generate a surface mesh.
		void init_surface_mesh();
		// Update the surface_vertex data based on the updated values of the integrator.
		void update_surface_mesh();
	};
}