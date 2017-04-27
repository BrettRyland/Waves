#pragma once
///@file

#include <vector>
#include <array>

#include "vertex.h"

namespace Waves {
	/// Surface mesh class (Qt version) for creating a surface mesh to pass to OpenGL
	class Mesh {
	public:
		std::vector<Vertex> vertices; ///< Vertices in the mesh
		std::vector<unsigned int> indices; ///< The index into vertices of each vertex
		std::vector<std::array<unsigned int, 3>> adjacency_information; ///< Adjacency information, i.e. indices of the cells to the right, above and above right of the current cell.

		/// Generate a surface mesh based on the current Waves::g_waves.
		void init_surface_mesh();
		/// Update the surface_vertex data based on the updated values of the integrator.
		void update_surface_mesh();
	};
}