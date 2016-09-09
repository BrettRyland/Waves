#include <cassert>
#include <limits>

#include <QVector3D>

#include "mesh.h"
#include "integrator.h"

namespace Waves {
	const unsigned int missing_index = std::numeric_limits<unsigned int>::max();

	// Generate a surface mesh.
	void Mesh::init_surface_mesh()
	{
		unsigned int vertex_count = g_waves.cells.size();
		vertices.resize(vertex_count); // An under estimate, but it will grow as needed (most likely only once).
		indices.clear();
		indices.reserve(6 * g_waves.cells.size());
		adjacency_information.resize(vertex_count);
		float scale{ 5.0f };

		// Assign the components of the position, normal, shininess and specular values, adding in extra vertices as necessary.
		// The non-static components (position[2] and normal[0-2]) are reassigned in update_surface_mesh().
		// FIXME: We should not draw cells for Dirichlet_right and Neumann_right.
		for (unsigned int c = 0; c < g_waves.cells.size(); ++c) {
			vertices[c] = { static_cast<float>(g_waves.position_information[c][0]), static_cast<float>(g_waves.position_information[c][1]) };
			switch (g_waves.cells[c].cell_type.first) {
				case Cell_Type::Normal:
				case Cell_Type::Dirichlet_left:
				case Cell_Type::Neumann_left:
					// Right vertex exists
					switch (g_waves.cells[c].cell_type.second) {
						case Cell_Type::Normal:
						case Cell_Type::Dirichlet_left:
						case Cell_Type::Neumann_left:
							// The top vertex exists
							switch (g_waves.cells[g_waves.adjacency_information[c][1]].cell_type.second) {
								case Cell_Type::Normal:
								case Cell_Type::Dirichlet_left:
								case Cell_Type::Neumann_left:
									// All vertices exist
									adjacency_information[c] = { g_waves.adjacency_information[c][1], g_waves.adjacency_information[c][3], g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3] };
									indices.insert(indices.end(), { c, g_waves.adjacency_information[c][1], g_waves.adjacency_information[c][3], g_waves.adjacency_information[c][1], g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3], g_waves.adjacency_information[c][3] });
									break;
								case Cell_Type::Periodic:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									vertices.emplace_back(Vertex{ static_cast<float>(g_waves.position_information[c][0] + g_waves.step_size_x),static_cast<float>(g_waves.position_information[c][1] + g_waves.step_size_y) });
									adjacency_information[c] = { g_waves.adjacency_information[c][1], g_waves.adjacency_information[c][3], vertex_count };
									indices.insert(indices.end(), { c, g_waves.adjacency_information[c][1], g_waves.adjacency_information[c][3], g_waves.adjacency_information[c][1], vertex_count, g_waves.adjacency_information[c][3] });
									++vertex_count;
									break;
								default:
									throw std::runtime_error("Invalid cell type.");
							}
							break;
						case Cell_Type::Periodic:
						case Cell_Type::Dirichlet_right:
						case Cell_Type::Neumann_right:
							// The top vertex is missing
							vertices.emplace_back(Vertex{ static_cast<float>(g_waves.position_information[c][0]),static_cast<float>(g_waves.position_information[c][1] + g_waves.step_size_y) });
							switch (g_waves.cells[g_waves.adjacency_information[c][1]].cell_type.second) {
								case Cell_Type::Normal:
								case Cell_Type::Dirichlet_left:
								case Cell_Type::Neumann_left:
									// Top right vertex exists
									adjacency_information[c] = { g_waves.adjacency_information[c][1], vertex_count, g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3] };
									indices.insert(indices.end(), { c, g_waves.adjacency_information[c][1], vertex_count, g_waves.adjacency_information[c][1], g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3], vertex_count });
									++vertex_count;
									break;
								case Cell_Type::Periodic:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									vertices.emplace_back(Vertex{ static_cast<float>(g_waves.position_information[c][0] + g_waves.step_size_x),static_cast<float>(g_waves.position_information[c][1] + g_waves.step_size_y) });
									adjacency_information[c] = { g_waves.adjacency_information[c][1], vertex_count, vertex_count + 1 };
									indices.insert(indices.end(), { c, g_waves.adjacency_information[c][1], vertex_count, g_waves.adjacency_information[c][1], vertex_count + 1, vertex_count });
									vertex_count += 2;
									break;
								default:
									throw std::runtime_error("Invalid cell type.");
							}
							break;
						default:
							throw std::runtime_error("Invalid cell type.");
					}
					break;
				case Cell_Type::Periodic:
				case Cell_Type::Dirichlet_right:
				case Cell_Type::Neumann_right:
					// Right vertex is missing
					vertices.emplace_back(Vertex{ static_cast<float>(g_waves.position_information[c][0] + g_waves.step_size_x), static_cast<float>(g_waves.position_information[c][1]) });
					switch (g_waves.cells[c].cell_type.second) {
						case Cell_Type::Normal:
						case Cell_Type::Dirichlet_left:
						case Cell_Type::Neumann_left:
							// Top vertex exists
							switch (g_waves.cells[g_waves.adjacency_information[c][3]].cell_type.first) {
								case Cell_Type::Normal:
								case Cell_Type::Dirichlet_left:
								case Cell_Type::Neumann_left:
									// Top right vertex exists
									adjacency_information[c] = { vertex_count, g_waves.adjacency_information[c][3], g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1] };
									indices.insert(indices.end(), { c, vertex_count, g_waves.adjacency_information[c][3], vertex_count, g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1], g_waves.adjacency_information[c][3] });
									++vertex_count;
									break;
								case Cell_Type::Periodic:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									vertices.emplace_back(Vertex{ static_cast<float>(g_waves.position_information[c][0] + g_waves.step_size_x),static_cast<float>(g_waves.position_information[c][1] + g_waves.step_size_y) });
									adjacency_information[c] = { vertex_count, g_waves.adjacency_information[c][3], vertex_count + 1 };
									indices.insert(indices.end(), { c, vertex_count, g_waves.adjacency_information[c][3], vertex_count, vertex_count + 1, g_waves.adjacency_information[c][3] });
									vertex_count += 2;
									break;
								default:
									throw std::runtime_error("Invalid cell type.");
							}
							break;
						case Cell_Type::Periodic:
						case Cell_Type::Dirichlet_right:
						case Cell_Type::Neumann_right:
							// Top vertex is missing
							vertices.emplace_back(Vertex{ static_cast<float>(g_waves.position_information[c][0]), static_cast<float>(g_waves.position_information[c][1] + g_waves.step_size_y) });
							// And hence, so is the top right vertex
							vertices.emplace_back(Vertex{ static_cast<float>(g_waves.position_information[c][0] + g_waves.step_size_x), static_cast<float>(g_waves.position_information[c][1] + g_waves.step_size_y) });
							adjacency_information[c] = { vertex_count, vertex_count + 1, vertex_count + 2 };
							indices.insert(indices.end(), { c, vertex_count, vertex_count + 1, vertex_count, vertex_count + 2, vertex_count + 1 });
							vertex_count += 3;
							break;
						default:
							throw std::runtime_error("Invalid cell type.");
					}
					break;
				default:
					throw std::runtime_error("Invalid cell type.");
			}
		}

		// Calculate surface heights and normals.
		update_surface_mesh();
	}

	// Update the Vertex data based on the updated values of the integrator g_waves.
	void Mesh::update_surface_mesh()
	{
		// Copy heights
#pragma omp parallel for
		for (int c = 0; c < g_waves.cells.size(); ++c)
			vertices[c].position[2] = g_waves.cells[c].U[0];

		// Copy periodic cells
#pragma omp parallel for
		for (int c = 0; c < g_waves.cells.size(); ++c) {
			if (g_waves.cells[c].cell_type.first == Cell_Type::Periodic) {
				vertices[adjacency_information[c][0]].position[2] = vertices[g_waves.adjacency_information[c][1]].position[2];
				assert(g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3] != missing_index); // assert that the top right vertex exists
				vertices[adjacency_information[c][2]].position[2] = vertices[g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3]].position[2];
			}
			if (g_waves.cells[c].cell_type.second == Cell_Type::Periodic) {
				vertices[adjacency_information[c][1]].position[2] = vertices[g_waves.adjacency_information[c][3]].position[2];
				assert(g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1] != missing_index); // assert that the top right vertex exists
				vertices[adjacency_information[c][2]].position[2] = vertices[g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1]].position[2];
			}
		}

		// Calculate normals
#pragma omp parallel for
		for (int c = 0; c < g_waves.cells.size(); ++c) {
			QVector3D u = { vertices[adjacency_information[c][0]].position[0] - vertices[c].position[0],vertices[adjacency_information[c][0]].position[1] - vertices[c].position[1], vertices[adjacency_information[c][0]].position[2] - vertices[c].position[2] };
			QVector3D v = { vertices[adjacency_information[c][1]].position[0] - vertices[c].position[0],vertices[adjacency_information[c][1]].position[1] - vertices[c].position[1], vertices[adjacency_information[c][1]].position[2] - vertices[c].position[2] };
			vertices[c].normal = QVector3D::normal(u, v);
		}

		// Copy periodic cells
#pragma omp parallel for
		for (int c = 0; c < g_waves.cells.size(); ++c) {
			if (g_waves.cells[c].cell_type.first == Cell_Type::Periodic) {
				vertices[adjacency_information[c][0]].normal = vertices[g_waves.adjacency_information[c][1]].normal;
				vertices[adjacency_information[c][2]].normal = vertices[g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3]].normal;
			}
			if (g_waves.cells[c].cell_type.second == Cell_Type::Periodic) {
				vertices[adjacency_information[c][1]].normal = vertices[g_waves.adjacency_information[c][3]].normal;
				vertices[adjacency_information[c][2]].normal = vertices[g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1]].normal;
			}
		}
	}
}
