//FIXME: all the pointers and horrid mallocs here should be redone in terms of std::unique_ptr

#include <gl/glew.h>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <memory>
#include <limits>
#include <algorithm>
#include <vector>
#include <array>
#include "mesh.h"
#include "integrator.h"

using glm::vec3;
using glm::vec4;
using glm::mat4;
const unsigned int missing_index = std::numeric_limits<unsigned int>::max();

namespace Waves::Mesh {
	// Pass the mesh generated by init_surface_mesh to the shader program.
	void init_mesh(surface_mesh& out_mesh, const std::vector<surface_vertex>& vertex_data, GLsizei vertex_count, const std::vector<GLuint>& element_data, GLsizei element_count, GLenum hint)
	{
		glGenBuffers(1, &(out_mesh.vertex_buffer));
		glGenBuffers(1, &(out_mesh.element_buffer));
		out_mesh.vertex_count = vertex_count;
		out_mesh.element_count = element_count;

		glBindBuffer(GL_ARRAY_BUFFER, out_mesh.vertex_buffer);
		glBufferData(
			GL_ARRAY_BUFFER,
			vertex_count * sizeof(Mesh::surface_vertex),
			vertex_data.data(),
			hint
		);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, out_mesh.element_buffer);
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			element_count * sizeof(GLuint),
			element_data.data(),
			GL_STATIC_DRAW // The element_data only gets updated when we change the boundary conditions, so it can be static.
		);
	}

	// Generate a surface mesh (later assigned to g_resoures.surface_vertex_array).
	std::vector<surface_vertex> init_surface_mesh(surface_mesh& out_mesh)
	{
		// Calculate vertex and element data and then call init_mesh with these.
		// The element data can then be thrown away, but the vertex data is needed later so we can update it in update_surface_mesh.
		unsigned int vertex_count = g_waves.cells.size();
		std::vector<surface_vertex> vertex_data(vertex_count); // under estimate, but it will grow as needed.
		std::vector<GLuint> element_data;
		element_data.reserve(6 * g_waves.cells.size());
		out_mesh.adjacency_information.resize(vertex_count);
		int s, t, i;
		GLuint index;
		GLfloat scale{ 5.0f };

		// Assign the components of the position, normal, shininess and specular values, adding in extra vertices as necessary.
		// The non-static components (position[2] and normal[0-2]) are reassigned in calculate_surface_vertex().
		// FIXME: We should not draw cells for Dirichlet_right and Neumann_right.
		for (unsigned int c = 0; c < g_waves.cells.size(); ++c) {
			vertex_data[c] = { static_cast<GLfloat>(g_waves.position_information[c][0]), static_cast<GLfloat>(g_waves.position_information[c][1]) };
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
									out_mesh.adjacency_information[c] = { g_waves.adjacency_information[c][1], g_waves.adjacency_information[c][3], g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3] };
									element_data.insert(element_data.end(), { c, g_waves.adjacency_information[c][1], g_waves.adjacency_information[c][3], g_waves.adjacency_information[c][1], g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3], g_waves.adjacency_information[c][3] });
									break;
								case Cell_Type::Periodic:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									vertex_data.emplace_back(surface_vertex{ static_cast<GLfloat>(g_waves.position_information[c][0] + g_waves.step_size_x),static_cast<GLfloat>(g_waves.position_information[c][1] + g_waves.step_size_y) });
									out_mesh.adjacency_information[c] = { g_waves.adjacency_information[c][1], g_waves.adjacency_information[c][3], vertex_count };
									element_data.insert(element_data.end(), { c, g_waves.adjacency_information[c][1], g_waves.adjacency_information[c][3], g_waves.adjacency_information[c][1], vertex_count, g_waves.adjacency_information[c][3] });
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
							vertex_data.emplace_back(surface_vertex{ static_cast<GLfloat>(g_waves.position_information[c][0]),static_cast<GLfloat>(g_waves.position_information[c][1] + g_waves.step_size_y) });
							switch (g_waves.cells[g_waves.adjacency_information[c][1]].cell_type.second) {
								case Cell_Type::Normal:
								case Cell_Type::Dirichlet_left:
								case Cell_Type::Neumann_left:
									// Top right vertex exists
									out_mesh.adjacency_information[c] = { g_waves.adjacency_information[c][1], vertex_count, g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3] };
									element_data.insert(element_data.end(), { c, g_waves.adjacency_information[c][1], vertex_count, g_waves.adjacency_information[c][1], g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3], vertex_count });
									++vertex_count;
									break;
								case Cell_Type::Periodic:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									vertex_data.emplace_back(surface_vertex{ static_cast<GLfloat>(g_waves.position_information[c][0] + g_waves.step_size_x),static_cast<GLfloat>(g_waves.position_information[c][1] + g_waves.step_size_y) });
									out_mesh.adjacency_information[c] = { g_waves.adjacency_information[c][1], vertex_count, vertex_count + 1 };
									element_data.insert(element_data.end(), { c, g_waves.adjacency_information[c][1], vertex_count, g_waves.adjacency_information[c][1], vertex_count + 1, vertex_count });
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
					vertex_data.emplace_back(surface_vertex{ static_cast<GLfloat>(g_waves.position_information[c][0] + g_waves.step_size_x), static_cast<GLfloat>(g_waves.position_information[c][1]) });
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
									out_mesh.adjacency_information[c] = { vertex_count, g_waves.adjacency_information[c][3], g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1] };
									element_data.insert(element_data.end(), { c, vertex_count, g_waves.adjacency_information[c][3], vertex_count, g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1], g_waves.adjacency_information[c][3] });
									++vertex_count;
									break;
								case Cell_Type::Periodic:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									vertex_data.emplace_back(surface_vertex{ static_cast<GLfloat>(g_waves.position_information[c][0] + g_waves.step_size_x),static_cast<GLfloat>(g_waves.position_information[c][1] + g_waves.step_size_y) });
									out_mesh.adjacency_information[c] = { vertex_count, g_waves.adjacency_information[c][3], vertex_count + 1 };
									element_data.insert(element_data.end(), { c, vertex_count, g_waves.adjacency_information[c][3], vertex_count, vertex_count + 1, g_waves.adjacency_information[c][3] });
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
							vertex_data.emplace_back(surface_vertex{ static_cast<GLfloat>(g_waves.position_information[c][0]), static_cast<GLfloat>(g_waves.position_information[c][1] + g_waves.step_size_y) });
							// And hence, so is the top right vertex
							vertex_data.emplace_back(surface_vertex{ static_cast<GLfloat>(g_waves.position_information[c][0] + g_waves.step_size_x), static_cast<GLfloat>(g_waves.position_information[c][1] + g_waves.step_size_y) });
							out_mesh.adjacency_information[c] = { vertex_count, vertex_count + 1, vertex_count + 2 };
							element_data.insert(element_data.end(), { c, vertex_count, vertex_count + 1, vertex_count, vertex_count + 2, vertex_count + 1 });
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
		calculate_surface_vertex(vertex_data, out_mesh.adjacency_information);

		// Pass the vertex and element data to init_mesh with the DYNAMIC flag as this data will change frequently.
		init_mesh(out_mesh, vertex_data, vertex_count, element_data, element_data.size(), GL_DYNAMIC_DRAW);

		return vertex_data;
	}

	// Update the surface_vertex data based on the updated values of the integrator g_waves.
	void calculate_surface_vertex(std::vector<surface_vertex>& vertex_data, const std::vector<std::array<unsigned int, 3>>& adjacency_information)
	{
		// Copy heights
#pragma omp parallel for
		for (int c = 0; c < g_waves.cells.size(); ++c)
			vertex_data[c].position[2] = g_waves.cells[c].U[0];

		// Copy periodic cells
#pragma omp parallel for
		for (int c = 0; c < g_waves.cells.size(); ++c) {
			if (g_waves.cells[c].cell_type.first == Cell_Type::Periodic) {
				vertex_data[adjacency_information[c][0]].position[2] = vertex_data[g_waves.adjacency_information[c][1]].position[2];
				assert(g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3] != missing_index); // assert that the top right vertex exists
				vertex_data[adjacency_information[c][2]].position[2] = vertex_data[g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3]].position[2];
			}
			if (g_waves.cells[c].cell_type.second == Cell_Type::Periodic) {
				vertex_data[adjacency_information[c][1]].position[2] = vertex_data[g_waves.adjacency_information[c][3]].position[2];
				assert(g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1] != missing_index); // assert that the top right vertex exists
				vertex_data[adjacency_information[c][2]].position[2] = vertex_data[g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1]].position[2];
			}
		}

		// Calculate normals
#pragma omp parallel for
		for (int c = 0; c < g_waves.cells.size(); ++c) {
			vec3 u = { vertex_data[adjacency_information[c][0]].position[0] - vertex_data[c].position[0],vertex_data[adjacency_information[c][0]].position[1] - vertex_data[c].position[1], vertex_data[adjacency_information[c][0]].position[2] - vertex_data[c].position[2] };
			vec3 v = { vertex_data[adjacency_information[c][1]].position[0] - vertex_data[c].position[0],vertex_data[adjacency_information[c][1]].position[1] - vertex_data[c].position[1], vertex_data[adjacency_information[c][1]].position[2] - vertex_data[c].position[2] };
			std::memcpy(&vertex_data[c].normal[0], glm::value_ptr(glm::normalize(glm::cross(u, v))), 3 * sizeof(GLfloat));
		}

		// Copy periodic cells
#pragma omp parallel for
		for (int c = 0; c < g_waves.cells.size(); ++c) {
			if (g_waves.cells[c].cell_type.first == Cell_Type::Periodic) {
				std::memcpy(&vertex_data[adjacency_information[c][0]].normal[0], &vertex_data[g_waves.adjacency_information[c][1]].normal[0], 3 * sizeof(GLfloat));
				std::memcpy(&vertex_data[adjacency_information[c][2]].normal[0], &vertex_data[g_waves.adjacency_information[g_waves.adjacency_information[c][1]][3]].normal[0], 3 * sizeof(GLfloat));
			}
			if (g_waves.cells[c].cell_type.second == Cell_Type::Periodic) {
				std::memcpy(&vertex_data[adjacency_information[c][1]].normal[0], &vertex_data[g_waves.adjacency_information[c][3]].normal[0], 3 * sizeof(GLfloat));
				std::memcpy(&vertex_data[adjacency_information[c][2]].normal[0], &vertex_data[g_waves.adjacency_information[g_waves.adjacency_information[c][3]][1]].normal[0], 3 * sizeof(GLfloat));
			}
		}
	}

	// Pass the updated surface_vertex data to the shader program.
	void update_surface_mesh(const surface_mesh& mesh, std::vector<surface_vertex>& vertex_data)
	{
		// Update the surface_vertex buffer.
		calculate_surface_vertex(vertex_data, mesh.adjacency_information);
		glBindBuffer(GL_ARRAY_BUFFER, mesh.vertex_buffer);
		glBufferData(GL_ARRAY_BUFFER, mesh.vertex_count * sizeof(Mesh::surface_vertex), vertex_data.data(), GL_DYNAMIC_DRAW);
	}
}