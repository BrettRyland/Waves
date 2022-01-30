///@file

#include <cassert>
#include <limits>
#include <QVector3D>
#include <QOpenGLFunctions>
#include "mesh.h"
#include <boost/compute/interop/opengl.hpp>

namespace Waves {
	// Generate a surface mesh.
	// Note: This only assigns "static" data that doesn't change each step. The surface heights and normals are assigned in update_surface_mesh, but that has to wait until the vertex buffer object is allocated on the rendering device.
	void Mesh::generate_surface_mesh(Integrator & integrator)
	{
		cl_uint vertex_count{ integrator.num_cells };
		m_vertices.resize(vertex_count); // An under estimate, but it will grow as needed (most likely only once).
		extra_vertices_duplicate_infomation.clear();
		m_indices.clear();
		m_indices.reserve(6 * integrator.num_cells);
		mesh_adjacency_information.resize(vertex_count);
		auto integrator_position_information = integrator.host_position_information;// integrator.get_position_information();
		auto integrator_adjacency_information = integrator.host_adjacency_information;// integrator.get_adjacency_information();
		auto cell_type_x = integrator.host_cell_type_x;// integrator.get_cell_type_x();
		auto cell_type_y = integrator.host_cell_type_y;// integrator.get_cell_type_y();
		assert(integrator_position_information.size() == 2 * vertex_count);
		assert(integrator_adjacency_information.size() == 4 * vertex_count);
		assert(cell_type_x.size() == vertex_count);
		assert(cell_type_y.size() == vertex_count);

		// Assign the components of the position, normal, shininess and specular values, adding in extra vertices as necessary.
		// The non-static components (position[2] and normal[0-2]) are reassigned in update_surface_mesh().
		// Note: of the extra vertices that are created, some of them may be duplicates, but the percentage of these should be low as they only occur along the boundaries.
		// Note: no cell can be both a left and right boundary at the same time.
		for (cl_uint c = 0; c < integrator.num_cells; ++c) {
			m_vertices[c] = { static_cast<float>(integrator_position_information[2 * c + 0]), static_cast<float>(integrator_position_information[2 * c + 1]) };
			switch (cell_type_x[c]) {
				case Cell_Type::Normal:
				case Cell_Type::Periodic_left:
				case Cell_Type::Dirichlet_left:
				case Cell_Type::Neumann_left:
					// Right vertex exists
					switch (cell_type_y[c]) {
						case Cell_Type::Normal:
						case Cell_Type::Periodic_left:
						case Cell_Type::Dirichlet_left:
						case Cell_Type::Neumann_left:
							// The top vertex exists
							mesh_adjacency_information[c] = { integrator_adjacency_information[4 * c + 1], integrator_adjacency_information[4 * c + 3] };
							switch (cell_type_y[integrator_adjacency_information[4 * c + 1]]) {
								case Cell_Type::Normal:
								case Cell_Type::Periodic_left:
								case Cell_Type::Dirichlet_left:
								case Cell_Type::Neumann_left:
									// All vertices exist
									m_indices.insert(m_indices.end(), { c, integrator_adjacency_information[4 * c + 1], integrator_adjacency_information[4 * c + 3], integrator_adjacency_information[4 * c + 1], integrator_adjacency_information[4 * integrator_adjacency_information[4 * c + 1] + 3], integrator_adjacency_information[4 * c + 3] });
									break;
								case Cell_Type::Periodic_right:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									m_vertices.emplace_back(Vertex{ static_cast<float>(integrator_position_information[2 * c + 0] + integrator.step_size_x),static_cast<float>(integrator_position_information[2 * c + 1] + integrator.step_size_y) });
									switch (cell_type_y[integrator_adjacency_information[4 * c + 1]]) {
										case Cell_Type::Periodic_right:
											extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * integrator_adjacency_information[4 * c + 1] + 3]); // Copy infomation from periodic node of right cell.
											break;
										default:
											extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * c + 1]); // Copy information from main node of right cell.
											break;
									}
									m_indices.insert(m_indices.end(), { c, integrator_adjacency_information[4 * c + 1], integrator_adjacency_information[4 * c + 3], integrator_adjacency_information[4 * c + 1], vertex_count, integrator_adjacency_information[4 * c + 3] });
									++vertex_count;
									break;
								default:
									throw std::runtime_error("Invalid cell type.");
							}
							break;
						case Cell_Type::Periodic_right:
							// The top vertex is missing
							m_vertices.emplace_back(Vertex{ static_cast<float>(integrator_position_information[2 * c + 0]),static_cast<float>(integrator_position_information[2 * c + 1] + integrator.step_size_y) });
							extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * c + 3]); // Copy infomation from periodic node of this cell.
							mesh_adjacency_information[c] = { integrator_adjacency_information[4 * c + 1], vertex_count };
							switch (cell_type_y[integrator_adjacency_information[4 * c + 1]]) {
								case Cell_Type::Normal:
								case Cell_Type::Periodic_left:
								case Cell_Type::Dirichlet_left:
								case Cell_Type::Neumann_left:
									// Top right vertex exists
									m_indices.insert(m_indices.end(), { c, integrator_adjacency_information[4 * c + 1], vertex_count, integrator_adjacency_information[4 * c + 1], integrator_adjacency_information[4 * integrator_adjacency_information[4 * c + 1] + 3], vertex_count });
									++vertex_count;
									break;
								case Cell_Type::Periodic_right:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									m_vertices.emplace_back(Vertex{ static_cast<float>(integrator_position_information[2 * c + 0] + integrator.step_size_x),static_cast<float>(integrator_position_information[2 * c + 1] + integrator.step_size_y) });
									switch (cell_type_y[integrator_adjacency_information[4 * c + 1]]) {
										case Cell_Type::Periodic_right:
											extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * integrator_adjacency_information[4 * c + 1] + 3]); // Copy infomation from periodic node of right cell.
											break;
										default:
											extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * c + 1]); // Copy information from main node of right cell.
											break;
									}
									m_indices.insert(m_indices.end(), { c, integrator_adjacency_information[4 * c + 1], vertex_count, integrator_adjacency_information[4 * c + 1], vertex_count + 1, vertex_count });
									vertex_count += 2;
									break;
								default:
									throw std::runtime_error("Invalid cell type.");
							}
							break;
						case Cell_Type::Dirichlet_right:
						case Cell_Type::Neumann_right:
							// This is a top boundary, so don't draw beyond it. We still need mesh_adjacency_information to be able to calculate the vertex normal though.
							mesh_adjacency_information[c] = { integrator_adjacency_information[4 * c + 2], integrator_adjacency_information[4 * c + 1] }; // use cross(-b,a) instead of cross(a,b)
							break;
						default:
							throw std::runtime_error("Invalid cell type.");
					}
					break;
				case Cell_Type::Periodic_right:
					// Right vertex is missing
					m_vertices.emplace_back(Vertex{ static_cast<float>(integrator_position_information[2 * c + 0] + integrator.step_size_x), static_cast<float>(integrator_position_information[2 * c + 1]) });
					extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * c + 1]); // Copy infomation from periodic node of this cell.
					switch (cell_type_y[c]) {
						case Cell_Type::Normal:
						case Cell_Type::Periodic_left:
						case Cell_Type::Dirichlet_left:
						case Cell_Type::Neumann_left:
							// Top vertex exists
							mesh_adjacency_information[c] = { vertex_count, integrator_adjacency_information[4 * c + 3] };
							switch (cell_type_x[integrator_adjacency_information[4 * c + 3]]) {
								case Cell_Type::Normal:
								case Cell_Type::Periodic_left:
								case Cell_Type::Dirichlet_left:
								case Cell_Type::Neumann_left:
									// Top right vertex exists
									m_indices.insert(m_indices.end(), { c, vertex_count, integrator_adjacency_information[4 * c + 3], vertex_count, integrator_adjacency_information[4 * integrator_adjacency_information[4 * c + 3] + 1], integrator_adjacency_information[4 * c + 3] });
									++vertex_count;
									break;
								case Cell_Type::Periodic_right:
								case Cell_Type::Dirichlet_right:
								case Cell_Type::Neumann_right:
									// Top right vertex is missing
									m_vertices.emplace_back(Vertex{ static_cast<float>(integrator_position_information[2 * c + 0] + integrator.step_size_x),static_cast<float>(integrator_position_information[2 * c + 1] + integrator.step_size_y) });
									switch (cell_type_x[integrator_adjacency_information[4 * c + 3]]) {
										case Cell_Type::Periodic_right:
											extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * integrator_adjacency_information[4 * c + 3] + 1]); // Copy infomation from periodic node of top cell.
											break;
										default:
											extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * c + 3]); // Copy information from main node of top cell.
											break;
									}
									m_indices.insert(m_indices.end(), { c, vertex_count, integrator_adjacency_information[4 * c + 3], vertex_count, vertex_count + 1, integrator_adjacency_information[4 * c + 3] });
									vertex_count += 2;
									break;
								default:
									throw std::runtime_error("Invalid cell type.");
							}
							break;
						case Cell_Type::Periodic_right:
						case Cell_Type::Dirichlet_right:
						case Cell_Type::Neumann_right:
							// Top vertex is missing
							m_vertices.emplace_back(Vertex{ static_cast<float>(integrator_position_information[2 * c + 0]), static_cast<float>(integrator_position_information[2 * c + 1] + integrator.step_size_y) });
							switch (cell_type_y[c]) {
								case Cell_Type::Periodic_right:
									extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * c + 3]); // Copy infomation from main node of the periodic cell.
									break;
								default:
									extra_vertices_duplicate_infomation.emplace_back(c); // Copy information from main node of this cell.
									break;
							}
							mesh_adjacency_information[c] = { vertex_count, vertex_count + 1 };
							// And hence, so is the top right vertex
							m_vertices.emplace_back(Vertex{ static_cast<float>(integrator_position_information[2 * c + 0] + integrator.step_size_x), static_cast<float>(integrator_position_information[2 * c + 1] + integrator.step_size_y) });
							switch (((cell_type_x[c] == Cell_Type::Periodic_right) ? (1 << 0) : 0) + ((cell_type_y[c] == Cell_Type::Periodic_right) ? (1 << 1) : 0)) {
								case 0: // Copy information from main node of this cell.
									extra_vertices_duplicate_infomation.emplace_back(c);
									break;
								case (1 << 0): // Copy infomation from main node of the x-periodic cell.
									extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * c + 1]);
									break;
								case (1 << 1): // Copy infomation from main node of the y-periodic cell.
									extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * c + 3]);
									break;
								case (1 << 0) + (1 << 1) : // Copy infomation from main node of the periodic [x2] cell.
									extra_vertices_duplicate_infomation.emplace_back(integrator_adjacency_information[4 * integrator_adjacency_information[4 * c + 1] + 3]);
									break;
							}
							m_indices.insert(m_indices.end(), { c, vertex_count, vertex_count + 1, vertex_count, vertex_count + 2, vertex_count + 1 });
							vertex_count += 3;
							break;
						default:
							throw std::runtime_error("Invalid cell type.");
					}
					break;
				case Cell_Type::Dirichlet_right:
				case Cell_Type::Neumann_right:
					switch (cell_type_y[c]) { // This is a right boundary, so don't draw beyond it. We still need mesh_adjacency_information to be able to calculate the vertex normal though.
						case Cell_Type::Normal:
						case Cell_Type::Periodic_left:
						case Cell_Type::Periodic_right:
						case Cell_Type::Dirichlet_left:
						case Cell_Type::Neumann_left:
							mesh_adjacency_information[c] = { integrator_adjacency_information[4 * c + 3], integrator_adjacency_information[4 * c + 0] };
							break;
						default:
							mesh_adjacency_information[c] = { integrator_adjacency_information[4 * c + 0], integrator_adjacency_information[4 * c + 2] };
							break;
					}
					break;
				default:
					throw std::runtime_error("Invalid cell type.");
			}
		}

		// Copy duplicate and adjacency information to the compute device
		m_duplicate_information = boost::compute::vector<cl_uint>(extra_vertices_duplicate_infomation.size(), integrator.m_context);
		boost::compute::copy(extra_vertices_duplicate_infomation.begin(), extra_vertices_duplicate_infomation.end(), m_duplicate_information.begin(), integrator.m_queue);
		[&]() {
			std::vector<cl_uint> tmp; // boost.compute doesn't like copying from std::array
			for (auto a : mesh_adjacency_information)
				tmp.insert(tmp.end(), a.begin(), a.end());
			m_adjacency_information = boost::compute::vector<cl_uint>(tmp.size(), integrator.m_context);
			boost::compute::copy(tmp.begin(), tmp.end(), m_adjacency_information.begin(), integrator.m_queue);
		}();
	}

	// Update the Vertex data based on the updated values of the integrator integrator.
	void Mesh::update_surface_mesh(Integrator & integrator)
	{
		// Copy from device back to host
		std::vector<integrator_precision> U(integrator.num_cells);
		boost::compute::copy(integrator.U.begin(), integrator.U.end(), U.begin(), integrator.m_queue);

		// Copy heights
		for (cl_uint c = 0; c < integrator.num_cells; ++c)
			m_vertices[c].position[2] = static_cast<float>(U[c*integrator.nodes_per_cell]);
		for (cl_uint c = 0; c < extra_vertices_duplicate_infomation.size(); ++c)
			m_vertices[integrator.num_cells + c].position[2] = m_vertices[extra_vertices_duplicate_infomation[c]].position[2];

		// Calculate normals
		for (cl_uint c = 0; c < integrator.num_cells; ++c) {
			auto u = m_vertices[mesh_adjacency_information[c][0]].position - m_vertices[c].position;
			auto v = m_vertices[mesh_adjacency_information[c][1]].position - m_vertices[c].position;
			m_vertices[c].normal = QVector3D::normal(u, v);
		}
		for (cl_uint c = 0; c < extra_vertices_duplicate_infomation.size(); ++c)
			m_vertices[integrator.num_cells + c].normal = m_vertices[extra_vertices_duplicate_infomation[c]].normal;
	}

	// Generate OpenCL kernel for copying the vertex data from the main nodes of the integrator into the position[2] component of the vertices on the device, then calculate normals and account for periodicity.
	void Mesh::generate_mesh_update_kernel(Integrator & integrator) {
		auto program = boost::compute::program::create_with_source_file("kernels/mesh.cl", integrator.m_context);
		program.build();
		m_kernel_copy_heights = program.create_kernel("copy_heights");
		m_kernel_calculate_normals = program.create_kernel("calculate_normals");
		m_kernel_copy_duplicate_heights = program.create_kernel("copy_duplicate_heights");
		m_kernel_copy_duplicate_normals = program.create_kernel("copy_duplicate_normals");
	}

	// Set the kernel arguments. This should be called whenever the vertex_buffer and integrator.U are structurally changed (resized, moved, etc.)
	void Mesh::configure_mesh_update_kernel(Integrator & integrator, QOpenGLBuffer & vertex_buffer) {
		m_vertex_buffer = boost::compute::opengl_buffer(integrator.m_context, vertex_buffer.bufferId(), CL_MEM_WRITE_ONLY);
		m_kernel_copy_heights.set_args(m_vertex_buffer, integrator.U.get_buffer(), integrator.nodes_per_cell);
		m_kernel_calculate_normals.set_args(m_vertex_buffer, m_adjacency_information.get_buffer());
		m_kernel_copy_duplicate_heights.set_args(m_vertex_buffer, m_duplicate_information.get_buffer(), integrator.num_cells);
		m_kernel_copy_duplicate_normals.set_args(m_vertex_buffer, m_duplicate_information.get_buffer(), integrator.num_cells);
	}

	// Update the surface mesh directly on the device
	void Mesh::update_surface_mesh_CL(Integrator &integrator) {
		// acquire the buffer so that it is accessible to OpenCL
		glFinish();
		boost::compute::opengl_enqueue_acquire_buffer(m_vertex_buffer, integrator.m_queue);

		// Update the vertex buffer
		integrator.m_queue.enqueue_1d_range_kernel(m_kernel_copy_heights, 0, integrator.num_cells, 0);
		if (extra_vertices_duplicate_infomation.size() != 0) // Cannot enqueue a kernel with 0 global work size
			integrator.m_queue.enqueue_1d_range_kernel(m_kernel_copy_duplicate_heights, 0, extra_vertices_duplicate_infomation.size(), 0);
		integrator.m_queue.enqueue_1d_range_kernel(m_kernel_calculate_normals, 0, integrator.num_cells, 0);
		if (extra_vertices_duplicate_infomation.size() != 0) // Cannot enqueue a kernel with 0 global work size
			integrator.m_queue.enqueue_1d_range_kernel(m_kernel_copy_duplicate_normals, 0, extra_vertices_duplicate_infomation.size(), 0);

		// release the buffer so that it is accessible to OpenGL
		clFinish(integrator.m_queue);
		boost::compute::opengl_enqueue_release_buffer(m_vertex_buffer, integrator.m_queue);
	}
}
