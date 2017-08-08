#pragma once
///@file

#include <vector>
#include <array>
#include <QOpenGLBuffer>

#include "vertex.h"
#include "integratorCL.h"
#include <boost/compute/interop/opengl.hpp>

namespace Waves {
	/** Surface mesh class (Qt version) for creating a surface mesh to pass to OpenGL.
	The Mesh class is friends with the Integrator class in order to access its raw surface height data (U) and to use its opencl context and queue.
	*/
	class Mesh {
	public:
		// Note: we can't use a const ref for integrator as the act of copying data from the integrator requires a non-const opencl queue
		/// Generate a surface mesh based on the current Waves::g_waves. The adjacency and duplicate information is copied to the compute device.
		void generate_surface_mesh(Integrator & integrator);
		/// Update the surface_vertex data based on the updated values of the integrator.
		void update_surface_mesh(Integrator & integrator);
		/// Update the surface mesh directly on the device
		void update_surface_mesh_CL(Integrator & integrator);

		/// Generate the OpenCL kernel responsible for copying the data from the integrator to the OpenGL vertex buffer
		void generate_mesh_update_kernel(Integrator & integrator);
		/// Configure the kernel (i.e. acquire the OpenCL buffer that accesses the OpenGL vertex buffer and set the kernel arguments)
		void configure_mesh_update_kernel(Integrator & integrator, QOpenGLBuffer & vertex_buffer);

		/** Accessors
		@{*/
		auto const get_vertex_data() const { return m_vertices.data(); } ///< Get the mesh vertex data
		auto const get_vertex_data_size() const { return m_vertices.size(); } ///< Get the size of the mesh vertex data
		auto const get_index_data() const { return m_indices.data(); } ///< Get the mesh index data
		auto const get_index_data_size() const { return m_indices.size(); } ///< Get the size of the mesh index data
		///@}

	private:
		std::vector<Vertex> m_vertices; ///< Vertices in the mesh
		std::vector<cl_uint> m_indices; ///< The index into vertices of each vertex
		std::vector<std::array<cl_uint, 2>> mesh_adjacency_information; ///< Adjacency information, i.e. indices of the cells to the right and above of the current cell.
		std::vector<cl_uint> extra_vertices_duplicate_infomation; ///< Indices of the real vertices that the extra indices should copy

		boost::compute::opengl_buffer m_vertex_buffer; ///< The OpenCL buffer for accessing the OpenGL buffer
		boost::compute::kernel m_kernel_copy_heights, ///< The OpenCL kernel for copying the surface heights from the integrator to the vertex buffer
			m_kernel_calculate_normals, ///< The OpenCL kernel for calculating the vertex normals
			m_kernel_copy_duplicate_heights, ///< The OpenCL kernel for copying the position of duplicate elements of the mesh required for rendering the boundary elements
			m_kernel_copy_duplicate_normals; ///< The OpenCL kernel for copying the normals of duplicate elements of the mesh required for rendering the boundary elements
		boost::compute::vector<cl_uint> m_duplicate_information, ///< Device side copy of duplicate information
			m_adjacency_information; ///< Device side copy of adjacency information
	};
}