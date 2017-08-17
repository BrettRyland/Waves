// The precision used by the integrator. This must be set in all .cl files to mirror that in integratorCL.h.
typedef float integrator_precision;

// This should be the C99 equivalent of the relevant components of the Vertex class in vertex.h
// Note: we can't use float3's or float4's here as the compiler promotes *everything* to float4's, which doesn't line up with the contents of the vertex_buffer
typedef struct tag_Vertex {
	float position[3];
	float normal[3];
	float shininess;
	float specular[4];
} Vertex;

// Copy the surface heights from the main nodes of the integrator to the vertex buffer
__kernel void copy_heights(__global Vertex *vertex_buffer, __global const integrator_precision *U, unsigned int nodes_per_cell) {
	size_t cell = get_global_id(0);

	// Copy surface heights
	vertex_buffer[cell].position[2] = (float)U[cell*nodes_per_cell];
}

// Calculate the normals of the vertices in the vertex buffer
__kernel void calculate_normals(__global Vertex *vertex_buffer, __global const unsigned int *adjacency_information) {
	size_t vertex = get_global_id(0);

	// Calculate normals
	float3 v1 = (float3)(vertex_buffer[adjacency_information[2 * vertex + 0]].position[0] - vertex_buffer[vertex].position[0], vertex_buffer[adjacency_information[2 * vertex + 0]].position[1] - vertex_buffer[vertex].position[1], vertex_buffer[adjacency_information[2 * vertex + 0]].position[2] - vertex_buffer[vertex].position[2]);
	float3 v2 = (float3)(vertex_buffer[adjacency_information[2 * vertex + 1]].position[0] - vertex_buffer[vertex].position[0], vertex_buffer[adjacency_information[2 * vertex + 1]].position[1] - vertex_buffer[vertex].position[1], vertex_buffer[adjacency_information[2 * vertex + 1]].position[2] - vertex_buffer[vertex].position[2]);
	float3 n = fast_normalize(cross(v1, v2));
	vertex_buffer[vertex].normal[0] = n.x;
	vertex_buffer[vertex].normal[1] = n.y;
	vertex_buffer[vertex].normal[2] = n.z;
}

// Copy the duplicate heights
__kernel void copy_duplicate_heights(__global Vertex *vertex_buffer, __global const unsigned int *duplicate_information, unsigned int offset) {
	size_t vertex = get_global_id(0);

	// Copy heights
	vertex_buffer[vertex + offset].position[2] = vertex_buffer[duplicate_information[vertex]].position[2];

}

// Copy the duplicate normals
__kernel void copy_duplicate_normals(__global Vertex *vertex_buffer, __global const unsigned int *duplicate_information, unsigned int offset) {
	size_t vertex = get_global_id(0);

	// Copy normals
	vertex_buffer[vertex + offset].normal[0] = vertex_buffer[duplicate_information[vertex]].normal[0];
	vertex_buffer[vertex + offset].normal[1] = vertex_buffer[duplicate_information[vertex]].normal[1];
	vertex_buffer[vertex + offset].normal[2] = vertex_buffer[duplicate_information[vertex]].normal[2];
}
