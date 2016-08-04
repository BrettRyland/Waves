#version 120

uniform mat4 p_matrix, mv_matrix;

attribute vec3 position, normal;
attribute float shininess;
attribute vec4 specular;

out vec3 frag_position, frag_normal;
out float frag_shininess;
out vec4 frag_specular;

void main()
{
    vec4 eye_position = mv_matrix * vec4(position, 1.0);
    gl_Position = p_matrix * eye_position;
    frag_position = eye_position.xyz;
    frag_normal   = (mv_matrix * vec4(normal, 0.0)).xyz;
    frag_shininess = shininess;
    frag_specular = specular;
}
