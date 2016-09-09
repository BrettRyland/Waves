#version 330

uniform mat4 projection, modelview;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in float shininess;
layout(location = 3) in vec4 specular;

out vec3 frag_position, frag_normal;
out float frag_shininess;
out vec4 frag_specular;
 
void main()
{
    vec4 eye_position = modelview * vec4(position, 1.0);
    gl_Position = projection * eye_position;
    frag_position = eye_position.xyz;
    frag_normal   = (modelview * vec4(normal, 0.0)).xyz;
    frag_shininess = shininess;
    frag_specular = specular;
}