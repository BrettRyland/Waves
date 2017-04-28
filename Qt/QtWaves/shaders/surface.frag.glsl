#version 330

uniform mat4 modelview;

in vec3 frag_position, frag_normal;
in float frag_shininess;
in vec4 frag_specular;

const vec3 light_direction = normalize(vec3(2, 1, -4));
const vec4 light_diffuse = vec4(0.8, 0.8, 0.8, 0.0);
const vec4 light_ambient = vec4(0.1, 0.1, 0.1, 1.0);
const vec4 light_specular = vec4(1.0, 1.0, 1.0, 1.0);

void main()
{
    vec3 mv_light_direction = (modelview * vec4(light_direction, 0.0)).xyz,
         normal = normalize(frag_normal),
         eye = normalize(frag_position),
         reflection = reflect(mv_light_direction, normal);

    vec4 diffuse_factor = max(-dot(normal, mv_light_direction), 0.0) * light_diffuse;
	vec4 frag_diffuse = frag_specular; // use the specular component as the diffuse component too
    vec4 ambient_diffuse_factor = diffuse_factor + light_ambient;
    vec4 specular_factor = max(pow(-dot(reflection, eye), frag_shininess), 0.0) * light_specular;
    
    gl_FragColor = specular_factor * frag_specular + ambient_diffuse_factor * frag_diffuse;
}
