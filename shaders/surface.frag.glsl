# version 120 

in vec4 color;
in vec3 mynormal; 
in vec4 myvertex; 

const int numLights = 10; 
uniform bool enablelighting; // are we lighting at all (global).
uniform vec4 lightposn[numLights]; // positions of lights 
uniform vec4 lightcolor[numLights]; // colors of lights
uniform int numused;               // number of lights used
uniform vec4 ambient; 
uniform vec4 diffuse; 
uniform vec4 specular; 
uniform vec4 emission; 
uniform float shininess; 

vec4 ComputeLight (const in vec3 direction, const in vec4 lightcolor, const in vec3 normal, const in vec3 halfvec, const in vec4 mydiffuse, const in vec4 myspecular, const in float myshininess)
{
	float nDotL = dot(normal, direction)  ;
	vec4 lambert = mydiffuse * lightcolor * max (nDotL, 0.0) ;
	float nDotH = dot(normal, halfvec) ;
	vec4 phong = myspecular * lightcolor * pow (max(nDotH, 0.0), myshininess) ;
	vec4 retval = lambert + phong ;
	return retval ;
}

void main (void) 
{       
	if (enablelighting) {       
		vec4 finalcolor = ambient + emission; 
		const vec3 eyepos = vec3(0,0,0);
		vec4 _mypos = gl_ModelViewMatrix * myvertex;
		vec3 mypos = _mypos.xyz / _mypos.w; // Dehomogenize current location
		vec3 eyedirn = normalize(eyepos - mypos);
		vec3 normal = normalize(gl_NormalMatrix * mynormal);
		for (int i = 0; i < numused; ++i) {
			vec3 lightpos, lightdirn;
			if (lightposn[i].w == 0) { // directional light
				lightdirn = normalize(lightposn[i].xyz);
			}
			else { // positional light
				lightpos = lightposn[i].xyz / lightposn[i].w; // Dehomogenize light position
				lightdirn = normalize (lightpos - mypos);
			}
			vec3 halfvec = normalize (lightdirn + eyedirn);
			finalcolor += ComputeLight(lightdirn, lightcolor[i], normal, halfvec, diffuse, specular, shininess);
		}
        gl_FragColor = finalcolor; 
    }
	else {
		gl_FragColor = color; 
	}
}
