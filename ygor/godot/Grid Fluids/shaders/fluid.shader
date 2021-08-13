shader_type canvas_item;

uniform sampler2D noise_img;
uniform float speed = 1.0;
uniform vec4 smoke_color : hint_color;

float interpolate(){
	
}

vec2 advect(vec2 v){
	// TODO
}

vec2 diffuse(vec2 v){
	// TODO
}

vec2 addForces(vec2 v){
	// pegar input do mouse forca repelindo
	// TODO
}

vec2 computePressure(vec2 v){
	// TODO
}

vec2 subtractPressureGradiant(vec2 v, float p){
	// TODO
}

void vertex(){
	vec2 v = (EXTRA_MATRIX * (WORLD_MATRIX * vec4(VERTEX, 0.0, 1.0))).xy;
	float p;
	
	// calculos
	v = addForces(v);
	v = advect(v);
	v = diffuse(v);
	p = computePressure(v);
	v = subtractPressureGradiant(v, p);
}

void fragment() {
	float u = texture(noise_img, UV).r; 
	float v = texture(noise_img, UV).g; 
	float p = texture(noise_img, UV).b; 
	
	float x_speed = u * TIME * speed;
	float y_speed = v * TIME * speed;
	
	
	// exibir resultado
	
}
