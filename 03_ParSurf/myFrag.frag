#version 130

in vec2 vs_out_pos;  //Munkamhoz legnagyobb segitseget innen vettem: http://jamie-wong.com/2016/07/15/ray-marching-signed-distance-functions/
out vec4 fs_out_col;

uniform vec2 mywin;
uniform vec3 cameye;
uniform vec3 camat;
uniform vec3 camup;
uniform float mytime;

const int MAX_MARCHING_STEPS = 140;
const float MIN_DIST = 0.0;
const float MAX_DIST = 40.0;
const float EPSILON = 0.005;


vec3 rotate( vec3 pos, float x, float y, float z )
{
	mat3 rotX = mat3( 1.0, 0.0, 0.0, 0.0, cos( x ), -sin( x ), 0.0, sin( x ), cos( x ) );
	mat3 rotY = mat3( cos( y ), 0.0, sin( y ), 0.0, 1.0, 0.0, -sin(y), 0.0, cos(y) );
	mat3 rotZ = mat3( cos( z ), -sin( z ), 0.0, sin( z ), cos( z ), 0.0, 0.0, 0.0, 1.0 );

	return rotX * rotY * rotZ * pos;
}

float sceneSDF( vec3 p ) //http://celarek.at/wp/wp-content/uploads/2014/05/realTimeFractalsReport.pdf
{						 //a fraktal tavolsagfuggvenyet ennek segitsegevel allitottam elo (algoritmus a konyvben)
	vec3 z = p.xyz;
	float dr = 1.0;
	float r = 0.0;
	for (int i=0; i<MAX_MARCHING_STEPS; ++i){
		r=length(z);
		if(r>MAX_DIST){
			break;
		}

		float theta = acos(z.z/r);
		float phi = atan(z.y, z.x);
		dr = pow(r, 8.0-1.0)* 8.0 * dr +1.0;

		float zr = pow(r, 8.0);
		theta= theta*8,0;
		phi= phi *8.0;
		z= zr *vec3(sin(theta)*cos(phi), sin(phi)*sin(theta),cos(theta));
		z= z+p;
		}
	
	return 0.5*log(r)*r/dr;
}

float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    float depth = start;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        float dist = sceneSDF(eye + depth * marchingDirection);
        if (dist < EPSILON) {
			return depth;
        }
        depth += dist;
        if (depth >= end) {
            return end;
        }
    }
    return end;
}

vec3 vecmul(vec3 v1, vec3 v2) { //vektorialis szorzas
  vec3 ret;
 
  ret.x = v1.y * v2.z - v1.z * v2.y;
  ret.y = v1.z * v2.x - v1.x * v2.z;
  ret.z = v1.x * v2.y - v1.y * v2.x;

  return ret;
}

vec3 rayDirection(float Degree, vec3 eye, vec3 at, vec3 up) { //iranyvektor
    vec3 forward= (eye-at);
	forward= normalize(forward);
	vec3 right= vecmul(forward,up);
	right= normalize(right);
	vec3 myup= vecmul(forward,right);
	float alpha=mywin.x*tan(radians(Degree));
	float beta=mywin.y*tan(radians(Degree));
	return (vs_out_pos.x * right * alpha + vs_out_pos.y *myup * beta - forward);
}

vec3 estimateNormal(vec3 p) { //feluleti normalis becslese
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity) { //feny: kodokat hasznaltam shadertoy 
    vec3 N = estimateNormal(p);									//oldalrol, elsosorban: Jamie Wong
    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));
    
    float dotLN = dot(L, N);
    float dotRV = dot(R, V);
    
    if (dotLN < 0.0) {
        // errol a pontrol nem latszik a feny
        return vec3(0.0, 0.0, 0.0);
    } 
    
    if (dotRV < 0.0) {
        // csak diffuz komponens, visszaverodes nem a szem iranyaba
        return lightIntensity * (k_d * dotLN);
    }
    return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    
    vec3 light1Pos = vec3(4.0, 2.0, 4.0); //*sin(mytime/1000.0)

    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light1Pos,
                                  light1Intensity);
    
	vec3 light2Pos = vec3(0.0, 0.0, 5.0);
    
    vec3 light2Intensity = vec3(0.8, 0.8, 0.8);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light2Pos,
                                  light2Intensity);    
    
	vec3 light3Pos = vec3(2.0, 6.0, 2.0);
    
    vec3 light3Intensity = vec3(0.7, 0.7, 0.7);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light3Pos,
                                  light3Intensity); 
    return color;
}

void main()
{
	vec3 dir = rayDirection(45.0, cameye, camat, camup);
	dir = normalize(dir);
    float dist = shortestDistanceToSurface(cameye, dir, MIN_DIST, MAX_DIST);
    
    if (dist > MAX_DIST - EPSILON) {
        // Didn't hit anything
        fs_out_col = vec4(0.0, 0.0, 0.0, 0.0);
		return;
    }

	vec3 p = cameye + dist * dir;

    vec3 K_a = vec3(0.2, 0.2, 0.2);
    vec3 K_d = vec3(0.0, 0.0, 1.0);
    vec3 K_s = vec3(0.0, 1.0, 0.0);
    float shininess = 27.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, cameye);
    
    fs_out_col = vec4(color, 1.0);
	
    //fs_out_col = vec4(1.0, 0.0, 0.0, 1.0);

}
/*---------------------------------------------------------------------------
ITT PRIMITIV MOZGAS NELKUL
#version 130

in vec2 vs_out_pos;
out vec4 fs_out_col;

#define M_PI 3.1415926535897932384626433832795

uniform vec2 mywin;
uniform vec3 cameye;
uniform vec3 camat;
uniform vec3 camup;
uniform float mytime;

const int MAX_MARCHING_STEPS = 300;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;

float intersectSDF(float distA, float distB) {
    return max(distA, distB);
}

float unionSDF(float distA, float distB) {
    return min(distA, distB);
}

float differenceSDF(float distA, float distB) {
    return max(distA, -distB);
}

float sphereSDF(vec3 samplePoint) {
    return length(samplePoint) - 1.0;
}

float cubeSDF(vec3 p) {
    vec3 d = abs(p) - vec3(1.0, 1.0, 1.0);
    float insideDistance = min(max(d.x, max(d.y, d.z)), 0.0);
    float outsideDistance = length(max(d, 0.0));
    
    return insideDistance + outsideDistance;
}

 float sdCappedCylinder( vec3 p, vec2 h )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

mat4 rotateZ(float theta) {
    float c = cos(theta);
    float s = sin(theta);

	return mat4(
        vec4(c, s, 0, 0),
        vec4(-s, c, 0, 0),
        vec4(0, 0, 1, 0),
        vec4(0, 0, 0, 1)
    );
}

mat4 rotateX(float theta) {
    float c = cos(theta);
    float s = sin(theta);

	return mat4(
        vec4(1, 0, 0, 0),
        vec4(0, c, s, 0),
        vec4(0, -s, c, 0),
        vec4(0, 0, 0, 1)
    );
}

float Cilinders(vec3 samplePoint) {
	vec3 CiliPoint2 = ((rotateZ(M_PI/2)) * vec4(samplePoint, 1.0)).xyz;
	vec3 CiliPoint3 = ((rotateX(M_PI/2)) * vec4(samplePoint, 1.0)).xyz;
	float CiliDist = sdCappedCylinder( samplePoint, vec2(1.6,3));
	float CiliDist2 = sdCappedCylinder( CiliPoint2, vec2(1.6,3));
	float CiliDist3 = sdCappedCylinder( CiliPoint3, vec2(1.6,3));
	return unionSDF(unionSDF(CiliDist, CiliDist2), CiliDist3);
}

float mysphere(vec3 samplePoint){
	return sphereSDF(samplePoint/3.3)*3.3;
}

float mycube(vec3 samplePoint) {
	return cubeSDF(samplePoint/2.5)*2.5;
}

float cubesphere(vec3 samplePoint) {
	float sphereDist= mysphere(samplePoint);
	float cubeDist= mycube(samplePoint);
	return intersectSDF(sphereDist, cubeDist);
}

float sceneSDF(vec3 samplePoint){
	float cs= cubesphere(samplePoint);
	float cili= Cilinders(samplePoint);
	return differenceSDF(cs, cili);
}

float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    float depth = start;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        float dist = sceneSDF(eye + depth * marchingDirection);
        if (dist < EPSILON) {
			return depth;
        }
        depth += dist;
        if (depth >= end) {
            return end;
        }
    }
    return end;
}

vec3 vecmul(vec3 v1, vec3 v2) {
  vec3 ret;
 
  ret.x = v1.y * v2.z - v1.z * v2.y;
  ret.y = v1.z * v2.x - v1.x * v2.z;
  ret.z = v1.x * v2.y - v1.y * v2.x;

  return ret;
}

vec3 rayDirection(float Degree, vec3 eye, vec3 at, vec3 up) {
    vec3 forward= (eye-at);
	forward= normalize(forward);
	vec3 right= vecmul(forward,up);
	right= normalize(right);
	vec3 myup= vecmul(forward,right);
	float alpha=mywin.x*tan(radians(Degree));
	float beta=mywin.y*tan(radians(Degree));
	return (vs_out_pos.x * right * alpha + vs_out_pos.y *myup * beta - forward);
}

vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity) {
    vec3 N = estimateNormal(p);
    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));
    
    float dotLN = dot(L, N);
    float dotRV = dot(R, V);
    
    if (dotLN < 0.0) {
        return vec3(0.0, 0.0, 0.0);
    } 
    
    if (dotRV < 0.0) {
        return lightIntensity * (k_d * dotLN);
    }
    return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    
    vec3 light1Pos = vec3(8.0 * sin(mytime/1000.0),
                          4.0,
                          8.0 * cos(mytime/1000.0));
    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light1Pos,
                                  light1Intensity);
    
    vec3 light2Pos = vec3(5.0 * sin(0.37 * mytime/1000.0),
                          5.0 * cos(0.37 * mytime/1000.0),
                          5.0);
    vec3 light2Intensity = vec3(0.8, 0.8, 0.8);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light2Pos,
                                  light2Intensity);    
    
    vec3 light3Pos = vec3(4.0,
                          10.0,
                          4.0);
    vec3 light3Intensity = vec3(0.7, 0.7, 0.7);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light3Pos,
                                  light3Intensity); 
    return color;
}

void main()
{
	vec3 dir = rayDirection(45.0, cameye, camat, camup);
	dir = normalize(dir);
    float dist = shortestDistanceToSurface(cameye, dir, MIN_DIST, MAX_DIST);
    
    if (dist > MAX_DIST - EPSILON) {
        // Didn't hit anything
        fs_out_col = vec4(0.0, 0.0, 0.0, 0.0);
		return;
    }

	vec3 p = cameye + dist * dir;

	if(mysphere(p)<EPSILON && mysphere(p)>-EPSILON){

    vec3 K_a = vec3(0.2, 0.2, 0.2);
    vec3 K_d = vec3(0.0, 1.0, 0.0);
    vec3 K_s = vec3(1.0, 0.0, 0.0);
    float shininess = 77.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, cameye);
    
    fs_out_col = vec4(color, 1.0);
	}

	if(mycube(p)<EPSILON && mycube(p)>-EPSILON){

    vec3 K_a = vec3(0.2, 0.2, 0.2);
    vec3 K_d = vec3(1.0, 0.0, 0.0);
    vec3 K_s = vec3(0.0, 0.0, 1.0);
    float shininess = 77.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, cameye);
    
    fs_out_col = vec4(color, 1.0);
	}

	if(Cilinders(p)<EPSILON && Cilinders(p)>-EPSILON){

    vec3 K_a = vec3(0.2, 0.2, 0.2);
    vec3 K_d = vec3(0.0, 0.0, 1.0);
    vec3 K_s = vec3(0.0, 1.0, 1.0);
    float shininess = 77.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, cameye);
    
    fs_out_col = vec4(color, 1.0);
	}

}
-----------------------------------------------------------------------------
ITT PRIMITIV MOZGASSAL
#version 130

in vec2 vs_out_pos;
out vec4 fs_out_col;

#define M_PI 3.1415926535897932384626433832795

uniform vec2 mywin;
uniform vec3 cameye;
uniform vec3 camat;
uniform vec3 camup;
uniform float mytime;

const int MAX_MARCHING_STEPS = 300;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;

float intersectSDF(float distA, float distB) {
    return max(distA, distB);
}

float unionSDF(float distA, float distB) {
    return min(distA, distB);
}

float differenceSDF(float distA, float distB) {
    return max(distA, -distB);
}

float sphereSDF(vec3 samplePoint) {
    return length(samplePoint) - 1.0;
}

float cubeSDF(vec3 p) {
    vec3 d = abs(p) - vec3(1.0, 1.0, 1.0);
    float insideDistance = min(max(d.x, max(d.y, d.z)), 0.0);
    float outsideDistance = length(max(d, 0.0));
    
    return insideDistance + outsideDistance;
}

 float sdCappedCylinder( vec3 p, vec2 h )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

mat4 rotateZ(float theta) {
    float c = cos(theta);
    float s = sin(theta);

	return mat4(
        vec4(c, s, 0, 0),
        vec4(-s, c, 0, 0),
        vec4(0, 0, 1, 0),
        vec4(0, 0, 0, 1)
    );
}

mat4 rotateX(float theta) {
    float c = cos(theta);
    float s = sin(theta);

	return mat4(
        vec4(1, 0, 0, 0),
        vec4(0, c, s, 0),
        vec4(0, -s, c, 0),
        vec4(0, 0, 0, 1)
    );
}

float Cilinders(vec3 samplePoint) {
	vec3 CiliPoint2 = ((rotateZ(M_PI/2)) * vec4(samplePoint, 1.0)).xyz;
	vec3 CiliPoint3 = ((rotateX(M_PI/2)) * vec4(samplePoint, 1.0)).xyz;
	float CiliDist = sdCappedCylinder( samplePoint, vec2(1.6,3));
	float CiliDist2 = sdCappedCylinder( CiliPoint2, vec2(1.6,3));
	float CiliDist3 = sdCappedCylinder( CiliPoint3, vec2(1.6,3));
	return unionSDF(unionSDF(CiliDist, CiliDist2), CiliDist3);
}

float mysphere(vec3 samplePoint){
	return sphereSDF(samplePoint/3.3)*3.3;
}

float mycube(vec3 samplePoint) {
	return cubeSDF(samplePoint/2.5)*2.5;
}

float cubesphere(vec3 samplePoint) {
	float sphereDist= mysphere(samplePoint);
	float cubeDist= mycube(samplePoint);
	return intersectSDF(sphereDist, cubeDist);
}

float sceneSDF(vec3 samplePoint){
	float cs= cubesphere(samplePoint + vec3(sin(mytime/1000.0), 0.0, cos(mytime/1000.0)));
	float cili= Cilinders(samplePoint + vec3(sin(mytime/1000.0), 0.0, cos(mytime/1000.0)));
	return differenceSDF(cs, cili);
}

float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    float depth = start;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        float dist = sceneSDF(eye + depth * marchingDirection);
        if (dist < EPSILON) {
			return depth;
        }
        depth += dist;
        if (depth >= end) {
            return end;
        }
    }
    return end;
}

vec3 vecmul(vec3 v1, vec3 v2) {
  vec3 ret;
 
  ret.x = v1.y * v2.z - v1.z * v2.y;
  ret.y = v1.z * v2.x - v1.x * v2.z;
  ret.z = v1.x * v2.y - v1.y * v2.x;

  return ret;
}

vec3 rayDirection(float Degree, vec3 eye, vec3 at, vec3 up) {
    vec3 forward= (eye-at);
	forward= normalize(forward);
	vec3 right= vecmul(forward,up);
	right= normalize(right);
	vec3 myup= vecmul(forward,right);
	float alpha=mywin.x*tan(radians(Degree));
	float beta=mywin.y*tan(radians(Degree));
	return (vs_out_pos.x * right * alpha + vs_out_pos.y *myup * beta - forward);
}

vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity) {
    vec3 N = estimateNormal(p);
    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));
    
    float dotLN = dot(L, N);
    float dotRV = dot(R, V);
    
    if (dotLN < 0.0) {
        return vec3(0.0, 0.0, 0.0);
    } 
    
    if (dotRV < 0.0) {
        return lightIntensity * (k_d * dotLN);
    }
    return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    
    vec3 light1Pos = vec3(8.0 * sin(mytime/1000.0),
                          4.0,
                          8.0 * cos(mytime/1000.0));
    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light1Pos,
                                  light1Intensity);
    
    vec3 light2Pos = vec3(5.0 * sin(0.37 * mytime/1000.0),
                          5.0 * cos(0.37 * mytime/1000.0),
                          5.0);
    vec3 light2Intensity = vec3(0.8, 0.8, 0.8);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light2Pos,
                                  light2Intensity);    
    
    vec3 light3Pos = vec3(4.0,
                          10.0,
                          4.0);
    vec3 light3Intensity = vec3(0.7, 0.7, 0.7);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light3Pos,
                                  light3Intensity); 
    return color;
}

void main()
{
	vec3 dir = rayDirection(45.0, cameye, camat, camup);
	dir = normalize(dir);
    float dist = shortestDistanceToSurface(cameye, dir, MIN_DIST, MAX_DIST);
    
    if (dist > MAX_DIST - EPSILON) {
        fs_out_col = vec4(0.0, 0.0, 0.0, 0.0);
		return;
    }

	vec3 p = cameye + dist * dir;

	if(mysphere(p + vec3(sin(mytime/1000.0), 0.0, cos(mytime/1000.0)))<EPSILON && mysphere(p + vec3(sin(mytime/1000.0), 0.0, cos(mytime/1000.0)))>-EPSILON){

    vec3 K_a = vec3(0.2, 0.2, 0.2);
    vec3 K_d = vec3(0.0, 1.0, 0.0);
    vec3 K_s = vec3(1.0, 0.0, 0.0);
    float shininess = 77.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, cameye);
    
    fs_out_col = vec4(color, 1.0);
	}

	if(mycube(p + vec3(sin(mytime/1000.0), 0.0, cos(mytime/1000.0)))<EPSILON && mycube(p + vec3(sin(mytime/1000.0), 0.0, cos(mytime/1000.0)) )>-EPSILON){

    vec3 K_a = vec3(0.2, 0.2, 0.2);
    vec3 K_d = vec3(1.0, 0.0, 0.0);
    vec3 K_s = vec3(0.0, 0.0, 1.0);
    float shininess = 77.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, cameye);
    
    fs_out_col = vec4(color, 1.0);
	}

	if(Cilinders(p + vec3(sin(mytime/1000.0), 0.0, cos(mytime/1000.0)))<EPSILON && Cilinders(p + vec3(sin(mytime/1000.0), 0.0, cos(mytime/1000.0)))>-EPSILON){

    vec3 K_a = vec3(0.2, 0.2, 0.2);
    vec3 K_d = vec3(0.0, 0.0, 1.0);
    vec3 K_s = vec3(0.0, 1.0, 1.0);
    float shininess = 77.0;
    
    vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, cameye);
    
    fs_out_col = vec4(color, 1.0);
	}

}


*/

