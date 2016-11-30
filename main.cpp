
// *********************************************************************************************************************************
// PA-3: Path Tracing
//
// Name: Brian Yee________________________
//
// ID: 00993104___________________________
//
// * Individual Work
// * Write the missing code -- i.e.,  // TO DO Task # ...
// * Compile and test it.
//
// TO DO Task #1: write the code for the procedure GenerateRandomFloat ()
// TO DO Task #2: in HemisphereSampling (), write the code for creating a random 3D ray (hemisphere sampling)
// TO DO Task #3: in m_PathTracer (), write the code for the core of the Path Tracing MC
// TO DO Task #4: in m_PathTracer (), write the code informing which procedure to call to find the direction of the diffuse ray
// TO DO Task #5: in main (), write the code to jitter (perturb) the camera's (x,y) position
//
// *********************************************************************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <ctime>
#include <vector>
#include <string>
#include <unordered_map>
#include <random>
#include <cstdint>
#include <algorithm>

// -----------------------------------------
//    RANDOM NUMBER GENERATION
// -----------------------------------------

// Mersenne Twister 19937 generator
// A Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
// http://www.cplusplus.com/reference/random/mt19937/

std::mt19937 MersenneTwisterPRNG;

// Uniform real distribution
// Random number distribution that produces floating-point values according to a uniform distribution, 
// which is described by the following probability density function (...) - more details here:
// http://www.cplusplus.com/reference/random/uniform_real_distribution/

std::uniform_real_distribution<double> m_URD;

// Hemisphere sampling (i.e., for diffuse reflection) 
// function "m_Vector HemisphereSampling(m_Vector m_normal)" below
// calls both m_RND_1 and m_RND_2
#define m_RND_1 (2.0*m_URD(MersenneTwisterPRNG)-1.0)
#define m_RND_2 (m_URD(MersenneTwisterPRNG))


// -----------------------------------------
//    PI & tolerances
// -----------------------------------------
#define m_PI 3.1415926536

const double m_inf=1e9;
const double m_eps=1e-6;

// -----------------------------------------
//    Final Path Traced Rendered Scene
//    Number of Samples per Pixel
// -----------------------------------------
const int m_pixmap_width = 255, m_pixmap_height = 255;
const double m_samples_per_pixel = 10;

using namespace std;

// =============================================
// This program has seven classes:
//
// m_Vector
// m_Ray, 
// m_Object
// m_Plane
// m_Sphere
// m_Intersection
// m_Scene
//
// =============================================

// -----------------------------------------
// VECTOR Class
// -----------------------------------------
struct m_Vector {
	double x, y, z;

	m_Vector (double x0=0, double y0=0, double z0=0){ x=x0; y=y0; z=z0; }
	m_Vector operator+(const m_Vector &b) const { return m_Vector(x+b.x,y+b.y,z+b.z); }
	m_Vector operator-(const m_Vector &b) const { return m_Vector(x-b.x,y-b.y,z-b.z); }
	m_Vector operator*(double b) const { return m_Vector(x*b,y*b,z*b); }
	m_Vector operator/(double b) const { return m_Vector(x/b,y/b,z/b); }
	m_Vector operator%(const m_Vector &b) const {return m_Vector(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}

	m_Vector Multiply(const m_Vector &b) const { return m_Vector(x*b.x,y*b.y,z*b.z); }
	m_Vector& Normalize(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
	double DotProduct(const m_Vector &b) const { return x*b.x+y*b.y+z*b.z; }
};

// -----------------------------------------
// RAY Class
// -----------------------------------------
struct m_Ray 
{
	m_Vector origin, direction; 

	m_Ray(m_Vector o0 = 0, m_Vector d0 = 0) 
	{ 
		origin = o0, 
		direction = d0.Normalize(); 
	}
};

// -----------------------------------------
// 3D OBJECT Class
// -----------------------------------------
class m_Object 
{
	public:
	m_Vector color;
	double emission;
	double kd, ks, kr;

	void setMaterial(m_Vector cl_ = 0, double emission_ = 0, double kd_ = 0, double ks_ = 0, double kr_ = 0) 
	{ color=cl_; emission=emission_; kd=kd_; ks=ks_; kr=kr_; }
	virtual double intersect(const m_Ray&) const = 0;
	virtual m_Vector normal(const m_Vector&) const = 0;
};

// -----------------------------------------
// PLANE Class
// -----------------------------------------
class m_Plane : public m_Object 
{
	public:
	double d;
	m_Vector n;
	
	m_Plane(double d_ = 0, m_Vector n_= 0) 
	{
		d=d_;
		n=n_;
	}

	double intersect(const m_Ray& ray) const 
	{
		double d0 = n.DotProduct(ray.direction);
		if(d0 != 0) 
		{
			double t = -1 * (((n.DotProduct(ray.origin))+d) / d0);
			return (t > m_eps) ? t : 0.0;
		}
		else return 0.0;
	}
	m_Vector normal(const m_Vector& p0) const { return n; }
};

// -----------------------------------------
// SPHERE Class
// -----------------------------------------
class m_Sphere : public m_Object 
{
	public:
	m_Vector center;
	double radius;

	m_Sphere(double r_= 0, m_Vector c_=0) 
	{ 
		center = c_; 
		radius = r_; 
	}

	double intersect(const m_Ray& ray) const 
	{
		double b = ((ray.origin-center)*2).DotProduct(ray.direction);
		double r2 = (radius*radius);
		double c = (ray.origin-center).DotProduct((ray.origin-center)) - r2;

		double disc = b*b - 4*c;
		if (disc < 0) return 0;
		else disc = sqrt(disc);

		double solution_1 = -b + disc;
		double solution_2 = -b - disc;
		return (solution_2 > m_eps) ? solution_2/2 : ((solution_1 > m_eps) ? solution_1/2 : 0);
	}

	m_Vector normal(const m_Vector& p0) const 
	{
		return (p0 - center).Normalize();
	}
};

// -----------------------------------------
// INTERSECTION Class
// -----------------------------------------
class m_Intersection 
{
	public:
		m_Object* object;

		// parametric distance (between 0.0 and 1.0) along the ray 
		// the value of t indicates how close or far the intersection was deteccted.
		double t;

	m_Intersection() { t = m_inf; object = nullptr; }
	m_Intersection(double t_, m_Object* object_) { t = t_; object = object_; }
	operator bool() { return object != nullptr; }
};

// -----------------------------------------
// SCENE Class
// -----------------------------------------
class m_Scene {
	// list of objects in the scene
	vector<m_Object*> objects;

	public:
	void add(m_Object* object) 
	{
		objects.push_back(object);
	}

	m_Intersection intersect(const m_Ray& ray) const 
	{
		m_Intersection closestIntersection;
		// test intersection with all objects in the scene
		for (auto iter = objects.begin(); iter != objects.end(); ++iter) 
		{
			double t = (*iter)->intersect(ray);

			if (t > m_eps && t < closestIntersection.t) 
			{
				closestIntersection.t = t;
				closestIntersection.object = *iter;
			}
		}
		return closestIntersection;
	}
};


// ==========================================================================
// Other Procedures:
//
// GenerateRandomFloat ()
// BuildOrthonormalSystem ()
// HemisphereSampling ()
//
// ProcessTransmissionRay ()
// m_PathTracer ()
// AddPlanes()
// AddSpheres()
// AddSphericalLightSource()
// ImagePlaneCoordinates ()
//
// main ()
//
// ==========================================================================

// ====================================================
// TO DO Task #1
//
// Write the code for the procedure
// GenerateRandomFloat ()
//
// Generate float random number between the
// float range min and max
//
// ====================================================
float GenerateRandomFloat(float min, float max)
{
    float r = static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX);  // random number in [0,1]
    r *= max - min;     // scale by delta between min/max
    r += min;           // offset by min
    return r;
}

// ----------------------------------------------------------------------------
// BuildOrthonormalSystem ()
//
// Generating outgoing ray directions by uniform sampling on a hemisphere 
//
// Input: vectors v1 and v2
// Output: vector v3 so (v1, v2, v3) form an orthonormal system
//         (assuming v1 is already normalized)
//
// ----------------------------------------------------------------------------
void BuildOrthonormalSystem(const m_Vector& v1, m_Vector& v2, m_Vector& v3) 
{
	float inverse_length, den;

    if (std::abs(v1.x) > std::abs(v1.y)) 
	{
		// project to the plane y = 0
		// construct a normalized orthogonal vector on this plane
		den = sqrtf(v1.x * v1.x + v1.z * v1.z);
		inverse_length = 1.f / den;
		v2 = m_Vector(-v1.z * inverse_length, 0.0f, v1.x * inverse_length);
    } 
	else 
	{
		// project to the plane x = 0 
		// construct a normalized orthogonal vector on this plane
		den = sqrtf(v1.y * v1.y + v1.z * v1.z);
		inverse_length = 1.0f / den; 
		v2 = m_Vector(0.0f, v1.z * inverse_length, -v1.y * inverse_length);
    }

	// construct v3 as the cross-product between v1 and v2
    v3 = v1 % v2;
}

// ----------------------------------------------------------------------------
// HemisphereSampling ()
//
// Generating outgoing ray directions by uniform sampling on a hemisphere 
//
// Input: normal vector at the intersection point
// Output: outgoing ray direction from uniform sampling on a hemisphere 
//
// ----------------------------------------------------------------------------

m_Vector HemisphereSampling(m_Vector m_normal)
{
// =================================================================
// TO DO Task #2
//
// Write the code for creating a random 3D ray (hemisphere sampling)
//
// vx = (...)?
// vy = (...)?
// vz = (...)?
//
// Refer to lecture slide file on Path Tracing:
// http://www.cpsc.ucalgary.ca/~mario/teaching/591-691/F16/topics/global/reading/slides/path%20tracing.pdf
//
// Slides #26 and #27
// Tip: use m_RND_2 for Xi_1 and Xi_2
// 
// =================================================================
    float r_1 = m_RND_2;
    float r_2 = m_RND_2;

    float r = sqrt(1 - r_1 * r_1);
    float phi = 2 * M_PI * r_2;

    double vx = cos(phi) * r;
    double vy = sin(phi) * r;
    double vz = r_1;

	m_Vector sampled_ray_direction = m_Vector(vx, vy, vz);

	// Now we build an otrhotnormal frame system 
	// to "rotate" the sampled ray to this system
	m_Vector v2, v3;
	BuildOrthonormalSystem (m_normal, v2, v3);

    // Construct the normalized rotated vector
	double vecx = m_Vector(v2.x, v3.x, m_normal.x).DotProduct(sampled_ray_direction);
    double vecy = m_Vector(v2.y, v3.y, m_normal.y).DotProduct(sampled_ray_direction);
    double vecz = m_Vector(v2.z, v3.z, m_normal.z).DotProduct(sampled_ray_direction);
	m_Vector m_rotated_ray_direction = m_Vector(vecx,vecy,vecz); 

	return m_rotated_ray_direction;
}

	
// Main scene
m_Scene scene;

// ----------------------------------------------------------------------------
// ProcessTransmissionRay ()
//
// Process a transmission ray (i.e. refraction, transparency)
//
// Input: normal vector at the intersection point and current ray
// Output: ray direction
//
// ----------------------------------------------------------------------------
m_Vector ProcessTransmissionRay (m_Vector normal, m_Vector ray_dir)
{
	m_Vector final_ray_direction;
	
	double n = 1.5;
	double a = 1.0 - n; 
	double b = 1.0 + n;
	double m_r0 = (a / b) * (a / b);

	if (normal.DotProduct(ray_dir) > 0) // ray is inside the medium
	{
		// invert the normal 
		normal = normal * -1.0;
		n = 1  / n;
	}

    // cosines of theta 1 and theta 2
	double cosine_t1 = (normal.DotProduct(ray_dir))*-1;
	double cosine_t2 = 1.0 - n*n*(1.0-cosine_t1*cosine_t1); 

	// Schlick's approximation:
	// approximates the contribution of the Fresnel factor in the specular reflection of light from a non-conducting interface (surface) between two media
	// http://www.cpsc.ucalgary.ca/~mario/teaching/591-691/F16/grading/PA-3/Schlick-1994.pdf
	// http://www.cpsc.ucalgary.ca/~mario/teaching/591-691/F16/grading/PA-3/Greve-2006.pdf
	//
	a = 1.0 - m_r0;
	b = pow(1.0 - cosine_t1, 5.0);
	double r_prob = m_r0 + a*b;
	m_Vector v1, v2;

	if (cosine_t2 > 0 && m_RND_2 > r_prob) 
	{
		// follow the refraction direction
		v1 = ray_dir * n;
		a = n * cosine_t1 - sqrt(cosine_t2);
		v2 = v1 + (normal * a);
		final_ray_direction = v2.Normalize();
	}
	else 
	{
		// follow the reflection direction
		a = 2*cosine_t1;
		v1 = normal * a;
		v2 = ray_dir + v1;
		final_ray_direction = v2.Normalize();
	}

	return final_ray_direction;
}

// -----------------------------------------
// m_PathTracer ()
//
// The main Path Tracing algorithm
//
// Input: ray and depth
// Output: color
//
// -----------------------------------------
void m_PathTracer (m_Ray &ray, int depth, m_Vector& color) 
{
	// Russian Roulette (rr) strategy
	// If depth >= 5
	// then Each recursive step will stop w/ a probability of 0.1
	double rr_factor = 1.0;
	if ((depth >= 5) & (m_RND_2 <= 0.1)) return;
	rr_factor = 1.0 / (1.0 - 0.1);
	
	// Find ray intsersection with the scene
	m_Intersection intersection = scene.intersect(ray);
	if (!intersection) return;

	// Compute intersection hit point position and its normal
	m_Vector hit_point = ray.origin + ray.direction * intersection.t;
	m_Vector normal_at_hit_point = intersection.object->normal(hit_point);
	ray.origin = hit_point;

	// Add the emission to the color, scaled e/ the Russian Roulette probability weight.
	const double emission = intersection.object->emission;
	color = color + m_Vector(emission, emission, emission) * rr_factor;

	// =================================================================
	// TO DO Task #3
	//
	// Write the code for the core of the Path Tracing MC
	//
	// ktot = (...)?
	// m_random_float = (...)?
	// 
	// Refer to lecture slide file on Path Tracing:
	// http://www.cpsc.ucalgary.ca/~mario/teaching/591-691/F16/topics/global/reading/slides/path%20tracing.pdf
	//
	// Slide #25
	//
	// =================================================================

    double ktot = intersection.object->kd + intersection.object->ks + intersection.object->kr;
    double m_random_float = GenerateRandomFloat(0, ktot);

	if (m_random_float < intersection.object->kd) // send a diffuse ray
	{
	
		// =========================================================================
		// TO DO Task #4
		//
		// Which procedure will you call to find the direction of the diffuse ray?
		//
	    // ray.direction = (...)?
		//
		// Refer to lecture slide file on Path Tracing:
		// http://www.cpsc.ucalgary.ca/~mario/teaching/591-691/F16/topics/global/reading/slides/path%20tracing.pdf
		//
		// Slides #26 and #27
		// =========================================================================

        ray.direction = HemisphereSampling(normal_at_hit_point);

		double cosine_t = ray.direction.DotProduct(normal_at_hit_point);
		m_Vector tmp;
		m_PathTracer(ray, depth+1, tmp);
		color = color + (tmp.Multiply(intersection.object->color)) * cosine_t * 0.1 * rr_factor;
	}
	else 
		if (m_random_float < (intersection.object->kd + intersection.object->ks)) // send a specular ray
		{
			
			double cosine_t = ray.direction.DotProduct(normal_at_hit_point);
			ray.direction = (ray.direction - normal_at_hit_point*(2*cosine_t)).Normalize();
			m_Vector tmp_color = m_Vector(0,0,0);
			m_PathTracer(ray, depth+1, tmp_color);
			color = color + tmp_color * rr_factor;
		}
		else // send a transmission (refraction) ray
		{
			ray.direction = ProcessTransmissionRay (normal_at_hit_point, ray.direction);
			m_Vector tmp_color;
			m_PathTracer(ray, depth+1, tmp_color);
			color = color + tmp_color * 1.15 * rr_factor;
		}
}

// -----------------------------------------
// Add planes to the scene
// -----------------------------------------
void AddPlanes()
{
	double m_d, m_emission;
	m_Vector m_color, m_normal;
	double m_kd, m_ks, m_kr;

	// -- BACK PLANE
	// Construct a plane
	m_d = 5.5; m_normal = m_Vector(0,0,1);
	m_Plane *m_back_plane = new m_Plane(m_d, m_normal);
	// Specify color, emission and material coefficients (kd, ks, kr)
	m_color = m_Vector(6, 6, 6); m_emission = 0; 
	m_kd = 1.0; m_ks = 0.0; m_kr = 0.0;
	m_back_plane->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
	// Add plane to the scene
	scene.add(m_back_plane);

	// -- LEFT PLANE
	m_d = 5.5; m_normal = m_Vector(1,0,0);
	m_Plane *m_left_plane = new m_Plane(m_d, m_normal);
	m_color = m_Vector(10, 2, 2); m_emission = 0;  
	m_kd = 1.0; m_ks = 0.0; m_kr = 0.0;
	m_left_plane->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
	scene.add(m_left_plane);

	// -- RIGHT PLANE
	m_d = 2.75; m_normal = m_Vector(-1, 0, 0);
	m_Plane *m_right_plane = new m_Plane(m_d, m_normal);
	m_color = m_Vector(2,10,2); m_emission = 0; 
	m_kd = 1.0; m_ks = 0.0; m_kr = 0.0;
	m_right_plane->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
	scene.add(m_right_plane);

	// -- CEILING PLANE
	m_d = 3.0; m_normal = m_Vector(0, -1, 0);
	m_Plane *m_ceiling_plane = new m_Plane(m_d, m_normal);
	m_color = m_Vector(6,6,6); m_emission = 0; 
	m_kd = 1.0; m_ks = 0.0; m_kr = 0.0;
	m_ceiling_plane->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
	scene.add(m_ceiling_plane);

	// -- FLOOR PLANE
	m_d = 2.5; m_normal = m_Vector(0, 1, 0);
	m_Plane *m_floor_plane = new m_Plane(m_d, m_normal);
	m_color = m_Vector(6,6,6); m_emission = 0; 
	m_kd = 0.0; m_ks = 1.0; m_kr = 0.0;
	m_floor_plane->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
	scene.add(m_floor_plane);
}

// -----------------------------------------
// Add spheres to the scene
// -----------------------------------------
void AddSpheres()
{
	double m_radius, m_emission;
	m_Vector m_position, m_color, m_normal;
	double m_kd, m_ks, m_kr;

	// -- SPHERE #1
	// Construct a sphere
	m_radius = 1.05; m_position = m_Vector(-0.75,-1.45,-4.4);
	m_Sphere *m_sphere1 = new m_Sphere(m_radius, m_position);
	// Specify color, emission and material coefficients (kd, ks, kr)
	m_color = m_Vector(4,8,4); m_emission = 0; 
	m_kd = 0.0; m_ks = 1.0; m_kr = 0.0;
	m_sphere1->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
    // Add sphere to the scene
	scene.add(m_sphere1);

	// -- SPHERE #2
	m_radius = 0.5; m_position = m_Vector(2.0,-2.05,-3.7);
	m_color = m_Vector(10,10,1); m_emission = 0; 
	m_kd = 0.7; m_ks = 0.0; m_kr = 0.3;
	m_Sphere *m_sphere2 = new m_Sphere(m_radius, m_position);
	m_sphere2->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
	//scene.add(m_sphere2);

	// -- SPHERE #3
	m_radius = 0.6; m_position = m_Vector(-1.75,-1.95,-3.1);
	m_color = m_Vector(10,10,1); m_emission = 0; 
	m_kd = 0.3; m_ks = 0.7; m_kr = 0.0;
	m_Sphere *m_sphere3 = new m_Sphere(m_radius, m_position);
	m_sphere3->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
	//scene.add(m_sphere3);
}

// -----------------------------------------
// Add spherical light source to the scene
// -----------------------------------------
void AddSphericalLightSource()
{
	double m_radius, m_emission;
	m_Vector m_position, m_color, m_normal;
	double m_kd, m_ks, m_kr;

	// Construct a sphere
	m_radius = 0.5; m_position = m_Vector(0,1.9,-3);
	m_Sphere *m_spherical_light_src = new m_Sphere(m_radius, m_position);
	// Specify color, emission and material coefficients (kd, ks, kr)
	m_color = m_Vector(4,8,4); m_emission = 10000; 
	m_kd = 1.0; m_ks = 0.0; m_kr = 0.0;
	m_spherical_light_src->setMaterial(m_color, m_emission, m_kd, m_ks, m_kr);
    // Add sphere to the scene
	scene.add(m_spherical_light_src);
}

// -----------------------------------------
// ImagePlaneCoordinates ()
//
// Input: pixel offset
// Output coordinate on the image plane
// -----------------------------------------

m_Vector ImagePlaneCoordinates(const double x, const double y) 
{
	double w = m_pixmap_width;
	double h = m_pixmap_height;

	float fovx = m_PI/4.0;
	float fovy = (h/w) * fovx;

	double vx = ((2*x-w)/w) * tan(fovx);
	double vy = -((2*y-h)/h) * tan(fovy);

	return m_Vector(vx, vy, -1.0);
}

// -----------------------------------------
// Main program
// -----------------------------------------
int main() {
	srand(time(NULL));

	// Buil d scene objects
	AddPlanes();
	AddSpheres();
	AddSphericalLightSource();

	// Initialize pixelmap
	m_Vector **m_pix = new m_Vector*[m_pixmap_width];
	for(int i = 0; i < m_pixmap_width; i++) 
	{
		m_pix[i] = new m_Vector[m_pixmap_height];
	}

	// start timing...
	clock_t start = clock();

    //#pragma omp parallel for schedule(dynamic)

	// for every pixel...
	for (int col = 0; col < m_pixmap_width; col++) 
	{
		float perc = (float) col/m_pixmap_width*100;
		fprintf(stdout,"\rRendering: %1.0f Samples per Pixel %8.2f%%",m_samples_per_pixel,(double)perc);

		for(int row = 0; row < m_pixmap_height; row++) 
		{
			// for every sample (i.e. ray) ...
			for(int s = 0; s < m_samples_per_pixel; s++) 
			{
				m_Vector color;
				m_Ray ray;

				ray.origin = (m_Vector(0,0,0)); 

				// construct the image plane coordinates
				m_Vector camera = ImagePlaneCoordinates(col,row); 
				
				// =========================================================================
				// TO DO Task #5
				//
				// Jitter (pertub) camera's (x,y) position
				//
				// camera.x = (...)?
				// camera.y = (...)?
				//
				// Refer to lecture slide file on Path Tracing:
				// http://www.cpsc.ucalgary.ca/~mario/teaching/591-691/F16/topics/global/reading/slides/path%20tracing.pdf
				//
				// Slide #20
				//
				// =================================================================
                float r_1 = m_RND_1;
                float r_2 = m_RND_1;
                float c = 10000;
                camera.x = camera.x + (r_1 / c);
                camera.y = camera.y + (r_2 / c);

				// Ray Direction: point from the origin to the camera plane
				ray.direction = (camera - ray.origin).Normalize(); 

				m_PathTracer(ray,0,color);

				// Update pixmap with color contribution
				m_pix[col][row] = m_pix[col][row] + color / m_samples_per_pixel;
			}
		}
	}

	// Saving final rendered scence in image file path-traced.ppm
	FILE *f = fopen("path_traced_scene.ppm", "w");
	fprintf(f, "P3\n%d %d\n%d\n ", m_pixmap_width, m_pixmap_height, 255);

	for (int row = 0; row < m_pixmap_height; row++) 
	{
		for (int col = 0; col < m_pixmap_width; col++) 
		{
			int x = (int)m_pix[col][row].x;
			int y = (int)m_pix[col][row].y;
			int z = (int)m_pix[col][row].z;

			fprintf(f,"%d %d %d ", min(x,255), min(y,255), min(z,255));
		}
		fprintf(f, "\n");
	}
	fclose(f);

	// stop timing
	clock_t end = clock();

	double t = (double)(end-start)/CLOCKS_PER_SEC;
	printf("\nRendering Time: %fs.\n",t);

	return 0;
}
