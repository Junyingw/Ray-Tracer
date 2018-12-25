/*
CSCI 480
Assignment 3 Raytracer

Name: <Junying Wang>
*/
/* Standard CPP Libraries */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
using namespace std;

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

// For recursive reflection, you can change it!
#define REFLECTION 3

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define FOV 60.0

/* Color buffer to hold the rgb values each 8-bit */
unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct point
{
	double x;
	double y;
	double z;
};

// Ray data structure
struct Ray 
{
  point org;
  point dir;
};

// Barycentric data structure
struct Bary 
{
  double x;
	double y;
	double z;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// Camera variables
point camera;

// Location of 4 corners of image plane
point Upleft;
point Upright;
point Downleft;
point Downright;
  

// Use to store screen pixel points
point** pix;

// Width and height of the image plane
double width, height;

// define the pi value
const double PI = 3.14159265;

bool isReflection = false;

bool isAntialiasing = false;

//Corss product of two points
point CrossProduct(point A, point B)
{
	point res;
	res.x = A.y*B.z - B.y*A.z;
	res.y = B.x*A.z - A.x*B.z;
	res.z = A.x*B.y - A.y*B.x;
	return res;
}

//Normalize point
point NormalizePoint(point p)
{
  point res;
  double base;
	base = sqrt(pow(p.x,2)+pow(p.y,2)+pow(p.z,2));
  if(base)
  {
    res.x = p.x / base;
		res.y = p.y / base;
		res.z = p.z / base;
  }
  return res;
}

// Calculate the dot product of two points 
double DotProduct(point A,point B)
{
	double res;
	res = A.x*B.x + A.y*B.y + A.z*B.z;
	return res;
}

//Get ray position, according to p(t) = p0 + dt
point GenerateRays(Ray ray, double t)
{
	point res;
	res.x = ray.org.x + t*(ray.dir.x);
	res.y = ray.org.y + t*(ray.dir.y);
	res.z = ray.org.z + t*(ray.dir.z);
	return res;
}

//Get reflection vector, given by norma,r=2(l.n)n-l
point GenerateReflections(point l, point n) 
{
  double d = DotProduct(l, n);
  point r;
  r.x = 2 * d * n.x - l.x;
  r.y = 2 * d * n.y - l.y;
  r.z = 2 * d * n.z - l.z;
  return r;
}

//Calculate distance between A and B
double CalculateDistance(point A, point B) 
{
	double res;
	res = sqrt(pow(A.x - B.x, 2) + pow(A.y - B.y, 2) + pow(A.z - B.z, 2));
  return res;
}

// Return the line segment between two vertices 
point VertexMinus(Vertex A, Vertex B)
{
	point res;
  res.x = A.position[0] - B.position[0];
	res.y = A.position[1] - B.position[1];
	res.z = A.position[2] - B.position[2];
	return res;
}

// Return the line segment between two points
point PointMinus(point A, point B)
{
	point res;
  res.x = A.x - B.x;
	res.y = A.y - B.y;
	res.z = A.z - B.z;
	return res;
}


// Get shpere normal at point p
point GetSphereNormal(point p, Sphere sphere)
{
	point n;
	n.x = (p.x-sphere.position[0])/sphere.radius;
	n.y = (p.y-sphere.position[1])/sphere.radius;
	n.z = (p.z-sphere.position[2])/sphere.radius;
	return n;
}


// Find magnitude of a point 
double CalculateMagnitude(point p)
{
	double res;
	res = sqrt(pow(p.x,2)+pow(p.y,2)+pow(p.z,2));
	return res;
}


// Get the divided value of a point 
point ScalarDivided(point A, double a)
{
	point B;

	B.x=0.0;
	B.y=0.0;
	B.z=0.0;

	if (abs(a)>1e-10) //a>0
	{
		B.x=A.x/a;
		B.y=A.y/a;
		B.z=A.z/a;
	}

	return B;
}

// Check if two points are equal
bool IsEqual(point A,point B)
{
	bool res;
	
	if ((abs(A.x-B.x)<1e-10) && (abs(A.y-B.y)<1e-10) && (abs(A.z-B.z)<1e-10))
		res=1;
	else res=0;

	return res;
}

// Change vertex to point 
point VertexToPoint(Vertex v)
{
	point res;
	res.x = v.position[0];
	res.y = v.position[1];
	res.z = v.position[2];
	return res;
}

//Get the parameter t, once ray intersect sphere
double SphereIntersection(Sphere sphere, Ray ray)
{
	double b,c,t1,t2;
	double t;
	double dis;

	t1 = 0;
	t2 = 0;

	// Sphere intersection equation	
	c = pow((ray.org.x-sphere.position[0]),2) 
			+ pow((ray.org.y-sphere.position[1]),2) 
			+ pow((ray.org.z-sphere.position[2]),2) 
			- pow(sphere.radius,2);
	b = 2*(ray.dir.x*(ray.org.x-sphere.position[0]) 
			+ ray.dir.y*(ray.org.y-sphere.position[1]) 
			+ ray.dir.z*(ray.org.z-sphere.position[2]));

  // Drop it if it is negtive
	if ( (pow(b,2)-4*c)<0 ) t = -1;

  // Get valid t
	if ( (pow(b, 2)-4*c)>=0 )
	{
		t1 = (((-1)*b)+sqrt(pow(b,2)-4*c))/2;
		t2 = (((-1)*b)-sqrt(pow(b,2)-4*c))/2;

		// if t1, t2 both < 0, drop them
  	if(t1 <= 0 && t2 <= 0) t = -1;

  	// Get the minimal value
 		if(t1 > 0 && t2 > 0) 
			t = min(t1, t2); 
  	else 
		{
			// if ray originates inside the sphere, drop it
			dis = max(t1, t2); 
  	  if(dis < 0.0001) 
				t = -1;
			else 
				t = dis;
		}
	}
	return t;
}

// Check if the ray is intersected with the sphere
bool CheckPoint_Sphere(Sphere sphere, Ray ray, double &distance)
{
	double b,c,t1,t2;
	double t;

	t1 = 0;
	t2 = 0;

	// Sphere intersection equation	
	b = 2*(ray.dir.x*(ray.org.x-sphere.position[0]) 
			+ray.dir.y*(ray.org.y-sphere.position[1]) 
			+ ray.dir.z*(ray.org.z-sphere.position[2]));

	c = pow((ray.org.x-sphere.position[0]),2) 
			+ pow((ray.org.y-sphere.position[1]),2) 
			+ pow((ray.org.z-sphere.position[2]),2) 
			- pow(sphere.radius,2);
	
	if ( (pow(b,2)-4*c)<0 ) 
		return false;

		// Points of intersection
		t1 = (((-1)*b)+sqrt(pow(b,2)-4*c))/2;
		t2 = (((-1)*b)-sqrt(pow(b,2)-4*c))/2;

	// If both of them < 0, return false
  if(t1 <= 0 && t2 <= 0) 
		return false;

  if(t1 > 0 && t2 > 0) 
		distance = min(t1, t2); 
  else 
		distance = max(t1, t2); 

  // If ray originates inside the sphere, return false
  if(distance < 0.0001) 
		return false;
  else
    return true;
}


//Get the parameter t, once ray intersect trangle plane
double TriangleIntersection(Triangle triangle, Ray ray)
{
	point AB, AC;
	double d = 0, t = 0;
	double res;
  point normal;

	// Get sides of the triangle to find the normal	
	AB = VertexMinus(triangle.v[0],triangle.v[1]);
	AC = VertexMinus(triangle.v[0],triangle.v[2]);

	// Get the normal
	normal = NormalizePoint(CrossProduct(AB, AC));

	// Calculate d, implicit form ax + by + cz +d = 0 
	d = -(normal.x)*triangle.v[0].position[0]
			-(normal.y)*triangle.v[0].position[1]
			- (normal.z)*triangle.v[0].position[2];
			
	// If ray parallel to plane, there is no intersection
  res = DotProduct(ray.dir, normal);
  if(res == 0.0) t=-1;
	else t=(-1)*(DotProduct(normal,ray.org)+d)/(DotProduct(normal,ray.dir));
  if(t <= 0.01)
    t = -1;
	return t;
}


// Check if the point lies inside the triangle
bool CheckPoint_Triangle(Triangle triangle, point P, Bary &bary)
{
	point A, B, C;
	Vertex a, b, c;
	point BA, AC, CB, AB;
	point PA, PC, PB;
	double whole, area1, area2, area3;
	point PAxPB, PBxPC, PCxPA;

	A = VertexToPoint(triangle.v[0]);
  B = VertexToPoint(triangle.v[1]);
	C = VertexToPoint(triangle.v[2]);

	PA = PointMinus(P, A);
	PC = PointMinus(P, C);
	PB = PointMinus(P, B);

	area3 = CalculateMagnitude(CrossProduct(PA,PB));
	area1 = CalculateMagnitude(CrossProduct(PB,PC));
	area2 = CalculateMagnitude(CrossProduct(PC,PA));

	whole = (area1 + area2 +area3);
	if(whole > 0)
	{
		bary.x = area1 / whole;
		bary.y = area2 / whole;
		bary.z = area3 / whole;
	}
	
	PAxPB = ScalarDivided(CrossProduct(PA,PB),area3);	
	PBxPC = ScalarDivided(CrossProduct(PB,PC),area1);	
	PCxPA = ScalarDivided(CrossProduct(PC,PA),area2);	

	// Check if the point is a valid point
	if(IsEqual(PAxPB,PBxPC) && IsEqual(PAxPB,PCxPA))
		return true;
	else 
		return false;
}

// Check if ray intersectes with light source 
bool LightIntersection(Ray ray, Light light, double &dis) 
{
  
  if(light.position[0] == ray.org.x && 
      light.position[1] == ray.org.y && 
      light.position[2] == ray.org.z)
    return false;

  dis = (light.position[0] - ray.org.x) / ray.dir.x;
  if(dis != (light.position[1] - ray.org.y) / ray.dir.y)
    return false;
  if(dis != (light.position[2] - ray.org.z)/ ray.dir.z)
    return false;
  return true;
}

// Interpolate specific color in triangle, once find the intersection point
point TriangleInterpolation(Triangle triangle, Bary &bary, int k)
{
	point P;
	
	if (k == 0)
	{
		P.x = bary.x*triangle.v[0].normal[0]
				 +bary.y*triangle.v[1].normal[0]
				 +bary.z*triangle.v[2].normal[0];

		P.y = bary.x*triangle.v[0].normal[1]
				 +bary.y*triangle.v[1].normal[1]
				 +bary.z*triangle.v[2].normal[1];

		P.z = bary.x*triangle.v[0].normal[2]
				 +bary.y*triangle.v[1].normal[2]
				 +bary.z*triangle.v[2].normal[2];
	}

	else if (k==1)
	{
		P.x = bary.x*triangle.v[0].color_diffuse[0]
				 +bary.y*triangle.v[1].color_diffuse[0]
				 +bary.z*triangle.v[2].color_diffuse[0];

		P.y = bary.x*triangle.v[0].color_diffuse[1]
				 +bary.y*triangle.v[1].color_diffuse[1]
				 +bary.z*triangle.v[2].color_diffuse[1];

		P.z = bary.x*triangle.v[0].color_diffuse[2]
				 +bary.y*triangle.v[1].color_diffuse[2]
				 +bary.z*triangle.v[2].color_diffuse[2];
	}

	else if (k==2)
	{
		P.x = bary.x*triangle.v[0].color_specular[0]
				 +bary.y*triangle.v[1].color_specular[0]
				 +bary.z*triangle.v[2].color_specular[0];

		P.y = bary.x*triangle.v[0].color_specular[1]
				 +bary.y*triangle.v[1].color_specular[1]
				 +bary.z*triangle.v[2].color_specular[1];

		P.z = bary.x*triangle.v[0].color_specular[2]
				 +bary.y*triangle.v[1].color_specular[2]
				 +bary.z*triangle.v[2].color_specular[2];
	}
	return P;
}


// If the intersect point is not in shadow, apply Phong shading
void PhongShading(point intersection, point normal, point &I, 
point kd, point ks, double shininess, point v) 
{
  for(int i = 0; i< num_lights; i++) 
	{
    bool isShadow = false;

    // Generate shadow ray to each of light
    point light, origin, direction;
  
    light.x = lights[i].position[0];
    light.y = lights[i].position[1];
    light.z = lights[i].position[2];

    origin = intersection;
    direction = NormalizePoint(PointMinus(light, origin));

    Ray shadowray; 
		shadowray.org = origin; 
		shadowray.dir = direction;

    // Calculate the distance to light source
    double LightDistance;
    LightDistance = CalculateDistance(light, origin);

    // Check if the showray intersects with other spheres before reaching the light source
    for(int j = 0; j < num_spheres; j++)
		{ 
      double SphereDistance;
      SphereDistance = SphereIntersection(spheres[j], shadowray);

      if(SphereDistance > 0)
      {
        if(CheckPoint_Sphere(spheres[j], shadowray, SphereDistance)) 
			  {
          point p = GenerateRays(shadowray, SphereDistance);
          SphereDistance = CalculateDistance(p, origin);
          if(SphereDistance <= LightDistance)
				  {
            isShadow = true;
          }
        }
      }      
    }

    // Check if the showray intersects with other triangles before reaching the light source
    for(int k = 0; k < num_triangles; k++) 
		{
			point p;
      Bary bary;
      double t;

			t = TriangleIntersection(triangles[k], shadowray);
			p = GenerateRays(shadowray,t);
      if(CheckPoint_Triangle(triangles[k], p, bary) && t > 0) 
			{
        double TriangleDistance;
        TriangleDistance = CalculateDistance(p, origin);
        if(TriangleDistance <= LightDistance)
				{
          isShadow = true;
        }
      }
    }

	  // Apply Phong shading
    if(!isShadow) 
		{
      point r;
      double LdotN, RdotV;  
      r = NormalizePoint(GenerateReflections(direction, normal));

      LdotN = DotProduct(direction, normal);
      if(LdotN < 0.0) LdotN = 0.0;
      if(LdotN > 1.0) LdotN = 1.0;
    
      RdotV = DotProduct(r, v);
      if(RdotV < 0.0) RdotV = 0.0;
      if(RdotV > 1.0) RdotV = 1.0;

      //Calculate intensity 
      I.x += lights[i].color[0] * (kd.x * LdotN + ks.x * pow(RdotV, shininess));
      I.y += lights[i].color[1] * (kd.y * LdotN + ks.y * pow(RdotV, shininess));
      I.z += lights[i].color[2] * (kd.z * LdotN + ks.z * pow(RdotV, shininess));
    }

  }
}

// Apply ray tracing
void RayTracing(Ray ray, point &color, int n) 
{
  // For recursive reflection
  if(n > REFLECTION) 
      return;
  bool isIntersect = false;
  double far = DBL_MAX;

  // When intersect sphere
  for(int i = 0; i < num_spheres; i++) 
	{
    double SphereDistance = 0.0;
    point normal;
    SphereDistance = SphereIntersection(spheres[i],ray);
    if(SphereDistance > 0)
    {
      if(CheckPoint_Sphere(spheres[i],ray,SphereDistance)) 
		  {
        if(SphereDistance < far) 
		  	{
          isIntersect = true;
          far = SphereDistance;
          
          point intersection = GenerateRays(ray, SphereDistance);
          normal = GetSphereNormal(intersection, spheres[i]);

          //Initial Phong color
          point phong;
          phong.x = 0.0; 
          phong.y = 0.0; 
          phong.z = 0.0;

          // Get the diffuse and specular colors;
          point kd, ks;

          kd.x = spheres[i].color_diffuse[0];
          kd.y = spheres[i].color_diffuse[1];
          kd.z = spheres[i].color_diffuse[2];

          ks.x = spheres[i].color_specular[0];
          ks.y = spheres[i].color_specular[1];
          ks.z = spheres[i].color_specular[2];

         
          double shininess = spheres[i].shininess;

          point v;
          v = NormalizePoint(PointMinus(camera, intersection)); 

          PhongShading(intersection, normal, phong, kd, ks, shininess, v);

          // Set the color to Phong color
          if(!isReflection)
          {
            color.x = phong.x * 255.0;
            color.y = phong.y * 255.0;
            color.z = phong.z * 255.0;
          } 
          
          // For recursive reflection
          if(isReflection) 
          {
            point reflectionColor; 
            reflectionColor.x = 0.0; 
            reflectionColor.y = 0.0; 
            reflectionColor.z = 0.0;
            color.x = pow(1 - ks.x, n+1) * phong.x * 255.0;
            color.y = pow(1 - ks.y, n+1) * phong.y * 255.0;
            color.z = pow(1 - ks.z, n+1) * phong.z * 255.0;
            RayTracing(ray, reflectionColor, n+1);
            color.x += ks.x * reflectionColor.x;
            color.y += ks.y * reflectionColor.y;
            color.z += ks.z * reflectionColor.z;
          }
        }
      } 
    }    
  }

  // When intersect triangle
  for(int j = 0; j< num_triangles; j++) 
	{
    point normal;
    Bary bary;
    double TriangleDistance = DBL_MAX;
		point P;

		TriangleDistance = TriangleIntersection(triangles[j], ray);
		P = GenerateRays(ray,TriangleDistance);

    if(CheckPoint_Triangle(triangles[j], P, bary) && TriangleDistance>0)
		{
      if(TriangleDistance < far) 
			{ 
				// Initiate phong color model
        point phong;
				point kd, ks;
				double shi;

        phong.x = 0.0; 
        phong.y = 0.0; 
        phong.z = 0.0;

        isIntersect = true;
        far = TriangleDistance;

        // Get the normal vector
				normal = TriangleInterpolation(triangles[j], bary, 0);
        
        // Get the diffuse and specular colors;
        kd = TriangleInterpolation(triangles[j], bary, 1);
        ks = TriangleInterpolation(triangles[j], bary, 2);

        // Get the shiniess, which is the shininess of the triangle
        shi = triangles[j].v[0].shininess * bary.x + 
            triangles[j].v[0].shininess * bary.y + 
            triangles[j].v[0].shininess * bary.z;

        // Get the v vector
        point v; 
        v = NormalizePoint(PointMinus(camera, P)); 

        // calculate the phong color
        PhongShading(P, normal, phong, kd, ks, shi, v);

         // Set the color to phong color
        if(!isReflection)
        {
          color.x = phong.x * 255.0;
          color.y = phong.y * 255.0;
          color.z = phong.z * 255.0;
        } 
         // For recursive reflection
        if(isReflection) 
        {
          point reflectionColor; 
          reflectionColor.x = 0.0; 
          reflectionColor.y = 0.0; 
          reflectionColor.z = 0.0; 
          color.x = pow(1 - ks.x, n+1) * phong.x * 255.0;
          color.y = pow(1 - ks.y, n+1) * phong.y * 255.0;
          color.z = pow(1 - ks.z, n+1) * phong.z * 255.0;
          RayTracing(ray, reflectionColor, n+1);
          color.x += ks.x * reflectionColor.x;
          color.y += ks.y * reflectionColor.y;
          color.z += ks.z * reflectionColor.z;
        }
      }
    }
  }

  // When intersect with lightsource
  for(int k  = 0; k < num_lights; k++) 
  {
    double LightDistance;
    if(LightIntersection(ray, lights[k], LightDistance)) 
    {
      if(LightDistance < far) 
      {
        isIntersect = true;
        far = LightDistance;
        color.x = lights[k].color[0] * 255.0;
        color.y = lights[k].color[1] * 255.0;
        color.z = lights[k].color[2] * 255.0;
      }
    }
  }

  if(!isIntersect) 
	{
    color.x = 255.0; 
    color.y = 255.0; 
    color.z = 255.0;
  }
  // Add ambient light once
  if(isIntersect) 
	{
    color.x += ambient_light[0] * 255;
    color.y += ambient_light[1] * 255;
    color.z += ambient_light[2] * 255;
  }
  
  color.x = max(min(color.x, 255.0), 0.0);
  color.y = max(min(color.y, 255.0), 0.0);
  color.z = max(min(color.z, 255.0), 0.0);
}


void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void GetScreenPixel()
{
  // Shoot the ray and do ray tracing
  double y = Downleft.y;

  for (int h = 0; h < HEIGHT; h++) 
	{ 
    double x = Downleft.x;
    for(int w = 0; w < WIDTH; w++) 
		{ 
      point pixelPos; 
      Ray newRay;
      point color; 
      point direction;

      pixelPos.x = x; 
			pixelPos.y = y; 
			pixelPos.z = -1.0; 

      direction = PointMinus(pixelPos, camera);
      direction = NormalizePoint(direction);

      newRay.org = camera;
      newRay.dir = direction;

			color.x = 0.0; 
			color.y = 0.0; 
			color.z = 0.0;
 
      RayTracing(newRay, color,0);

      pix[w][h].x = color.x;
      pix[w][h].y = color.y;
      pix[w][h].z = color.z;

      x += width / (WIDTH );
    }
    y += height / (HEIGHT);
  }
}

//Using average filtering to reduce aliasing
void AverageFiltering()
{
  int col = 0;
	for(int y=0; y<HEIGHT; y++) 
  {
    glPointSize(2.0);
    glBegin(GL_POINTS);
    int row = 0;
		for(int x=0; x<WIDTH; x++) 
    {
      double r = 0.0, g = 0.0, b = 0.0;
      if(x!=0 && x!=WIDTH-1 && y!=0 && y!=HEIGHT-1) 
      {
			  r = ((pix[x-1][y-1].x) + (pix[x][y-1].x)+ (pix[x+1][y-1].x)
						+ (pix[x-1][y].x) + (pix[x][y].x) + (pix[x+1][y].x)
            + (pix[x-1][y+1].x) + (pix[x][y+1].x) + (pix[x+1][y+1].x))/9;

        g = ((pix[x-1][y-1].y) + (pix[x][y-1].y)+ (pix[x+1][y-1].y)
						+ (pix[x-1][y].y) + (pix[x][y].y) + (pix[x+1][y].y)
            + (pix[x-1][y+1].y) + (pix[x][y+1].y) + (pix[x+1][y+1].y))/9;

        b = ((pix[x-1][y-1].z) + (pix[x][y-1].z)+ (pix[x+1][y-1].z)
           + (pix[x-1][y].z) + (pix[x][y].z) + (pix[x+1][y].z)
           + (pix[x-1][y+1].z) + (pix[x][y+1].z) + (pix[x+1][y+1].z))/9;
      }
      else
      {
        r = pix[x][y].x;
        g = pix[x][y].y;
        b = pix[x][y].z;
      }
      plot_pixel(row, col, r, g, b);
      row++;
    }
      glEnd();
      glFlush();
      col++;
  }
  printf("Done!\n");
  fflush(stdout);
}

void draw_scene()
{
  if(isAntialiasing)
  {
    AverageFiltering();
  }
  else
  {
    int col = 0;
	  for(int y=0; y<HEIGHT; y++) 
    {
      glPointSize(2.0);
      glBegin(GL_POINTS);
      int row = 0;
		  for(int x=0; x<WIDTH; x++) 
      {
        double r = 0.0, g = 0.0,b = 0.0;
        r = pix[x][y].x;
        g = pix[x][y].y;
        b = pix[x][y].z;
        plot_pixel(row, col, r, g, b);
        row++;
      }
      glEnd();
      glFlush();
      col++;
    }
    printf("Done!\n");
    fflush(stdout);
  }
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}


void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void SetupImagePlane() 
{
  double x = ((double)WIDTH/(double)HEIGHT) * tan(FOV / 2 * (PI / 180));
  double y = tan(FOV / 2 * (PI / 180));
  double z = -1.0;

  width = 0.0;
  height = 0.0;

  // Set up four corners of image plane
  Upleft.x = -x;    
  Upleft.y = y;       
  Upleft.z = z;

  Upright.x = x;    
  Upright.y = y;      
  Upright.z = z;

  Downleft.x = -x; 
  Downleft.y = -y;   
  Downleft.z = z;

  Downright.x = x; 
  Downright.y = -y;  
  Downright.z = z;

  // Get width and height of image plane
  width = 2.0 * x;
  height = 2.0 * y;
}

// Initialize the screen pixel
void initScreenPixel() 
{
  pix = new point* [WIDTH];
  for(int i = 0; i < WIDTH; i++) 
  {
    pix[i] = new point[HEIGHT];
  }
}

void display()
{
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
	  save_jpg();
  }
  once=1;
}

// Use keyboard to quit
void keyboardFunc( unsigned char key, int x, int y )
{
    if(key == 'q'|| key == 'Q')
    exit(0);    
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);

  // Generate empty image plane
  initScreenPixel();

  // Set up the image plane, according to aspect ratio
  SetupImagePlane();

  //Initialize camera position
  camera.x = 0; 
  camera.y = 0; 
  camera.z = 0;

	glLoadIdentity();
	GetScreenPixel();
}

void idle()
{

}

int main (int argc, char ** argv)
{
	
  if ((argc<2 || argc > 5)&& argc !=3)
  {  
    printf ("usage: %s <screenfile> [1 for recursive reflection][1 for antialiasing][jpegname]\n", argv[0]);
  }
  if(argc == 5)
  {
    if(atoi(argv[2]) == 1) isReflection = true;
    if(atoi(argv[3])== 1) isAntialiasing = true;
    mode = MODE_JPEG;
    filename = argv[4];
  }
  if(argc == 4)
  {
    if(atoi(argv[2]) == 1) isReflection = true;
    if(atoi(argv[3])== 1) isAntialiasing = true;
    mode = MODE_DISPLAY;
  }
  if(argc == 2)
  mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);
  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glutKeyboardFunc(keyboardFunc);
  init();
  glutMainLoop();
}