#include <bits/stdc++.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream> // header in standard library

#include "bitmap_image.hpp"
#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif

#define pi (2 * acos(0.0))
#define epsilon .000001

using namespace std;

class Object;
// class PointLight;
// class SpotLight;
// vector<Object *> objects;
// vector<PointLight> pointLights;
// vector<SpotLight> spotLights;
// int recursionLevel;




class Vector3D
{
public:
    double x;
    double y;
    double z;

    Vector3D()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }

    Vector3D(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vector3D normalize()
    {
        double length = sqrt(x * x + y * y + z * z);
        x = x / length;
        y = y / length;
        z = z / length;
        return *this;
    }

    double dot(Vector3D v2)
    {
        return x * v2.x + y * v2.y + z * v2.z;
    }

    Vector3D cross(Vector3D b)
    {
        Vector3D res;
        res.x = y * b.z - z * b.y;
        res.y = (-1) * (x * b.z - z * b.x);
        res.z = x * b.y - y * b.x;
        return res;
    }
        Vector3D operator-(Vector3D b) {
        Vector3D res;
        res.x = x - b.x;
        res.y = y - b.y;
        res.z = z - b.z;
        return res;
    }
    Vector3D operator+(Vector3D b) {
        Vector3D res;
        res.x = x + b.x;
        res.y = y + b.y;
        res.z = z + b.z;
        return res;
    }
    Vector3D operator*(double b) {
        Vector3D res;
        res.x = x * b;
        res.y = y * b;
        res.z = z * b;
        return res;
    }

    double distance(Vector3D b)
    {
        return sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y) + (z - b.z) * (z - b.z));
    }
};

class Ray
{
public:
    Vector3D start;
    Vector3D dir; // normalized direction vector

    Ray(Vector3D start, Vector3D dir)
    {
        this->start = start;
        this->dir = dir.normalize();
    }
};


// A point light emits light uniformly in all directions from a single point in space.
class PointLight
{
public:
    Vector3D light_pos;
    double color[3];
    // draw the pointlight
    void draw()
    {
        glPointSize(3);
        glBegin(GL_POINTS);
        glColor3f(color[0], color[1], color[2]);
        glVertex3f(light_pos.x, light_pos.y, light_pos.z);
        glEnd();
    }
};
// A spotlight emits light in a specific direction, typically defined by a cone.
class SpotLight
{
public:
    PointLight point_light;
    Vector3D light_direction;
    double cutoff_angle;
    void draw()
    {
        glPointSize(7);
        glBegin(GL_POINTS);
        glColor3f(point_light.color[0], point_light.color[1], point_light.color[2]); // Set the color of the point using the RGB values from the 'point_light' color array
        glVertex3f(point_light.light_pos.x, point_light.light_pos.y, point_light.light_pos.z);
        glEnd();
    }
};

vector<Object *> objects;
vector<PointLight> pointLights;
vector<SpotLight> spotLights;
int recursionLevel;

class Object
{
    string type;
    double generalCoefficients[10]; // general coefficients A,B,C,D,E,F,G,H,I,J
protected:
    double height, width, length;
    double color[3];
    double color_bmp[3]; // image texture color
    double coEfficients[4]; // ambient, diffuse, specular, reflection coefficients
    int shine;              // exponent term of specular component
public:
    Vector3D reference_point; // should have x, y, z
    Object(){

    };
    //for general object
    Object(double *generalCoefficients, Vector3D reference_point, double length, double width, double height)
    {
        for (int i = 0; i < 10; i++)
        {
            this->generalCoefficients[i] = generalCoefficients[i];
        }
        this->reference_point = reference_point;
        this->length = length;
        this->width = width;
        this->height = height;
        this->type = "general";
    }
    void setType(string type)
    {
        this->type = type;
    }
    string getType()
    {
        return type;
    }
    virtual void draw()
    {
   
    }
    void setColor(double color1, double color2, double color3)
    {
        this->color[0] = this->color_bmp[0] = color1;
        this->color[1] = this->color_bmp[1] = color2;
        this->color[2] = this->color_bmp[2] = color3;
    }
    void setShine(int shine)
    {
        this->shine = shine;
    }
    void setCoEfficients(double ambient, double diffuse, double specular, double recursive_reflection_coefficient)
    {
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = recursive_reflection_coefficient;
    }
    double getLength()
    {
        return length;
    }
    double *getColor()
    {
        return color;
    }

    void setColorBmp(double *colors)
    {
        this->color_bmp[0] = colors[0];
        this->color_bmp[1] = colors[1];
        this->color_bmp[2] = colors[2];
    }
    double *getColorBmp()
    {
        return color_bmp;
    }
   
   // this function handles the diffuse and specular reflection of a general object
    virtual void lightingWithPhongModel(Vector3D normal, Ray *v, Vector3D intersection_point, int level)
    {
    double *initColor;
    double *finalColor=new double[3];
    initColor = getColor();
   
    for (int i = 0; i < 3; i++)
    {
        finalColor[i] = initColor[i] * coEfficients[0]; // take ambient light into account
    }

    vector<PointLight> tempLights = pointLights;
     
     //iterates over spotlight sources and determines whether each spotlight contributes to the illumination of the current point on the surface. 

    for (auto &light : spotLights)
    {
        light.light_direction = light.light_direction.normalize(); 

        Vector3D dir;                                           ///Is, the direction from the light source to the point on the surface
        dir.x = intersection_point.x - light.point_light.light_pos.x;
        dir.y = intersection_point.y - light.point_light.light_pos.y;
        dir.z = intersection_point.z - light.point_light.light_pos.z;

        dir = dir.normalize();
     
      //claculate theta of spotlight
        double angle = acos(light.light_direction.dot(dir));
        angle = angle*(180 / pi);


        //This angle is compared against the cutoff angle of the spotlight to determine whether the light
        ///// from that spotlight contributes to the illumination of the current point on the surface.


        if (angle <= light.cutoff_angle)
        {
            tempLights.push_back(light.point_light);
        }
    }
   ////after converting spotlight to pointlight, we can use the same code for pointlight. Now all the lights are pointlights

    for (auto &light : tempLights)
    {
        Vector3D incidentRayDirection; //(L;the direction from the light source to the point on the surface)
        incidentRayDirection.x = intersection_point.x - light.light_pos.x;
        incidentRayDirection.y = intersection_point.y - light.light_pos.y;
        incidentRayDirection.z = intersection_point.z - light.light_pos.z;

        Ray *incidentRay = new Ray(light.light_pos, incidentRayDirection); //(L)
        Vector3D reflectedRayDirection;

        double l_dot_n = incidentRay->dir.dot(normal); //dir is the vector
        
        //reflected ray=incident ray−2×(incident ray dot normal)×normal.
        reflectedRayDirection.x = incidentRay->dir.x - 2 * l_dot_n * normal.x;
        reflectedRayDirection.y = incidentRay->dir.y - 2 * l_dot_n * normal.y;
        reflectedRayDirection.z = incidentRay->dir.z - 2 * l_dot_n * normal.z;

        Ray *reflectedRay = new Ray(intersection_point, reflectedRayDirection);  //R
        double *dummycolor;
    

        double dist = light.light_pos.distance(intersection_point);
        double t, tMin = INT_MAX;
       
       //shadow computation
        for (auto &o : objects)
        {
            // for each object, o in objects
            t = o->intersect(incidentRay, dummycolor, 0);

            //if any object is found between the light source and the intersection point,
            // the point is in shadow, and there's no need to consider further objects
            if (t < (dist - epsilon) && t > 0)
            {
                tMin = t;

                break;
            }
        }

        if (tMin != INT_MAX)
            continue;

        /**********************************************    diffuse and speccular reflection  *******************************************************************************/
        double diffuseComponent = max(0.0, -(incidentRay->dir.dot(normal)));  //diffuse reflection(Lambert value)

        //(R.V)^k  where V is the direction of the viewer and R is the direction of the reflected light,k is the shininess of the material
        double specularComponent = pow(max(0.0, -(v->dir.dot(reflectedRay->dir))), shine); //specular reflection (Phong value)

      // diffuse and specular reflection color calculation
        for (int i = 0; i < 3; i++)
        {
            finalColor[i] += light.color[i] * (coEfficients[1] * diffuseComponent + coEfficients[2] * specularComponent) * initColor[i];
         
            if (finalColor[i] < 0)
                finalColor[i] = 0;

               if (finalColor[i] > 1)
                finalColor[i] = 1;
        }
    }

    setColorBmp(finalColor);

    if (level >= recursionLevel)
        return;

     // After code for finding color components,now it is time for recursive reflection!!!!!!!!!!!!!!!!!!!

    //recursive reflection  (work with reflection coefficient)    **********************************************************
    Vector3D reflectedRayDirection;
   //v will now be the new L
    double l_dot_n = v->dir.dot(normal);
    reflectedRayDirection.x = v->dir.x - 2 * l_dot_n * normal.x;
    reflectedRayDirection.y = v->dir.y - 2 * l_dot_n * normal.y;
    reflectedRayDirection.z = v->dir.z - 2 * l_dot_n * normal.z;
     //construct reflected ray from intersection point
     // actually slightly forward from the point (by moving the
    // start a little bit towards the reflection direction)
    // to avoid self intersection

    Vector3D rr_start;
    rr_start.x = intersection_point.x + reflectedRayDirection.x * epsilon;
    rr_start.y = intersection_point.y + reflectedRayDirection.y * epsilon;
    rr_start.z = intersection_point.z + reflectedRayDirection.z * epsilon;
    Ray *rr = new Ray(rr_start, reflectedRayDirection);
    

    /**
    find tmin from the nearest intersecting object, using
    intersect() method, as done in the capture() method
    if found, call intersect(rreflected, colorreflected, level+1)
    // colorreflected will be updated while in the subsequent call
    // update color using the impact of reflection
    color += colorreflected * coEfficient[REC_REFFLECTION]
     * 
    */

    double t=INT_MAX;
    double tMin = INT_MAX;
    int nearest = -1;
    int idx = 0;
    double *color_reflected;
    double *dummycolor;
    for (auto &o : objects)
    {
        // for each object, o in objects
        t = o->intersect(rr, dummycolor, 0);
        if (t < tMin && t > 0)
        {
            tMin = t;
            nearest = idx;
        }
        idx++;
    }
    if (tMin != INT_MAX)
    {
        Object *n = objects[nearest];
        tMin = n->intersect(rr, color_reflected, level + 1);
        color_reflected = n->getColorBmp();
        for (int i = 0; i < 3; i++)
        {
            finalColor[i] += color_reflected[i] * coEfficients[3];
            if (finalColor[i] > 1)
                finalColor[i] = 1;
            if (finalColor[i] < 0)
                finalColor[i] = 0;
        }
    }
    setColorBmp(finalColor);
}

    // calculates the intersection of a ray with a general object defined by coefficients
    virtual double intersect(Ray *r, double *color, int level)
    {
        Vector3D translated;
        translated.x = r->start.x - reference_point.x;
        translated.y = r->start.y - reference_point.y;
        translated.z = r->start.z - reference_point.z;
      //Equation: F(x,y,z) = Ax2+By2+Cz2+Dxy+Exz+Fyz+Gx+Hy+Iz+J = 0

      // mapping this to the equation of a quadric surface at2+bt+c=0
        double a = generalCoefficients[0] * r->dir.x * r->dir.x +  // This term corresponds to the coefficient of x^2 in the general equation
                    generalCoefficients[1] * r->dir.y * r->dir.y + // This term corresponds to the coefficient of y^2 in the general equation
                    generalCoefficients[2] * r->dir.z * r->dir.z + // This term corresponds to the coefficient of z^2 in the general equation
                   generalCoefficients[3] * r->dir.x * r->dir.y +  // This term corresponds to the coefficient of xy in the general equation
                   generalCoefficients[4] * r->dir.y * r->dir.z + // This term corresponds to the coefficient of yz in the general equation
                   generalCoefficients[5] * r->dir.x * r->dir.z; // This term corresponds to the coefficient of xz in the general equation

        double b = 2 * generalCoefficients[0] * translated.x * r->dir.x + 
                   2 * generalCoefficients[1] * translated.y * r->dir.y +
                   2 * generalCoefficients[2] * translated.z * r->dir.z + 
                    generalCoefficients[3] * (r->dir.x * translated.y + translated.x * r->dir.y) + 
                    generalCoefficients[4] * (translated.y * r->dir.z + translated.z * r->dir.y) + 
                    generalCoefficients[5] * (r->dir.z * translated.x + translated.z * r->dir.x) + 
                    generalCoefficients[6] * r->dir.x +
                    generalCoefficients[7] * r->dir.y +
                    generalCoefficients[8] * r->dir.z;
                    
        double c = generalCoefficients[0] * translated.x * translated.x + 
                generalCoefficients[1] * translated.y * translated.y + 
                generalCoefficients[2] * translated.z * translated.z +
                 generalCoefficients[3] * translated.x * translated.y + 
                 generalCoefficients[4] * translated.y * translated.z + 
                 generalCoefficients[5] * translated.x * translated.z +
                  generalCoefficients[6] * translated.x + 
                  generalCoefficients[7] * translated.y + 
                  generalCoefficients[8] * translated.z + 
                  generalCoefficients[9];

        double determinant = sqrt(b * b - 4 * a * c);
        //no solution
        if (determinant < 0)
            return -1;

        double t1 = (-b - determinant) / (2 * a);
        double t2 = (-b + determinant) / (2 * a);
      /*
      The ray is represented by the equation r(t) = r0 + t * d, where r0 is the starting point of the ray and d is the direction of the ray.
      */
        Vector3D intersection_point1;
        intersection_point1.x = r->start.x + t1 * r->dir.x;
        intersection_point1.y = r->start.y + t1 * r->dir.y;
        intersection_point1.z = r->start.z + t1 * r->dir.z;

        Vector3D intersection_point2;
        intersection_point2.x = r->start.x + t2 * r->dir.x;
        intersection_point2.y = r->start.y + t2 * r->dir.y;
        intersection_point2.z = r->start.z + t2 * r->dir.z;

        if (length != 0.0)
        {  // check if the intersection point is within the length of the object
            if (intersection_point1.x < reference_point.x || intersection_point1.x > reference_point.x + length)
                t1 = -1;
            if (intersection_point2.x < reference_point.x || intersection_point2.x > reference_point.x + length)
                t2 = -1;
        }

        if (width != 0.0)
        {   // check if the intersection point is within the width of the object
            if (intersection_point1.y < reference_point.y || intersection_point1.y > reference_point.y + width)
                t1 = -1;
            if (intersection_point2.y < reference_point.y || intersection_point2.y > reference_point.y + width)
                t2 = -1;
        }

        if (height != 0.0)
        {   // check if the intersection point is within the height of the object
            if (intersection_point1.z < reference_point.z || intersection_point1.z > reference_point.z + height)
                t1 = -1;
            if (intersection_point2.z < reference_point.z || intersection_point2.z > reference_point.z + height)
                t2 = -1;
        }
       //choose the t
        double t;
        if (t1 > 0)
            t = t1;
        else
            t = t2;
        //When the level is 0, the purpose of the intersect() method is to determine the
        //nearest object only. No color computation is required
        if (level == 0)
            return t;
        
        //when the level is greater than 0, lighting codes according to the Phong model i.e. compute ambient,
        //diffuse, and specular components for each of the light sources and combine them.


        Vector3D intersection_point, normal;
        //based on which t is chosen, the intersection point is calculated accordingly
        intersection_point.x = r->start.x + t * r->dir.x;
        intersection_point.y = r->start.y + t * r->dir.y;
        intersection_point.z = r->start.z + t * r->dir.z;


    //The normal at the intersection point is calculated using the general equation of the object, equation is F(x,y,z) = Ax2+By2+Cz2+Dxy+Exz+Fyz+Gx+Hy+Iz+J = 0
    //The equation of normal is given by N = ∇F(x,y,z) = (F_x, F_y, F_z) where F_x, F_y, F_z are the partial derivatives of F(x,y,z) with respect to x, y, z respectively.
       double F_x = 2 * generalCoefficients[0] * intersection_point.x + generalCoefficients[3] * intersection_point.y + generalCoefficients[4] * intersection_point.z + generalCoefficients[6];
        double F_y = 2 * generalCoefficients[1] * intersection_point.y + generalCoefficients[3] * intersection_point.x + generalCoefficients[5] * intersection_point.z + generalCoefficients[7];
        double F_z = 2 * generalCoefficients[2] * intersection_point.z + generalCoefficients[4] * intersection_point.x + generalCoefficients[5] * intersection_point.y + generalCoefficients[8];

        normal.x = F_x;
        normal.y = F_y;
        normal.z = F_z;
        normal = normal.normalize();

        // diffuse and specular
        lightingWithPhongModel(normal, r, intersection_point, level);

        return t;
    }
};












class Sphere : public Object
{
public:
    Sphere(Vector3D center, double radius)
    {
        reference_point = center;
        length = radius;
    }
    void draw();
    double intersect(Ray *r, double *color, int level);
};

void Sphere::draw()
{

    double radius = length;
    glColor3f(color[0], color[1], color[2]);
    Vector3D points[100][100];
    int slices=24;
	int stacks=30;
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}
            glEnd();
		}
	}   

}

double Sphere::intersect(Ray *r, double *clr, int level)
{   
    //r is the castedRay initially
    // start-center
    Vector3D R_o; //origin of the ray
    R_o.x = r->start.x - reference_point.x;
    R_o.y = r->start.y - reference_point.y;
    R_o.z = r->start.z - reference_point.z;
    //at^2+bt+c=0
    //a=1,(because r.dir is normalized)
    //b=2R_d.R_o
    //c=R_o.R_o-r^2
    //algebraic solution
    double a = 1;
    double b = 2 * R_o.dot(r->dir);
    double c = R_o.dot(R_o) - length * length; //here length is the radius of the sphere

    double determinant = sqrt(b * b - 4 * a * c);
    double t1 = (-b - determinant) / (2 * a);
    double t2 = (-b + determinant) / (2 * a);
    double t;
    if (determinant < 0)
        return -1;
    else
    {
        if (t1 > 0)
            t = t1;
        else
            t = t2;
    }
    if (level == 0)
        return t;

    setColorBmp(getColor());

    Vector3D intersection_point, normal;
    //The ray is represented by the equation r(t) = r0 + t * d, where r0 is the starting point of the ray and d is the direction of the ray.
    intersection_point.x = r->start.x + t * r->dir.x;
    intersection_point.y = r->start.y + t * r->dir.y;
    intersection_point.z = r->start.z + t * r->dir.z;

    normal.x = intersection_point.x - reference_point.x;
    normal.y = intersection_point.y - reference_point.y;
    normal.z = intersection_point.z - reference_point.z;

    normal = normal.normalize();

    // diffuse and specular
    lightingWithPhongModel(normal, r, intersection_point, level);

    return t;
}

double calculateDet(double mat[3][3])
{
    double ans;
    ans = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    return ans;
}

class Triangle : public Object
{

    Vector3D point1;
    Vector3D point2;
    Vector3D point3;

public:
    Triangle(Vector3D point1, Vector3D point2, Vector3D point3)
    {
        reference_point.x = reference_point.y = reference_point.z = 0;
        this->point1 = point1;
        this->point2 = point2;
        this->point3 = point3;
    }
    void draw()
    {

    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_TRIANGLES);
    {
        glVertex3f(point1.x, point1.y, point1.z);
        glVertex3f(point2.x, point2.y, point2.z);
        glVertex3f(point3.x, point3.y, point3.z);
    }
    glEnd();
}
    double intersect(Ray *r, double *color, int level);
};


//using Barycentric Coordinates
double Triangle::intersect(Ray *r, double *color, int level)
{
    //firstly construct tthe left   side matrix
    //formula is leftSideMatrix = [point1-point2, point1-point3, r.dir]
    double leftSideMatrix[3][3];

    leftSideMatrix[0][0] = point1.x - point2.x;
    leftSideMatrix[0][1] = point1.x - point3.x;
    leftSideMatrix[0][2] = r->dir.x;
    leftSideMatrix[1][0] = point1.y - point2.y;
    leftSideMatrix[1][1] = point1.y - point3.y;
    leftSideMatrix[1][2] = r->dir.y;
    leftSideMatrix[2][0] = point1.z - point2.z;
    leftSideMatrix[2][1] = point1.z - point3.z;
    leftSideMatrix[2][2] = r->dir.z;

    double leftMatDet = calculateDet(leftSideMatrix);
    //if determinant is zero, then the matrix is singular and the ray is parallel to the plane of the triangle,so no intersection
    if (leftMatDet == 0)
        return -1;

    double matBeta[3][3];
    //formula of matBeta is [point1-start, point1-point3, r.dir]
   

    matBeta[0][0] = point1.x - r->start.x;
    matBeta[0][1] = point1.x - point3.x;
    matBeta[0][2] = r->dir.x;
    matBeta[1][0] = point1.y - r->start.y;
    matBeta[1][1] = point1.y - point3.y;
    matBeta[1][2] = r->dir.y;
    matBeta[2][0] = point1.z - r->start.z;
    matBeta[2][1] = point1.z - point3.z;
    matBeta[2][2] = r->dir.z;

    double detBeta = calculateDet(matBeta);
    double beta = detBeta / leftMatDet;

   ///////if beta is less than 0 or greater than 1, then the intersection point is outside the triangle
    if (beta <= 0)
        return -1;

    double matGamma[3][3];
    //formulaa of matGamma is [point1-point2, point1-start, r.dir]           
    matGamma[0][0] = point1.x - point2.x;    
    matGamma[0][1] = point1.x - r->start.x;
    matGamma[0][2] = r->dir.x;
    matGamma[1][0] = point1.y - point2.y;
    matGamma[1][1] = point1.y - r->start.y;
    matGamma[1][2] = r->dir.y;
    matGamma[2][0] = point1.z - point2.z;
    matGamma[2][1] = point1.z - r->start.z;
    matGamma[2][2] = r->dir.z;

    double detGamma = calculateDet(matGamma);
    double gamma = detGamma / leftMatDet;
//if gamma is less than 0 or greater than 1, then the intersection point is outside the triangle
    if (gamma <= 0)
        return -1;
//another condition is beta+gamma should be less than 1
    if (beta + gamma >= 1)
        return -1;

    double tMatrix[3][3];
    //formula of tMatrix is [point1-point2, point1-point3, point1-start]
    tMatrix[0][0] = point1.x - point2.x;
    tMatrix[0][1] = point1.x - point3.x;
    tMatrix[0][2] = point1.x - r->start.x;
    tMatrix[1][0] = point1.y - point2.y;
    tMatrix[1][1] = point1.y - point3.y;
    tMatrix[1][2] = point1.y - r->start.y;
    tMatrix[2][0] = point1.z - point2.z;
    tMatrix[2][1] = point1.z - point3.z;
    tMatrix[2][2] = point1.z - r->start.z;

    double detT = calculateDet(tMatrix);
    double t = detT / leftMatDet;
   
   ////if t is less than 0, then the intersection point is behind the ray's starting point
    if (t <= 0)
        return -1;
    ///If the recursion level is 0, the method returns the intersection parameter t without further calculations.
    if (level == 0)
        return t;

    Vector3D intersection_point, normal;
    //calculate the intersection point in 3D space using the parametric equation of the ray: intersection_point = r->start + t * r->dir
    //The ray is represented by the equation r(t) = r0 + t * d, where r0 is the starting point of the ray and d is the direction of the ray. (A point on the rayyy)

    intersection_point.x = r->start.x + t * r->dir.x;
    intersection_point.y = r->start.y + t * r->dir.y;
    intersection_point.z = r->start.z + t * r->dir.z;

    Vector3D AB, AC;
    AB.x = point2.x - point1.x;
    AB.y = point2.y - point1.y;
    AB.z = point2.z - point1.z;
    AC.x = point3.x - point1.x;
    AC.y = point3.y - point1.y;
    AC.z = point3.z - point1.z;
    normal = AB.cross(AC);

    normal = normal.normalize();

    // diffuse and specular
    lightingWithPhongModel(normal, r, intersection_point, level);

    return t;
}



class Prism : public Object
{
private:
    Vector3D vertices[6]; // Six vertices for a prism
public:
    Prism(Vector3D p1, Vector3D p2, Vector3D p3, Vector3D p4, Vector3D p5, Vector3D p6)
    {
        reference_point.x = reference_point.y = reference_point.z = 0;
        vertices[0] = p1;
        vertices[1] = p2;
        vertices[2] = p3;
        vertices[3] = p4;
        vertices[4] = p5;
        vertices[5] = p6;
    }

    void draw() override;
    double intersect(Ray *r, double *color, int level) override;
};

void Prism::draw()
{
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_TRIANGLES);
    {
        // Draw the triangular faces of the prism
        glVertex3f(vertices[0].x, vertices[0].y, vertices[0].z);
        glVertex3f(vertices[1].x, vertices[1].y, vertices[1].z);
        glVertex3f(vertices[2].x, vertices[2].y, vertices[2].z);

        glVertex3f(vertices[3].x, vertices[3].y, vertices[3].z);
        glVertex3f(vertices[4].x, vertices[4].y, vertices[4].z);
        glVertex3f(vertices[5].x, vertices[5].y, vertices[5].z);
    }
    glEnd();
}

double Prism::intersect(Ray *r, double *color, int level)
{
 
       double tMin = -1; // Initialize tMin to an invalid value

    // Check intersection with each face of the prism
    for (int i = 0; i < 5; ++i)
    {
        Vector3D p1 = vertices[i % 3];
        Vector3D p2 = vertices[(i + 1) % 3];
        Vector3D p3 = vertices[i >= 3 ? i - 2 : i + 3];

        Vector3D normal = (p2 - p1).cross(p3 - p1).normalize();
        double denom = normal.dot(r->dir);

        // Check if the ray is parallel or nearly parallel to the face
        if (std::abs(denom) < epsilon)
            continue;

        double t = normal.dot(p1 - r->start) / denom;

        // Check if the intersection point is behind the ray's starting point
        if (t <= 0)
            continue;

        // Check if the intersection point is inside the triangular region
        Vector3D intersection_point = r->start +  r->dir*t;
        Vector3D v1 = p2 - p1;
        Vector3D v2 = p3 - p1;
        Vector3D v3 = intersection_point - p1;

        double dot00 = v1.dot(v1);
        double dot01 = v1.dot(v2);
        double dot02 = v1.dot(v3);
        double dot11 = v2.dot(v2);
        double dot12 = v2.dot(v3);

        double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        if (u >= 0 && v >= 0 && u + v <= 1)
        {
            // Update tMin with the smallest positive intersection
            if (tMin == -1 || t < tMin)
                tMin = t;

            // Apply lighting model if necessary (similar to Triangle::intersect)
            if (level > 0)
            {
                // Calculate normal at the intersection point
                normal = normal.normalize();

                // Apply lighting model
                lightingWithPhongModel(normal, r, intersection_point, level);
            }

            // Set color and return tMin
            setColorBmp(getColor());
            return tMin;
        }
    }

    // No intersection with any face
    return -1;
}


class Floor : public Object
{
    int floorWidth;
    int tileWidth;

public:
    Floor(double floorWidth, double tileWidth)
    {
        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;
        reference_point.x = -floorWidth / 2;
        reference_point.y = -floorWidth / 2;
        reference_point.z = 0;
        length = tileWidth;
        this->setType("floor");
        setCoEfficients(0.4, 0.2, 0.2, 0.2);
        shine = 5;
    }
    void draw()
    {

        glPushMatrix();

        int row = (int)(floorWidth / tileWidth);
        int column = row;

        double init_x = reference_point.x;
        double init_y = reference_point.y;

        int color = 0;

        for (int i = 0; i < row; i++)
        {
            color = 1 - color;
            glColor3f(color, color, color);
            init_x = i * tileWidth;
            for (int j = 0; j < column; j++)
            {
                color = 1 - color;
                glColor3f(color, color, color);

                init_y = j * tileWidth;
                glBegin(GL_QUADS);
                {

                    // lines parallel to Y-axis
                    glVertex3f(init_x, init_y, 0);
                    glVertex3f(init_x + tileWidth, init_y, 0);

                    // lines parallel to X-axis
                    glVertex3f(init_x + tileWidth, init_y + tileWidth, 0);
                    glVertex3f(init_x, init_y + tileWidth, 0);
                }
                glEnd();
            }
        }
        glPopMatrix();
    }
    double intersect(Ray *r, double *color, int level);
};

double Floor::intersect(Ray *r, double *color, int level)
{
    Vector3D p0 = reference_point;
    Vector3D normal;
    normal.x = 0;
    normal.y = 0;
    normal.z = 1;

    double denom = r->dir.dot(normal);
    if (denom == 0)
        return -1;
    double nom = p0.dot(normal) - r->start.dot(normal);
    //
    double t = nom / denom;
    //
    Vector3D intersection_point;
    intersection_point.x = r->start.x + t * r->dir.x;
    intersection_point.y = r->start.y + t * r->dir.y;
    intersection_point.z = r->start.z + t * r->dir.z;

    if (intersection_point.x < -floorWidth / 2 || intersection_point.x > floorWidth / 2 || intersection_point.y < -floorWidth / 2 || intersection_point.y > floorWidth / 2)
        return -1;
    //
    //
    int row_index = (intersection_point.x - reference_point.x) / tileWidth;
    int col_index = (intersection_point.y - reference_point.y) / tileWidth;

    if (row_index % 2 == 0)
    {
        if (col_index % 2 == 0)
        {
            setColor(1, 1, 1);
        }
        else
        {
            setColor(0, 0, 0);
        }
    }
    else
    {
        if (col_index % 2 == 0)
        {
            setColor(0, 0, 0);
        }
        else
        {
            setColor(1, 1, 1);
        }
    }

    if (level == 0)
        return t;

    //    //diffuse and specular
    lightingWithPhongModel(normal, r, intersection_point, level);

    return t;
}
