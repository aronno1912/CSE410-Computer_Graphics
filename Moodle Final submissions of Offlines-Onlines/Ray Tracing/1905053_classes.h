#include <bits/stdc++.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream> 
#include "bitmap_image.hpp"
#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif


#define SMALLDELTA .000001
#define pi (2 * acos(0.0))

using namespace std;

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
        double val = sqrt(x * x + y * y + z * z);
        x = x / val;
        y = y / val;
        z = z / val;
        return *this;
    }

    double dot(Vector3D v2)
    {
        double res = x * v2.x + y * v2.y + z * v2.z;
        return res;
    }

    Vector3D cross(Vector3D b)
    {
        Vector3D res;
        res.x = y * b.z - z * b.y;
        res.y = (-1) * (x * b.z - z * b.x);
        res.z = x * b.y - y * b.x;
        return res;
    }
    Vector3D operator-(Vector3D b)
    {
        Vector3D res;
        res.x = x - b.x;
        res.y = y - b.y;
        res.z = z - b.z;
        return res;
    }
    Vector3D operator+(Vector3D b)
    {
        Vector3D res;
        res.x = x + b.x;
        res.y = y + b.y;
        res.z = z + b.z;
        return res;
    }
    Vector3D operator*(double b)
    {
        Vector3D res;
        res.x = x * b;
        res.y = y * b;
        res.z = z * b;
        return res;
    }

    Vector3D operator/(double b)
    {
        Vector3D res;
        res.x = x / b;
        res.y = y / b;
        res.z = z / b;
        return res;
    }

    double distance(Vector3D b)
    {
        double dis= sqrt((x - b.x) * (x - b.x) + (y - b.y) * (y - b.y) + (z - b.z) * (z - b.z));
        return dis;
    }

    void print()
    {
        cout << x << " " << y << " " << z << endl;
    }

    void setValue(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    void setValue(Vector3D v)
    {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
    }

    Vector3D getVector()
    {
        return *this;
    }

    ~Vector3D()
    {
    }

 
};

struct point
{
    GLfloat x, y, z;
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

    Ray()
    {
        start = Vector3D(0, 0, 0);
        dir = Vector3D(0, 0, 0);
    }
    ~Ray()
    {
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

    void print()
    {
        cout << "Light Position: ";
        light_pos.print();
        cout << "Color: " << color[0] << " " << color[1] << " " << color[2] << endl;
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
    void print()
    {
        cout << "Light Position: ";
        point_light.light_pos.print();
        cout << "Color: " << point_light.color[0] << " " << point_light.color[1] << " " << point_light.color[2] << endl;
        cout << "Light Direction: ";
        light_direction.print();
        cout << "Cutoff Angle: " << cutoff_angle << endl;
    }
};
class Object;
vector<Object *> objects;
vector<PointLight> pointLights;
vector<SpotLight> spotLights;
int recursionLevel;

class Object
{   
    protected:
    string type;
    double generalCoefficients[10]; // general coefficients A,B,C,D,E,F,G,H,I,J
    double height, width, length;
    double color[3];
    double photoColor[3];    // image texture color
    double coEfficients[4]; // ambient, diffuse, specular, reflection coefficients
    int shine;              // exponent term of specular component
    double refractiveIndex[3]; //for prism******************************************
public:
    Vector3D reference_point; // should have x, y, z
    Object()
    {

    };
    // for general object
    Object(double *generalCoefficients, Vector3D reference_point, double length, double width, double height)
    {
        this->type = "general";
        for (int i = 0; i < 10; i++)
        {
            this->generalCoefficients[i] = generalCoefficients[i];
        }
        this->reference_point = reference_point;
        this->length = length;
        this->width = width;
        this->height = height;
        
    }
    virtual void showObjectDetails()
    {
        cout << "Type: " << type << endl;
        cout << "Reference Point: ";
        reference_point.print();
        cout << "Color: " << color[0] << " " << color[1] << " " << color[2] << endl;
        cout << "Shine: " << shine << endl;
        cout << "Ambient Coefficient: " << coEfficients[0] << " Diffuse Coefficient: " << coEfficients[1] << " Specular Coefficient: " << coEfficients[2] << " Reflection Coefficient: " << coEfficients[3] << endl;
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
        this->color[0] = this->photoColor[0] = color1;
        this->color[1] = this->photoColor[1] = color2;
        this->color[2] = this->photoColor[2] = color3;
    }
    void setShine(int shine)
    {
        this->shine = shine;
    }
    double getShine()
    {
        return shine;
    }
    void setCoEfficients(double ambient, double diffuse, double specular, double recursive_reflection_coefficient)
    {
        this->coEfficients[0] = ambient;
        this->coEfficients[1] = diffuse;
        this->coEfficients[2] = specular;
        this->coEfficients[3] = recursive_reflection_coefficient;
    }
    void setLength(double length)
    {
        this->length = length;
    }
    void setWidth(double width)
    {
        this->width = width;
    }
    void setHeight(double height)
    {
        this->height = height;
    }
    double getLength()
    {
        return length;
    }
    double getWidth()
    {
        return width;
    }
    double getHeight()
    {
        return height;
    }
    double *getColor()
    {
        return color;
    }

    void setPhotoColor(double *colors)
    {
        this->photoColor[0] = colors[0];
        this->photoColor[1] = colors[1];
        this->photoColor[2] = colors[2];
    }
    double *getPhotoColor()
    {
        return photoColor;
    }

    virtual void setRefractiveIndex(double r, double g, double b)  //****************************for prism
    {
     
    }
    

    //after ray casting is done,here  v is the casted ray
    //now intersection point is already obtained, now we need to calculate the color of the intersection point
   //all reflections handled heeeeeeeere to determine color of the intersection point
    virtual void colorComputation(Vector3D intersection_point, Vector3D normal, Ray *v,  int level)
    {
        double *intersectionPointColor; //initial color of the object
        double *finalColor = new double[3];
        intersectionPointColor = getColor(); 
        
        //**************************************** Ambient Reflection ****************************************
        //color = intersectionPointColor*coEfficient[AMB]

        for (int i = 0; i < 3; i++)
        {
            finalColor[i] = intersectionPointColor[i] * coEfficients[0]; // take ambient light into account
        }

        vector<PointLight> lightSources = pointLights;

        // iterates over spotlight sources and determines whether each spotlight contributes to the illumination of the current point on the surface.

        for (auto &light : spotLights)
        {
            light.light_direction = light.light_direction.normalize(); //the constant direction of the spotlight

            Vector3D lightToIpVector; /// the vector from the light source point to the point of intersection on the surface
            lightToIpVector.x = intersection_point.x - light.point_light.light_pos.x;
            lightToIpVector.y = intersection_point.y - light.point_light.light_pos.y;
            lightToIpVector.z = intersection_point.z - light.point_light.light_pos.z;

            lightToIpVector = lightToIpVector.normalize();

            // claculate theta of spotlight
            double angle = acos(light.light_direction.dot(lightToIpVector));
            angle = angle * (180 / pi);

            // This angle is compared against the cutoff angle of the spotlight to determine whether the light
            ///// from that spotlight contributes to the illumination of the current point on the surface.

            if (angle <= light.cutoff_angle)
            {
                lightSources.push_back(light.point_light);
            }
        }
        ////after converting spotlight to pointlight, we can use the same code for pointlight. Now all the lights are pointlights
        //Phong Model calculation starts here ************************************************************************************************
        //Consider a point Light source p and viewpoiint v. wwwhat should be the color of light refleecttttedd into vieweer's eye from point Q on the surface?
        for (auto &light : lightSources)
        {
            Vector3D incidentRayDirection; //(L;the direction from the light source to the point on the surface)
            incidentRayDirection.x = intersection_point.x - light.light_pos.x;
            incidentRayDirection.y = intersection_point.y - light.light_pos.y;
            incidentRayDirection.z = intersection_point.z - light.light_pos.z;

            Ray *incidentRay = new Ray(light.light_pos, incidentRayDirection); //(L)
            Vector3D reflectedRayDirection;

            double l_dot_n = incidentRay->dir.dot(normal); // dir is the vector

            // reflected ray=incident ray−2×(incident ray dot normal)×normal.
            reflectedRayDirection.x = incidentRay->dir.x - 2 * l_dot_n * normal.x;
            reflectedRayDirection.y = incidentRay->dir.y -2 * l_dot_n * normal.y;
            reflectedRayDirection.z = incidentRay->dir.z - 2 * l_dot_n * normal.z;

            Ray *reflectedRay = new Ray(intersection_point, reflectedRayDirection); // R
            
            double dist = light.light_pos.distance(intersection_point);
            double t, tMin = INT_MAX;
            double *any;


            // shadow computation
            //// if intersectionPoint is in shadow, the diffuse and specular components need not be calculated

            for (auto &obj : objects)
            {
                // for each object, o in objects
                t = obj->intersect(incidentRay, any, 0);

                // if any object is found between the light source and the intersection point,
                //  the point is in shadow, and there's no need to consider further objects
                if (t < (dist - SMALLDELTA) && t > 0)
                {
                    tMin = t;

                    break;
                }
            }

            if (tMin != INT_MAX)
                continue;
            
            //if ray1 is not obscured by any object,then calculate diffuse and specular reflection

            /**********************************************    diffuse and speccular reflection  *******************************************************************************/
           // Calculate diffuse reflection (Lambert value)
            double diffuseComponent = -(incidentRay->dir.dot(normal));
            diffuseComponent = (diffuseComponent > 0) ? diffuseComponent : 0;

            //(R.V)^k  where V is the direction of the viewer and R is the direction of the reflected light,k is the shininess of the material
           
           // Calculate specular reflection (Phong value)
            double specularComponent = -(v->dir.dot(reflectedRay->dir));
            specularComponent = (specularComponent > 0) ? pow(specularComponent, shine) : 0;

            // diffuse and specular reflection color calculation
            for (int i = 0; i < 3; i++)
            {   
                //consider light color also
                finalColor[i] += light.color[i] * coEfficients[1] * diffuseComponent  * intersectionPointColor[i];
                finalColor[i] += light.color[i] * coEfficients[2] * specularComponent* intersectionPointColor[i];

               // Clamp to [0, 1]
                finalColor[i] = (finalColor[i] < 0) ? 0 : ((finalColor[i] > 1) ? 1 : finalColor[i]);
            }
        }

        setPhotoColor(finalColor);

        if (level >= recursionLevel)
            return;

        // After code for finding color components,now it is time for recursive reflection!!!!!!!!!!!!!!!!!!!

        //*****************************************recursive reflection  (work with reflection coefficient)    **********************************************************
        Vector3D recursiveReflectedRayDirection;
        // v will now be the new L
        double l_dot_n = v->dir.dot(normal);
        recursiveReflectedRayDirection.x = v->dir.x - 2 * l_dot_n * normal.x;
        recursiveReflectedRayDirection.y = v->dir.y - 2 * l_dot_n * normal.y;
        recursiveReflectedRayDirection.z = v->dir.z - 2 * l_dot_n * normal.z;
        // construct reflected ray from intersection point
        //  actually slightly forward from the point (by moving the
        // start a little bit towards the reflection direction)
        // to avoid self intersection

        Vector3D rr_start;
        rr_start.x = intersection_point.x + recursiveReflectedRayDirection.x * SMALLDELTA;
        rr_start.y = intersection_point.y + recursiveReflectedRayDirection.y * SMALLDELTA;
        rr_start.z = intersection_point.z + recursiveReflectedRayDirection.z * SMALLDELTA;
        Ray *rr = new Ray(rr_start, recursiveReflectedRayDirection);

        /**
        find tmin from the nearest intersecting object, using
        intersect() method, as done in the capture() method
        if found, call intersect(rreflected, colorreflected, level+1)
        // colorreflected will be updated while in the subsequent call
        // update color using the impact of reflection
        color += colorreflected * coEfficient[REC_REFFLECTION]
         *
        */

        double t = INT_MAX;
        double tMin = INT_MAX;
        int nearest = -1;
        int idx = 0;
        double *color_reflected;
        double *any;
        for (auto &o : objects)
        {
            // for each object, o in objects
            t = o->intersect(rr, any, 0);
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
            tMin = n->intersect(rr, color_reflected, level + 1); //level++
            color_reflected = n->getPhotoColor();
            for (int i = 0; i < 3; i++)
            {
                finalColor[i] += color_reflected[i] * coEfficients[3];
                finalColor[i] = (finalColor[i] < 0) ? 0 : ((finalColor[i] > 1) ? 1 : finalColor[i]);
            }
        }
        setPhotoColor(finalColor);
    }







    // calculates the intersection of a ray with a general object defined by coefficients
virtual double intersect(Ray *r, double *color, int level)
    {   
        //starting point of the ray relative to the reference point
        Vector3D R_0;
        R_0.x = r->start.x - reference_point.x;
        R_0.y = r->start.y - reference_point.y;
        R_0.z = r->start.z - reference_point.z;
        // general Quadratic Surface Equation: F(x,y,z) = Ax2+By2+Cz2+Dxy+Exz+Fyz+Gx+Hy+Iz+J = 0
        //Ray equation: P(t)=R0+t*R_d
        //So, P_x=Rx0+t*R_dx, P_y=Ry0+t*R_dy, P_z=Rz0+t*R_dz
        //put the values of P_x, P_y, P_z as x,y,z in the general equation of the object 


        /*
        
        Quadratic equation: (A*R_dx^2 + B*R_dy^2 + C*R_dz^2+ D*R_dx*R_dy +E*R_dz*R_dx + F*R_dy*R_dz) * t^2 +

                   (2*A*R_0x*R_dx + 2*B*R_0y*R_dy + 2*C*R_0z*R_dz +
                    D*R_dx*R_dy + E*R_dz*R_dx + F*R_dy*R_dz) * t +
                    G*R_0x + H*R_0y + I*R_0z 

                   (A*R_0x^2 + B*R_0y^2 + C*R_0z^2 +
                    D*R_0x*R_0y + E*R_0z*R_0x + F*R_0y*R_0z +
                    G*R_0x + H*R_0y + I*R_0z + J) = 0
        */

        // mapping this to the equation of a quadric surface at2+bt+c=0
        double a = generalCoefficients[0] * r->dir.x * r->dir.x + 
                   generalCoefficients[1] * r->dir.y * r->dir.y + 
                   generalCoefficients[2] * r->dir.z * r->dir.z +
                   generalCoefficients[3] * r->dir.x * r->dir.y + 
                   generalCoefficients[4] * r->dir.y * r->dir.z + 
                   generalCoefficients[5] * r->dir.x * r->dir.z;  

        double b = 2 * generalCoefficients[0] * R_0.x * r->dir.x +
                   2 * generalCoefficients[1] * R_0.y * r->dir.y +
                   2 * generalCoefficients[2] * R_0.z * r->dir.z +
                   generalCoefficients[3] * (r->dir.x * R_0.y + R_0.x * r->dir.y) +
                   generalCoefficients[4] * (R_0.y * r->dir.z + R_0.z * r->dir.y) +
                   generalCoefficients[5] * (r->dir.z * R_0.x + R_0.z * r->dir.x) +
                   generalCoefficients[6] * r->dir.x +
                   generalCoefficients[7] * r->dir.y +
                   generalCoefficients[8] * r->dir.z;

        double c = generalCoefficients[0] * R_0.x * R_0.x +
                   generalCoefficients[1] * R_0.y * R_0.y +
                   generalCoefficients[2] * R_0.z * R_0.z +
                   generalCoefficients[3] * R_0.x * R_0.y +
                   generalCoefficients[4] * R_0.y * R_0.z +
                   generalCoefficients[5] * R_0.x * R_0.z +
                   generalCoefficients[6] * R_0.x +
                   generalCoefficients[7] * R_0.y +
                   generalCoefficients[8] * R_0.z +
                   generalCoefficients[9];

        double determinant = sqrt(b * b - 4 * a * c);
        // no solution
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
        { // check if the intersection point is within the length of the object
            if (intersection_point1.x < reference_point.x || intersection_point1.x > reference_point.x + length)
                t1 = -1;
            if (intersection_point2.x < reference_point.x || intersection_point2.x > reference_point.x + length)
                t2 = -1;
        }

        if (width != 0.0)
        { // check if the intersection point is within the width of the object
            if (intersection_point1.y < reference_point.y || intersection_point1.y > reference_point.y + width)
                t1 = -1;
            if (intersection_point2.y < reference_point.y || intersection_point2.y > reference_point.y + width)
                t2 = -1;
        }

        if (height != 0.0)
        { // check if the intersection point is within the height of the object
            if (intersection_point1.z < reference_point.z || intersection_point1.z > reference_point.z + height)
                t1 = -1;
            if (intersection_point2.z < reference_point.z || intersection_point2.z > reference_point.z + height)
                t2 = -1;
        }
        // choose the t(the smallest non-negative availabe)
        double t;
        if (t1 > 0 && t2<0) 
            t = t1;
        else if (t2 > 0 && t1<0)
            t = t2;
        else if (t1 > 0 && t2 > 0)
        {
            t = min(t1, t2);
        }
        else
        {
            return -1;
        }


        // When the level is 0, the purpose of the intersect() method is to determine the
        // nearest object only. No color computation is required
        if (level == 0)
            return t;

        // when the level is greater than 0, lighting codes according to the Phong model i.e. compute ambient,
        // diffuse, and specular components for each of the light sources and combine them.

        Vector3D intersection_point, normal;
        // based on which t is chosen, the intersection point is calculated accordingly
        intersection_point.x = r->start.x + t * r->dir.x;
        intersection_point.y = r->start.y + t * r->dir.y;
        intersection_point.z = r->start.z + t * r->dir.z;

        // The normal at the intersection point is calculated using the general equation of the object, equation is F(x,y,z) = Ax2+By2+Cz2+Dxy+Exz+Fyz+Gx+Hy+Iz+J = 0
        // The equation of normal is given by N = ∇F(x,y,z) = (F_x, F_y, F_z) where F_x, F_y, F_z are the partial derivatives of F(x,y,z) with respect to x, y, z respectively.
        double F_x = 2 * generalCoefficients[0] * intersection_point.x + generalCoefficients[3] * intersection_point.y + generalCoefficients[4] * intersection_point.z + generalCoefficients[6];
        double F_y = 2 * generalCoefficients[1] * intersection_point.y + generalCoefficients[3] * intersection_point.x + generalCoefficients[5] * intersection_point.z + generalCoefficients[7];
        double F_z = 2 * generalCoefficients[2] * intersection_point.z + generalCoefficients[4] * intersection_point.x + generalCoefficients[5] * intersection_point.y + generalCoefficients[8];

        normal.x = F_x;
        normal.y = F_y;
        normal.z = F_z;
        normal = normal.normalize();

        //apply lighting model
        colorComputation(intersection_point,normal, r, level);

        return t;
    }
};


double calculateDet(double mat[3][3])
{
    double res;
    res = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) 
            - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) 
            + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    return res;
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

    void showObjectDetails()
    {
        cout << "Type: " << type << endl;
        cout << "Point 1: ";
        point1.print();
        cout << "Point 2: ";
        point2.print();
        cout << "Point 3: ";
        point3.print();
        cout << "Color: " << color[0] << " " << color[1] << " " << color[2] << endl;
        cout << "Shine: " << shine << endl;
        cout << "Ambient Coefficient: " << coEfficients[0] << " Diffuse Coefficient: " << coEfficients[1] << " Specular Coefficient: " << coEfficients[2] << " Reflection Coefficient: " << coEfficients[3] << endl;
    }

    double intersect(Ray *r, double *color, int level); //overriding the intersect method of the Object class
};

// Ray-triangle intersection
// Normal = cross product of two vectors along the edges e.g. (b-a)X(c-a)

// using Barycentric Coordinates
double Triangle::intersect(Ray *r, double *color, int level)
{
    // firstly construct tthe left   side matrix
    // formula is leftSideMatrix = [point1-point2, point1-point3, r.dir]
    double leftSideMatrix[3][3];  //A

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
    // if determinant is zero, then the matrix is singular and the ray is parallel to the plane of the triangle,so no intersection
    if (round(leftMatDet) == 0)
        return -1;

    double matBeta[3][3];
    // formula of matBeta is [point1-start, point1-point3, r.dir]

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
    // formulaa of matGamma is [point1-point2, point1-start, r.dir]
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
    // if gamma is less than 0 or greater than 1, then the intersection point is outside the triangle
    if (gamma <= 0)
        return -1;
    // another condition is beta+gamma should be less than 1
    if (beta + gamma >= 1)
        return -1;

    double tMatrix[3][3];
    // formula of tMatrix is [point1-point2, point1-point3, point1-start]
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
    /// If the recursion level is 0, the method returns the intersection parameter t without further calculations.
    if (level == 0)
        return t;

    Vector3D intersection_point, normal;
    // calculate the intersection point in 3D space using the parametric equation of the ray: intersection_point = r->start + t * r->dir
    // The ray is represented by the equation r(t) = r0 + t * d, where r0 is the starting point of the ray and d is the direction of the ray. (A point on the rayyy)

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

    // apply lighting model
    colorComputation(intersection_point,normal, r, level);

    return t;
}



class Sphere : public Object
{
public:
    Sphere(Vector3D center, double radius)
    {
        reference_point = center;
        length = radius;
    }

    void showObjectDetails()
    {
        cout << "Type: " << type << endl;
        cout << "Reference Point: ";
        reference_point.print();
        cout << "Radius: " << length << endl;
        cout << "Color: " << color[0] << " " << color[1] << " " << color[2] << endl;
        cout << "Shine: " << shine << endl;
        cout << "Ambient Coefficient: " << coEfficients[0] << " Diffuse Coefficient: " << coEfficients[1] << " Specular Coefficient: " << coEfficients[2] << " Reflection Coefficient: " << coEfficients[3] << endl;
    }


    void draw()
    {
        double radius = length;
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix();
        glutSolidSphere(radius, 24, 30);
        glPopMatrix();
    }
    
    
    
    double intersect(Ray *r, double *color, int level);  //overriding the intersect method of the Object class
};




double Sphere::intersect(Ray *r, double *clr, int level)
{
    // r is the castedRay initially
    //  start-center
    //Ray-Sphere intersection
    //H(P)=P^2-r^2=0
    //P(t)=R_o+tR_d
    //putting the value of P(t) in H(P) we get a quadratic equation in t
    Vector3D R_o; // origin of the ray
    R_o.x = r->start.x - reference_point.x;
    R_o.y = r->start.y - reference_point.y;
    R_o.z = r->start.z - reference_point.z;
    // at^2+bt+c=0
    // a=1,(because r.dir is normalized)
    // b=2R_d.R_o
    // c=R_o.R_o-r^2
    // algebraic solution
    double a = 1;
    double b = 2 * R_o.dot(r->dir);
    double c = R_o.dot(R_o) - length * length; // here length is the radius of the sphere

    double determinant = sqrt(b * b - 4 * a * c);
    double t1 = (-b - determinant) / (2 * a);
    double t2 = (-b + determinant) / (2 * a);
    double t;
    if (determinant < 0)
        return -1;
    else 
    ////take the smallest positive t
    {
        if (t1 > 0)
            t = t1;

        else if (t2 > 0)
            t = t2;

        else if(t1>0 && t2>0)
            t = min(t1, t2);

        else
            return -1;

    }
    if (level == 0)
        return t;

    setPhotoColor(getColor());

    Vector3D intersection_point, normal;
    // The ray is represented by the equation r(t) = r0 + t * d, where r0 is the starting point of the ray and d is the direction of the ray.
    intersection_point.x = r->start.x + t * r->dir.x;
    intersection_point.y = r->start.y + t * r->dir.y;
    intersection_point.z = r->start.z + t * r->dir.z;

    normal.x = intersection_point.x - reference_point.x;
    normal.y = intersection_point.y - reference_point.y;
    normal.z = intersection_point.z - reference_point.z;

    normal = normal.normalize();

    //apply lighting model
    colorComputation( intersection_point,normal, r, level);

    return t;
}

//////////////////////////////////////// Bonus Part //////////////////////////////////////////////////////  
void transposeMatrix(double** A, double** B, int r, int c)
{
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            B[j][i] = A[i][j];
}  

class Prism : public Object
{
private:
    Vector3D vertices[6]; // Six vertices for a prism
    double refractive_index[3]; // Refractive indices for the three colors
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

    void showObjectDetails()
    {
        cout << "Type: " << type << endl;
        cout << "Vertices: ";
        for (int i = 0; i < 6; i++)
        {
            vertices[i].print();
        }
        cout << "Color: " << color[0] << " " << color[1] << " " << color[2] << endl;
        cout << "Shine: " << shine << endl;
        cout << "Ambient Coefficient: " << coEfficients[0] << " Diffuse Coefficient: " << coEfficients[1] << " Specular Coefficient: " << coEfficients[2] << " Reflection Coefficient: " << coEfficients[3] << endl;
    }

    void setRefractiveIndex(double r, double g, double b)  //overide
    {
        refractive_index[0] = r;
        refractive_index[1] = g;
        refractive_index[2] = b;
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
        if (std::abs(denom) < SMALLDELTA)
            continue;

        double t = normal.dot(p1 - r->start) / denom;

        // Check if the intersection point is behind the ray's starting point
        if (t <= 0)
            continue;

        // Check if the intersection point is inside the triangular region
        Vector3D intersection_point = r->start + r->dir * t;
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
                colorComputation( intersection_point,normal, r, level);
            }

            // Set color and return tMin
            setPhotoColor(getColor());
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
        setType("floor");
        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;
        Vector3D *rp = new Vector3D(-floorWidth / 2, -floorWidth / 2, 0);
        reference_point = *rp;
        length = tileWidth;
        setCoEfficients(0.4, 0.2, 0.1, 0.3);
        shine = 5;
    }

    void showObjectDetails()
    {
        cout << "Type: " << type << endl;
        cout << "Floor Width: " << floorWidth << " Tile Width: " << tileWidth << endl;
        cout << "Shine: " << shine << endl;
        cout << "Ambient Coefficient: " << coEfficients[0] << " Diffuse Coefficient: " << coEfficients[1] << " Specular Coefficient: " << coEfficients[2] << " Reflection Coefficient: " << coEfficients[3] << endl;
    }

    void draw()
    {
        glPushMatrix();

        int rows = static_cast<int>(floorWidth / tileWidth);
        int columns = rows;

        double startX = reference_point.x;
        double startY = reference_point.y;

        for (int i = 0; i < rows; i++)
        {
            // startX is updated based on the current row.
            startX = i * tileWidth;
            for (int j = 0; j < columns; j++)
            {
                double color = ((i + j) % 2 == 0) ? 1.0 : 0.0;
                glColor3f(color, color, color);
                // startY is updated based on the current column.
                startY = j * tileWidth;

                glBegin(GL_QUADS);
                {
                    // Vertices of the quad
                    glVertex3f(startX, startY, 0);
                    glVertex3f(startX + tileWidth, startY, 0);
                    glVertex3f(startX + tileWidth, startY + tileWidth, 0);
                    glVertex3f(startX, startY + tileWidth, 0);
                }
                glEnd();
            }
        }

        glPopMatrix();
    }
    double intersect(Ray *r, double *objectColor, int level);
};

/**
 *
○ Ray-plane intersection followed by checking if the intersection point lies within
the span of the floor
○ Normal = (0,0,1)
○ Obtaining color will depend on which tile the intersection point lies on
*/

double Floor::intersect(Ray *r, double *color, int level)
{
    // P(t)=R0+t*Rd , where P(t) is the point on the ray(as well as plane), R0 is the starting point of the ray, Rd is the direction of the ray
    // H(P)=n.P+D=0, where H(P) is the plane equation, n is the normal to the plane, P is the point on the plane, D is the distance from the origin to the plane
    // n.P+D=0 => n.(R0+t*Rd)+D=0 => t=-(n.(R0)+D)/n.Rd

    // so equation for finding t is t=-(n.(R0)+D)/n.Rd

    Vector3D normal;
    normal.x = 0;
    normal.y = 0;
    normal.z = 1; // as the floor is on xy plane

    double denominator = r->dir.dot(normal);

    if (denominator == 0) // the ray is parallel to the plane, and no intersection occurs.
        return -1;
    Vector3D p0 = reference_point;
    double numerator = -(p0.dot(normal) + r->start.dot(normal)); 
    //
    double t = numerator / denominator;
    //
    Vector3D intersection_point;
    intersection_point.x = r->start.x + t * r->dir.x;
    intersection_point.y = r->start.y + t * r->dir.y;
    intersection_point.z = r->start.z + t * r->dir.z;

    if (intersection_point.x < -floorWidth / 2) // if the x-coordinate of the intersection point is to the left of the leftmost boundary of the floor.
        return -1;

    if (intersection_point.x > floorWidth / 2) // if the x-coordinate of the intersection point is to the right of the rightmost boundary of the floor.
        return -1;

    if (intersection_point.y < -floorWidth / 2) // if the y-coordinate of the intersection point is below the bottom boundary of the floor.
        return -1;

    if (intersection_point.y > floorWidth / 2) // if the y-coordinate of the intersection point is above the top boundary of the floor.
        return -1;
   
    int rowIdx = (intersection_point.x - reference_point.x) / tileWidth;
    int colIdx = (intersection_point.y - reference_point.y) / tileWidth;

    bool isEvenRow = (rowIdx % 2 == 0);
    bool isEvenCol = (colIdx % 2 == 0);

    if ((isEvenRow && isEvenCol) || (!isEvenRow && !isEvenCol))
    {
        setColor(1, 1, 1); // White
    }
    else
    {
        setColor(0, 0, 0); // Black
    }

    if (level == 0)
        return t;

    // Apply lighting model
    colorComputation( intersection_point,normal, r, level);

    return t;
}
