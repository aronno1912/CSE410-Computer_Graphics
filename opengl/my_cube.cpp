#define _USE_MATH_DEFINES
#include <bits/stdc++.h>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif

#define RADIAN M_PI / 180.0f


float eyex = 2, eyey = 2, eyez = 2;
float centerx = -2, centery = -2, centerz = -2;  //forward direction
float upx = 0, upy = 1, upz = 0;          // Up vector
float rightx = 1, righty = 0, rightz = 0; // right direction
float movementFactor = 0.1;
float changeinRadius = 1.0f / (16.0f * sqrt(2)) ;// 16 means in 16 steps there will be a complete sphere or no sphere



void init()
{
    // glClearColor(0.1f, .0f, 0.0f, 1.0f); // Set background color to black and opaque

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 1, 100);
}
struct point
{
   GLfloat x, y, z;
};


GLfloat radius = 0.0f;
double rate = 0.01f;   // Rate of change of movement
double revolve = 0.0f; // Rate of change of rotation

void drawCylinderSegment(float angle, float radius, float height)
{
   int numSegments = 100;
   float segmentAngle = angle * RADIAN;

   glPushMatrix();
   glTranslatef(0, -height / 2, 0); // Translate to the base of the cylinder.This translation moves the base of the cylinder to the origin (0, 0, 0) in the local coordinate system
   glRotatef(angle / 2, 0.0f, 1.0f, 0.0f); // Rotate around the y-axis

   glBegin(GL_TRIANGLE_STRIP); //an efficient way to render a cylindrical shape.
   for (int i = 0; i <= numSegments; ++i)
   {
      float theta = i * segmentAngle / numSegments;
      float x = radius * cos(theta);
      float z = radius * sin(theta);

      glVertex3f(x, 0.0f, z);  // Vertex on the bottom circle
      glVertex3f(x, height, z); // Corresponding vertex on the top circle
   }
   glEnd();
   glPopMatrix();
}


void drawAllYellowBoundaries()
{

   //need to draw total 12 cylinder like segments
   double angle = 70.53;
   double height = sqrt(2) - 2.0f * radius; //initially full height(equal to side of the triangle) as radius is 0.. when full radius is 1/sqrt(2) height is vanished. height is decreased from both side
   double tr = (1.0f - sqrt(2) * radius) / 2.0f;  // ????

   glColor3f(1, 1, 0);

   //4 from top center of the octahedron

   for (int i = 0; i < 4; i++)
   {
      glPushMatrix();
      glRotatef(90 * i, 0, 1, 0);
      glTranslatef(tr, tr, 0);
      glRotatef(45, 0, 0, 1);
      drawCylinderSegment(angle, radius, height);
      glPopMatrix();
   }

   // //4 from bottom center of the octahedron
   for (int i = 0; i < 4; i++)
   {
      glPushMatrix();
      glRotatef(180, 1, 0, 0);
      glRotatef(90 * i, 0, 1, 0);
      glTranslatef(tr, tr, 0);
      glRotatef(45, 0, 0, 1);
      drawCylinderSegment(angle, radius, height);
      glPopMatrix();


   }


   // // 4 for the base of the octahedron
   for (int i = 0; i < 4; i++)
   {
      glPushMatrix();
      glRotatef(90 * i, 0, 1, 0);
      glRotatef(90, 1, 0, 0);
      glTranslatef(tr, tr, 0);
      glRotatef(45, 0, 0, 1);
      drawCylinderSegment(angle, radius, height);
      glPopMatrix();
   }


}

void drawTriangle()
{
   GLfloat sc = 1.0f - sqrt(2) * radius;  //initially fullscale as radius is 0.. when full radius is 1/sqrt(2) trangle is vanished
   glPushMatrix();
   glTranslatef(radius / sqrt(3), radius / sqrt(3), radius / sqrt(3));    ////////////////////////////////////////////////////////
   glScalef(sc, sc, sc);
   glBegin(GL_TRIANGLES);
   //triangle with each side root(2) unit
   glVertex3f(0,1,0);
   glVertex3f(0,0,1);
   glVertex3f(1,0,0);
   glEnd(); 
   glPopMatrix();
}

void drawOctahedron()
{
   glPushMatrix();              // Create a new scope
   glColor3f(0.0f, 1.0f, 1.0f); // skyblue
   drawTriangle();
   glPopMatrix();
   
   //the right one
   glPushMatrix(); // Create a new scope
   glColor3f(1.0f, 0.0f, 1.0f);//pink
   glRotatef(90, 0, 1, 0);
   drawTriangle();
   glPopMatrix();
   
   //the left one
   glPushMatrix(); // Create a new scope
   glColor3f(1.0f, 0.0f, 1.0f); //pink
   glRotatef(-90, 0, 1, 0);
   drawTriangle();
   glPopMatrix();
   
   //the back one
   glPushMatrix(); // Create a new scope
   glColor3f(0.0f, 1.0f, 1.0f); //skyblue
   glRotatef(180, 0, 1, 0);
   drawTriangle();
   glPopMatrix();

   //now just mirror wrt x axis

   glPushMatrix(); // Create a new scope
   glRotatef(180, 1, 0, 0); //for mirror
   glColor3f(0.0f, 1.0f, 1.0f); //skyblue
   drawTriangle();
   glPopMatrix();

   glPushMatrix(); // Create a new scope
   glColor3f(1.0f, 0.0f, 1.0f);
   glRotatef(180, 1, 0, 0); //for mirror
   glRotatef(90, 0, 1, 0);
   drawTriangle();
   glPopMatrix();

   glPushMatrix(); // Create a new scope
   glColor3f(1.0f, 0.0f, 1.0f);
   glRotatef(180, 1, 0, 0);
   glRotatef(-90, 0, 1, 0);
   drawTriangle();
   glPopMatrix();

   glPushMatrix(); // Create a new scope
   glColor3f(0.0f, 1.0f, 1.0f);
   glRotatef(180, 1, 0, 0);
   glRotatef(180, 0, 1, 0);
   drawTriangle();
   glPopMatrix();
}


void keyboardListener(unsigned char key, int x, int y)
{

   switch (key)
   {
   // Control eye (location of the eye)
case '1':

        // change where it is looking now...means facing forward is changes
        centerx = centerx * cos(movementFactor) - rightx * sin(movementFactor);
        centery = centery * cos(movementFactor) - righty * sin(movementFactor);
        centerz = centerz * cos(movementFactor) - rightz * sin(movementFactor);

        break;
    case '2':

        // change where it is looking now...means facing forward is changes
        centerx = centerx * cos(-movementFactor) - rightx * sin(-movementFactor);
        centery = centery * cos(-movementFactor) - righty * sin(-movementFactor);
        centerz = centerz * cos(-movementFactor) - rightz * sin(-movementFactor);
        break;
    case '3':

        // center and up are perpendicular to each other

        centerx = centerx * cos(movementFactor) + upx * sin(movementFactor);
        centery = centery * cos(movementFactor) + upy * sin(movementFactor);
        centerz = centerz * cos(movementFactor) + upz * sin(movementFactor);

        break;
    case '4':
        centerx = centerx * cos(-movementFactor) + upx * sin(-movementFactor);
        centery = centery * cos(-movementFactor) + upy * sin(-movementFactor);
        centerz = centerz * cos(-movementFactor) + upz * sin(-movementFactor);
        break;

    case '5':
        // update the up vector by rotating it around the right vector

        upx = upx * cos(movementFactor) + rightx * sin(movementFactor);
        upy = upy * cos(movementFactor) + righty * sin(movementFactor);
        upz = upz * cos(movementFactor) + rightz * sin(movementFactor);

        //now adjust right vector accordingly
        rightx = rightx * cos(movementFactor) - upx * sin(movementFactor);
        righty = righty * cos(movementFactor) - upy * sin(movementFactor);
        rightz = rightz * cos(movementFactor) - upz * sin(movementFactor);


        break;

    case '6':
        upx = upx * cos(-movementFactor) + rightx * sin(-movementFactor);
        upy = upy * cos(-movementFactor) + righty * sin(-movementFactor);
        upz = upz * cos(-movementFactor) + rightz * sin(-movementFactor);

        // now adjust right vector accordingly
        rightx = rightx * cos(-movementFactor) - upx * sin(-movementFactor);
        righty = righty * cos(-movementFactor) - upy * sin(-movementFactor);
        rightz = rightz * cos(-movementFactor) - upz * sin(-movementFactor);

        break;


   case ',':
      // cout << "radius: " << radius << endl;
       //circle inscribed in a square whose each side of the sqare is root(2) unit,so radius is root(2)/2=1/sqrt(2)
      radius = radius + changeinRadius;
      if (radius >= 1/sqrt(2))
      {
         radius =1/sqrt(2);
      }

      break;


   case '.':
      radius = radius - changeinRadius;
      if (radius < 0)
      {
         radius = 0;
      }
      break;

    
   case 'a':
      revolve = revolve - 1;
      break;
   case 'd':
      revolve = revolve + 1;
      break;
   }
   glutPostRedisplay(); // Post a paint request to activate display()
}


/* Callback handler for special-key event */
void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_UP:
        eyex += centerx * movementFactor;
        eyey += centery * movementFactor;
        eyez += centerz * movementFactor;

        break;
    case GLUT_KEY_DOWN:
        eyex -= centerx * movementFactor;
        eyey -= centery * movementFactor;
        eyez -= centerz * movementFactor;
        break;

    case GLUT_KEY_RIGHT:
        eyex += rightx * movementFactor;
        eyey += righty * movementFactor;
        eyez += rightz * movementFactor;
        break;
    case GLUT_KEY_LEFT:
        eyex -= rightx * movementFactor;
        eyey -= righty * movementFactor;
        eyez -= rightz * movementFactor;

        break;

    case GLUT_KEY_PAGE_UP:
        eyex += upx * movementFactor;
        eyey += upy * movementFactor;
        eyez += upz * movementFactor;
        break;

    case GLUT_KEY_PAGE_DOWN:
        eyex -= upx * movementFactor;
        eyey -= upy * movementFactor;
        eyez -= upz * movementFactor;
        break;

    default:
        break;
    }
    glutPostRedisplay();
}

// from this link this portion is implemented https://www.songho.ca/opengl/gl_sphere.html#example_cubesphere
vector<vector<point>> buildUnitPositiveX(int subdivision)
{
   const float DEG2RAD = acos(-1) / 180.0f;


   float n1[3]; // normal of longitudinal plane rotating along Y-axis
   float n2[3]; // normal of latitudinal plane rotating along Z-axis
   float v[3];  // direction vector intersecting 2 planes, n1 x n2
   float a1;    // longitudinal angle along Y-axis
   float a2;    // latitudinal angle along Z-axis

   // compute the number of vertices per row, 2^n + 1
   int pointsPerRow = (int)pow(2, subdivision) + 1;

   vector<vector<point>> vertices(pointsPerRow);

   // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
   for (unsigned int i = 0; i < pointsPerRow; ++i)
   {
      // normal for latitudinal plane
      // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
      // therefore, it is rotating (0,1,0) vector by latitude angle a2
      a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
      n2[0] = -sin(a2);
      n2[1] = cos(a2);
      n2[2] = 0;

      // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
      for (unsigned int j = 0; j < pointsPerRow; ++j)
      {
         // normal for longitudinal plane
         // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
         // therefore, it is rotating (0,0,-1) vector by longitude angle a1
         a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
         n1[0] = -sin(a1);
         n1[1] = 0;
         n1[2] = -cos(a1);

         // find direction vector of intersected line, n1 x n2
         v[0] = n1[1] * n2[2] - n1[2] * n2[1];
         v[1] = n1[2] * n2[0] - n1[0] * n2[2];
         v[2] = n1[0] * n2[1] - n1[1] * n2[0];

         // normalize direction vector
         float scale = 1 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
         v[0] *= scale;
         v[1] *= scale;
         v[2] *= scale;

         // add a vertex into array
         vertices[i].push_back({v[0], v[1], v[2]});
      }
   }

   return vertices;
}

//build quads with generated points

void drawOneFace(int subdivision)
{  
    //generate those points on that surface
   vector<vector<point>> verSet = buildUnitPositiveX(subdivision);

   // draw +X face of the cube using quadrilaterals
   glBegin(GL_QUADS);
   for (unsigned int i = 0; i < verSet.size() - 1; ++i)
   {
      for (unsigned int j = 0; j < verSet[i].size() - 1; ++j)
      {
         // draw a quadrilateral using 4 vertices of 2 adjacent rows
         glVertex3f(verSet[i][j].x, verSet[i][j].y, verSet[i][j].z);
         glVertex3f(verSet[i][j + 1].x, verSet[i][j + 1].y, verSet[i][j + 1].z);
         glVertex3f(verSet[i + 1][j + 1].x, verSet[i + 1][j + 1].y, verSet[i + 1][j + 1].z);
         glVertex3f(verSet[i + 1][j].x, verSet[i + 1][j].y, verSet[i + 1][j].z);
      }
   }
   glEnd();
}

void drawSphere(int subdivision)
{
   glPushMatrix(); // Create a new scope     for +X face
   glColor3f(1, 0, 0);
   //when full radius,translate to the position of +X triangle which was at (1,0,0)
   glTranslatef((1 - sqrt(2) * radius), 0, 0);
   glScalef(radius, radius, radius);
   drawOneFace(subdivision);
   glPopMatrix();


   glPushMatrix(); // Create a new scope     for +Y face
   glColor3f(0, 0, 1);
   glTranslatef(0, (1 - sqrt(2) * radius), 0);
   glRotatef(90, 0, 0, 1);
   glScalef(radius, radius, radius);
   drawOneFace(subdivision);
   glPopMatrix();



   glPushMatrix(); // Create a new scope     for +Z face
   glColor3f(0, 1, 0);
   glTranslatef(0, 0, (1 - sqrt(2) * radius));
   glRotatef(270, 0, 1, 0);
   glScalef(radius, radius, radius);
   drawOneFace(subdivision);
   glPopMatrix();


   glPushMatrix(); // Create a new scope     for -X face
   glColor3f(1, 0, 0);
   glTranslatef(-(1 - sqrt(2) * radius), 0, 0);
   glRotatef(180, 0, 1, 0); //mirror
   glScalef(radius, radius, radius);
   drawOneFace(subdivision);
   glPopMatrix();



   glPushMatrix(); // Create a new scope     for -Y face
   glColor3f(0, 0, 1);
   glTranslatef(0, -(1 - sqrt(2) * radius), 0);
   glRotatef(-90, 0, 0, 1);
   glScalef(radius, radius, radius);
   drawOneFace(subdivision);
   glPopMatrix();


   glPushMatrix(); // Create a new scope     for -Z face
   glColor3f(0, 1, 0);
   glTranslatef(0, 0, -(1 - sqrt(2) * radius));
   glRotatef(90, 0, 1, 0); // Rotate 90 degree around Y-axis
   glScalef(radius, radius, radius);
   drawOneFace(subdivision);
   glPopMatrix();





}

/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes()
{
  glLineWidth(3);
   glBegin(GL_LINES);
   glColor3f(1, 0, 0); // Red
   // X axis
   glVertex3f(0, 0, 0);
   glVertex3f(1.5, 0, 0);

   glColor3f(0, 1, 0); // Green
   // Y axis
   glVertex3f(0, 0, 0);
   glVertex3f(0, 1.5, 0);

   glColor3f(0, 0, 1); // Blue
   // Z axis
   glVertex3f(0, 0, 0);
   glVertex3f(0, 0, 1.5);
   glEnd();
}

void drawOneFace(float r, float g, float b, float x, float y, float z) {
    glColor3f(r, g, b);
    glPushMatrix();
    glTranslatef(x, y, z);
    glutSolidSphere(radius, 30, 30);
    glPopMatrix();
}

void drawSolidSphere(int subdivision)
{

    drawOneFace(1.0, 0.0, 0.0, 1 - sqrt(2) * radius, 0, 0);   // +X face 
    drawOneFace(0.0, 1.0, 0.0, 0, 0, -(1 - sqrt(2) * radius)); // -Z face 
    drawOneFace(1.0, 0.0, 0.0, -(1 - sqrt(2) * radius), 0, 0);  // -X face 
    drawOneFace(0.0, 1.0, 0.0, 0, 0, 1 - sqrt(2) * radius);    // +Z face 
    drawOneFace(0.0, 0.0, 1.0, 0, 1 - sqrt(2) * radius, 0);    // +Y face 
    drawOneFace(0.0, 0.0, 1.0, 0, -(1 - sqrt(2) * radius), 0); // -Y face 
   
}

void display()
{
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // look at e eye use hoy...eye.x,eye.y,eye.z ,center o
    gluLookAt(
        eyex, eyey, eyez,
        centerx+eyex, centery+eyey, centerz+eyez,
        upx, upy, upz);
    glPushMatrix();
    glRotatef(revolve, 0, 1, 0);
    //drawAxes();
    //glutSolidSphere(0.5f, 30, 30);
    drawSphere(6); //more subdivion more smooth
    drawOctahedron();
    drawAllYellowBoundaries();
    glPopMatrix();
    
    glutSwapBuffers(); // If you are using double buffering
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Hello");
    init();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboardListener);                       // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);                      // Register callback handler for special-key event
    
    glEnable(GL_DEPTH_TEST);
    glutMainLoop();
    return 0;
}