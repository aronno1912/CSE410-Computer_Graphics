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
#define SCALING_FACTOR 1.0f / (16.0f * sqrt(2))

void init(){
    printf("Do your initialization here\n");
    //drawaxes = 1;
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
    
    // glMatrixMode(GL_PROJECTION);
    // glLoadIdentity();

   glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
   glLoadIdentity();            // Reset the projection matrix
   gluPerspective(45.0f, 1, 0.1f, 100.0f);
    
   glEnable(GL_DEPTH_TEST);  
    

}
struct point
{
   GLfloat x, y, z;
};
// struct point eye = {2, 2, 2};       // Camera position
struct point eye = {4, 4, 4};       // Camera position
struct point center = {-2, -2, -2}; // Look-at point object
//struct point center = {0, 0, 0}; // Look-at point
struct point up = {0, 1, 0};        // Up vector
struct point r = {1, 0, 0};         // right direction

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


void drawCylinder()
{
   double angle = 70.53;
   double height = sqrt(2) - 2.0f * radius; //????????????????????????????????????????????????????////
   double translate = (1.0f - sqrt(2) * radius) / 2.0f;    //????????????????????????????????????????????????/

   glColor3f(1, 1, 0);
   //3 for each axis

   for (int i = 0; i < 4; i++)
   {
      glPushMatrix();
      glRotatef(90 * i, 0, 1, 0);
      glTranslatef(translate, translate, 0);
      glRotatef(45, 0, 0, 1);

      drawCylinderSegment(angle, radius, height);
      glPopMatrix();
   }

   for (int i = 0; i < 4; i++)
   {
      glPushMatrix();
      glRotatef(180, 1, 0, 0);
      glRotatef(90 * i, 0, 1, 0);
      glTranslatef(translate, translate, 0);
      glRotatef(45, 0, 0, 1);

      drawCylinderSegment(angle, radius, height);
      glPopMatrix();
   }

   for (int i = 0; i < 4; i++)
   {
      glPushMatrix();
      glRotatef(90 * i, 0, 1, 0);
      glRotatef(90, 1, 0, 0);
      glTranslatef(translate, translate, 0);
      glRotatef(45, 0, 0, 1);

      drawCylinderSegment(angle, radius, height);
      glPopMatrix();
   }


}

void drawTriangle()
{
   GLfloat factor = 1.0f - sqrt(2) * radius;  //??????????????????????????????
   glPushMatrix();
   glTranslatef(radius / sqrt(3), radius / sqrt(3), radius / sqrt(3));   // ???????????????????????????????
   glScalef(factor, factor, factor);
   glBegin(GL_TRIANGLES);
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
      r.x = r.x * cos(rate) + center.x * sin(rate);
      r.y = r.y * cos(rate) + center.y * sin(rate);
      r.z = r.z * cos(rate) + center.z * sin(rate);

      center.x = center.x * cos(rate) - r.x * sin(rate);
      center.y = center.y * cos(rate) - r.y * sin(rate);
      center.z = center.z * cos(rate) - r.z * sin(rate);
      break;

   case '2':
      r.x = r.x * cos(-rate) + center.x * sin(-rate);
      r.y = r.y * cos(-rate) + center.y * sin(-rate);
      r.z = r.z * cos(-rate) + center.z * sin(-rate);

      center.x = center.x * cos(-rate) - r.x * sin(-rate);
      center.y = center.y * cos(-rate) - r.y * sin(-rate);
      center.z = center.z * cos(-rate) - r.z * sin(-rate);
      break;

   case '3':
      center.x = center.x * cos(rate) + up.x * sin(rate);
      center.y = center.y * cos(rate) + up.y * sin(rate);
      center.z = center.z * cos(rate) + up.z * sin(rate);

      up.x = up.x * cos(rate) - center.x * sin(rate);
      up.y = up.y * cos(rate) - center.y * sin(rate);
      up.z = up.z * cos(rate) - center.z * sin(rate);
      break;

   case '4':
      center.x = center.x * cos(-rate) + up.x * sin(-rate);
      center.y = center.y * cos(-rate) + up.y * sin(-rate);
      center.z = center.z * cos(-rate) + up.z * sin(-rate);

      up.x = up.x * cos(-rate) - center.x * sin(-rate);
      up.y = up.y * cos(-rate) - center.y * sin(-rate);
      up.z = up.z * cos(-rate) - center.z * sin(-rate);
      break;

   case '5':
      up.x = up.x * cos(rate) + r.x * sin(rate);
      up.y = up.y * cos(rate) + r.y * sin(rate);
      up.z = up.z * cos(rate) + r.z * sin(rate);

      r.x = r.x * cos(rate) - up.x * sin(rate);
      r.y = r.y * cos(rate) - up.y * sin(rate);
      r.z = r.z * cos(rate) - up.z * sin(rate);
      break;

   case '6':
      up.x = up.x * cos(-rate) + r.x * sin(-rate);
      up.y = up.y * cos(-rate) + r.y * sin(-rate);
      up.z = up.z * cos(-rate) + r.z * sin(-rate);

      r.x = r.x * cos(-rate) - up.x * sin(-rate);
      r.y = r.y * cos(-rate) - up.y * sin(-rate);
      r.z = r.z * cos(-rate) - up.z * sin(-rate);
      break;

   case ',':
      // cout << "radius: " << radius << endl;

      radius = radius + SCALING_FACTOR;
      if (radius >= 1.0f / sqrt(2))
      {
         radius = 1.0f / sqrt(2);
      }

      // cout << "After radius: " << radius << endl;
      break;


   case '.':
      radius = radius - SCALING_FACTOR;
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
   case GLUT_KEY_UP: // down arrow key
      eye.x += center.x * rate;
      eye.y += center.y * rate;
      eye.z += center.z * rate;
      break;
   case GLUT_KEY_DOWN: // up arrow key
      eye.x -= center.x * rate;
      eye.y -= center.y * rate;
      eye.z -= center.z * rate;
      break;

   case GLUT_KEY_RIGHT:
      eye.x += r.x * rate;
      eye.y += r.y * rate;
      eye.z += r.z * rate;
      break;
   case GLUT_KEY_LEFT:
      eye.x -= r.x * rate;
      eye.y -= r.y * rate;
      eye.z -= r.z * rate;
      break;

   case GLUT_KEY_PAGE_UP:

      eye.x -= up.x * rate;
      eye.y -= up.y * rate;
      eye.z -= up.z * rate;
      break;

   case GLUT_KEY_PAGE_DOWN:
      eye.x += up.x * rate;
      eye.y += up.y * rate;
      eye.z += up.z * rate;
      break;

   case GLUT_KEY_INSERT:
      break;

   case GLUT_KEY_HOME:
      break;
   case GLUT_KEY_END:
      break;

   default:
      break;
   }
   glutPostRedisplay();
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

void drawSphere(int subdivision)
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
    // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // glMatrixMode(GL_MODELVIEW);
    // glLoadIdentity();
    // gluLookAt(100,100,100,	0,0,0,	0,0,1);
    // drawAxes();
    // glutSwapBuffers();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    gluLookAt(eye.x, eye.y, eye.z,
             center.x+eye.x,center.y+eye.y, center.z+eye.z,
             up.x, up.y, up.z);
    glPushMatrix();
    glRotatef(revolve, 0, 1, 0);
    //drawAxes();
    //glutSolidSphere(0.5f, 30, 30);
    drawSphere(30);
    drawOctahedron();
    drawCylinder();
    glPopMatrix();
    
    glutSwapBuffers(); // If you are using double buffering
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowPosition(50, 50);
    glutInitWindowSize(640, 640);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Hello");
    glutDisplayFunc(display);
    //glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);                       // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);                      // Register callback handler for special-key event
    init();
    glutMainLoop();
    return 0;
}