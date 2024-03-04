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
#define pi 3.14159



float eyex = 2, eyey = 2, eyez = 2;
float forwardX = -2, forwardY = -2, forwardZ = -2;  //forward direction.where camera is looking at
float upx = 0, upy = 1, upz = 0;          // Up vector
float rightx = 1, righty = 0, rightz = 0; // right direction
float movementFactor = 0.1;
float changeinRadius = 1.0f / (16.0f * sqrt(2)) ;// 16 means in 16 steps there will be a complete sphere or no sphere
float radius = 0.0f;
double revolve = 0.0f; // Rate of change of rotation


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
void drawPyramid() {
    glBegin(GL_TRIANGLES);

    // Front face
    glColor3f(1.0f, 0.0f, 0.0f);  // Red
    glVertex3f(0.0f, 1.0f, 0.0f);  // Top vertex
    glVertex3f(-1.0f, -1.0f, 1.0f);  // Bottom-left vertex
    glVertex3f(1.0f, -1.0f, 1.0f);  // Bottom-right vertex

    // Right face
    glColor3f(0.0f, 1.0f, 0.0f);  // Green
    glVertex3f(0.0f, 1.0f, 0.0f);  // Top vertex
    glVertex3f(1.0f, -1.0f, 1.0f);  // Bottom-left vertex
    glVertex3f(1.0f, -1.0f, -1.0f);  // Bottom-right vertex

    // Back face
    glColor3f(0.0f, 0.0f, 1.0f);  // Blue
    glVertex3f(0.0f, 1.0f, 0.0f);  // Top vertex
    glVertex3f(1.0f, -1.0f, -1.0f);  // Bottom-left vertex
    glVertex3f(-1.0f, -1.0f, -1.0f);  // Bottom-right vertex

    // Left face
    glColor3f(1.0f, 1.0f, 0.0f);  // Yellow
    glVertex3f(0.0f, 1.0f, 0.0f);  // Top vertex
    glVertex3f(-1.0f, -1.0f, -1.0f);  // Bottom-left vertex
    glVertex3f(-1.0f, -1.0f, 1.0f);  // Bottom-right vertex

    glEnd();
}


void drawCylinderSegmentWithQuads(float angle, float radius, float height)
{
    int segments = 20;
    glPushMatrix();
    glTranslatef(0, -height / 2, 0);
    glRotatef(angle / 2, 0.0f, 1.0f, 0.0f);
    
    float segmentAngle = angle * pi / 180.0f;

    glBegin(GL_QUADS);
    //with these 20 quads a cylinder segment is drawn
    for (int i = 0; i < segments; ++i)
    {
        float p1 = i * segmentAngle / segments; //angular span
        float p2 = (i + 1) * segmentAngle / segments;

        float x1 = radius * cos(p1);
        float z1 = radius * sin(p1);

        float x2 = radius * cos(p2);
        float z2 = radius * sin(p2);
       // draw a quad with 4 vertices
        glVertex3f(x1, 0.0f, z1);      // Bottom left vertex
        glVertex3f(x1, height, z1);    // Top left vertex
        glVertex3f(x2, height, z2);    // Top right vertex
        glVertex3f(x2, 0.0f, z2);      // Bottom right vertex
    }
    glEnd();
    glPopMatrix();
}

void drawAllYellowBoundaries()
{

   //need to draw total 12 cylinder like segments
   glColor3f(1, 1, 0);
   double angle = 70.53;
   double height = sqrt(2) - 2.0f * radius; //initially full height(equal to side of the triangle) as radius is 0.. when full radius is 1/sqrt(2) height is vanished. height is decreased from both side
   double tr = (1.0f - sqrt(2) * radius) / 2.0f; //to draw it at boundary of triangles

   //4 from top center of the octahedron

   for (int i = 0; i < 4; i++)
   {
      glPushMatrix();
      glRotatef(90 * i, 0, 1, 0);
      glTranslatef(tr, tr, 0);
      glRotatef(45, 0, 0, 1);
      drawCylinderSegmentWithQuads(angle, radius, height);
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
      drawCylinderSegmentWithQuads(angle, radius, height);
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
      drawCylinderSegmentWithQuads(angle, radius, height);
      glPopMatrix();


   }

}
void drawCircle(float centerX, float centerY, float radius, int segments) 
{
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(centerX, centerY); // Center of the circle
    for (int i = 0; i <= segments; ++i) {
        float theta = 2.0f * M_PI * static_cast<float>(i) / static_cast<float>(segments);
        float x = centerX + radius * std::cos(theta);
        float y = centerY + radius * std::sin(theta);
        glVertex2f(x, y);
    }
    glEnd();
}

void drawTriangle()
{
   GLfloat sc = 1.0f - sqrt(2) * radius;  //initially fullscale as radius is 0.. when full radius is 1/sqrt(2) trangle is vanished
   glPushMatrix();
   glTranslatef(radius / sqrt(3), radius / sqrt(3), radius / sqrt(3));    //translate to the center of the octahedron
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
   glPushMatrix();              
   glColor3f(0.0f, 1.0f, 1.0f); // skyblue
   drawTriangle();
   glPopMatrix();
   
   //the right one
   glPushMatrix(); 
   glColor3f(1.0f, 0.0f, 1.0f);//pink
   glRotatef(90, 0, 1, 0);
   drawTriangle();
   glPopMatrix();
   
   //the left one
   glPushMatrix(); 
   glColor3f(1.0f, 0.0f, 1.0f); //pink
   glRotatef(-90, 0, 1, 0);
   drawTriangle();
   glPopMatrix();
   
   //the back one
   glPushMatrix(); 
   glColor3f(0.0f, 1.0f, 1.0f); //skyblue
   glRotatef(180, 0, 1, 0);
   drawTriangle();
   glPopMatrix();

   //now just mirror wrt x axis

   glPushMatrix(); 
   glRotatef(180, 1, 0, 0); //for mirror
   glColor3f(0.0f, 1.0f, 1.0f); //skyblue
   drawTriangle();
   glPopMatrix();

   glPushMatrix(); 
   glColor3f(1.0f, 0.0f, 1.0f);
   glRotatef(180, 1, 0, 0); //for mirror
   glRotatef(90, 0, 1, 0);
   drawTriangle();
   glPopMatrix();

   glPushMatrix(); 
   glColor3f(1.0f, 0.0f, 1.0f);
   glRotatef(180, 1, 0, 0);
   glRotatef(-90, 0, 1, 0);
   drawTriangle();
   glPopMatrix();

   glPushMatrix(); 
   glColor3f(0.0f, 1.0f, 1.0f);
   glRotatef(180, 1, 0, 0);
   glRotatef(180, 0, 1, 0);
   drawTriangle();
   glPopMatrix();
}

void draw3DQuadrilateral(float x1, float y1, float z1,
                          float x2, float y2, float z2,
                          float x3, float y3, float z3,
                          float x4, float y4, float z4) 
    {
    glBegin(GL_QUADS);
    {
        glVertex3f(x1, y1, z1);
        glVertex3f(x2, y2, z2);
        glVertex3f(x3, y3, z3);
        glVertex3f(x4, y4, z4);
    }
    glEnd();
}

void keyboardListener(unsigned char key, int x, int y)
{

   switch (key)
   {
   // Control eye (location of the eye)
case '1':

        // change where it is looking now...means facing forward is changes
        forwardX = forwardX * cos(movementFactor) - rightx * sin(movementFactor);
        forwardY = forwardY * cos(movementFactor) - righty * sin(movementFactor);
        forwardZ = forwardZ * cos(movementFactor) - rightz * sin(movementFactor);

        break;
    case '2':

        // change where it is looking now...means facing forward is changes
        forwardX = forwardX * cos(-movementFactor) - rightx * sin(-movementFactor);
        forwardY = forwardY * cos(-movementFactor) - righty * sin(-movementFactor);
        forwardZ = forwardZ * cos(-movementFactor) - rightz * sin(-movementFactor);
        break;
    case '3':

        // center and up are perpendicular to each other

        forwardX = forwardX * cos(movementFactor) + upx * sin(movementFactor);
        forwardY = forwardY * cos(movementFactor) + upy * sin(movementFactor);
        forwardZ = forwardZ * cos(movementFactor) + upz * sin(movementFactor);

        break;
    case '4':
        forwardX = forwardX * cos(-movementFactor) + upx * sin(-movementFactor);
        forwardY = forwardY * cos(-movementFactor) + upy * sin(-movementFactor);
        forwardZ = forwardZ * cos(-movementFactor) + upz * sin(-movementFactor);
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

    case GLUT_KEY_UP:
        eyex += forwardX * movementFactor;
        eyey += forwardY * movementFactor;
        eyez += forwardZ * movementFactor;

        break;
    case GLUT_KEY_DOWN:
        eyex -= forwardX * movementFactor;
        eyey -= forwardY * movementFactor;
        eyez -= forwardZ * movementFactor;
        break;

    default:
        break;
    }
    glutPostRedisplay();
}

// from this link this portion is implemented https://www.songho.ca/opengl/gl_sphere.html#example_cubesphere
// generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
vector<vector<point>> generatePoints(int subdivision)
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
//The other 5 faces can be generated by repeating the above procedure, or swapping and/or negating axis of the vertices of +X face to optimize redundant sine/cosine computations.
// For example, the vertices of -X face are only negating x and z coordinates of +X face, and +Y face requires swapping x → y, y → -z, and z → -x.
void drawOneFace(int subdivision)
{  
    //generate those points on that surface
   vector<vector<point>> verSet = generatePoints(subdivision);

   // draw +X face of the cube using quadrilaterals
   glBegin(GL_QUADS);

   for (int i = 0; i < verSet.size() - 1; ++i)
   {
      for (int j = 0; j < verSet[i].size() - 1; ++j)
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
        forwardX+eyex, forwardY+eyey, forwardZ+eyez,
        upx, upy, upz);
    glPushMatrix();
    //in how many degree should it revolve
    glRotatef(revolve, 0, 1, 0);
   
    //glutSolidSphere(0.5f, 30, 30);
    drawSphere(8); //more subdivion more smooth
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

/**
 * g++ my_cube.cpp -o my_cube -lglut -lGLU -lGL
./my_cube
 * 
*/
