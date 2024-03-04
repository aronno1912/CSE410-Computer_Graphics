#include <bits/stdc++.h>
using namespace std;

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif

#define pi 3.14159
float eyex = 7, eyey = 7, eyez = 7;
float forwardX = -2, forwardY = -2, forwardZ = -2;
float upx = 0, upy = 1, upz = 0;          // Up vector
float rightx = 1, righty = 0, rightz = 0; // right directi

double rotation=0;
void init()
{
    //glClearColor(0.1f, .0f, 0.0f, 1.0f); // Set background color to black and opaque

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 1, 100);
}


void axes()
{
    glLineWidth(3);
    glBegin(GL_LINES);
    {
        glColor3f(1.0, 0.0, 0.0); // x red
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);
        glColor3f(0.0, 1.0, 0.0); // y green
        glVertex3f(0, -100, 0);
        glVertex3f(0, 100, 0);
        glColor3f(0.0, 0.0, 1.0); // z blue
        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    }
    glEnd();
}


void drawCircle(float radius) {
    glBegin(GL_TRIANGLE_FAN);
    for (int i = 0; i <= 360; i++) {
        float angle = i * 3.14159265359 / 180.0;
        glVertex3f(radius * cos(angle),0, radius * sin(angle));
    }
    glEnd();
}



void drawRectangle(float width, float height) {
    glBegin(GL_QUADS);
    glVertex3f(-width / 2,0, -height / 2);
    glVertex3f(width / 2, 0,-height / 2);
    glVertex3f(width / 2, 0,height / 2);
    glVertex3f(-width / 2, 0,height / 2);
    glEnd();
}


void drawSwing()

{
    //basically a rectangle tied with two straight line hanging from the top circle
 
    //draw the rectangle
    glColor3f(1, 0, 0);
    glPushMatrix();
    glTranslatef(0, -2, 0);
    drawRectangle(0.5, 0.5);
    glPopMatrix();

    //draw two lines from the two corners of the rectangle to tie it with the top of circle
    glColor3f(1, 1, 1);
    glBegin(GL_LINES);
    glVertex3f(0.25, -1.75, 0);
    glVertex3f(0, 0, 0);
    glEnd();

    glBegin(GL_LINES);
    glVertex3f(-0.25, -1.75, 0);
    glVertex3f(0, 0, 0);
    glEnd();




}

void drawMerryGoRound(int numHorses) {
    // Draw central pole
    glColor3f(0.5, 0.5, 0.5); 
    drawRectangle(0.05, 1.0);

    // Draw circular platform
    glColor3f(0.7, 0.7, 0.7); 
    drawCircle(1.0);

    //draw 5 swings with translation and rotations
    for (int i = 0; i < numHorses; i++) {
        float angle = i * (360.0 / numHorses);
        glPushMatrix();
        glRotatef(angle, 0, 1, 0);
        glTranslatef(1, 0, 0);
        drawSwing();
        glPopMatrix();
    }

   //revolve the whole merry go round
  
  glRotatef(rotation, 0, 1, 0);
    


    
    
     
    
    // glColor3f(1, 1, 1); 
    // for (int i = 0; i < numHorses; i++) {
    //     float angle = i * (360.0 / numHorses);
    //     float x = 1.0 * cos(angle * 3.14159265359 / 180.0);
    //     float y = 1.0 * sin(angle * 3.14159265359 / 180.0);
        
    //     glPushMatrix();
    //     glTranslatef(x, y-2, 0.0);
    //     drawRectangle(0.5, 0.5); 
    //     glPopMatrix();
    // }
}
void Timer(int value)
{
    //it will revolve the marry go round
     glutPostRedisplay();
    glutTimerFunc(80,Timer,0);
    rotation+=0.5;
 
    
    //    glutPostRedisplay();
    // glutTimerFunc(1000 / 60, Timer, 0);

}

void display()
{
     glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);// 3d
    glLoadIdentity(); //3d

    // look at e eye use hoy...eye.x,eye.y,eye.z ,center o
    //3d
    gluLookAt(
        eyex, eyey, eyez,
        forwardX+eyex, forwardY+eyey, forwardZ+eyez,
        upx, upy, upz);
    glPushMatrix();
    axes();
    drawMerryGoRound(5);
    glFlush();
    glPopMatrix();
    glutSwapBuffers();

}



int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Practice 1");
    init();
    glutDisplayFunc(display);
    glutTimerFunc(0,Timer,0);
    //glutIdleFunc(doYourOwnMove);
    // glutTimerFunc(0, eventDrivenSchedule, 0);
    // glutSpecialFunc(specialKeyListener);
    // glutKeyboardFunc(keyboardListener);
    glEnable(GL_DEPTH_TEST);
    glutMainLoop();

    return 0;
}





