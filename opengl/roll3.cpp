#include <bits/stdc++.h>
using namespace std;

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif

#define pi 3.14159
// float eyex = 10, eyey = 10, eyez = 10;
float eyex = 7, eyey = 7, eyez = 7;
float centerx = -2, centery = -2, centerz = -2;
float upx = 0, upy = 1, upz = 0;          // Up vector
float rightx = 1, righty = 0, rightz = 0; // right direction
float movementFactor = 0.1;
bool autoMove = false; // Whether the ball should move automatically or not
float stacks = 10, sectors = 10;
float ballLocation[3] = {0, 0.5, 0};
float ballDirection[3] = {0, 0, 0};
float radius = 0.5;
float forwardAngle = 0;
float ballRotation = 0;

float directionX = 0.0f; // Initial direction along x-axis
float directionY = 0.0f; // Initial direction along y-axis

bool hit = false;
int currDir = 45;

struct point
{
    GLfloat x, y, z;
};

void init()
{
    // glClearColor(0.1f, .0f, 0.0f, 1.0f); // Set background color to black and opaque

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

    //  glBegin(GL_LINES);
    // {
    //     glVertex3f(ballLocation[0], ballLocation[1], ballLocation[2]);
    //     glVertex3f(ballLocation[0]+cos(currDir * pi / (double)180)*2, 0,ballLocation[2]+sin(currDir * pi / (double)180)*2 );
        
    //     // glVertex3f(0,0,0);
    //     // glVertex3f(10,0,0);
    // }
    // glEnd();
    // Draw the arrowhead
// glBegin(GL_TRIANGLES);
// {
//     // Calculate the tip of the arrowhead
//     float arrowTipX = ballLocation[0] + cos(currDir * pi / 180) * 2;
//     float arrowTipY = 0;
//     float arrowTipZ = ballLocation[2] + sin(currDir * pi / 180) * 2;

//     // Calculate two points for the base of the arrowhead
//     float base1X = arrowTipX - 0.2 * cos((currDir + 150) * pi / 180);
//     float base1Y = arrowTipY;
//     float base1Z = arrowTipZ - 0.2 * sin((currDir + 150) * pi / 180);

//     float base2X = arrowTipX - 0.2 * cos((currDir - 150) * pi / 180);
//     float base2Y = arrowTipY;
//     float base2Z = arrowTipZ - 0.2 * sin((currDir - 150) * pi / 180);

//     // Draw the triangle
//     glVertex3f(arrowTipX, arrowTipY, arrowTipZ);
//     glVertex3f(base1X, base1Y, base1Z);
//     glVertex3f(base2X, base2Y, base2Z);
// }
// glEnd();
}

void square(double a)
{

    glBegin(GL_QUADS);
    {
        glVertex3f(a, a, 0);
        glVertex3f(a, -a, 0);
        glVertex3f(-a, -a, 0);
        glVertex3f(-a, a, 0);
    }
    glEnd();
}

void drawCheckerboardXZ(int size, float squareSize)
{
    for (int i = -size; i < size; ++i)
    {
        for (int j = -size; j < size; ++j)
        {
            // Calculate coordinates for each square on the XZ plane
            float x = i * squareSize;
            float z = j * squareSize;

            // Alternate between black and white squares
            if ((i + j) % 2 == 0)
            {
                // Black square
                glColor3f(0.0f, 0.0f, 0.0f);
            }
            else
            {
                // White square
                glColor3f(1.0f, 1.0f, 1.0f);
            }

            // Draw the square on the XZ plane
            glBegin(GL_QUADS);
            glVertex3f(x, 0.0f, z);
            glVertex3f(x + squareSize, 0.0f, z);
            glVertex3f(x + squareSize, 0.0f, z + squareSize);
            glVertex3f(x, 0.0f, z + squareSize);
            glEnd();
        }
    }
}

void oneBorder(int length, int width)
{
    glVertex3f(0.0f, width / 2, length / 2);
    glVertex3f(0.0f, -width / 2, length / 2);
    glVertex3f(0.0f, -width / 2, -length / 2);
    glVertex3f(0.0f, width / 2, -length / 2);
}
void drawBoundary(float length, float width)
{
    glPushMatrix();
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_QUADS);
    glVertex3f(0.0f, width / 2, length / 2);
    glVertex3f(0.0f, -width / 2, length / 2);
    glVertex3f(0.0f, -width / 2, -length / 2);
    glVertex3f(0.0f, width / 2, -length / 2);
    glEnd();
    glPopMatrix();
}

void drawSegmentedSphere()
{

    float stackAngle = pi / stacks, phi = -pi / 2;
    float sectorAngle = 2 * pi / sectors, theta;

    glPushMatrix();
    glTranslatef(ballLocation[0], ballLocation[1], ballLocation[2]);
    glRotatef(forwardAngle, 0, 0, 1);

    // glRotatef(ballRotation, 0, 0, 1);
    for (int i = 0; i <= stacks; i++)
    {
        theta = 0;
        glBegin(GL_QUAD_STRIP);
        {
            for (int j = 0; j <= sectors; j++)
            {
                if ((i - 1) * 2 >= stacks)
                    // if ((i + j) % 2 == 0)
                    glColor3f(j % 2, (j + 1) % 2, 0);
                // glColor3f(1, 0, 0);
                else
                    glColor3f((j + 1) % 2, j % 2, 0);
                // glColor3f(0, 1, 0);
                glVertex3f(radius * cos(phi - stackAngle) * cos(theta), radius * cos(phi - stackAngle) * sin(theta), radius * sin(phi - stackAngle));
                glVertex3f(radius * cos(phi) * cos(theta), radius * cos(phi) * sin(theta), radius * sin(phi));
                theta += sectorAngle;
            }
        }
        glEnd();
        phi += stackAngle;
    }

    // draw arrow
    glColor3f(0, 0, 1);
    glPushMatrix();
    // glTranslatef(ballLocation[0], ballLocation[1], ballLocation[2]);
    //glRotatef(currDir, 1, 0, 0);
    glLineWidth(5);
    // glBegin(GL_LINES);
    // {
        
    //     glVertex3f(ballLocation[0]+sin(currDir * pi / (double)180), 0,ballLocation[2]+cos(currDir * pi / (double)180) );
    //     glVertex3f(ballLocation[0], ballLocation[1], ballLocation[2]);
    //     // glVertex3f(0,0,0);
    //     // glVertex3f(10,0,0);
    // }
    // glEnd();
    // glBegin(GL_TRIANGLES);{
    // 	glVertex3f(10,0,0);
    // 	glVertex3f(8, 2, 0);
    // 	glVertex3f(8, -2, 0);
    // }glEnd();
    //glLineWidth(1);
    glPopMatrix();

    glPopMatrix();
}
float thres=0.5;
float border=4.5f;
void checkReflect()
{
    // Check if the ball is hitting the boundary
    int counter = 0;
    float dx= cos(currDir * pi / (double)180) *2;
    float dz= sin(currDir * pi / (double)180) *2;
    // if (ballLocation[0] - radius-thres < -5.0f || ballLocation[0] + radius +thres> 5.0f ||
    //     ballLocation[2] - radius -thres< -5.0f || ballLocation[2] + radius+thres > 5.0f)
    // {
    //     printf("hit  %d\n", counter++);
    //     hit = !hit;
    //     currDir = (currDir*2) % 360;
    // }

    // printf("ball location %f %f %f\n", ballLocation[0], ballLocation[1], ballLocation[2]);
    //  printf("\n");

    if (ballLocation[0] - radius+dx < -border || ballLocation[0] + radius +dx> border)
        
    {
        printf("hit  %d\n", counter++);
        hit = !hit;
        currDir = 180-currDir;
    }
    else if(ballLocation[2] - radius -dz< -border || ballLocation[2] + radius+dz > border)
    {
        printf("hit  %d\n", counter++);
        hit = !hit;
        currDir = 360-currDir;
    }
}

void moveBackward()
{   
    //if(hit)currDir=(currDir+180)%360;
    // double distanceTrav;
    // forwardAngle = (forwardAngle + 10);
    // distanceTrav = 2 * pi * 5.0 * (double)(10) / (double)(360);
    // ballLocation[0] -= cos(currDir * pi / (double)180) * distanceTrav;
    // ballLocation[2] -= sin(currDir * pi / (double)180) * distanceTrav;

    // if(hit)currDir=(currDir+180)%360;
    // //double distanceTrav;
    // //forwardAngle = (forwardAngle - 10);
    // //distanceTrav = 2 * pi * 5.0 * (double)(-10) / (double)(360);
    // double dx= cos(currDir * pi / (double)180) *2;
    // double dz= sin(currDir * pi / (double)180) *2;
    // ballLocation[0] += dx;
    // ballLocation[2] += dz;
    // forwardAngle = (forwardAngle - 10);

    
}
void moveForward()
{   
    // //if(hit)currDir=(currDir+180)%360;
    // double distanceTrav;
    // forwardAngle = (forwardAngle - 10);
    // distanceTrav = 2 * pi * 5.0 * (double)(-10) / (double)(360);
    // ballLocation[0] -= cos(currDir * pi / (double)180) * distanceTrav;
    // ballLocation[2] -= sin(currDir * pi / (double)180) * distanceTrav;

    //if(hit)currDir=(currDir+90)%360;
    //double distanceTrav;
    //forwardAngle = (forwardAngle - 10);
    //distanceTrav = 2 * pi * 5.0 * (double)(-10) / (double)(360);
    // double dx= cos(currDir * pi / (double)180) *2;
    // double dz= sin(currDir * pi / (double)180) *2;
    // ballLocation[0] += dx;
    // ballLocation[2] += dz;
    // forwardAngle = (forwardAngle + 10);


    //hreeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee


      float dx= cos(currDir * pi / (double)180) *2;
    float dz= sin(currDir * pi / (double)180) *2;
 
    if( ballLocation[0] + radius >= 4.5f)
    {
        //dx=border-radius-ballLocation[0];
         currDir = 180-currDir;
    }

    else if (ballLocation[0] - radius <= -4.5f  )
        
    {
        // printf("hit  \n");
        // hit = !hit;
        dx=-border+radius-ballLocation[0];
        currDir = 180-currDir;
    }
    else if( ballLocation[2] + radius >= 4.5f)
    {
        //dz=border-radius-ballLocation[2];
        currDir = 180-currDir;
    }
    else if(ballLocation[2] - radius <=-4.5f )
    {
        // printf("hit  \n" );
        // hit = !hit;
        //dz=-border+radius-ballLocation[2];
        currDir = 180-currDir;
    }
    ballLocation[0] += dx;
    ballLocation[2] += dz;
    forwardAngle +=180*2/(pi*radius);
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

    drawCheckerboardXZ(100, 1);
    axes();
    glPushMatrix();
    glTranslatef(-5.0f, 0.0f, 0.0f);
    drawBoundary(10, 1);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(5.0f, 0.0f, 0.0f);
    drawBoundary(10, 1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, 0, 5);
    glRotatef(90, 0, 1, 0);
    drawBoundary(10, 1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0, 0, -5);
    glRotatef(90, 0, 1, 0);
    drawBoundary(10, 1);
    glPopMatrix();

    glPushMatrix();
    drawSegmentedSphere();
    glPopMatrix();
    //checkReflect();
    glPopMatrix();
    glutSwapBuffers();
}
void doYourOwnMove()
{
    if (autoMove)
    {
       moveForward();   
        //printf("called");
        glutPostRedisplay();
    }
}
void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
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

    case 'i': // Move the ball forward
        
            moveForward();
        
        // if (hit)
        // {   currDir = (currDir*2) % 360;
        //     moveBackward();
            
        // }

        printf("ball location %f %f %f\n", ballLocation[0], ballLocation[1], ballLocation[2]);
        printf("\n");

        break;
    case 'k': // Move the ball backward
        
            moveBackward();
        
        

        printf("ball location %f %f %f\n", ballLocation[0], ballLocation[1], ballLocation[2]);
        printf("\n");
        break;

    case 'j': // direction will rotate counterclockwise

        currDir = (currDir - 1) % 360;
        break;

    case 'l': // direction will rotate clockwise

        currDir = (currDir + 1) % 360;
        break;

    case ' ':
        autoMove = !autoMove; // Toggle the autoMove variable

        break;

    default:
        break;
    }
    glutPostRedisplay();
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

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Practice 1");
    init();
    glutDisplayFunc(display);
    glutIdleFunc(doYourOwnMove);
    glutSpecialFunc(specialKeyListener);
    glutKeyboardFunc(keyboardListener);
    glEnable(GL_DEPTH_TEST);
    glutMainLoop();

    return 0;
}