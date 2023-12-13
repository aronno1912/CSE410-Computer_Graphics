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
float upx =0, upy=1, upz=0;        // Up vector
float rightx = 1,righty= 0,rightz= 0;         // right direction
float movementFactor=0.1;
bool autoMove = false; // Whether the ball should move automatically or not

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
		glBegin(GL_LINES);{
            glColor3f(1.0, 0.0, 0.0); //x red
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);
             glColor3f(0.0, 1.0, 0.0); // y green
			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);
            glColor3f(0.0, 0.0, 1.0); //z blue
			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
    
}


void square(double a){
    
    glBegin(GL_QUADS);
    {
        glVertex3f(a, a, 0);
        glVertex3f(a, -a, 0);
        glVertex3f(-a, -a, 0);
        glVertex3f(-a, a, 0);
    }glEnd();
}



void drawCheckerboard(int size, float squareSize) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            // Calculate coordinates for each square
            float x = i * squareSize;
            float y = j * squareSize;

            // Alternate between black and white squares
            if ((i + j) % 2 == 0) {
                // Black square
                glColor3f(0.0f, 0.0f, 0.0f);
            } else {
                // White square
                glColor3f(1.0f, 1.0f, 1.0f);
            }

            // Draw the square
            glBegin(GL_QUADS);
            glVertex2f(x, y);
            glVertex2f(x + squareSize, y);
            glVertex2f(x + squareSize, y + squareSize);
            glVertex2f(x, y + squareSize);
            glEnd();
        }
    }
}

void drawCheckerboardXZ(int size, float squareSize) {
    for (int i = -size; i < size; ++i) {
        for (int j = -size; j < size; ++j) {
            // Calculate coordinates for each square on the XZ plane
            float x = i * squareSize;
            float z = j * squareSize;

            // Alternate between black and white squares
            if ((i + j) % 2 == 0) {
                // Black square
                glColor3f(0.0f, 0.0f, 0.0f);
            } else {
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

 void oneBorder(int length,int width)
 {
    glVertex3f( 0.0f, width/2,length/2);
    glVertex3f( 0.0f, -width/2,length/2);
    glVertex3f( 0.0f, -width/2,-length/2);
    glVertex3f( 0.0f, width/2,-length/2);
 }
void drawBoundary(float length,float width)
{   
    glPushMatrix();
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_QUADS);
    glVertex3f( 0.0f, width/2,length/2);
    glVertex3f( 0.0f, -width/2,length/2);
    glVertex3f( 0.0f, -width/2,-length/2);
    glVertex3f( 0.0f, width/2,-length/2);
    glEnd();
    glPopMatrix();


   
}
float stacks=10,sectors=10;
float ballLocation[3]={0,0.5,0};
float ballDirection[3]={0,0,0};
float radius=0.5;
float forwardAngle=0;
float ballRotation=0;

float directionX = 0.0f; // Initial direction along x-axis
float directionY = 0.0f; // Initial direction along y-axis
void drawSegmentedSphere()
{
    
    float stackAngle = pi / stacks, phi = -pi / 2;
    float sectorAngle = 2 * pi / sectors, theta;
 
    glPushMatrix();
    glTranslatef(ballLocation[0], ballLocation[1], ballLocation[2]);
    //glRotatef(forwardAngle, 0, 0, 1);

     glRotatef(ballRotation, 0, 0, 1);
    for (int i = 0; i <= stacks; i++)
    {
        theta = 0;
        glBegin(GL_QUAD_STRIP);
        {
            for (int j = 0; j <= sectors; j++)
            {
                if ((i - 1) * 2 >= stacks)
                //if ((i + j) % 2 == 0)
                    glColor3f(j % 2, (j + 1) % 2, 0);
                   // glColor3f(1, 0, 0);
                else
                    glColor3f((j + 1) % 2, j % 2, 0);
                    //glColor3f(0, 1, 0);
                glVertex3f(radius * cos(phi - stackAngle) * cos(theta), radius * cos(phi - stackAngle) * sin(theta), radius * sin(phi - stackAngle));
                glVertex3f(radius * cos(phi) * cos(theta), radius * cos(phi) * sin(theta), radius * sin(phi));
                theta += sectorAngle;
            }
        }
        glEnd();
        phi += stackAngle;
    }
    
    // glEnd();
     //glPushMatrix();
    //  //Draw the arrow
    // float arrowLength = 1.5;
    // glBegin(GL_LINES);
    // glColor3f(0, 0, 1); // Red color for the arrow
    // glVertex3f(0, 0, 0);
    // glVertex3f(arrowLength * cos(forwardAngle * pi / 180.0f),0, arrowLength * sin(forwardAngle * pi / 180.0f));
    // glEnd();
    //glPopMatrix();
    
    glPopMatrix();
}

// void checkReflect()
// {   //printf("check reflect\n");
//     // Check if the ball is hitting the boundary
//     int counter=0;
//     if (ballLocation[0] - radius < -4.5f || ballLocation[0] + radius > 5.5f ||
//         ballLocation[1] - radius < -4.5f || ballLocation[1] + radius > 5.5f)
//     {

//         printf("hit  %d\n",counter++);
//         // Calculate the normal vector of the boundary (assumed to be a square boundary)
//         float normalX = 0, normalZ = 0;
//         if (ballLocation[0] - radius < -4.5f || ballLocation[0] + radius > 5.5f)
//             normalX = (ballLocation[0] - radius < -4.5f) ? 1.0f : -1.0f; // Reflect along the x-axis

//         if (ballLocation[1] - radius < -4.5f || ballLocation[1] + radius > 5.5f)
//             normalZ = (ballLocation[1] - radius < -4.5f) ? 1.0f : -1.0f; // Reflect along the z-axis

//         // Calculate the reflection direction
//         float dotProduct = normalX * directionX + normalZ * directionY;
//         float reflectionX = directionX - 2 * dotProduct * normalX;
//         float reflectionY = directionY - 2 * dotProduct * normalZ;

//         // Update the direction vector with the reflection direction
//         directionX = reflectionX;
//         directionY = reflectionY;

//          ballRotation -= 15;
//     forwardAngle -= .5;
//     ballLocation[0] -= 0.1f * cos(forwardAngle * pi / 180.0f);
//     ballLocation[1] -= 0.1f * sin(forwardAngle * pi / 180.0f);
//     }
// }
bool hit=false;
float currDir=45;

void checkReflect()
{
    // Check if the ball is hitting the boundary
    int counter = 0;
    if (ballLocation[0] - radius < -4.5f || ballLocation[0] + radius > 5.5f ||
        ballLocation[1] - radius < -4.5f || ballLocation[1] + radius > 5.5f)
    {
        printf("hit  %d\n", counter++);
        // Calculate the normal vector of the boundary (assumed to be a square boundary)
        // float normalX = 0, normalZ = 0;
        // if (ballLocation[0] - radius < -4.5f || ballLocation[0] + radius > 5.5f)
        //     normalX = (ballLocation[0] - radius < -4.5f) ? 1.0f : -1.0f; // Reflect along the x-axis

        // if (ballLocation[1] - radius < -4.5f || ballLocation[1] + radius > 5.5f)
        //     normalZ = (ballLocation[1] - radius < -4.5f) ? 1.0f : -1.0f; // Reflect along the z-axis

        // // Calculate the reflection direction
        // float dotProduct = normalX * directionX + normalZ * directionY;
        // float reflectionX = directionX - 2 * dotProduct * normalX;
        // float reflectionY = directionY - 2 * dotProduct * normalZ;

        // // Update the direction vector with the reflection direction
        // directionX = reflectionX;
        // directionY = reflectionY;
        hit=!hit;
    }

    // Update the ball's position based on the direction vector
    ballLocation[0] += 0.1f * directionX;
    ballLocation[1] += 0.1f * directionY;

    printf("ball location %f %f %f\n", ballLocation[0], ballLocation[1], ballLocation[2]);
    printf("\n");
}
void display()
{
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // look at e eye use hoy...eye.x,eye.y,eye.z ,center o
    gluLookAt(
        eyex, eyey, eyez, 
        centerx, centery, centerz, 
        0, 1, 0
        );
    glPushMatrix();
   
    drawCheckerboardXZ(100, 1);
     axes();
     glPushMatrix();
    glTranslatef(-5.0f, 0.0f, 0.0f);
    drawBoundary(10,1);
    glPopMatrix();
    glPushMatrix();
    glTranslatef(5.0f, 0.0f, 0.0f); 
    drawBoundary(10,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,5);
    glRotatef(90,0,1,0);
    drawBoundary(10,1);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(0,0,-5);
    glRotatef(90,0,1,0);
    drawBoundary(10,1);
    glPopMatrix();
    

    glPushMatrix();
    //sphere drawn on xy plane. own axis is z axis.. 
    //glRotatef(90,1,0,0);
    drawSegmentedSphere();
    glPopMatrix();
    checkReflect();
    glPopMatrix();
    glutSwapBuffers();

}
void doYourOwnMove()
{
    if(!hit)
    {
        ballRotation += 15;
    forwardAngle += .5;
    ballLocation[0] += 0.1f * cos(forwardAngle * pi / 180.0f);
    ballLocation[1] += 0.1f * sin(forwardAngle * pi / 180.0f);
    }
    if(hit)
    {     ballRotation -= 15;
    forwardAngle -= .5;
        ballLocation[0] -= 0.1f * cos(forwardAngle * pi / 180.0f);
        ballLocation[1] -= 0.1f * sin(forwardAngle * pi / 180.0f);
    }


}
void keyboardListener(unsigned char key, int x, int y)
{
   switch (key)
   {
   case '1':
    // rightx=rightx*cos(movementFactor)+centerx*sin(movementFactor);
    // righty=righty*cos(movementFactor)+centery*sin(movementFactor);
    // rightz=rightz*cos(movementFactor)+centerz*sin(movementFactor);
    
    //change where it is looking now...means facing forward is changes
    centerx=centerx*cos(movementFactor)-rightx*sin(movementFactor);
    centery=centery*cos(movementFactor)-righty*sin(movementFactor);
    centerz=centerz*cos(movementFactor)-rightz*sin(movementFactor);
        
      break;
     case '2': 
        // rightx=rightx*cos(-movementFactor)+centerx*sin(-movementFactor);
        // righty=righty*cos(-movementFactor)+centery*sin(-movementFactor);
        // rightz=rightz*cos(-movementFactor)+centerz*sin(-movementFactor);

        //change where it is looking now...means facing forward is changes
        centerx=centerx*cos(-movementFactor)-rightx*sin(-movementFactor);
        centery=centery*cos(-movementFactor)-righty*sin(-movementFactor);
        centerz=centerz*cos(-movementFactor)-rightz*sin(-movementFactor);
        break;
   case '3':

  // center and up are perpendicular to each other

    
    centerx=centerx*cos(movementFactor)+upx*sin(movementFactor);
    centery=centery*cos(movementFactor)+upy*sin(movementFactor);
    centerz=centerz*cos(movementFactor)+upz*sin(movementFactor);

    // upx=upx*cos(movementFactor)-centerx*sin(movementFactor);
    // upy=upy*cos(movementFactor)-centery*sin(movementFactor);
    // upz=upz*cos(movementFactor)-centerz*sin(movementFactor);
        
      break;
     case '4': 
        centerx=centerx*cos(-movementFactor)+upx*sin(-movementFactor);
        centery=centery*cos(-movementFactor)+upy*sin(-movementFactor);
        centerz=centerz*cos(-movementFactor)+upz*sin(-movementFactor);

        // upx=upx*cos(-movementFactor)*centerx*sin(-movementFactor);
        // upy=upy*cos(-movementFactor)-centery*sin(-movementFactor);
        // upz=upz*cos(-movementFactor)-centerz*sin(-movementFactor);
        break;

   case '5':
      // update the up vector by rotating it around the right vector

    upx = upx * cos(movementFactor) + rightx * sin(movementFactor);
    upy = upy * cos(movementFactor) + righty * sin(movementFactor);
    upz = upz * cos(movementFactor) + rightz * sin(movementFactor);

    //now adjust right vector accordingly
    rightx=rightx*cos(movementFactor)-upx*sin(movementFactor);
    righty=righty*cos(movementFactor)-upy*sin(movementFactor);
    rightz=rightz*cos(movementFactor)-upz*sin(movementFactor);

    //  // Recompute the center vector to be perpendicular to the new up and right vectors
    // centerx = upy * rightz - upz * righty;
    // centery = upz * rightx - upx * rightz;
    // centerz = upx * righty - upy * rightx;

  

      break;

    case '6':   
    upx = upx * cos(-movementFactor) + rightx * sin(-movementFactor);
    upy = upy * cos(-movementFactor) + righty * sin(-movementFactor);
    upz = upz * cos(-movementFactor) + rightz * sin(-movementFactor);

    //now adjust right vector accordingly
    rightx=rightx*cos(-movementFactor)-upx*sin(-movementFactor);
    righty=righty*cos(-movementFactor)-upy*sin(-movementFactor);
    rightz=rightz*cos(-movementFactor)-upz*sin(-movementFactor);

    //  // Recompute the center vector to be perpendicular to the new up and right vectors
    // centerx = upy * rightz - upz * righty;
    // centery = upz * rightx - upx * rightz;
    // centerz = upx * righty - upy * rightx;


        break;

    case 'i': // Move the ball forward
    if(!hit)
    {
        ballRotation += 15;
    forwardAngle += .5;
    ballLocation[0] += 0.1f * cos(forwardAngle * pi / 180.0f);
    ballLocation[1] += 0.1f * sin(forwardAngle * pi / 180.0f);
    }
    if(hit)
    {     ballRotation -= 15;
    forwardAngle -= .5;
        ballLocation[0] -= 0.1f * cos(forwardAngle * pi / 180.0f);
        ballLocation[1] -= 0.1f * sin(forwardAngle * pi / 180.0f);
    }

    printf("ball location %f %f %f\n",ballLocation[0],ballLocation[1],ballLocation[2]);
    printf("\n");


    break;
    case 'k': // Move the ball backward
    if(!hit)
    {
    ballRotation -= 15;
    forwardAngle -= .5;
    ballLocation[0] -= 0.1f * cos(forwardAngle * pi / 180.0f);
    ballLocation[1] -= 0.1f * sin(forwardAngle * pi / 180.0f);}
    if(hit)
    {    
    ballRotation += 15;
    forwardAngle += .5;
    ballLocation[0] += 0.1f * cos(forwardAngle * pi / 180.0f);
    ballLocation[1] += 0.1f * sin(forwardAngle * pi / 180.0f);
    }

    printf("ball location %f %f %f\n",ballLocation[0],ballLocation[1],ballLocation[2]);
    printf("\n");
    break;

    case 'j': //direction will rotate counterclockwise

    forwardAngle +=.5;
    break;

    case 'l': //direction will rotate clockwise

    forwardAngle -=.5;
    break;

    case ' ':
        autoMove = !autoMove; // Toggle the autoMove variable
        printf("space pressed\n");
        if(autoMove)
        {
            glutIdleFunc(doYourOwnMove);
        }
        else
        {
            glutIdleFunc(NULL);
        }
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
      eyey += righty* movementFactor;
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


int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(100, 100);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutCreateWindow("Practice 1");
    init();
    glutDisplayFunc(display);
    glutSpecialFunc(specialKeyListener); 
    glutKeyboardFunc(keyboardListener);
    //glutIdleFunc(animate);
     glEnable(GL_DEPTH_TEST);
    glutMainLoop();

    return 0;
}