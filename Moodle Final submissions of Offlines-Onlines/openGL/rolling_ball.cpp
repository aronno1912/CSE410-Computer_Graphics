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
float forwardX = -2, forwardY = -2, forwardZ = -2;
float upx = 0, upy = 1, upz = 0;          // Up vector
float rightx = 1, righty = 0, rightz = 0; // right direction
float movementFactor = 0.1;
bool autoMove = false; // Whether the ball should move automatically or not
float stacks = 8, sectors = 8;
float ballLocation[3] = {0, 0.5, 0};
float ballDirection[3] = {0, 0, 0};
float radius = 0.5;
float revolvingAngle = 0;
float ballRotation = 0;

float directionX = 0.0f; // Initial direction along x-axis
float directionY = 0.0f; // Initial direction along y-axis

bool hit = false;
int currDir = 45;

struct point
{
    GLfloat x, y, z;
};

void computeDirection()
{
    float angleInRadians = currDir * pi / 180.0f;
    directionX = cos(angleInRadians);
    directionY = sin(angleInRadians);
}
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

void drawCheckerboard(int size, float squareSize)
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
    glRotatef(revolvingAngle, 0, 0, 1);

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


    glPopMatrix();
}

void drawSphere(double radius,int sectors,int stacks)
	{

    glPushMatrix();
    glTranslatef(ballLocation[0], ballLocation[1], ballLocation[2]);
    //glRotatef(revolvingAngle, 1, 0, 1);
    glRotatef(revolvingAngle, sin(revolvingAngle* pi/180),0,cos(revolvingAngle*pi/180));
		struct point points[50][50];
		int i,j;
		double h,r;
		//generate points
		for(i=0;i<=stacks;i++)
		{
			h=radius*sin(((double)i/(double)stacks)*(pi/2)); //the height at the current latitude.
			r=radius*cos(((double)i/(double)stacks)*(pi/2));  //the radius at the current latitude.
			for(j=0;j<=sectors;j++)
			{
				points[i][j].x=r*cos(((double)j/(double)sectors)*2*pi);
				points[i][j].y=r*sin(((double)j/(double)sectors)*2*pi);
				points[i][j].z=h;
			}
		}
		//draw quads using generated points
		for(i=0;i<stacks;i++)
		{
			int nowColor = 0;
			// glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
			for(j=0;j<sectors;j++)
			{
				if(nowColor){
					glColor3f(1,0,0);
				}else{
					glColor3f(0,1,0);
				}
				glBegin(GL_QUADS);{
					//upper hemisphere
					glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
					glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
					glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
					glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
				}glEnd();


				if(nowColor==1)
                {
					glColor3f(0,1,0);
				}
                else
                {

					glColor3f(1,0,0);
				}

				glBegin(GL_QUADS);
                {
					//lower hemisphere
					glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
					glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
					glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
					glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
				}glEnd();

				if(j % (2) == 0) nowColor = 1-nowColor;
			}
		}
        //glRotatef(-90,1,0,0);
        glPopMatrix();
    }
float thres=0.5;
float border=4.9f;
void checkReflect()
{

float angleInRadians = currDir * pi / 180.0f;

// Calculate the direction components
computeDirection();

// Check if the ball hits the boundary
if (ballLocation[0] - radius < -border || ballLocation[0] + radius > border ||
    ballLocation[2] - radius < -border || ballLocation[2] + radius > border)
{

 // Reflect along the x-axis
    if (ballLocation[0] - radius < -border || ballLocation[0] + radius > border)
        directionX = -directionX;

    // Reflect along the z-axis
    if (ballLocation[2] - radius < -border || ballLocation[2] + radius > border)
        directionY = -directionY;

    // Update the ball's direction
    currDir = atan2(directionY, directionX) * 180.0f / pi;
    hit=!hit;

    printf("hit \n");
}
 
  
}

void moveBackward()
{   
    checkReflect();
    ballLocation[0] -= 0.5f * directionX;
    ballLocation[2] -= 0.5f * directionY;
    revolvingAngle = (revolvingAngle - 15);

    
}
void moveForward()
{  
    //hreeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
     checkReflect();

   ballLocation[0] += 0.5f * directionX;
    ballLocation[2] += 0.5f * directionY;
    revolvingAngle = (revolvingAngle + 15);


}

void showAlwaysDirection()
{
    glLineWidth(5);
  glColor3f(0.0, 0.0, 1.0); // z blue
 glBegin(GL_LINES);
    {
        glVertex3f(ballLocation[0], ballLocation[1], ballLocation[2]);
        glVertex3f(ballLocation[0]+cos(currDir * pi / (double)180)*1.5, 0,ballLocation[2]+sin(currDir * pi / (double)180)*1.5 );
        
        // glVertex3f(0,0,0);
        // glVertex3f(10,0,0);
    }
    glEnd();


     // Draw the arrow tip (triangle)
    glBegin(GL_TRIANGLES);
    {
        // Calculate the points of the triangle
        float arrowTipX = ballLocation[0] + cos(currDir * pi / (double)180) * 1.5;
        float arrowTipY = 0;
        float arrowTipZ = ballLocation[2] + sin(currDir * pi / (double)180) * 1.5;

        float angle1 = currDir + 135; // Angle for one side of the triangle
        float angle2 = currDir - 135; // Angle for the other side of the triangle

        // Calculate the points of the triangle based on the angles
        float tipX1 = arrowTipX + 0.2 * cos(angle1 * pi / 180);
        float tipY1 = arrowTipY;
        float tipZ1 = arrowTipZ + 0.2 * sin(angle1 * pi / 180);

        float tipX2 = arrowTipX + 0.2 * cos(angle2 * pi / 180);
        float tipY2 = arrowTipY;
        float tipZ2 = arrowTipZ + 0.2 * sin(angle2 * pi / 180);

        // Draw the triangle
        glVertex3f(arrowTipX, arrowTipY, arrowTipZ);
        glVertex3f(tipX1, tipY1, tipZ1);
        glVertex3f(tipX2, tipY2, tipZ2);
    }
    glEnd();
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
    drawCheckerboard(100, 1);
    //axes();
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
    //glRotatef(-90, 0, 0, 1);
    //drawSegmentedSphere();
    drawSphere(0.5, 25, 25);
    glPopMatrix();
    //checkReflect();
    showAlwaysDirection();
    glPopMatrix();
    glutSwapBuffers();
}

float calculateMinTimeToHitBoundary()
{
    // Initialize minimum time to a large value
    float minTime = std::numeric_limits<float>::infinity();

    // Check the four sides of the boundary (left, right, front, back)
    float times[4];

    // Check left boundary
    if (directionX < 0)
    {
        times[0] = (ballLocation[0] - radius + border) / -directionX;
        if (times[0] < minTime)
        {
            minTime = times[0];
        }
    }

    // Check right boundary
    if (directionX > 0)
    {
        times[1] = (border - ballLocation[0] - radius) / directionX;
        if (times[1] < minTime)
        {
            minTime = times[1];
        }
    }

    // Check front boundary
    if (directionY < 0)
    {
        times[2] = (ballLocation[2] - radius + border) / -directionY;
        if (times[2] < minTime)
        {
            minTime = times[2];
        }
    }

    // Check back boundary
    if (directionY > 0)
    {
        times[3] = (border - ballLocation[2] - radius) / directionY;
        if (times[3] < minTime)
        {
            minTime = times[3];
        }
    }

    return minTime;
}
void eventDrivenSchedule(int val)
{
    int minTime=5;
    if (autoMove)
    {
         moveForward();
         minTime=calculateMinTimeToHitBoundary();
         printf("%d\n",minTime);
        //  glutTimerFunc(minTime, eventDriven, 0);
        //  glutPostRedisplay();

        
    }
   
    glutTimerFunc(minTime, eventDrivenSchedule, 0);
     glutPostRedisplay();
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

    case 'i': // Move the ball forward
        
            moveForward();

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

    case 'w': //move camera up without changing reference point
         eyex += upx * movementFactor;
        eyey += upy * movementFactor;
        eyez += upz * movementFactor;


        break;

    case 's': //move camera down without changing reference point
         eyex -= upx * movementFactor;
        eyey -= upy * movementFactor;
        eyez -= upz * movementFactor;

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
        eyex += forwardX * movementFactor;
        eyey += forwardY * movementFactor;
        eyez += forwardZ * movementFactor;

        break;
    case GLUT_KEY_DOWN:
        eyex -= forwardX * movementFactor;
        eyey -= forwardY * movementFactor;
        eyez -= forwardZ * movementFactor;
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


        forwardX += upx * movementFactor;
        forwardY += upy * movementFactor;
        forwardZ += upz * movementFactor;
        break;

    case GLUT_KEY_PAGE_DOWN:
        eyex -= upx * movementFactor;
        eyey -= upy * movementFactor;
        eyez -= upz * movementFactor;


        forwardX -= upx * movementFactor;
        forwardY -= upy * movementFactor;
        forwardZ -= upz * movementFactor;
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
    //glutIdleFunc(doYourOwnMove);
    glutTimerFunc(0, eventDrivenSchedule, 0);
    glutSpecialFunc(specialKeyListener);
    glutKeyboardFunc(keyboardListener);
    glEnable(GL_DEPTH_TEST);
    glutMainLoop();

    return 0;
}

