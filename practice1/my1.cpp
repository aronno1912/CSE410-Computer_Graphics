#include <bits/stdc++.h>
using namespace std;

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif
int counter = 0;
void display()
{
   printf("It is display function %d\n",++counter);
}
int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(500, 500);
    glutCreateWindow("Hello");
    glutDisplayFunc(display);
    glutMainLoop();


return 0;
}