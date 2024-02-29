#include "myC.h"


extern vector <Object*> objects;
extern vector <PointLight> pointLights;
extern vector <SpotLight> spotLights;

extern int recursionLevel;
int pixels=0;
int num_of_objects=0;
int num_point_lights=0;
int num_spotlights=0;
int img_cap=1;

double cameraHeight;
double cameraAngle;
double theta;
int drawgrid;
int drawaxes;




Vector3D eye(200, 100, 100);
Vector3D u(0, 0, 1);
Vector3D r(-0.707,0.707, 0);
Vector3D l(-0.707, -0.707, 0);
int movementFactor=3;
float movementAngle=0.03;



double radius=20;
double height;

double windowHeight = 500;
double windowWidth = 500;

void clear_memory()
{
    cout << "clearing memory" << endl;
    for (auto p : objects)
    {
        delete p;
    }
    objects.clear();
    pointLights.clear();
    spotLights.clear();
}


void capture()
{
   cout << "capturing image " << img_cap << endl;
   bitmap_image image(pixels,pixels);
/**
 * planeDistance = (windowHeight/2.0) /
tan(viewAngle/2.0)
topleft = eye + l*planeDistance - r*windowWidth/2 +
u*windowHeight/2
du = windowWidth/imageWidth
dv = windowHeight/imageHeight
// Choose middle of the grid cell
topleft = topleft + r*(0.5*du) - u*(0.5*dv)

*/


   // Calculate the distance from the eye to the image plane
   double planeDistance = (windowHeight/2.0)/tan((cameraAngle/2.0)*(pi/180));
   // Calculate the top-left corner of the image plane
   Vector3D topleft;
   topleft.x = eye.x + (l.x)*planeDistance - (r.x)*windowWidth/2 + (u.x)*windowHeight/2;
   topleft.y = eye.y  + (l.y )*planeDistance - (r.y )*windowWidth/2 + (u.y )*windowHeight/2;
   topleft.z = eye.z + (l.z)*planeDistance - (r.z)*windowWidth/2 + (u.z)*windowHeight/2;
   double imageWidth,imageHeight;
   imageHeight = imageWidth = pixels;

   double du = (double) ((windowWidth*1.0)/(imageWidth*1.0));  // step size in the horizontal direction
   double dv = (double) ((windowHeight*1.0)/(imageHeight*1.0)); // step size in the vertical  direction
   topleft.x = topleft.x + r.x*(0.5*du) - u.x*(0.5*dv); //center of the pixel
   topleft.y = topleft.y + r.y*(0.5*du) - u.y*(0.5*dv);
   topleft.z = topleft.z + r.z*(0.5*du) - u.z*(0.5*dv);
   int nearest;
   double t,tMin;
   for(int i=0;i<pixels;i++)
   {
       for(int j=0;j<pixels;j++)
       {
           //calculate curPixel using topleft,r,u,i,j,du,dv
           //topLeft+r*(column*du)-u*(row*dv)
           image.set_pixel(i,j,0,0,0);
           Vector3D curPixel;
           curPixel.x = topleft.x + r.x*i*du - u.x*j*dv;
           curPixel.y = topleft.y + r.y*i*du - u.y*j*dv;
           curPixel.z = topleft.z + r.z*i*du - u.z*j*dv;
           //cast ray from eye to (curPixel-eye) direction
           Vector3D rd;
           rd.x = curPixel.x-eye.x;
           rd.y = curPixel.y-eye.y;
           rd.z = curPixel.z-eye.z;

           Ray *castedRay = new Ray(eye,rd);

           double *color;

           tMin = INT_MAX;
           nearest = 0;
           int idx=0;
           for (auto & o : objects)
           {
               //for each object, o in objects,find the nearest intersection point

               t = o->intersect(castedRay, color, 0);

//                update t so that it stores min +ve value
               if(t<tMin && t>0)
               {
                   tMin=t;
                   nearest = idx; //which index object is nearest
               }
               idx++;
           }
           // If there is an intersection, calculate the color and set the pixel in the image
           if(tMin!=INT_MAX)
           {
                Object* n = objects[nearest];
                tMin = n->intersect(castedRay,color,1);
                color = n->getColorBmp();
                image.set_pixel(i,j,color[0]*255.0,color[1]*255.0,color[2]*255.0);
           }

        delete castedRay;

       }
   }

cout<<"output"<<img_cap<<".bmp"<<endl;
image.save_image("output"+std::to_string(img_cap)+".bmp");
   cout << "captured image " << img_cap << endl;
   img_cap++;
}



void drawAllComponents()
{
    for (auto & individual : objects)
    {
       glPushMatrix();
       glTranslatef(individual->reference_point.x,individual->reference_point.y,individual->reference_point.z);
       individual->draw();
       cout<<individual->reference_point.x<<" "<<individual->reference_point.y<<" "<<individual->reference_point.z<<endl;
       cout<<"drawn"<<endl;
       glPopMatrix();
    }
    for (auto & individual : pointLights)
    {
       individual.draw();
    }
    for (auto & individual : spotLights)
    {
       individual.draw();
    }

}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
  

        case '0':
            capture();
            break;
        case '1':
           
			l.x = l.x*cos(movementAngle)-r.x*sin(movementAngle);
			l.y = l.y*cos(movementAngle)-r.y*sin(movementAngle);
			l.z = l.z*cos(movementAngle)-r.z*sin(movementAngle);

			r.x = l.x*sin(movementAngle)+r.x*cos(movementAngle);
			r.y = l.y*sin(movementAngle)+r.y*cos(movementAngle);
			r.z = l.z*sin(movementAngle)+r.z*cos(movementAngle);
			break;
    
        case '2':
           
			l.x = l.x*cos(-movementAngle)-r.x*sin(-movementAngle);
			l.y = l.y*cos(-movementAngle)-r.y*sin(-movementAngle);
			l.z = l.z*cos(-movementAngle)-r.z*sin(-movementAngle);

			r.x = l.x*sin(-movementAngle)+r.x*cos(-movementAngle);
			r.y = l.y*sin(-movementAngle)+r.y*cos(-movementAngle);
			r.z = l.z*sin(-movementAngle)+r.z*cos(-movementAngle);
			break;

 
       case '3':
            
			l.x = l.x*cos(movementAngle)-u.x*sin(movementAngle);
			l.y = l.y*cos(movementAngle)-u.y*sin(movementAngle);
			l.z = l.z*cos(movementAngle)-u.z*sin(movementAngle);

			u.x = l.x*sin(movementAngle)+u.x*cos(movementAngle);
			u.y = l.y*sin(movementAngle)+u.y*cos(movementAngle);
			u.z = l.z*sin(movementAngle)+u.z*cos(movementAngle);
			break;
 
        case '4':
           
			l.x = l.x*cos(-movementAngle)-u.x*sin(-movementAngle);
			l.y = l.y*cos(-movementAngle)-u.y*sin(-movementAngle);
			l.z = l.z*cos(-movementAngle)-u.z*sin(-movementAngle);

			u.x = l.x*sin(-movementAngle)+u.x*cos(-movementAngle);
			u.y = l.y*sin(-movementAngle)+u.y*cos(-movementAngle);
			u.z = l.z*sin(-movementAngle)+u.z*cos(-movementAngle);
			break;
 
        case '5':
            
			u.x = u.x*cos(-movementAngle)-r.x*sin(-movementAngle);
			u.y = u.y*cos(-movementAngle)-r.y*sin(-movementAngle);
			u.z = u.z*cos(-movementAngle)-r.z*sin(-movementAngle);

			r.x = u.x*sin(-movementAngle)+r.x*cos(-movementAngle);
			r.y = u.y*sin(-movementAngle)+r.y*cos(-movementAngle);
			r.z = u.z*sin(-movementAngle)+r.z*cos(-movementAngle);
			break;
 
        case '6':
            
			u.x = u.x*cos(movementAngle)-r.x*sin(movementAngle);
			u.y = u.y*cos(movementAngle)-r.y*sin(movementAngle);
			u.z = u.z*cos(movementAngle)-r.z*sin(movementAngle);

			r.x = u.x*sin(movementAngle)+r.x*cos(movementAngle);
			r.y = u.y*sin(movementAngle)+r.y*cos(movementAngle);
			r.z = u.z*sin(movementAngle)+r.z*cos(movementAngle);
			break;



		default:
			break;
	}
    glutPostRedisplay();
}


void specialKeyListener(int key, int x,int y){
	switch(key){

		case GLUT_KEY_UP:		

			eye.x = eye.x + movementFactor*l.x;
			eye.y = eye.y + movementFactor*l.y;
			eye.z = eye.z + movementFactor*l.z;
			break;

		case GLUT_KEY_DOWN:		

			eye.x = eye.x - movementFactor*l.x;
			eye.y = eye.y - movementFactor*l.y;
			eye.z = eye.z - movementFactor*l.z;
			break;


		case GLUT_KEY_RIGHT:

            eye.x = eye.x + movementFactor*r.x;
			eye.y = eye.y + movementFactor*r.y;
			eye.z = eye.z + movementFactor*r.z;
			break;
		case GLUT_KEY_LEFT:

            eye.x = eye.x - movementFactor*r.x;
			eye.y = eye.y - movementFactor*r.y;
			eye.z = eye.z - movementFactor*r.z;
			break;

		case GLUT_KEY_PAGE_UP:

            eye.x = eye.x + movementFactor*u.x;
			eye.y = eye.y + movementFactor*u.y;
			eye.z = eye.z + movementFactor*u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:

            eye.x = eye.x - movementFactor*u.x;
			eye.y = eye.y - movementFactor*u.y;
			eye.z = eye.z - movementFactor*u.z;
			break;



		default:
			break;
	}

    glutPostRedisplay();
}





void init(){


	cameraAngle=80;
	//clear the screen
	glClearColor(0,0,0,0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(cameraAngle,	1,	1,	1000.0);

}


void loadData()
{
    std::ifstream in("scene.txt");
    std::streambuf *cinbuf = std::cin.rdbuf();
    std::cin.rdbuf(in.rdbuf()); 

    cin >> recursionLevel;
    cin >> pixels;
    cin >> num_of_objects;



    string object_type;
    double color[3];
    double ambient, diffuse, specular, recursive_reflection_coefficient;
    int shine;

    for(int i=0;i<num_of_objects;i++)
    {
        Object *obj;
        cin >> object_type;
     
        if(object_type=="triangle")
        {
            Vector3D* points[3];
            for(int j=0;j<3;j++)
            {
                points[j] = new Vector3D();
                cin >> points[j]->x >> points[j]->y >> points[j]->z;
            }
            cin >> color[0] >> color[1] >> color[2];
            cin >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
            cin >> shine;

            obj = new Triangle(*points[0],*points[1],*points[2]);
            obj->setColor(color[0],color[1],color[2]);
            obj->setCoEfficients(ambient, diffuse, specular, recursive_reflection_coefficient);
            obj->setShine(shine);
        }
        else if(object_type=="sphere")
        {
            Vector3D *center = new Vector3D();
            cin >> center->x >> center->y >> center->z;

            double radius;

            cin >> radius;
            cin >> color[0] >> color[1] >> color[2];
            cin >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
            cin >> shine;

            obj = new Sphere(*center,radius);
            obj->setColor(color[0],color[1],color[2]);
            obj->setCoEfficients(ambient, diffuse, specular, recursive_reflection_coefficient);
            obj->setShine(shine);

        }
        else
        {
            double generalCoefficients[10],length,width,height;
            Vector3D *ref_point = new Vector3D();
            for(int j=0;j<10;j++)
            {
                cin >> generalCoefficients[j];
            }
            cin >> ref_point->x >> ref_point->y >> ref_point->z >> length >> width >> height;
            cin >> color[0] >> color[1] >> color[2];
            cin >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
            cin >> shine;

            obj = new Object(generalCoefficients,*ref_point,length,width,height);
            obj->setColor(color[0],color[1],color[2]);
            obj->setCoEfficients(ambient, diffuse, specular, recursive_reflection_coefficient);
            obj->setShine(shine);
        }
        obj->setType(object_type);
        objects.push_back(obj);
    }

    cin >> num_point_lights;
    for(int i=0;i<num_point_lights;i++)
    {
        PointLight pointlight;
        Vector3D pos;
        double color[3];
        cin >> pos.x >> pos.y >> pos.z;
        cin >> color[0] >> color[1] >> color[2];
        for(int j=0;j<3;j++)
        {
            pointlight.color[j] = color[j];
        }
        pointlight.light_pos = pos;
        pointLights.push_back(pointlight);
    }

    cin >> num_spotlights;
    for(int i=0;i<num_spotlights;i++)
    {
        SpotLight spotlight;
        Vector3D pos,dir;
        double color[3],cutoff_angle;
        cin >> pos.x >> pos.y >> pos.z;
        cin >> color[0] >> color[1] >> color[2];
        cin >> dir.x >> dir.y >> dir.z;
        cin >> cutoff_angle;

        for(int j=0;j<3;j++)
        {
            spotlight.point_light.color[j] = color[j];
        }
        spotlight.point_light.light_pos = pos;
        spotlight.light_direction = dir;
        spotlight.cutoff_angle = cutoff_angle;

        spotLights.push_back(spotlight);
    }

    // lasssstly push the floor
    Object* floor = new Floor(1000,20);
    objects.push_back(floor);


}

void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    gluLookAt(eye.x,eye.y,eye.z,	eye.x+l.x,eye.y+l.y,eye.z+
              l.z,	u.x,u.y,u.z);
	glMatrixMode(GL_MODELVIEW);
    drawAllComponents();
	glutSwapBuffers();
}


int main(int argc, char **argv){

    loadData();

	glutInit(&argc,argv);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color
	glutCreateWindow("Practice 1");
	init();
    glEnable(GL_CLIP_PLANE0);
	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
    atexit(clear_memory);
	glutMainLoop();		//The main loop of OpenGL

	return 0;
}


