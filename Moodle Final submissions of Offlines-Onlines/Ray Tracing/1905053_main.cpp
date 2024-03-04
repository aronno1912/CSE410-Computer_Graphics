#include "1905053_classes.h"


// Vector3D eye(210, 120, 100);
// Vector3D forwardVector(-0.707, -0.707, 0);
// Vector3D upVector(0, 0, 1);
// Vector3D rightVector(-0.707,0.707, 0);

Vector3D eye(190, 140, 100);
Vector3D forwardVector(-0.707, -0.707, 0);
Vector3D upVector(0, 0, 1);
Vector3D rightVector(-0.707,0.707, 0);


extern vector <Object*> objects;
extern vector <PointLight> pointLights;
extern vector <SpotLight> spotLights;
extern int recursionLevel;

int movementFactor=3;
double movementAngle=0.03;
double cameraAngle;
int pixels=0;
double windowHeight = 500.0;
double windowWidth = 500.0;

int imgcount=1;

void capture()
{
   cout << "capturing image " << imgcount << endl;
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
   double imageWidth,imageHeight;
   double planeDistance = (windowHeight/2.0)/tan((cameraAngle/2.0)*(pi/180));
   // Calculate the top-left corner of the image plane
   Vector3D topleft;
   topleft.x = eye.x + (forwardVector.x)*planeDistance - (rightVector.x)*windowWidth/2 + (upVector.x)*windowHeight/2;
   topleft.y = eye.y  + (forwardVector.y )*planeDistance - (rightVector.y )*windowWidth/2 + (upVector.y )*windowHeight/2;
   topleft.z = eye.z + (forwardVector.z)*planeDistance - (rightVector.z)*windowWidth/2 + (upVector.z)*windowHeight/2;
   
   imageHeight =pixels;
   imageWidth = pixels;

   double du = (double) ((windowWidth)/(imageWidth*1.0));  // step size in the horizontal direction
   double dv = (double) ((windowHeight)/(imageHeight*1.0)); // step size in the vertical  direction
   topleft.x = topleft.x + rightVector.x*(0.5*du) - upVector.x*(0.5*dv); //center of the pixel
   topleft.y = topleft.y + rightVector.y*(0.5*du) - upVector.y*(0.5*dv);
   topleft.z = topleft.z + rightVector.z*(0.5*du) - upVector.z*(0.5*dv);
   double t,tMin;
   int nearestObject;
   for(int i=0;i<pixels;i++)
   {
       for(int j=0;j<pixels;j++)
       {
           //calculate currentPixel using topleft,r,u,i,j,du,dv
           //topLeft+r*(column*du)-u*(row*dv)
           image.set_pixel(i,j,0,0,0);
           double red;
           double green;
            double blue;
           Vector3D currentPixel;
           currentPixel.x = topleft.x + rightVector.x*i*du - upVector.x*j*dv;
           currentPixel.y = topleft.y + rightVector.y*i*du - upVector.y*j*dv;
           currentPixel.z = topleft.z + rightVector.z*i*du - upVector.z*j*dv;
           //cast ray from eye to (currentPixel-eye) direction
           Vector3D rd;
           rd.x = currentPixel.x-eye.x;
           rd.y = currentPixel.y-eye.y;
           rd.z = currentPixel.z-eye.z;

           Ray *castedRay = new Ray(eye,rd);

           tMin = INT_MAX;
           nearestObject = 0;
           int whichObject=0;
           double *color;
           for (auto & o : objects)
           {
               //for each object, o in objects,find the nearest intersection point

               t = o->intersect(castedRay, color, 0);

//                update t so that it stores min +ve value
               if(t<tMin && t>0)
               {
                   tMin=t;
                   nearestObject = whichObject; //which index object is nearest
               }
               whichObject++;
           }
           // If there is an intersection, calculate the color and set the pixel in the image
           if(tMin!=INT_MAX)
           {
                Object* n = objects[nearestObject];
                tMin = n->intersect(castedRay,color,1);
                color = n->getPhotoColor();
                red = color[0]*255.0;
                green = color[1]*255.0;
                blue = color[2]*255.0;
                image.set_pixel(i,j,red,green,blue);
           }

        delete castedRay;

       }
   }

cout<<"output"<<imgcount<<".bmp"<<endl;
image.save_image("output"+to_string(imgcount)+".bmp");
cout << "captured image " << imgcount << endl;
imgcount++;
}



void drawAllComponents()
{
    for (auto & individual : objects)
    {
       glPushMatrix();
       glTranslatef(individual->reference_point.x,individual->reference_point.y,individual->reference_point.z);
       individual->draw();
       //cout<<individual->reference_point.x<<" "<<individual->reference_point.y<<" "<<individual->reference_point.z<<endl;
       //cout<<"drawn"<<endl;
       glPopMatrix();
    }
    for (auto & individual : pointLights)
    {
       individual.draw();
       //cout<<"pointLight position: "<<individual.light_pos.x<<" "<<individual.light_pos.y<<" "<<individual.light_pos.z<<endl;
    }
    for (auto & individual : spotLights)
    {
       individual.draw();
      //cout<<"spotLight position: "<<individual.point_light.light_pos.x<<" "<<individual.point_light.light_pos.y<<" "<<individual.point_light.light_pos.z<<endl; 
    }

}

void freeMemory()
{
    cout << "freeing" << endl;

    pointLights.clear();
    spotLights.clear();


    for (auto ob : objects)
    {
        delete ob;
    }
   
    objects.clear();

}


void keyboardListener(unsigned char key, int x,int y){
	switch(key){
  

        case '0':
            capture();
            break;
        case '1':
           
			forwardVector.x = forwardVector.x*cos(movementAngle)-rightVector.x*sin(movementAngle);
			forwardVector.y = forwardVector.y*cos(movementAngle)-rightVector.y*sin(movementAngle);
			forwardVector.z = forwardVector.z*cos(movementAngle)-rightVector.z*sin(movementAngle);

            rightVector.x=forwardVector.x*sin(movementAngle)+rightVector.x*cos(movementAngle);
            rightVector.y=forwardVector.y*sin(movementAngle)+rightVector.y*cos(movementAngle);
            rightVector.z=forwardVector.z*sin(movementAngle)+rightVector.z*cos(movementAngle);

			break;
    
        case '2':
           
            forwardVector.x = forwardVector.x*cos(-movementAngle)-rightVector.x*sin(-movementAngle);
            forwardVector.y = forwardVector.y*cos(-movementAngle)-rightVector.y*sin(-movementAngle);
            forwardVector.z = forwardVector.z*cos(-movementAngle)-rightVector.z*sin(-movementAngle);

             rightVector.x=forwardVector.x*sin(-movementAngle)+rightVector.x*cos(-movementAngle);
            rightVector.y=forwardVector.y*sin(-movementAngle)+rightVector.y*cos(-movementAngle);
            rightVector.z=forwardVector.z*sin(-movementAngle)+rightVector.z*cos(-movementAngle);

			break;

 
       case '3':
            forwardVector.x = forwardVector.x*cos(-movementAngle)-upVector.x*sin(-movementAngle);
            forwardVector.y = forwardVector.y*cos(-movementAngle)-upVector.y*sin(-movementAngle);
            forwardVector.z = forwardVector.z*cos(-movementAngle)-upVector.z*sin(-movementAngle);

            upVector.x=forwardVector.x*sin(-movementAngle)+upVector.x*cos(-movementAngle);
            upVector.y=forwardVector.y*sin(-movementAngle)+upVector.y*cos(-movementAngle);
            upVector.z=forwardVector.z*sin(-movementAngle)+upVector.z*cos(-movementAngle);


			break;
 
        case '4':
           

            forwardVector.x = forwardVector.x*cos(movementAngle)-upVector.x*sin(movementAngle);
            forwardVector.y = forwardVector.y*cos(movementAngle)-upVector.y*sin(movementAngle);
            forwardVector.z = forwardVector.z*cos(movementAngle)-upVector.z*sin(movementAngle);

              upVector.x=forwardVector.x*sin(movementAngle)+upVector.x*cos(movementAngle);
            upVector.y=forwardVector.y*sin(movementAngle)+upVector.y*cos(movementAngle);
            upVector.z=forwardVector.z*sin(movementAngle)+upVector.z*cos(movementAngle);


			
			break;
 
        case '5':
            
			upVector.x = upVector.x*cos(-movementAngle)-rightVector.x*sin(-movementAngle);
            upVector.y = upVector.y*cos(-movementAngle)-rightVector.y*sin(-movementAngle);
            upVector.z = upVector.z*cos(-movementAngle)-rightVector.z*sin(-movementAngle);

			rightVector.x = upVector.x*sin(-movementAngle)+rightVector.x*cos(-movementAngle);
            rightVector.y = upVector.y*sin(-movementAngle)+rightVector.y*cos(-movementAngle);
            rightVector.z = upVector.z*sin(-movementAngle)+rightVector.z*cos(-movementAngle);
			break;
 
        case '6':
            
			upVector.x = upVector.x*cos(movementAngle)-rightVector.x*sin(movementAngle);
            upVector.y = upVector.y*cos(movementAngle)-rightVector.y*sin(movementAngle);
            upVector.z = upVector.z*cos(movementAngle)-rightVector.z*sin(movementAngle);

            rightVector.x = upVector.x*sin(movementAngle)+rightVector.x*cos(movementAngle); 
            rightVector.y = upVector.y*sin(movementAngle)+rightVector.y*cos(movementAngle);
            rightVector.z = upVector.z*sin(movementAngle)+rightVector.z*cos(movementAngle);
			break;



		default:
			break;
	}
    glutPostRedisplay();
}


void specialKeyListener(int key, int x,int y){
	switch(key){

		case GLUT_KEY_UP:		

			eye.x = eye.x + movementFactor*forwardVector.x;
			eye.y = eye.y + movementFactor*forwardVector.y;
			eye.z = eye.z + movementFactor*forwardVector.z;
			break;

		case GLUT_KEY_DOWN:		

			eye.x = eye.x - movementFactor*forwardVector.x;
			eye.y = eye.y - movementFactor*forwardVector.y;
			eye.z = eye.z - movementFactor*forwardVector.z;
			break;


		case GLUT_KEY_RIGHT:

            eye.x = eye.x + movementFactor*rightVector.x;
			eye.y = eye.y + movementFactor*rightVector.y;
			eye.z = eye.z + movementFactor*rightVector.z;
			break;
		case GLUT_KEY_LEFT:

            eye.x = eye.x - movementFactor*rightVector.x;
			eye.y = eye.y - movementFactor*rightVector.y;
			eye.z = eye.z - movementFactor*rightVector.z;
			break;

		case GLUT_KEY_PAGE_UP:

            eye.x = eye.x + movementFactor*upVector.x;
			eye.y = eye.y + movementFactor*upVector.y;
			eye.z = eye.z + movementFactor*upVector.z;
			break;
		case GLUT_KEY_PAGE_DOWN:

            eye.x = eye.x - movementFactor*upVector.x;
			eye.y = eye.y - movementFactor*upVector.y;
			eye.z = eye.z - movementFactor*upVector.z;
			break;



		default:
			break;
	}

    glutPostRedisplay();
}





void init(){


	cameraAngle=80.0;
	//clear the screen
	glClearColor(0,0,0,0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(cameraAngle,	1,	1,	1000.0);

}

int totalObjects=0;
int totalPointLightss=0;
int totalSpotLights=0;


void loadData()
{   
    double color[3];
    double ambient;
    double diffuse;
    double specular;
    double recursive_reflection_coefficient;
    int shine;
    string object_type;
    ifstream in("scene.txt");
    streambuf *cinbuf = cin.rdbuf();
    cin.rdbuf(in.rdbuf()); 
    //handle error in file opening

    if(!in.is_open())
    {
        cout << "Error in opening file" << endl;
        return;
    } 
    //if file not found
    if(in.fail())
    {
        cout << "File not found" << endl;
        return;
    }  

    cin >> recursionLevel;
    cin >> pixels;
    cin >> totalObjects;

    for(int i=0;i<totalObjects;i++)
    {
        Object *obj;
        cin >> object_type;
     
        if(object_type=="triangle")
        {
            Vector3D* vertices[3];
            for(int j=0;j<3;j++)
            {
                vertices[j] = new Vector3D();
                cin >> vertices[j]->x >> vertices[j]->y >> vertices[j]->z;
            }
            cin >> color[0] >> color[1] >> color[2];
            cin >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
            cin >> shine;

            obj = new Triangle(*vertices[0],*vertices[1],*vertices[2]);
            cout << "triangle" << endl;
            cout << vertices[0]->x << " " << vertices[0]->y << " " << vertices[0]->z << endl;
            
            obj->setCoEfficients(ambient, diffuse, specular, recursive_reflection_coefficient);
            obj->setColor(color[0],color[1],color[2]);
            obj->setShine(shine);
        }


        else if(object_type=="sphere")
        {
            Vector3D *center = new Vector3D();
            double radius;
            cin >> center->x >> center->y >> center->z;
            cin >> radius;
            cin >> color[0] >> color[1] >> color[2];
            cin >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
            cin >> shine;

            obj = new Sphere(*center,radius);
            cout << "sphere" << endl;
            cout << center->x << " " << center->y << " " << center->z << endl;
            obj->setCoEfficients(ambient, diffuse, specular, recursive_reflection_coefficient);
            obj->setColor(color[0],color[1],color[2]);
            obj->setShine(shine);

        }
        else if(object_type=="prism")
        {

            Vector3D *vertices[6];
            for(int j=0;j<6;j++)
            {
                vertices[j] = new Vector3D();
                cin >> vertices[j]->x >> vertices[j]->y >> vertices[j]->z;
            }
           

            cin >> color[0] >> color[1] >> color[2];
            cin >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
            cin >> shine;
            double refractive_index[3];
            for (int j = 0; j < 3; j++)
            {
                cin >> refractive_index[j];
            }
            obj->setRefractiveIndex(refractive_index[0],refractive_index[1],refractive_index[2]);
            cout<<"prism"<<endl;
            obj = new Prism(*vertices[0],*vertices[1],*vertices[2],*vertices[3],*vertices[4],*vertices[5]);
            obj->setColor(color[0],color[1],color[2]);
            obj->setCoEfficients(ambient, diffuse, specular, recursive_reflection_coefficient);
            obj->setShine(shine);
        }
        //general object
        else
        {
            double generalCoefficients[10];
            Vector3D *referencePoint = new Vector3D();
            double length;
            double width;
            double height;
            
            for(int j=0;j<10;j++)
            {
                cin >> generalCoefficients[j];
            }
            cout<<"general object"<<endl;
            
            cin >> referencePoint->x >> referencePoint->y >> referencePoint->z >> length >> width >> height;
            cin >> color[0] >> color[1] >> color[2];
            cin >> ambient >> diffuse >> specular >> recursive_reflection_coefficient;
            cin >> shine;
            cout<<"reference points"<<referencePoint->x<<" "<<referencePoint->y<<" "<<referencePoint->z<<endl;

            obj = new Object(generalCoefficients,*referencePoint,length,width,height);
            obj->setColor(color[0],color[1],color[2]);
            obj->setCoEfficients(ambient, diffuse, specular, recursive_reflection_coefficient);
            obj->setShine(shine);
        }
        obj->setType(object_type);
        objects.push_back(obj);
    }

    cin >> totalPointLightss;
    for(int i=0;i<totalPointLightss;i++)
    {
        PointLight pointlight;
        Vector3D pos;
        double color[3];
        cin >> pos.x >> pos.y >> pos.z;
        cin >> color[0] >> color[1] >> color[2];
        cout<<"point light found"<<endl;
       
        pointlight.light_pos = pos;
        pointlight.color[0] = color[0];
        pointlight.color[1] = color[1];
        pointlight.color[2] = color[2];
        pointLights.push_back(pointlight);
    }

    cin >> totalSpotLights;
    for(int i=0;i<totalSpotLights;i++)
    {
        SpotLight spotlight;
        Vector3D pos;
        Vector3D dir;
        double color[3];
        double cutoff_angle;
        cin >> pos.x >> pos.y >> pos.z;
        cin >> color[0] >> color[1] >> color[2];
        cin >> dir.x >> dir.y >> dir.z;
        cin >> cutoff_angle;

        spotlight.cutoff_angle = cutoff_angle;
        spotlight.point_light.light_pos = pos;
        spotlight.point_light.color[0] = color[0];
        spotlight.point_light.color[1] = color[1];
        spotlight.point_light.color[2] = color[2];
        spotlight.light_direction = dir;
        cout<<"spot light found"<<endl;

        spotLights.push_back(spotlight);
    }

    // lasssstly push the floor
    Object* checkBoard = new Floor(1000,20);
    objects.push_back(checkBoard);


}

void printData()
{
    cout << "recursionLevel: " << recursionLevel << endl;
    cout << "pixels: " << pixels << endl;
    cout << "totalObjects: " << totalObjects << endl;
    for (auto & individual : objects)
    {
        cout << "object type: " << individual->getType() << endl;
        cout << "object reference point: " << individual->reference_point.x << " " << individual->reference_point.y << " " << individual->reference_point.z << endl;
        //cout << "object color: " << individual->color[0] << " " << individual->color[1] << " " << individual->color[2] << endl;
        //cout << "object coEfficients: " << individual->coEfficients[0] << " " << individual->coEfficients[1] << " " << individual->coEfficients[2] << " " << individual->coEfficients[3] << endl;
        //cout << "object shine: " << individual->shine << endl;
    }
    cout << "totalPointLightss: " << totalPointLightss << endl;
    for (auto & individual : pointLights)
    {
        cout << "pointLight position: " << individual.light_pos.x << " " << individual.light_pos.y << " " << individual.light_pos.z << endl;
        cout << "pointLight color: " << individual.color[0] << " " << individual.color[1] << " " << individual.color[2] << endl;
    }
    cout << "totalSpotLights: " << totalSpotLights << endl;
    for (auto & individual : spotLights)
    {
        cout << "spotLight position: " << individual.point_light.light_pos.x << " " << individual.point_light.light_pos.y << " " << individual.point_light.light_pos.z << endl;
        cout << "spotLight color: " << individual.point_light.color[0] << " " << individual.point_light.color[1] << " " << individual.point_light.color[2] << endl;
        cout << "spotLight direction: " << individual.light_direction.x << " " << individual.light_direction.y << " " << individual.light_direction.z << endl;
        cout << "spotLight cutoff_angle: " << individual.cutoff_angle << endl;
    }
}


void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    gluLookAt(eye.x,eye.y,eye.z,eye.x+forwardVector.x,eye.y+forwardVector.y,eye.z+forwardVector.z,upVector.x,upVector.y,upVector.z);	
    glPushMatrix();
    drawAllComponents();
    glPopMatrix();
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
	glEnable(GL_DEPTH_TEST);	//enable Depth Testing
    

	glutDisplayFunc(display);	
	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
    atexit(freeMemory);
	glutMainLoop();	

	return 0;
}


/**************** Run This Script
g++ 1905053_main.cpp -o me3 -lGL -lGLU -lglut
./me3
*/