#ifndef CAMERA
#define CAMERA

#include "Vector3D.h"

class Camera
{
public:
    Vector3D *position = new Vector3D();    // where the camera is
    Vector3D *look = new Vector3D();        // where the camera is facing
    Vector3D *up = new Vector3D();          // the up vector
    Vector3D *right = new Vector3D();       // the right vector
    Vector3D *center  = new Vector3D();     // currently where the camera is looking at

    double linear_speed = .3, angular_speed = .01;

    Camera(Vector3D *position, Vector3D *center, Vector3D *up){
        this->position = position;
        this->center = center;
        this->up = up;

        *this->look = *this->position - *this->center;
        *this->right = (*this->up) * (*this->look);         // up, look and right are perpendicular to each other
        *this->right->normalize_vector();
    }

    // assume that the direction is normalized
    // need to update the look and the position
    void move_camera(Vector3D * direction) { 
        *this->position = (*this->position + *direction);
        *this->center = (*this->center + *direction);
    }

    //The logic involves \
    subtracting the current position (pos) \
    from the target, normalizing the result, \
    multiplying it by the magnitude of the \
    original direction, and then adding the \
    position back

    void rotate_camera(Vector3D* direction){
        // *this->center = (((*this->center - *this->position)
        // ->normalize_vector()) * *direction) + *this->position;  
    }

    void move_forward(){
        Vector3D direction = *this->look * -linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    void move_backward(){
        Vector3D direction = *this->look * linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    void move_right(){
        Vector3D direction = *this->right * (-linear_speed);
        direction.normalize_vector();
        move_camera(&direction);
    }

    void move_left(){
        Vector3D direction = *this->right * (linear_speed);
        direction.normalize_vector();
        move_camera(&direction);
    }

    void move_up(){
        Vector3D direction = *this->up * linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    void move_down(){
        Vector3D direction = *this->up * -linear_speed;
        direction.normalize_vector();
        move_camera(&direction);
    }

    void look_up(){

    }





    ~Camera(){
        delete position;
        delete look;
        delete up;
        delete right;
        delete center;
    }
    

};




#endif
