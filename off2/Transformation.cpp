#ifndef TRANSFORMATION_CPP
#define TRANSFORMATION_CPP

#include <vector>
#include <cmath>
#include "Point.cpp"

#define PI acos(-1.0)
#define degToRad(x) (x * PI / 180.0)

using namespace std;

class Transformation {
public:
    vector<vector<double>> matrix;

    Transformation(vector<vector<double>> matrix) : matrix(matrix) {}
    
    // by default transformation matrix is identity matrix
    Transformation() {
        matrix = vector<vector<double>>(4, vector<double>(4, 0));
        matrix[0][0] = 1;
        matrix[1][1] = 1;
        matrix[2][2] = 1;
        matrix[3][3] = 1;
    }

    static Transformation translationMatrix(double tx, double ty, double tz) {
        Transformation T;
        T.matrix[0][3] = tx;
        T.matrix[1][3] = ty;
        T.matrix[2][3] = tz;
        return T;
    }

    static Transformation scalingMatrix(double sx, double sy, double sz) {
        Transformation T;
        T.matrix[0][0] = sx;
        T.matrix[1][1] = sy;
        T.matrix[2][2] = sz;
        return T;
    }

    static Point normalize(Point p) {
        double norm = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        return Point(p.x / norm, p.y / norm, p.z / norm);
    }

    static Transformation rotationMatrix(double angle, double ax, double ay, double az) {
        double theta = degToRad(angle);
        double c = cos(theta);
        double s = sin(theta);
        double t = 1 - c;

        Point p = Point(ax, ay, az);
        p = normalize(p);
        ax = p.x;
        ay = p.y;
        az = p.z;


        Transformation T;
        T.matrix[0][0] = t * ax * ax + c;
        T.matrix[0][1] = t * ax * ay - s * az;
        T.matrix[0][2] = t * ax * az + s * ay;

        T.matrix[1][0] = t * ax * ay + s * az;
        T.matrix[1][1] = t * ay * ay + c;
        T.matrix[1][2] = t * ay * az - s * ax;

        T.matrix[2][0] = t * ax * az - s * ay;
        T.matrix[2][1] = t * ay * az + s * ax;
        T.matrix[2][2] = t * az * az + c;

        return T;
    }

    static Transformation viewMatrix(Point eye, Point look, Point up) {
        Point l = normalize(look - eye);
        Point r = normalize(l.cross(up));
        Point u = normalize(r.cross(l));

        Transformation T = translationMatrix(-eye.x, -eye.y, -eye.z);

        Transformation R = Transformation();
        R.matrix[0][0] = r.x;
        R.matrix[0][1] = r.y;
        R.matrix[0][2] = r.z;

        R.matrix[1][0] = u.x;
        R.matrix[1][1] = u.y;
        R.matrix[1][2] = u.z;

        R.matrix[2][0] = -l.x;
        R.matrix[2][1] = -l.y;
        R.matrix[2][2] = -l.z;

        Transformation V = R * T;
        return V;
    }

    static Transformation projectionMatrix(double fovY, double aspectRatio, double near, double far) {
        double fovX = fovY * aspectRatio;
        double t = near * tan(degToRad(fovY) / 2);
        double r = near * tan(degToRad(fovX) / 2);

        Transformation P = Transformation(vector<vector<double>>(4, vector<double>(4, 0)));
        P.matrix[0][0] = near / r;
        P.matrix[1][1] = near / t;
        P.matrix[2][2] = -(far + near) / (far - near);
        P.matrix[2][3] = -(2 * far * near) / (far - near);
        P.matrix[3][2] = -1;

        return P;
    }


    Point transform(Point p) {
        double x = matrix[0][0] * p.x + matrix[0][1] * p.y + matrix[0][2] * p.z + matrix[0][3] * p.w;
        double y = matrix[1][0] * p.x + matrix[1][1] * p.y + matrix[1][2] * p.z + matrix[1][3] * p.w;
        double z = matrix[2][0] * p.x + matrix[2][1] * p.y + matrix[2][2] * p.z + matrix[2][3] * p.w;
        double w = matrix[3][0] * p.x + matrix[3][1] * p.y + matrix[3][2] * p.z + matrix[3][3] * p.w;


        return (Point(x, y, z))/w;
    }

    Transformation operator*(const Transformation &t) const {
        Transformation M = Transformation(vector<vector<double>>(4, vector<double>(4, 0)));
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                for(int k=0; k<4; k++){
                    M.matrix[i][j] += matrix[i][k] * t.matrix[k][j];
                }
            }
        }
        return M;
    }

    void print() {
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                cout << matrix[i][j] << " ";
            }
            cout << endl;
        }
    }


};   

#endif