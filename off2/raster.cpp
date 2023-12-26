#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include<array>
#include <cmath>
using namespace std;

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif
#define PI acos(-1.0)
#define degToRad(x) (x * PI / 180.0)
constexpr int rows = 4;
constexpr int cols = 4;
using Matrix = array<array<double, cols>, rows>;
struct Point {
    double x, y, z;
};

struct Vector3D {
    double x, y, z;

    // Normalize the vector
    void normalize() {
        double length = std::sqrt(x * x + y * y + z * z);
        x /= length;
        y /= length;
        z /= length;
    }

       // Scalar multiplication
    Vector3D operator*(double scalar) const {
        Vector3D result;
        result.x = x * scalar;
        result.y = y * scalar;
        result.z = z * scalar;
        return result;
    }

      // Vector addition
    Vector3D operator+(const Vector3D& other) const {
        Vector3D result;
        result.x = x + other.x;
        result.y = y + other.y;
        result.z = z + other.z;
        return result;
    }
};

// Function to compute the dot product of two vectors
double dotProduct(const Vector3D& a, const Vector3D& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Function to compute the cross product of two vectors
Vector3D crossProduct(const Vector3D& a, const Vector3D& b) {
    Vector3D result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}
// Function to perform the vector rotation using Rodrigues' formula
//𝑅(𝑥⃗, 𝑎⃗, 𝜃) = 𝑐𝑜𝑠𝜃 𝑥⃗ + (1 − 𝑐𝑜𝑠𝜃) (𝑎⃗̇. 𝑥⃗)𝑎⃗ + 𝑠𝑖𝑛𝜃 (𝑎⃗ × 𝑥⃗)

Vector3D rotateVector(const Vector3D& x, const Vector3D& a, double angle) {
    double cosTheta = std::cos(degToRad(angle));
    double sinTheta = std::sin(degToRad(angle));

    Vector3D term1 = x*cosTheta;
    Vector3D term2 = a * (1 - cosTheta) * dotProduct(a, x);
    Vector3D term3 = crossProduct(a, x) * sinTheta;

    return term1 + term2 + term3;
}

// Function to transform a 3D point using a 4x4 matrix
Point transformPoint(const double matrix[4][4], const Point& p) {
    Point result;
    double w   =  matrix[3][0] * p.x + matrix[3][1] * p.y + matrix[3][2] * p.z + matrix[3][3] * 1;
    result.x = (matrix[0][0] * p.x + matrix[0][1] * p.y + matrix[0][2] * p.z + matrix[0][3])/w;
    result.y = (matrix[1][0] * p.x + matrix[1][1] * p.y + matrix[1][2] * p.z + matrix[1][3])/w;
    result.z = (matrix[2][0] * p.x + matrix[2][1] * p.y + matrix[2][2] * p.z + matrix[2][3])/w;
    return result;
}

// Function to multiply two 4x4 matrices
void matrixMultiply(double result[4][4], const double matrix1[4][4], const double matrix2[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = 0.0;
            for (int k = 0; k < 4; ++k) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}

void printMatrix(const double matrix[rows][cols]) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

int main()
{   
    ifstream input("scene.txt");
    ofstream output("stage1.txt");
    output << fixed << setprecision(7);
    //stack<double [4][4]> S;
    //stack declaration to store pointers to double[4][4]
    stack<double (*)[4]> S;
    // Initialize M as the identity matrix
    double M[4][4] = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0}
    };

    double temp[4][4];
    memcpy(temp, M, sizeof(temp));
    double eyeX, eyeY, eyeZ, lookX, lookY, lookZ, upX, upY, upZ, fovY, aspectRatio, near, far;
    input >> eyeX >> eyeY >> eyeZ >> lookX >> lookY >> lookZ >> upX >> upY >> upZ >> fovY >> aspectRatio >> near >> far;
    //double fovX = fovY * aspectRatio;
    cout << eyeX << " " << eyeY << " " << eyeZ << endl;
    cout << lookX << " " << lookY << " " << lookZ << endl;
    cout << upX << " " << upY << " " << upZ << endl;
    cout << fovY << " " << aspectRatio << " " << near << " " << far << endl;
    string word="";

    while(true)
    {
        input >> word;

        if(word == "triangle")
        {   
            // double x1, y1, z1, x2, y2, z2, x3, y3, z3;
            // input >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
            // cout << x1 << " " << y1 << " " << z1 << endl;
            // cout << x2 << " " << y2 << " " << z2 << endl;
            // cout << x3 << " " << y3 << " " << z3 << endl;
            Point p1, p2, p3;
            input >> p1.x >> p1.y >> p1.z >> p2.x >> p2.y >> p2.z >> p3.x >> p3.y >> p3.z;
   
                Point transformedPoint = transformPoint(M, p1);
                output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
                transformedPoint = transformPoint(M, p2);
                output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
                transformedPoint = transformPoint(M, p3);
                output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
                output<<endl;
                //output<<"ekta korsiiiiiiiiiiii"<<endl;

            
        }
        else if(word == "translate")
        {
            double tx, ty, tz;
            input >> tx >> ty >> tz;
            cout<<"translate er moddhe"<<endl;  
            cout << tx << " " << ty << " " << tz << endl;
             // Generate the corresponding translation matrix T
            double T[4][4] = {
                {1.0, 0.0, 0.0, tx},
                {0.0, 1.0, 0.0, ty},
                {0.0, 0.0, 1.0, tz},
                {0.0, 0.0, 0.0, 1.0}
            };
          // Update M by multiplying with T
            double result[4][4];
            matrixMultiply(result, M, T);
            // Copy the result to M
            memcpy(M, result, sizeof(result));
            cout<<"after translate"<<endl;
            printMatrix(M);
        }
        else if(word == "scale")
        {
            double sx, sy, sz;
            input >> sx >> sy >> sz;
            cout << sx << " " << sy << " " << sz << endl;
            // Generate the corresponding scaling matrix T
            double T[4][4] = {
                {sx, 0.0, 0.0, 0.0},
                {0.0, sy, 0.0, 0.0},
                {0.0, 0.0, sz, 0.0},
                {0.0, 0.0, 0.0, 1.0}
            };
            // Update M by multiplying with T
            double result[4][4];
            matrixMultiply(result, M, T);
            // Copy the result to M
            memcpy(M, result, sizeof(result));
            cout<<"after scale"<<endl;
            printMatrix(M);
        }
        else if(word == "rotate")
        {
            double angle, ax, ay, az;
            input >> angle >> ax >> ay >> az;
            cout<<"rotate er moddhe"<<endl;
            cout << angle << " " << ax << " " << ay << " " << az << endl;
            Vector3D a={ax,ay,az};
            a.normalize();
            cout<<"normalized vector"<<endl;
            cout<<a.x<<" "<<a.y<<" "<<a.z<<endl;
            // Generate the corresponding rotation matrix T
             
              // Calculate the rotation matrix columns
            Vector3D c1 = rotateVector({1, 0, 0}, a, angle);
            Vector3D c2 = rotateVector({0, 1, 0}, a, angle);
            Vector3D c3 = rotateVector({0, 0, 1}, a, angle);
            double T[4][4] = {
                {c1.x, c2.x, c3.x, 0.0},
                {c1.y, c2.y, c3.y, 0.0},
                {c1.z, c2.z, c3.z, 0.0},
                {0.0, 0.0, 0.0, 1.0}
            };
            cout<<"T matrix print kori"<<endl;
            printMatrix(T);
            // Update M by multiplying with T
            double result[4][4];
            matrixMultiply(result, M, T);
            // Copy the result to M
            cout<<"ekhane"<<endl;
            memcpy(M, result, sizeof(result));
            cout<<"matrix print kori"<<endl;
            printMatrix(M);



        }
        else if(word == "push")
        {
             cout << "push" << endl;
            double (*temp)[4] = new double[4][4];
             memcpy(temp, M, sizeof(double) * 16);
            S.push(temp);
        }
        else if(word == "pop")
        {
            cout << "pop" << endl;
            if (!S.empty())
            {
            double (*temp)[4] = S.top();
            memcpy(M, temp, sizeof(double) * 16);
            S.pop();
            delete[] temp;
         }
        
        }
        else if(word == "end")
        {
            cout << "end" << endl;
            break;
        }

    }

    return 0;
}