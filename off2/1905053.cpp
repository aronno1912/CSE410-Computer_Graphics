#include <bits/stdc++.h>
#include <iostream>
#include <stack>
#include <array>
#include <cmath>
#include "bitmap_image.hpp"
using namespace std;

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif
#define PI acos(-1.0)
#define degToRad(x) (x * PI / 180.0)

struct Point
{
    double x, y, z;
};
static unsigned long int g_seed = 1;

inline int randomNum()
{
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}

bitmap_image image;

struct Vector3D
{
    double x, y, z;

    // Normalize the vector
    void normalize()
    {
        double length = std::sqrt(x * x + y * y + z * z);
        x /= length;
        y /= length;
        z /= length;
    }

    // Scalar multiplication
    Vector3D operator*(double scalar) const
    {
        Vector3D result;
        result.x = x * scalar;
        result.y = y * scalar;
        result.z = z * scalar;
        return result;
    }

    // Vector addition
    Vector3D operator+(const Vector3D &other) const
    {
        Vector3D result;
        result.x = x + other.x;
        result.y = y + other.y;
        result.z = z + other.z;
        return result;
    }
    // Vector subtraction
    Vector3D operator-(const Vector3D &other) const
    {
        Vector3D result;
        result.x = x - other.x;
        result.y = y - other.y;
        result.z = z - other.z;
        return result;
    }
};

// Function to compute the dot product of two vectors
double dotProduct(const Vector3D &a, const Vector3D &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Function to compute the cross product of two vectors
Vector3D crossProduct(const Vector3D &a, const Vector3D &b)
{
    Vector3D result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}
// Function to perform the vector rotation using Rodrigues' formula
// 𝑅(𝑥⃗, 𝑎⃗, 𝜃) = 𝑐𝑜𝑠𝜃 𝑥⃗ + (1 − 𝑐𝑜𝑠𝜃) (𝑎⃗̇. 𝑥⃗)𝑎⃗ + 𝑠𝑖𝑛𝜃 (𝑎⃗ × 𝑥⃗)

Vector3D rotateVector(const Vector3D &x, const Vector3D &a, double angle)
{
    double cosTheta = std::cos(degToRad(angle));
    double sinTheta = std::sin(degToRad(angle));

    Vector3D term1 = x * cosTheta;
    Vector3D term2 = a * (1 - cosTheta) * dotProduct(a, x);
    Vector3D term3 = crossProduct(a, x) * sinTheta;

    return term1 + term2 + term3;
}

// Function to transform a 3D point using a 4x4 matrix
Point transformPoint(const double matrix[4][4], const Point &p)
{
    Point result;
    double w = matrix[3][0] * p.x + matrix[3][1] * p.y + matrix[3][2] * p.z + matrix[3][3] * 1;
    result.x = (matrix[0][0] * p.x + matrix[0][1] * p.y + matrix[0][2] * p.z + matrix[0][3]) / w;
    result.y = (matrix[1][0] * p.x + matrix[1][1] * p.y + matrix[1][2] * p.z + matrix[1][3]) / w;
    result.z = (matrix[2][0] * p.x + matrix[2][1] * p.y + matrix[2][2] * p.z + matrix[2][3]) / w;
    return result;
}

// Function to multiply two 4x4 matrices
void matrixMultiply(double result[4][4], const double matrix1[4][4], const double matrix2[4][4])
{
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            result[i][j] = 0.0;
            for (int k = 0; k < 4; ++k)
            {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}
class Triangle
{

public:
    Point p1, p2, p3;
    int R, G, B; // Color values
    double a, b, c, d;

    void normalToPlane()
    {
        Vector3D v1, v2;
        v1.x = p2.x - p1.x;
        v1.y = p2.y - p1.y;
        v1.z = p2.z - p1.z;

        v2.x = p3.x - p1.x;
        v2.y = p3.y - p1.y;
        v2.z = p3.z - p1.z;

        Vector3D normal = crossProduct(v1, v2);
        cout << "normal is " << normal.x << " " << normal.y << " " << normal.z << endl;
        normal.normalize();
        Point p;
        p.x = normal.x;
        p.y = normal.y;
        p.z = normal.z;

        a = p.x;
        b = p.y;
        c = p.z;
        d = -(a * p1.x + b * p1.y + c * p1.z);
    }

    Triangle(const Point &p1, const Point &p2, const Point &p3)
    {

        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
        R = randomNum() % 256;
        G = randomNum() % 256;
        B = randomNum() % 256;
        normalToPlane();
    }
    void printPoints() const
    {
        std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl;
        std::cout << p2.x << " " << p2.y << " " << p2.z << std::endl;
        std::cout << p3.x << " " << p3.y << " " << p3.z << std::endl;
    }
};

// void printMatrix(const double matrix[rows][cols]) {
//     for (int i = 0; i < rows; ++i) {
//         for (int j = 0; j < cols; ++j) {
//             std::cout << matrix[i][j] << "\t";
//         }
//         std::cout << std::endl;
//     }
// }

double top_y, left_x, bottom_y, right_x;
double dx, dy;
int screen_width, screen_height;

// Function to clip the top of a triangle
int clipTop(const Triangle &triangle)
{
    double max_y = std::max({triangle.p1.y, triangle.p2.y, triangle.p3.y});

    // cout<<"max y is "<<max_y<<endl;
    // cout<<"top y is "<<top_y<<endl;

    // If max_y is above Top_Y, clip the portion above Top_Y
    if (max_y > top_y)
    {
        cout << "clipped top" << endl;
        return ceil(top_y / dy);
    }

    // Otherwise, no clipping needed
    return ceil(max_y / dy);
}

// Function to clip the bottom of a triangle
int clipBottom(const Triangle &triangle)
{
    double min_y = std::min({triangle.p1.y, triangle.p2.y, triangle.p3.y});

    // If min_y is below Bottom_Y, clip the portion below Bottom_Y
    if (min_y < bottom_y)
    {
        cout << "clipped bottom" << endl;
        return floor(bottom_y / dy);
    }

    // Otherwise, no clipping needed
    return floor(min_y / dy);
}


int getLeftIntersectionCol(const Triangle &triangle, float y)
{
    double minX = 99999;
    // For each edge of the triangle, it checks if the horizontal line intersects that edge.If an intersection is found, it calculates the x-coordinate of the intersection using linear interpolation.

    // Edge 1
    double x1 = triangle.p1.x;
    double y1 = triangle.p1.y;
    double x2 = triangle.p2.x;
    double y2 = triangle.p2.y;
 
    // y1 has the smaller y-coordinate, swap if necessary
    if (y1 > y2)
    {
        swap(x1, x2);
        swap(y1, y2);
    }
    // check if the horizontal line intersects the edge. If yes, calculate the x-coordinate of the intersection using linear interpolation.
    if (y >= y1 && y <= y2)
    {
        double x = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
       
        minX = min(minX, x);
    }

    // Edge 2
    x1 = triangle.p2.x;
    y1 = triangle.p2.y;
    x2 = triangle.p3.x;
    y2 = triangle.p3.y;
    if (y1 > y2)
    {
        swap(x1, x2);
        swap(y1, y2);
    }
    if (y >= y1 && y <= y2)
    {
        double x = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
        minX = min(minX, x);
    }

    // Edge 3
    x1 = triangle.p3.x;
    y1 = triangle.p3.y;
    x2 = triangle.p1.x;
    y2 = triangle.p1.y;
    if (y1 > y2)
    {
        swap(x1, x2);
        swap(y1, y2);
    }
    if (y >= y1 && y <= y2)
    {
        double x = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
        minX = min(minX, x);
    }

    // clip if necessary
    int leftIntCol = round(max(minX / dx, left_x / dx));
    return leftIntCol;
}

int getRightIntersectionCol(const Triangle &triangle, float y)
{
    double maxX = -99999;
    // For each edge of the triangle, it checks if the horizontal line intersects that edge.If an intersection is found, it calculates the x-coordinate of the intersection using linear interpolation.

    // Edge 1
    double x1 = triangle.p1.x;
    double y1 = triangle.p1.y;
    double x2 = triangle.p2.x;
    double y2 = triangle.p2.y;
    // y1 has the smaller y-coordinate, swap if necessary
    if (y1 > y2)
    {
        swap(x1, x2);
        swap(y1, y2);
    }
    // check if the horizontal line intersects the edge. If yes, calculate the x-coordinate of the intersection using linear interpolation.
    if (y >= y1 && y <= y2)
    {
        double x = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
        maxX = max(maxX, x);
    }

    // Edge 2
    x1 = triangle.p2.x;
    y1 = triangle.p2.y;
    x2 = triangle.p3.x;
    y2 = triangle.p3.y;
    if (y1 > y2)
    {
        swap(x1, x2);
        swap(y1, y2);
    }
    if (y >= y1 && y <= y2)
    {
        double x = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
        maxX = max(maxX, x);
    }

    // Edge 3
    x1 = triangle.p3.x;
    y1 = triangle.p3.y;
    x2 = triangle.p1.x;
    y2 = triangle.p1.y;
    if (y1 > y2)
    {
        swap(x1, x2);
        swap(y1, y2);
    }
    if (y >= y1 && y <= y2)
    {
        double x = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
        maxX = max(maxX, x);
    }

    // clip if necessary
    int rightIntCol = round(min(maxX / dx, right_x / dx));
    return rightIntCol;
}



double calculateZValue(double a, double b, double c, double d, double x, double y)
{
    return -(a * x + b * y + d) / c;
}

// Function to update z-buffer 

void processTriangle(const Triangle &triangle, vector<vector<double>> &zBuffer, double z_front_limit)
{

    // Clip the triangle and find top and bottom scanlines
    int top_scanline = clipTop(triangle);
    cout << "top scanline is " << top_scanline << endl;
    int bottom_scanline = clipBottom(triangle);
    cout << "bottom scanline is " << bottom_scanline << endl;

    // Iterate over each row from top to bottom
    for (int row_no = top_scanline; row_no >= bottom_scanline; row_no--)
    {
        // Clip the triangle and find left and right intersecting columns

        int left_intersecting_column = getLeftIntersectionCol(triangle, row_no * dy);
        int right_intersecting_column = getRightIntersectionCol(triangle, row_no * dy);

        // Iterate over each column from left to right
        for (int col_no = left_intersecting_column; col_no <= right_intersecting_column; col_no++)
        {
            double x_val = col_no * dx;
            double y_val = row_no * dy;

            // Calculate z value for the current pixel
            double z_value = calculateZValue(triangle.a, triangle.b, triangle.c, triangle.d, x_val, y_val);

            // During scanning from top to bottom and left to right, check for the middle values of each cell.
            int row = top_y / dy - row_no;
            int col = col_no - left_x / dx;

            // Check if the point is in front of z_front_limit
            if (z_value < zBuffer[row][col] && z_value > z_front_limit)
            {
                // Compare with z-buffer and update if required
                zBuffer[row][col] = z_value;
                image.set_pixel(col, row, triangle.R, triangle.G, triangle.B); //instead of using frame buffer we can directly use image.set_pixel
            }
        }
    }
}

// Function to render all triangles
void renderTriangles(const vector<Triangle> &triangles, vector<vector<double>> &zBuffer, double z_front_limit)
{
    for (const Triangle &triangle : triangles)
    {
        processTriangle(triangle, zBuffer, z_front_limit);
        // output<<"one done"<<endl;
        // cout<<"a,b,c,d is "<<triangle.a<<" "<<triangle.b<<" "<<triangle.c<<" "<<triangle.d<<endl;
    }
}

int main()
{
    int count = 0;
    ifstream input("scene.txt");
    ofstream output("stage1.txt");

    output << fixed << setprecision(7);
    // stack<double [4][4]> S;
    // stack declaration to store pointers to double[4][4]
    stack<double(*)[4]> S;
    // Initialize M as the identity matrix
    double M[4][4] = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0}};

    double temp[4][4];
    memcpy(temp, M, sizeof(temp));
    double fovY, aspectRatio, near, far;
    Vector3D eye, look, up;

    input >> eye.x >> eye.y >> eye.z >> look.x >> look.y >> look.z >> up.x >> up.y >> up.z >> fovY >> aspectRatio >> near >> far;
    // double fovX = fovY * aspectRatio;
    cout << fovY << " " << aspectRatio << " " << near << " " << far << endl;
    string word = "";

    while (true)
    {
        input >> word;

        if (word == "triangle")
        {
            count++;
            Point p1, p2, p3;
            input >> p1.x >> p1.y >> p1.z >> p2.x >> p2.y >> p2.z >> p3.x >> p3.y >> p3.z;

            Point transformedPoint = transformPoint(M, p1);
            output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
            transformedPoint = transformPoint(M, p2);
            output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
            transformedPoint = transformPoint(M, p3);
            output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
            output << endl;
            // output<<"ekta korsiiiiiiiiiiii"<<endl;
        }
        else if (word == "translate")
        {
            double tx, ty, tz;
            input >> tx >> ty >> tz;

            // Generate the corresponding translation matrix T
            double T[4][4] = {
                {1.0, 0.0, 0.0, tx},
                {0.0, 1.0, 0.0, ty},
                {0.0, 0.0, 1.0, tz},
                {0.0, 0.0, 0.0, 1.0}};
            // Update M by multiplying with T
            double result[4][4];
            matrixMultiply(result, M, T);
            // Copy the result to M
            memcpy(M, result, sizeof(result));
        }
        else if (word == "scale")
        {
            double sx, sy, sz;
            input >> sx >> sy >> sz;

            // Generate the corresponding scaling matrix T
            double T[4][4] = {
                {sx, 0.0, 0.0, 0.0},
                {0.0, sy, 0.0, 0.0},
                {0.0, 0.0, sz, 0.0},
                {0.0, 0.0, 0.0, 1.0}};
            // Update M by multiplying with T
            double result[4][4];
            matrixMultiply(result, M, T);
            // Copy the result to M
            memcpy(M, result, sizeof(result));
        }
        else if (word == "rotate")
        {
            double angle, ax, ay, az;
            input >> angle >> ax >> ay >> az;
            Vector3D a = {ax, ay, az};
            a.normalize();

            // Generate the corresponding rotation matrix T

            // Calculate the rotation matrix columns
            Vector3D c1 = rotateVector({1, 0, 0}, a, angle);
            Vector3D c2 = rotateVector({0, 1, 0}, a, angle);
            Vector3D c3 = rotateVector({0, 0, 1}, a, angle);
            double T[4][4] = {
                {c1.x, c2.x, c3.x, 0.0},
                {c1.y, c2.y, c3.y, 0.0},
                {c1.z, c2.z, c3.z, 0.0},
                {0.0, 0.0, 0.0, 1.0}};

            // Update M by multiplying with T
            double result[4][4];
            matrixMultiply(result, M, T);
            // Copy the result to M

            memcpy(M, result, sizeof(result));
        }
        else if (word == "push")
        {

            double(*temp)[4] = new double[4][4];
            memcpy(temp, M, sizeof(double) * 16);
            S.push(temp);
        }
        else if (word == "pop")
        {

            if (!S.empty())
            {
                double(*temp)[4] = S.top();
                memcpy(M, temp, sizeof(double) * 16);
                S.pop();
                delete[] temp;
            }
        }
        else if (word == "end")
        {

            break;
        }
    }
    input.close();
    output.close();

    input.open("stage1.txt");
    output.open("stage2.txt");

    // stage 2
    Vector3D l, r, u;
    l = look - eye;
    l.normalize();
    r = crossProduct(l, up);
    r.normalize();
    u = crossProduct(r, l);
    u.normalize();
    double T[4][4] = {
        {1.0, 0.0, 0.0, -eye.x},
        {0.0, 1.0, 0.0, -eye.y},
        {0.0, 0.0, 1.0, -eye.z},
        {0.0, 0.0, 0.0, 1.0}};

    double R[4][4] = {
        {r.x, r.y, r.z, 0.0},
        {u.x, u.y, u.z, 0.0},
        {-l.x, -l.y, -l.z, 0.0},
        {0.0, 0.0, 0.0, 1.0}};
    double V[4][4];
    matrixMultiply(V, R, T);
    for (int i = 0; i < count; i++)
    {
        Point p1, p2, p3; // 3 points of a triangle
        input >> p1.x >> p1.y >> p1.z >> p2.x >> p2.y >> p2.z >> p3.x >> p3.y >> p3.z;
        Point transformedPoint = transformPoint(V, p1);
        output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
        transformedPoint = transformPoint(V, p2);
        output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
        transformedPoint = transformPoint(V, p3);
        output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
        output << endl;
    }
    input.close();
    output.close();

    // Stage 3: Projection Transformation
    input.open("stage2.txt");
    output.open("stage3.txt");
    double fovX = fovY * aspectRatio;
    double t = near * tan(degToRad(fovY / 2));
    double r1 = near * tan(degToRad(fovX / 2));

    double P[4][4] = {
        {near / r1, 0.0, 0.0, 0.0},
        {0.0, near / t, 0.0, 0.0},
        {0.0, 0.0, -(far + near) / (far - near), -(2 * far * near) / (far - near)},
        {0.0, 0.0, -1.0, 0.0}};

    for (int i = 0; i < count; i++)
    {
        Point p1, p2, p3; // 3 points of a triangle
        input >> p1.x >> p1.y >> p1.z >> p2.x >> p2.y >> p2.z >> p3.x >> p3.y >> p3.z;
        Point transformedPoint = transformPoint(P, p1);
        output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
        transformedPoint = transformPoint(P, p2);
        output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
        transformedPoint = transformPoint(P, p3);
        output << transformedPoint.x << " " << transformedPoint.y << " " << transformedPoint.z << std::endl;
        output << endl;
    }
    input.close();
    output.close();

    //************************Stage 4:  Clipping & scan conversion using Z-buffer algorithm*****************************
    input.open("config.txt");
    input >> screen_width >> screen_height;
    input.close();
    input.open("stage3.txt");
    output.open("z_buffer.txt");

    // read the triangles from the file
    vector<Triangle> triangles;
    int trk = 0;
    while (trk != count)
    {
        Point p1, p2, p3; // 3 points of a triangle
        input >> p1.x >> p1.y >> p1.z >> p2.x >> p2.y >> p2.z >> p3.x >> p3.y >> p3.z;
        Triangle triangle(p1, p2, p3); // a,b,c,d calculated already!!!!!!!!!!!!!!!!!!!!!!!!
        // triangle.printPoints();
        // cout<<"ekta korsi"<<endl;
        cout << endl;
        triangles.push_back(triangle);
        trk++;
    }

    cout << "size is " << triangles.size() << endl;

    // Initialize z-buffer and frame buffer

    // Create a z-buffer, a two dimensional array of Screen_Width X Screen_Height dimension.
    // Initialize all the values in z-buffer with z_max. In the aforementioned examples, z_max = 2.0.
    // The memory for z-buffer should be dynamically allocated (using STL is allowed).

    vector<vector<double>> zBuffer(screen_height, vector<double>(screen_width, 1.0));

    double z_front_limit = -1.0;
    double z_max = 1.0;

    // Create a pixel mapping between the x-y range values and the Screen_Width X Screen_height range.

    double x_min = -1.0, x_max = 1.0, y_min = -1.0, y_max = 1.0;
    dx = (x_max - x_min) / screen_width;
    dy = (y_max - y_min) / screen_height;
    top_y = y_max - dy / 2;
    left_x = x_min + dx / 2;
    bottom_y = y_min + dy / 2;
    right_x = x_max - dx / 2;

    // Create a bitmap_image object with Screen_Width X Screen_Height resolution and initialize its background color with black.
    image.setwidth_height(screen_width, screen_height);
    image.set_all_channels(0, 0, 0);

    renderTriangles(triangles, zBuffer, z_front_limit);
    // output zBuffer values
    for (int i = 0; i < screen_height; i++)
    {
        for (int j = 0; j < screen_width; j++)
        {
            if (zBuffer[i][j] < z_max)
                output << zBuffer[i][j] << " ";
        }
        output << endl;
    }
    output.close();
    image.save_image("out.bmp");

    // free the memory allocated for z-buffer
    zBuffer.clear();
    image.clear();

    return 0;
}

//******************** Run this script *********************************************************

/***
 * 
 * 
 * g++ 1905053.cpp -o demonew -lGL -lGLU -lglut
   ./demonew
*/