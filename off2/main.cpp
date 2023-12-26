#include <bits/stdc++.h>
#include "Point.cpp"
#include "Transformation.cpp"

using namespace std;

int main() {

    // Stage 1
    // Modeling Transformation

    ifstream input("scene.txt");
    ofstream output("stage1.txt");

    Point eye;
    Point look;
    Point up;
    double fovY;
    double aspectRatio;
    double near;
    double far;

    input >> eye.x >> eye.y >> eye.z;
    input >> look.x >> look.y >> look.z;
    input >> up.x >> up.y >> up.z;
    input >> fovY >> aspectRatio >> near >> far;

    output << fixed << setprecision(7);

    stack<Transformation> S;
    Transformation M;
    S.push(M);

    string command = "";
    int lineNumber = 4;
    int triangleCount = 0;

    while(true){
        input >> command;
        lineNumber++;

        if(command=="triangle"){
            triangleCount++;
            Point v1, v2, v3;
            input >> v1.x >> v1.y >> v1.z;
            input >> v2.x >> v2.y >> v2.z;
            input >> v3.x >> v3.y >> v3.z;
            lineNumber += 3;

            v1 = S.top().transform(v1);
            v2 = S.top().transform(v2);
            v3 = S.top().transform(v3);

            // stage1out << "triangle" << endl;
            output << v1.x << " " << v1.y << " " << v1.z << endl;
            output << v2.x << " " << v2.y << " " << v2.z << endl;
            output << v3.x << " " << v3.y << " " << v3.z << endl;
            output << endl;
        } else if(command=="translate"){
            double tx, ty, tz;
            input >> tx >> ty >> tz;
            lineNumber++;

            Transformation T = Transformation::translationMatrix(tx, ty, tz);
            M = S.top() * T;
            S.pop();
            S.push(M);
            cout<<"After translate M is"<<endl;
            M.print();
        } else if(command=="scale"){
            double sx, sy, sz;
            input >> sx >> sy >> sz;
            lineNumber++;

            Transformation T = Transformation::scalingMatrix(sx, sy, sz);
            M = S.top() * T;
            S.pop();
            S.push(M);
            cout<<"After scale M is"<<endl;
            M.print();
        } else if(command=="rotate"){
            double angle, ax, ay, az;
            input >> angle >> ax >> ay >> az;
            lineNumber++;

            Transformation T = Transformation::rotationMatrix(angle, ax, ay, az);
            cout << "Rotation Matrix:" << endl;
            T.print();

            M = S.top() * T;
            S.pop();
            cout<<"Now M is"<<endl;
            M.print();
            S.push(M);
        } else if(command=="push"){
            M = S.top();
            S.push(M);
            // S.top().print();
        } else if(command=="pop"){
            S.pop();
        } else if(command=="end"){
            break;
        } else {
            cout << "Invalid command at line " << lineNumber << endl;
            break;
        }
    }
    input.close();
    output.close();

    // Stage 2
    // View Transformation

    input.open("stage1.txt");
    output.open("stage2.txt");

    Transformation V = Transformation::viewMatrix(eye, look, up);
    // cout << "View Matrix:" << endl;
    // V.print();
    for (int i = 0; i < triangleCount; i++) {
        Point v1, v2, v3;
        input >> v1.x >> v1.y >> v1.z;
        input >> v2.x >> v2.y >> v2.z;
        input >> v3.x >> v3.y >> v3.z;

        v1 = V.transform(v1);
        v2 = V.transform(v2);
        v3 = V.transform(v3);

        output << v1.x << " " << v1.y << " " << v1.z << endl;
        output << v2.x << " " << v2.y << " " << v2.z << endl;
        output << v3.x << " " << v3.y << " " << v3.z << endl;
        output << endl;
    }

    input.close();
    output.close();

    // Stage 3
    // Projection Transformation

    input.open("stage2.txt");
    output.open("stage3.txt");

    Transformation P = Transformation::projectionMatrix(fovY, aspectRatio, near, far);
    // cout << "Projection Matrix:" << endl;
    // P.print();

    for (int i = 0; i < triangleCount; i++) {
        Point v1, v2, v3;
        input >> v1.x >> v1.y >> v1.z;
        input >> v2.x >> v2.y >> v2.z;
        input >> v3.x >> v3.y >> v3.z;

        v1 = P.transform(v1);
        v2 = P.transform(v2);
        v3 = P.transform(v3);

        output << v1.x << " " << v1.y << " " << v1.z << endl;
        output << v2.x << " " << v2.y << " " << v2.z << endl;
        output << v3.x << " " << v3.y << " " << v3.z << endl;
        output << endl;
    }

    input.close();
    output.close();

    // Stage 4

    return 0;
}


// class Triangle {
// public:
//     Point v1, v2, v3;

//     Triangle(Point v1, Point v2, Point v3) : v1(v1), v2(v2), v3(v3) {}
// };