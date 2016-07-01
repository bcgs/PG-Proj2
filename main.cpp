#include <iostream>
#include <GLUT/GLUT.h>

#include "math.h"

const int WINDOW_W = 800;
const int WINDOW_H = 600;


/* Class declarations */
class Vector;

/* Function declarations */
Vector crossProduct(Vector u, Vector v);
void normalize(Vector& v);
Vector orthogonalization(Vector u, Vector v);

/* Variable declarations */
double zBuffer[600][800];


/* CLASSES */

class Color {
public:
    int R, G, B;
    
    Color() {};
    Color(int R, int G, int B) : R(R), G(G), B(B) {};
    ~Color() {};
};

class Vector {
public:
    double x, y, z;
    
    Vector() {};
    Vector(double x, double y, double z) : x(x), y(y), z(z) {};
    ~Vector() {};
};

class Vertex {
public:
    double x, y, z;
    Vector normal;
    Vector luminousIntensity;
    
    Vertex() {};
    Vertex(double x, double y, double z) : x(x), y(y), z(z) {
        normal = Vector(0,0,0);
        luminousIntensity = Vector(0,0,0);
    };
    ~Vertex() {};
};

class Triangle {
public:
    Vertex a, b, c;
    Vector normal;
    
    Triangle() {};
    Triangle(Vertex v1, Vertex v2, Vertex v3) : a(v1), b(v2), c(v3) {
        setNormal();
    }
    ~Triangle() {};
    
    void setNormal() {
        Vector u(b.x - a.x, b.y - a.y, b.z - a.z);
        Vector v(c.x - a.x, c.y - a.y, c.z - a.z);
        normal = crossProduct(u, v);
        normalize(normal);
    }
};

class Camera {
public:
    double x, y, z; // C
    Vector N, U, V;
    double d, hx, hy;
    
    Camera() {}
    Camera(double x, double y, double z, Vector N, Vector V, double d,
           double hx, double hy) : x(x), y(y), z(z), N(N), d(d), hx(hx), hy(hy)
    {
        this->V = orthogonalization(N, V);   // V => V'
        normalize(this->N);
        normalize(this->V);
        this->U = crossProduct(this->V, this->N);
    }
    ~Camera() {};
    
    /* Global coordinate -> Local coordinates */
    Vertex globalTOlocal(Vector P) {
        Vertex Pv;
        Vector P_C(P.x - x, P.y - y, P.z - z);  // P-C
        
        Pv.x = U.x * P_C.x + U.y * P_C.y + U.z * P_C.z;
        Pv.y = V.x * P_C.x + V.y * P_C.y + V.z * P_C.z;
        Pv.z = N.x * P_C.x + N.y * P_C.y + N.z * P_C.z;
        
        return Pv;
    }
    
    Vertex pointOnUserScreen(Vector P) {
        Vertex Pv = globalTOlocal(P);
        
        double Xp = (Pv.x * d)/(hx * Pv.z);
        double Yp = (Pv.y * d)/(hy * Pv.z);
        
        double xScreen = ((Xp + 1)/2) * WINDOW_W;
        double yScreen = ((1 - Yp)/2) * WINDOW_H;
        
        return Vertex(xScreen, yScreen, 0);
    }
};


/* FUNCTIONS */

double scalarProduct(Vector u, Vector v) {
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

Vector crossProduct(Vector u, Vector v) {
    Vector s;
    s.x = u.y * v.z - u.z * v.y;
    s.y = u.z * v.x - u.x * v.z;
    s.z = u.x * v.y - u.y * v.x;
    
    return s;
}

void normalize(Vector& v) {
    double norm = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    v.x = v.x/norm;
    v.y = v.y/norm;
    v.z = v.z/norm;
}

Vector orthogonalization(Vector u, Vector v) {
    Vector s;
    double c = scalarProduct(u, v)/scalarProduct(u, u);
    s.x = v.x - c * u.x;
    s.y = v.y - c * u.y;
    s.z = v.z - c * u.z;
    
    return s;
}

void initializeBuffer() {
    for (int i = 0; i < 600; ++i) {
        for (int j = 0; j < 800; ++j) {
            zBuffer[i][j] = __DBL_MAX__;
        }
    }
}


/* GLUT */

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glFlush();
}

void reshape(int w, int h) {
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0f, WINDOW_W, WINDOW_H, 0.0f, 0.0f, 5.0f); //(LEFT|RIGHT|BOTTOM|TOP|NEAR|FAR)
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void handleKeypress(unsigned char key, int x, int y) {
    switch (key) {
        case 27:
            exit(0);
            break;
    }
}

int main(int argc, char ** argv) {
    glutInit(&argc, argv);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WINDOW_W, WINDOW_H);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    
    glutCreateWindow("PG-Proj2");
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glMatrixMode(GL_MODELVIEW); // scene model matrix
    glLoadIdentity();
    
    glutDisplayFunc(display);
    glutKeyboardFunc(handleKeypress);
    glutReshapeFunc(reshape);
    
    glutMainLoop();
    
    return 0;
}

