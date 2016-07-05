#include <GLUT/GLUT.h>
#include <iostream>
#include <fstream>
#include "math.h"
#include <vector>
#include <string>

const int WINDOW_W = 800;
const int WINDOW_H = 600;


/* Class declarations */
class Vector;

/* Function declarations */
Vector crossProduct(Vector u, Vector v);
void normalize(Vector& v);
Vector orthogonalization(Vector u, Vector v);
Vector vectorAddition(Vector u, Vector v);

/* Global variable declarations */
double zBuffer[600][800];
int vertexListSize;
int triangleListSize;


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
    Vertex(double x, double y, double z) : x(x), y(y), z(z) {};
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
        this->U = crossProduct(this->N, this->V);
    }
    ~Camera() {};
    
    /* Global coordinate -> Local coordinates */
    void globalTOlocal(Vertex& Pv) {
        Vector P_C(Pv.x - x, Pv.y - y, Pv.z - z);  // P-C
        
        Pv.x = U.x * P_C.x + U.y * P_C.y + U.z * P_C.z;
        Pv.y = V.x * P_C.x + V.y * P_C.y + V.z * P_C.z;
        Pv.z = N.x * P_C.x + N.y * P_C.y + N.z * P_C.z;
    }
    
    Vertex pointOnUserScreen(Vertex Pv) {
        double Xp = (Pv.x * d)/(hx * Pv.z);
        double Yp = (Pv.y * d)/(hy * Pv.z);
        
        double xScreen = ((Xp + 1)/2) * WINDOW_W;
        double yScreen = ((1 - Yp)/2) * WINDOW_H;
        
        return Vertex(xScreen, yScreen, 0);
    }
};

class Illumination {
public:
    double Ka, Kd, Ks;
    double n;
    Vector Od;
    
    Color Ia, Il;
    Vertex Pl;
    
    Illumination() {};
    Illumination(double Ka, double Kd, double Ks, double n, Vector Od,
                Color Ia, Color Il, Vertex Pl) : Ka(Ka), Kd(Kd), Ks(Ks),
                n(n), Od(Od), Ia(Ia), Il(Il), Pl(Pl) {};
    ~Illumination() {};
};

class Object {
public:
    std::vector<Vertex> vertices;
    std::vector<Vertex> vertices_local;
    std::vector<Triangle> triangles;
    
    Object() {};
    ~Object() {};
    
    void addVertice(Vertex vertex) {
        vertices.push_back(vertex);
    }
    
    void addVertice_local(Vertex vertex) {
        vertices_local.push_back(vertex);
    }
    
    void addTriangle(int v1, int v2, int v3) {
        triangles.push_back(Triangle(vertices[v1-1], vertices[v2-1], vertices[v3-1]));
        
        // Update vertices' normal
        vertices[v1-1].normal = vectorAddition(vertices[v1-1].normal, triangles.back().normal);
        vertices[v2-1].normal = vectorAddition(vertices[v2-1].normal, triangles.back().normal);
        vertices[v3-1].normal = vectorAddition(vertices[v3-1].normal, triangles.back().normal);
        
        vertices_local[v1-1].normal = vectorAddition(vertices_local[v1-1].normal, triangles.back().normal);
        vertices_local[v2-1].normal = vectorAddition(vertices_local[v2-1].normal, triangles.back().normal);
        vertices_local[v3-1].normal = vectorAddition(vertices_local[v3-1].normal, triangles.back().normal);
    }
    
    void normalizeVerticesNormal() {
        for (int i = 0; i < vertices.size(); ++i) {
            normalize(vertices[i].normal);
            normalize(vertices_local[i].normal);
        }
    }
};

/* READ FILES */

Camera camera;
void readCamera(std::string path) {
    std::ifstream file;
    file.open(path);
    
    if(file.is_open()) {
        double x, y, z;
        
        file >> x;
        file >> y;
        file >> z;
        
        Vector N, V;
        
        file >> N.x;
        file >> N.y;
        file >> N.z;
        
        file >> V.x;
        file >> V.y;
        file >> V.z;
        
        double d, hx, hy;
        
        file >> d;
        file >> hx;
        file >> hy;
        
        camera = Camera(x, y, z, N, V, d, hx, hy);
        
    } else {
        printf("File not found.");
    }
    file.close();
}

Illumination illumination;
void readIllumination(std::string path) {
    std::ifstream file;
    file.open(path);
    
    if(file.is_open()) {
        Vertex Pl;
        
        file >> Pl.x;
        file >> Pl.y;
        file >> Pl.z;
        camera.globalTOlocal(Pl);   //OBS: ler a camera antes da iluminação
        
        double Ka;
        
        file >> Ka;
        
        Color Ia;
        
        file >> Ia.R;
        file >> Ia.G;
        file >> Ia.B;
        
        double Kd;
        
        file >> Kd;
        
        Vector Od;
        
        file >> Od.x;
        file >> Od.y;
        file >> Od.z;
        
        double Ks;
        
        file >> Ks;
        
        Color Il;
        
        file >> Il.R;
        file >> Il.G;
        file >> Il.B;
        
        double n;
        
        file >> n;
        
        illumination = Illumination(Ka, Kd, Ks, n, Od, Ia, Il, Pl);
        
    } else {
        printf("File not found.");
    }
    file.close();
}

Object object;
void readObject(std::string path) {
    std::ifstream file;
    file.open(path);
    
    if(file.is_open()) {
        file >> vertexListSize;
        file >> triangleListSize;
        
        Vertex vertex;
        
        for (int i = 0; i < vertexListSize; ++i) {
            file >> vertex.x;
            file >> vertex.y;
            file >> vertex.z;
            
            object.addVertice(vertex);
            camera.globalTOlocal(vertex);
            object.addVertice_local(vertex);
        }
        
        double v1, v2, v3;
        
        for (int i = 0; i < triangleListSize; ++i) {
            file >> v1;
            file >> v2;
            file >> v3;
            
            object.addTriangle(v1, v2, v3);
        }
        object.normalizeVerticesNormal();
        
    } else {
         printf("File not found.");
    }
    file.close();
}

/* FUNCTIONS */

Vector vectorAddition(Vector u, Vector v) {
    return Vector(u.x + v.x, u.y + v.y, u.z + v.z);
}

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
//    readCamera("/Users/bcgs/Downloads/PG/Cameras/01_Camera.cfg.txt");
//    readIllumination("/Users/bcgs/Downloads/PG/Iluminacao.txt");
//    readObject("/Users/bcgs/Downloads/PG/Objetos/01_Objeto.byu");
//    std::cout << object.triangles[0].normal.x << "," <<
//    object.triangles[0].normal.y << "," <<
//    object.triangles[0].normal.z << std::endl;
//    std::cout << object.vertices[1304].normal.x << std::endl;
    
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

