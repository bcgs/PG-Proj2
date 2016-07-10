#include <GLUT/GLUT.h>
#include <iostream>
#include <fstream>
#include "math.h"
#include <vector>
#include <string>
#include <algorithm>

const int WINDOW_W = 800;
const int WINDOW_H = 600;


/* Class declarations */
class Vector;

/* Function declarations */
Vector crossProduct(Vector u, Vector v);
double scalarProduct(Vector u, Vector v);
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
    
    /* Pixel projection on user screen based on percentage */
    void pointOnUserScreen(Vertex& Pscreen) {
        double Xp = (Pscreen.x * d)/(hx * Pscreen.z);
        double Yp = (Pscreen.y * d)/(hy * Pscreen.z);
        
        Pscreen.x = ((Xp + 1)/2) * WINDOW_W;
        Pscreen.y = ((1 - Yp)/2) * WINDOW_H;
        Pscreen.z = 0;
    }
};

class Illumination {
public:
    double Ka, Kd, Ks, n;
    Vector Od;
    Color Ia, Il;
    Vertex Pl;
    
    Illumination() {};
    ~Illumination() {};
    Illumination(double Ka, double Kd, double Ks, double n, Vector Od,
                 Color Ia, Color Il, Vertex Pl) : Ka(Ka), Kd(Kd), Ks(Ks),
    n(n), Od(Od), Ia(Ia), Il(Il), Pl(Pl) {};
    
    Vector getL(Vertex _3DPoint) {
        Vector L(Pl.x - _3DPoint.x, Pl.y - _3DPoint.y, Pl.z - _3DPoint.z);
        normalize(L);
        
        return L;
    }
    
    Vertex getAmbComp() {
        return Vertex(Ka * Ia.R, Ka * Ia.G, Ka * Ia.B);
    }
    
    Vertex getDiffuseComp(Vertex _3DPoint) {
        Vector L = getL(_3DPoint);
        Vector N = _3DPoint.normal;
        Vertex diffuse_comp;
        
        diffuse_comp.x = scalarProduct(N, L) * Il.R * Od.x * Kd;
        diffuse_comp.y = scalarProduct(N, L) * Il.G * Od.y * Kd;
        diffuse_comp.z = scalarProduct(N, L) * Il.B * Od.z * Kd;
        
        return diffuse_comp;
    }
    
    Vector getR(Vertex _3DPoint) {
        Vector L = getL(_3DPoint);
        Vector N = _3DPoint.normal;
        Vector R;
        
        R.x = (2 * N.x * scalarProduct(N, L)) - L.x;
        R.y = (2 * N.y * scalarProduct(N, L)) - L.y;
        R.z = (2 * N.z * scalarProduct(N, L)) - L.z;
        normalize(R);
        
        return R;
    }
    
    Vertex getSpecularComp(Vertex _3DPoint, Vector V) {
        Vertex specularComp;
        Vector R = getR(_3DPoint);
        
        specularComp.x = pow(scalarProduct(V, R), n) * Ks * Il.R;
        specularComp.y = pow(scalarProduct(V, R), n) * Ks * Il.G;
        specularComp.z = pow(scalarProduct(V, R), n) * Ks * Il.B;
        
        return specularComp;
    }
};

class Object {
public:
    std::vector<Vertex> vertices_global;
    std::vector<Vertex> vertices_local;
    std::vector<Vertex> vertices2D;
    
    std::vector<Triangle> triangles2D;
    std::vector<Triangle> triangles3D;
    
    Object() {};
    ~Object() {};
    
    void addVertice_global(Vertex vertex) {
        vertices_global.push_back(vertex);
    }
    
    void addVertice_local(Vertex vertex) {
        vertices_local.push_back(vertex);
    }
    
    void addVertice2D(Vertex vertex) {
        vertices2D.push_back(vertex);
    }
    
    void addTriangle(int v1, int v2, int v3) {
        triangles3D.push_back(Triangle(vertices_local[v1-1], vertices_local[v2-1], vertices_local[v3-1]));
        triangles2D.push_back(Triangle(vertices2D[v1-1], vertices2D[v2-1], vertices2D[v3-1]));
        
        //update vertices_local's normal
        vertices_local[v1-1].normal = vectorAddition(vertices_local[v1-1].normal, triangles3D.back().normal);
        vertices_local[v2-1].normal = vectorAddition(vertices_local[v2-1].normal, triangles3D.back().normal);
        vertices_local[v3-1].normal = vectorAddition(vertices_local[v3-1].normal, triangles3D.back().normal);
    }
    
    void normalizeVerticesNormal() {
        for (int i = 0; i < vertices_global.size(); ++i) {
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
        
        std::cout << "Camera OK!" << std::endl;
        
    } else {
        printf("Camera not found.\n");
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
        camera.globalTOlocal(Pl);
        
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
        std::cout << "Illumination OK!" << std::endl;
        
    } else {
        printf("Illumination not found.\n");
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
            
            object.addVertice_global(vertex);
            
            camera.globalTOlocal(vertex);
            object.addVertice_local(vertex);
            
            camera.pointOnUserScreen(vertex);
            object.addVertice2D(vertex);
        }
        
        double v1, v2, v3;
        
        for (int i = 0; i < triangleListSize; ++i) {
            file >> v1;
            file >> v2;
            file >> v3;
            
            object.addTriangle(v1, v2, v3);
        }
        object.normalizeVerticesNormal();
        
        std::cout << "Object OK!" << std::endl;
    } else {
        printf("Object not found.\n");
    }
    file.close();
}


/* FUNCTIONS */

double triangleArea(Triangle triangle) {
    
    Vector AB(triangle.b.x - triangle.a.x, triangle.b.y - triangle.a.y, triangle.b.z - triangle.a.z);
    Vector AC(triangle.c.x - triangle.a.x, triangle.c.y - triangle.a.y, triangle.c.z - triangle.a.z);
    Vector BC(triangle.c.x - triangle.b.x, triangle.c.y - triangle.b.y, triangle.c.z - triangle.b.z);
    
    double sideAB = sqrt(AB.x * AB.x + AB.y * AB.y);
    double sideAC = sqrt(AC.x * AC.x + AC.y * AC.y);
    double sideBC = sqrt(BC.x * BC.x + BC.y * BC.y);
    
    double semiperimeter = (sideAB + sideAC + sideBC)/2;
    double area = sqrt(semiperimeter * (semiperimeter-sideAB) * (semiperimeter-sideAC) * (semiperimeter-sideBC));
    
    return area;
}

// Barycentric coordinate system

double alpha(Triangle triangle, Vertex p) {
    return (triangleArea(Triangle(p, triangle.b, triangle.c)) / triangleArea(triangle));
}

double beta(Triangle triangle, Vertex p) {
    return (triangleArea(Triangle(triangle.a, p, triangle.c)) / triangleArea(triangle));
}

void scanLine(double xmin, double xmax, int yScan, Triangle triangle2D, Triangle triangle3D) {
    
    for (int x = xmin; x <= xmax; ++x) {
        
        Vertex pixel(x, yScan, 0);  // current vertex
        
        double _alpha = alpha(triangle2D, pixel);
        double _beta = beta(triangle2D, pixel);
        double _gamma = 1 - (_alpha + _beta);
        
        Vertex _3DPoint;
        
        _3DPoint.x = _alpha * triangle3D.a.x + _beta * triangle3D.b.x + _gamma * triangle3D.c.x;
        _3DPoint.y = _alpha * triangle3D.a.y + _beta * triangle3D.b.y + _gamma * triangle3D.c.y;
        _3DPoint.z = _alpha * triangle3D.a.z + _beta * triangle3D.b.z + _gamma * triangle3D.c.z;
        
        // Check if the current pixel is within the user screen
        
        if(x <= WINDOW_W && x >= 0 && yScan <= WINDOW_H && yScan >= 0) {
            
            if(_3DPoint.z < zBuffer[x][yScan]) {
                if(_3DPoint.z > 1) {                    // Is not behind the camera
                    zBuffer[x][yScan] = _3DPoint.z;
                    
                    _3DPoint.normal.x = triangle3D.a.normal.x * _alpha + triangle3D.b.normal.x * _beta + triangle3D.c.normal.x * _gamma;
                    _3DPoint.normal.y = triangle3D.a.normal.y * _alpha + triangle3D.b.normal.y * _beta + triangle3D.c.normal.y * _gamma;
                    _3DPoint.normal.z = triangle3D.a.normal.z * _alpha + triangle3D.b.normal.z * _beta + triangle3D.c.normal.z * _gamma;
                    normalize(_3DPoint.normal);
                    
                    // Headed to the camera (user eyes)
                    Vector V(-_3DPoint.x, -_3DPoint.y, -_3DPoint.z);
                    normalize(V);
                    
                    // If negative then normal is pointing to the opposite side
                    if(scalarProduct(_3DPoint.normal, V) < 0) {
                        _3DPoint.normal.x = -_3DPoint.normal.x;
                        _3DPoint.normal.y = -_3DPoint.normal.y;
                        _3DPoint.normal.z = -_3DPoint.normal.z;
                    }
                    
                    Vertex luminousIntensity = illumination.getAmbComp();
                    
                    // Headed to the light
                    Vector L = illumination.getL(_3DPoint);
                    
                    // If negative neither diffuse nor specular
                    if(scalarProduct(_3DPoint.normal, L) > 0) {
                        Vertex diffuseComp = illumination.getDiffuseComp(_3DPoint);
                        
                       luminousIntensity.x += diffuseComp.x;
                       luminousIntensity.y += diffuseComp.y;
                       luminousIntensity.z += diffuseComp.z;
                        
                        Vector R = illumination.getR(_3DPoint);
                        
                        // If negative no specular
                        if(scalarProduct(V, R) > 0) {
                            
                            Vertex specularComp = illumination.getSpecularComp(_3DPoint, V);
                            
                            luminousIntensity.x += specularComp.x;
                            luminousIntensity.y += specularComp.y;
                            luminousIntensity.z += specularComp.z;
                            
                            if(luminousIntensity.x > 255) luminousIntensity.x = 255;
                            if(luminousIntensity.y > 255) luminousIntensity.y = 255;
                            if(luminousIntensity.z > 255) luminousIntensity.z = 255;
                            
                        }
                    }
                    glColor3f(luminousIntensity.x/255,luminousIntensity.y/255,luminousIntensity.z/255);
                    glBegin(GL_POINTS);
                    glVertex2f(x, yScan);
                    glEnd();
                }
            }
        }
    }
}

// Angular coefficient of a segment
double angularCoefficient(Vertex a, Vertex b) {
    if(b.x - a.x == 0) return 0;
    return (b.y - a.y)/(b.x - a.x);
}

// Compare function to sort Y values in order to set yScan
bool compare(Vertex a, Vertex b) {
    return a.y != b.y ? a.y < b.y : a.x < b.x;
}

void scanConversion(Triangle triangle2D, Triangle triangle3D) {
    
    Vertex v1 = triangle2D.a;
    v1.y = (int) v1.y;
    Vertex v2 = triangle2D.b;
    v2.y = (int) v2.y;
    Vertex v3 = triangle2D.c;
    v3.y = (int) v3.y;
    
    Vertex toSort[] = {v1, v2, v3};
    
    // C++ function to sort
    std::sort(toSort, toSort+3, compare);
    
    Vertex smaller, middle, taller;
    smaller = toSort[0];
    middle  = toSort[1];
    taller  = toSort[2];
    
    int yScan = smaller.y;
    
    // Angular coefficient of a segment AB
    double a, b, c;
    a = angularCoefficient(smaller, middle);
    b = angularCoefficient(smaller, taller);
    c = angularCoefficient(middle, taller);
    
    double xmin, xmax;
    xmin = smaller.x;
    xmax = smaller.x;
    
    if(middle.y == smaller.y) {
        xmin = smaller.x;
        xmax = middle.x;
    }
    
    // _1a from xmax' = 1/a + xmax
    double _1a, _1b, _1c;
    _1a = 1/a;
    _1b = 1/b;
    _1c = 1/c;
    
    // if 1/a, 1/b and/or tend/s to infinity... make it 0
    if(_1a == std::numeric_limits<double>::infinity()) _1a = 0;
    if(_1b == std::numeric_limits<double>::infinity()) _1b = 0;
    if(_1c == std::numeric_limits<double>::infinity()) _1c = 0;
    
    int set = 0;
    if(smaller.y != middle.y) {         // if smaller-middle segment is not parallel to axis x
        if(middle.y != taller.y) {      // if middle-taller segment is not parallel to axis x
            if(_1a <= _1b)
                set = 1;
        } else {
            set = 1;
        }
    }
    
    // xmin2 & xmax2 are used when a triangle is divided in two subtriangles
    int _xmin, _xmax, _xmin2, _xmax2;
    
    if(set == 1) {
        _xmin  = 0; _xmax  = 1;
        _xmin2 = 2; _xmax2 = 1;
    } else {
        _xmin  = 1; _xmax  = 0;
        _xmin2 = 1; _xmax2 = 2;
    }
    
    while(yScan <= taller.y) {
        scanLine(xmin, xmax, yScan, triangle2D, triangle3D);
        ++yScan;
        
        if(yScan <= middle.y) {
            switch (_xmin) {
                case 0:
                    xmin += _1a;
                    break;
                case 1:
                    xmin += _1b;
                    break;
            }
            switch (_xmax) {
                case 0:
                    xmax += _1a;
                    break;
                case 1:
                    xmax += _1b;
                    break;
            }
        } else {
            switch (_xmin2) {
                case 1:
                    xmin += _1b;
                    break;
                case 2:
                    xmin += _1c;
                    break;
            }
            switch (_xmax2) {
                case 1:
                    xmax += _1b;
                    break;
                case 2:
                    xmax += _1c;
                    break;
            }
        }
    }
}

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
    std::cout << "Buffer OK.\n" << std::endl;
}

void execute() {
    readCamera("/Users/bcgs/Downloads/PG/Cameras/vaso.cfg.txt");
    readObject("/Users/bcgs/Downloads/PG/Objetos/vaso.byu");
    readIllumination("/Users/bcgs/Downloads/PG/Iluminacao.txt");
    
    initializeBuffer();
}

/* GLUT */

void drawTriangle(std::vector<Triangle> triangles) {
    glBegin(GL_LINES);
    glColor3f(1.0f, 1.0f, 0.0f);
    for(auto p : triangles) {
        glVertex2f(p.a.x, p.a.y);
        glVertex2f(p.b.x, p.b.y);
        
        glVertex2f(p.b.x, p.b.y);
        glVertex2f(p.c.x, p.c.y);
        
        glVertex2f(p.c.x, p.c.y);
        glVertex2f(p.a.x, p.a.y);
    }
    glEnd();
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    
//    drawTriangle(object.triangles2D);
    for (int i = 0; i < triangleListSize; ++i) {
        scanConversion(object.triangles2D[i], object.triangles3D[i]);
    }
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
    execute();
    
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
