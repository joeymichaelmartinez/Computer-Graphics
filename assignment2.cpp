/***
 Assignment-2: Geometric Modeling of a Scene

 Name: Wong, Alex (Please write your name in Last Name, First Name format)

 Collaborators: Doe, John; Doe, Jane
 ** Note: although the assignment should be completed individually
 you may speak with classmates on high level algorithmic concepts. Please
 list their names in this section

 Project Summary: A short paragraph (3-4 sentences) describing the work you
 did for the project.
***/



#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

/**************************************************
 *  Rectangular Prisms via Hierarchical Modeling  *
 *                                                *
 *  using planes as building blocks, build a unit *
 *  cube with transformations that will serve as  *
 *  a primitive for modeling objects in the scene *
 *                                                *
 *************************************************/

// Initializes a square plane of unit lengths
vector<GLfloat> init_plane() {
    vector<GLfloat> vertices = {
        +0.5,   +0.5,   +0.0,
        -0.5,   +0.5,   +0.0,
        -0.5,   -0.5,   +0.0,
        +0.5,   -0.5,   +0.0
    };
    return vertices;
}

// Converts a vector to an array
GLfloat* vector2array(vector<GLfloat> vec) {
    GLfloat* arr = new GLfloat[vec.size()];
    for (int i = 0; i < vec.size(); i++) {
        arr[i] = vec[i];
    }
    return arr;
}

// Converts Cartesian coordinates to homogeneous coordinates
vector<GLfloat> to_homogenous_coord(vector<GLfloat> cartesian_coords) {
    // vector<GLfloat> result;

    for (int i = cartesian_coords.size(); i > 0; i=i-3) {
        cartesian_coords.insert(cartesian_coords.begin() + i, +1.0);
    }

    // cout << cartesian_coords << "\n"
    return cartesian_coords;
    // return result;
}

// Converts Cartesian coordinates to homogeneous coordinates
vector<GLfloat> to_cartesian_coord(vector<GLfloat> homogenous_coords) {
    // vector<GLfloat> result;
    
    for (int i = homogenous_coords.size() - 1; i > 0; i=i-4) {
        homogenous_coords.erase(homogenous_coords.begin() + i);
    }
    
    return homogenous_coords;
}

// Definition of a translation matrix
vector<GLfloat> translation_matrix (float dx, float dy, float dz) {
    vector<GLfloat> translate_mat = {
        +1.0f, +0.0f, +0.0f, dx,
        +0.0f, +1.0f, +0.0f, dy,
        +0.0f, +0.0f, +1.0f, dz,
        +0.0f, +0.0f, +0.0f, +1.0f
    };    
    
    return translate_mat;
}

// Definition of a scaling matrix
vector<GLfloat> scaling_matrix (float sx, float sy, float sz) {
    vector<GLfloat> scale_mat = {
        sx, +0.0f, +0.0f, +0.0f,
        +0.0f, sy, +0.0f, +0.0f,
        +0.0f, +0.0f, sz, +0.0f,
        +0.0f, +0.0f, +0.0f, +1.0f
    };
    
    return scale_mat;
}

// Definition of a rotation matrix along the x-axis theta degrees
vector<GLfloat> rotation_matrix_x (float theta) {
    vector<GLfloat> rotate_mat_x= {
        +1.0f, +0.0f, +0.0f, +0.0f,
        +0.0f, +1.0f, +0.0f, +0.0f,
        +0.0f, +0.0f, cos(theta), -sin(theta),
        +0.0f, +0.0f, sin(theta), cos(theta)
    };

    return rotate_mat_x;
}


// Definition of a rotation matrix along the y-axis by theta degrees
vector<GLfloat> rotation_matrix_y (float theta) {
    vector<GLfloat> rotate_mat_y= {
        cos(theta), +0.0f, sin(theta), +0.0f,
        +0.0f, +1.0f, +0.0f, +0.0f,
        -sin(theta), +0.0f, cos(theta), +0.0f,
        +0.0f, +0.0f, +0.0f, +1.0f
    };
    
    return rotate_mat_y;
}


// Definition of a rotation matrix along the z-axis by theta degrees
vector<GLfloat> rotation_matrix_z (float theta) {
    vector<GLfloat> rotate_mat_z= {
        cos(theta), -sin(theta), +0.0f, +0.0f,
        sin(theta), cos(theta), +0.0f, +0.0f,
        +0.0f, +0.0f, +1.0f, +0.0f,
        +0.0f, +0.0f, +0.0f, +1.0f
    };
    
    return rotate_mat_z;
}

// Perform matrix multiplication for A B
vector<GLfloat> mat_mult(vector<GLfloat> A, vector<GLfloat> B) {
    vector<GLfloat> result;
    GLfloat multiplied_answer;
    // cout << "Hello \n";
    // for(int i = 0; i < A.size(); i++) {
    //     cout << "A" << A[i] << "\n";
    // }

    // for(int i = 0; i < B.size(); i++) {
    //     cout << "B" << B[i] << "\n";
    // }


        for (int i = 0; i < A.size(); i++) {
            multiplied_answer = 0;
            vector<GLfloat> set_of_coords;
            set_of_coords.clear();
            for(int k = i*4; k < i*4+4; k++) {
                set_of_coords.push_back(B[k]);
                // cout << "set " << B[k] << " \n";
            }
            for(int j = 0; j < 4; j++){
                cout << "A " << A[i+j] << " \n" ;
                cout << "B " << set_of_coords[j] << " \n";
                multiplied_answer += A[i+j] * set_of_coords[j];
            }
                cout << "multiplied_answer " << multiplied_answer << " \n";
                result.push_back(multiplied_answer);  
        }
    
    return result;
}

// Builds a unit cube centered at the origin
vector<GLfloat> build_cube() {
    vector<GLfloat> result;
    
    vector<GLfloat> initial_plane = to_homogenous_coord(init_plane());
    vector<GLfloat> translated_matrix = translation_matrix(0.0f, 0.0f, 0.0f);
    vector<GLfloat> side_1 = mat_mult(translated_matrix, initial_plane);
    // for(int i = 0; i < initial_plane.size(); i++) {
    //     cout << "Side one vertices before " << initial_plane[i] << " \n";
    //     // result.push_back(side_1[i]);
    // }
    side_1 = to_cartesian_coord(side_1);
    for(int i = 0; i < side_1.size(); i++) {
        // cout << "Side one vertices after " << side_1[i] << " \n";
        result.push_back(side_1[i]);
    }
    
    return result;
}

/**************************************************
 *            Camera and World Modeling           *
 *                                                *
 *  create a scene by applying transformations to *
 *  the objects built from planes and position    *
 *  the camera to view the scene by using setting *
 *  the projection viewing matrices               *
 *                                                *
 *************************************************/

void setup() {
    // Enable the vertex array functionality
    glEnableClientState(GL_VERTEX_ARRAY);
    // Enable the color array functionality (so we can specify a color for each vertex)
    glEnableClientState(GL_COLOR_ARRAY);
    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);
    // Set up some default base color
    glColor3f(0.5, 0.5, 0.5);
    // Set up white background
    glClearColor(1.0, 1.0, 1.0, 0.0);
}

void init_camera() {
    // Camera parameters
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // Define a 50 degree field of view, 1:1 aspect ratio, near and far planes at 3 and 7
    gluPerspective(50.0, 1.0, 2.0, 10.0);
    // Position camera at (2, 3, 5), attention at (0, 0, 0), up at (0, 1, 0)
    gluLookAt(2.0, 6.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}

// Construct the scene using objects built from cubes/prisms
GLfloat* init_scene() {
    GLfloat* vertex_cube = vector2array(build_cube());
    return vertex_cube;
}

// Construct the color mapping of the scene
// GLfloat* init_color() {
//     return nullptr;
// }

void display_func() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    cout << "hello \n";

    //pass the color pointer
    // glColorPointer(3,     
    //                GL_FLOAT, 
    //                0,       
    //                colors);
    

    GLfloat* cube = init_scene();
    
    glVertexPointer(3,          
                    GL_FLOAT,  
                    0,         
                    cube);
    // Draw quad point planes: each 4 vertices
    glDrawArrays(GL_QUADS, 0, 4*6);
    
    glFlush();          //Finish rendering
    glutSwapBuffers();
    
}


int main (int argc, char **argv) {
    // vector<GLfloat> homogenous_coords = to_homogenous_coord(init_plane());

    // vector<GLfloat> cartesian_coords = to_cartesian_coord(homogenous_coords);

    // vector<GLfloat> test = {2,2,2,1};

    // vector<GLfloat> answer = mat_mult(homogenous_coords, test);    
    
    // for (int i = 0; i < answer.size(); i++) {
    //     cout << answer[i] << " ";
    //     if(i % 3 == 2) {
    //         cout << " \n";
    //     }
    // }
    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    // Create a window with rendering context and everything else we need
    glutCreateWindow("Assignment 2");
    
    setup();
    init_camera();
    
    // Set up our display function
    glutDisplayFunc(display_func);

    // Render our world
    glutMainLoop();


    // Remember to call "delete" on your dynmically allocated arrays
    // such that you don't suffer from memory leaks. e.g.
    // delete arr;
    
    return 0;
}

