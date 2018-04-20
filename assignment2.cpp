/***
 Assignment-2: Geometric Modeling of a Scene

 Name: Martinez, Joseph 
 Project Summary: I first found how to turn turn a vector into 
 having homogenous coordinates or cartesian coordinates. I then 
 designed the different matrices that were necessary for transformations. 
 Once I had these I tweaked the initial cube with matrix multiplication
 until I figured out how to get matrix multiplication just right. Once
 this was done, I then assembled the initial cube, and manipulated it 
 to form all the various objects that I needed.
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

int number_of_sides;

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

    for (int i = cartesian_coords.size(); i > 0; i=i-3) {
        cartesian_coords.insert(cartesian_coords.begin() + i, +1.0);
    }

    return cartesian_coords;
}

// Converts Cartesian coordinates to homogeneous coordinates
vector<GLfloat> to_cartesian_coord(vector<GLfloat> homogenous_coords) {
    
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

    double rad = (double)M_PI*theta/180.0f;
    vector<GLfloat> rotate_mat_x= {
        +1.0f, +0.0f, +0.0f, +0.0f,
        +0.0f, (float)::cos(rad), (float)-::sin(rad), +0.0f,
        +0.0f, (float)::sin(rad), (float)::cos(rad), +0.0f,
        +0.0f, +0.0f, +0.0f, 1
    };

    return rotate_mat_x;
}


// Definition of a rotation matrix along the y-axis by theta degrees
vector<GLfloat> rotation_matrix_y (float theta) {

    double rad = (double)M_PI*theta/180.0f;
    vector<GLfloat> rotate_mat_y= {
        (float)::cos(rad), +0.0f, (float)::sin(rad), +0.0f,
        +0.0f, +1.0f, +0.0f, +0.0f,
        (float)-::sin(rad), +0.0f, (float)::cos(rad), +0.0f,
        +0.0f, +0.0f, +0.0f, +1.0f
    };
    
    return rotate_mat_y;
}


// Definition of a rotation matrix along the z-axis by theta degrees
vector<GLfloat> rotation_matrix_z (float theta) {

    double rad = (double)M_PI*theta/180.0f;
    vector<GLfloat> rotate_mat_z= {
        (float)::cos(rad), (float)-::sin(rad), +0.0f, +0.0f,
        (float)::sin(rad), (float)::cos(rad), +0.0f, +0.0f,
        +0.0f, +0.0f, +1.0f, +0.0f,
        +0.0f, +0.0f, +0.0f, +1.0f
    };
    
    return rotate_mat_z;
}

// Perform matrix multiplication for A B
vector<GLfloat> mat_mult(vector<GLfloat> A, vector<GLfloat> B) {
    vector<GLfloat> result;
    float multiplied_answer;
    
    for(int k = 0; k < B.size(); k+=4) {
        for (int i = 0; i < A.size(); i+=4) {
                multiplied_answer = 0;
                int l = 0;
                for(int j = k; j < k + 4; j++) {
                    multiplied_answer += A[i+l] * B[j];
                    l++;
                }
                    result.push_back(multiplied_answer);  
        }
    }
    return result;
}

// Builds a unit cube centered at the origin
vector<GLfloat> build_cube() {
    vector<GLfloat> result;
    
    vector<GLfloat> initial_plane = to_homogenous_coord(init_plane());

    vector<vector<GLfloat>> sides_vector; 

    sides_vector.push_back(mat_mult(translation_matrix(0.0f, 0.0f, 0.5f), initial_plane));

    sides_vector.push_back(mat_mult(translation_matrix(-0.5f, 0.0f, 0.0f), mat_mult(rotation_matrix_y(-90), initial_plane)));

    sides_vector.push_back(mat_mult(translation_matrix(0.5f, 0.0f, 0.0f), mat_mult(rotation_matrix_y(90), initial_plane)));

    sides_vector.push_back(mat_mult(translation_matrix(0.0f, 0.0f, -0.5f), mat_mult(rotation_matrix_y(180), initial_plane)));

    sides_vector.push_back(mat_mult(translation_matrix(0.0f, 0.5f, 0.0f), mat_mult(rotation_matrix_x(-90), initial_plane)));

    sides_vector.push_back(mat_mult(translation_matrix(0.0f, -0.5f, 0.0f), mat_mult(rotation_matrix_x(90), initial_plane)));
    
    for(int i = 0; i < sides_vector.size(); i++){
        for(int j = 0; j < sides_vector[i].size(); j++) {
            result.push_back(sides_vector[i][j]);
        }
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
    gluPerspective(50.0, 1.0, 2.0, 50.0);
    gluLookAt(20.0, 15.0, -15.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}

// Construct the scene using objects built from cubes/prisms
GLfloat* init_scene() {

    vector<GLfloat> unit_cube = build_cube();

    vector<vector<GLfloat>> collection_of_chair_pieces;
    vector<vector<GLfloat>> collection_of_table_pieces;
    vector<GLfloat> objects_vector;
    vector<GLfloat> chair_1;
    vector<GLfloat> chair_2;
    vector<GLfloat> table;

    vector<GLfloat> chair_leg_1 = mat_mult(translation_matrix(-2.25f, 0.0f, 2.75f), mat_mult(scaling_matrix(0.5, 5.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_leg_2 = mat_mult(translation_matrix(2.25f, 0.0f, 2.75f), mat_mult(scaling_matrix(0.5f, 5.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_leg_3 = mat_mult(translation_matrix(-2.25f, 0.0f, -2.75f), mat_mult(scaling_matrix(0.5, 5.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_leg_4 = mat_mult(translation_matrix(2.25f, 0.0f, -2.75f), mat_mult(scaling_matrix(0.5f, 5.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_seat = mat_mult(translation_matrix(0.0f, 2.5f, 0.0f), mat_mult(scaling_matrix(5.0f, 1.0f, 6.0f), unit_cube));
    vector<GLfloat> chair_back_1 = mat_mult(translation_matrix(-2.25f, 6.0f, -2.75f), mat_mult(scaling_matrix(0.5f, 6.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_back_2 = mat_mult(translation_matrix(-2.25f, 6.0f, -1.375f), mat_mult(scaling_matrix(0.5f, 6.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_back_3 = mat_mult(translation_matrix(-2.25f, 6.0f, 0.0f), mat_mult(scaling_matrix(0.5f, 6.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_back_4 = mat_mult(translation_matrix(-2.25f, 6.0f, 1.375f), mat_mult(scaling_matrix(0.5f, 6.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_back_5 = mat_mult(translation_matrix(-2.25f, 6.0f, 2.75f), mat_mult(scaling_matrix(0.5f, 6.0f, 0.5f), unit_cube));
    vector<GLfloat> chair_back_6 = mat_mult(translation_matrix(-2.25f, 9.0f, 0.0f), mat_mult(scaling_matrix(0.5f, 1.0f, 6.0f), unit_cube));

    collection_of_chair_pieces.push_back(chair_leg_1);
    collection_of_chair_pieces.push_back(chair_leg_2);
    collection_of_chair_pieces.push_back(chair_leg_3);
    collection_of_chair_pieces.push_back(chair_leg_4);
    collection_of_chair_pieces.push_back(chair_seat);
    collection_of_chair_pieces.push_back(chair_back_1);
    collection_of_chair_pieces.push_back(chair_back_2);
    collection_of_chair_pieces.push_back(chair_back_3);
    collection_of_chair_pieces.push_back(chair_back_4);
    collection_of_chair_pieces.push_back(chair_back_5);
    collection_of_chair_pieces.push_back(chair_back_6);

    for(int i = 0; i < collection_of_chair_pieces.size(); i++){
        for(int j = 0; j < collection_of_chair_pieces[i].size(); j++) {
            chair_1.push_back(collection_of_chair_pieces[i][j]);
        }
    } 

    vector<GLfloat> chair_1_transformed = mat_mult(translation_matrix(0.0f, 0.0f, -8.0f), mat_mult(rotation_matrix_y(-20), chair_1));
    vector<GLfloat> chair_2_transformed = mat_mult(translation_matrix(0.0f, 0.0f, 6.0f), mat_mult(rotation_matrix_y(70), chair_1));

    for(int i = 0; i < chair_1_transformed.size(); i++){
        objects_vector.push_back(chair_1_transformed[i]);
    } 

    for(int i = 0; i < chair_2_transformed.size(); i++){
        objects_vector.push_back(chair_2_transformed[i]);
    }

    vector<GLfloat> table_leg_1 = mat_mult(translation_matrix(-3.0f, 2.0f, 4.25f), mat_mult(scaling_matrix(0.5, 9.0f, 0.5f), unit_cube));
    vector<GLfloat> table_leg_2 = mat_mult(translation_matrix(3.0f, 2.0f, 4.25f), mat_mult(scaling_matrix(0.5f, 9.0f, 0.5f), unit_cube));
    vector<GLfloat> table_leg_3 = mat_mult(translation_matrix(-3.0f, 2.0f, -4.25f), mat_mult(scaling_matrix(0.5, 9.0f, 0.5f), unit_cube));
    vector<GLfloat> table_leg_4 = mat_mult(translation_matrix(3.0f, 2.0f, -4.25f), mat_mult(scaling_matrix(0.5f, 9.0f, 0.5f), unit_cube));
    vector<GLfloat> table_top = mat_mult(translation_matrix(0.0f, 7.0f, 0.0f), mat_mult(scaling_matrix(8.0f, 0.5f, 10.0f), unit_cube));

    collection_of_table_pieces.push_back(table_leg_1);
    collection_of_table_pieces.push_back(table_leg_2);
    collection_of_table_pieces.push_back(table_leg_3);
    collection_of_table_pieces.push_back(table_leg_4);
    collection_of_table_pieces.push_back(table_top);

    for(int i = 0; i < collection_of_table_pieces.size(); i++){
        for(int j = 0; j < collection_of_table_pieces[i].size(); j++) {
            objects_vector.push_back(collection_of_table_pieces[i][j]);
        }
    }

    objects_vector = to_cartesian_coord(objects_vector);

    number_of_sides = collection_of_chair_pieces.size()*2 + collection_of_table_pieces.size();
    cout << number_of_sides << "\n";

    GLfloat* objects = vector2array(objects_vector);
    return objects;
}

float theta = 0.0;

void display_func() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    GLfloat* objects = init_scene();
    glVertexPointer(3,          
                    GL_FLOAT,  
                    0,         
                    objects);

    glDrawArrays(GL_QUADS, 0, number_of_sides*4*6);
    
    delete objects;
    glFlush();         
    glutSwapBuffers();
    
}

void idle_func() {
    theta = theta+0.3;
    display_func();
}

int main (int argc, char **argv) {
    vector<GLfloat> homogenous_coords = to_homogenous_coord(init_plane());
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Assignment 2");
    
    setup();
    init_camera();
    
    glutDisplayFunc(display_func);

    glutIdleFunc(idle_func);
    glutMainLoop();
    
    return 0;
}

