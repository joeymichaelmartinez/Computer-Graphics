/***
 Assignment-3: Shading via Lighting and Colors
 
 Name: Martinez, Joseph (Please write your name in Last Name, First Name format)
 
 Collaborators: Hardy, John
 
 Project Summary: /***
    Found how to calculate the normals and the h vector based on given points. Also found how to add color to my shapes, and 
    then update that based on the shading algorithm. Finally, I made all of the new objects and applied this shading.
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

#include <math.h>
#include <vector>
#include <iostream>
#include <cstdlib>
using namespace std;

GLfloat* objects;
GLfloat* colors;
vector<GLfloat> light_source = {0.0, 2.0, 5.0};
vector<GLfloat> camera = {35.0, 15.0, -30.0};
/**************************************************
 *              Object Model Class                *
 *                                                *
 *  Stores the points of the object as a vector   *
 *  along with the colors and normals for each    *
 *  point. Normals are computed from the points.  *
 *                                                *
 *************************************************/
class ObjectModel {
    vector<GLfloat> _points;
    vector<GLfloat> _normals;
    vector<GLfloat> _base_colors;
    vector<GLfloat> _colors;
public:
    ObjectModel() { };
    vector<GLfloat> get_points() { return _points; };
    vector<GLfloat> get_normals() { return _normals; };
    vector<GLfloat> get_base_colors() { return _base_colors; };
    vector<GLfloat> get_colors() { return _colors; };
    void set_points(vector<GLfloat> points) { _points = points; };
    void set_normals(vector<GLfloat> normals) { _normals = normals; };
    void set_base_colors(vector<GLfloat> base_colors) { _base_colors = base_colors; };
    void set_colors(vector<GLfloat> colors) { _colors = colors; };
};

ObjectModel chair_1;
ObjectModel chair_2;
ObjectModel table;
ObjectModel carpet;
ObjectModel lamp;
ObjectModel cup;

/**************************************************
 *  Rectangular Prisms via Hierarchical Modeling  *
 *                                                *
 *  Using planes as building blocks, build a unit *
 *  cube with transformations that will serve as  *
 *  a primitive for modeling objects in the scene.*
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
    result = to_cartesian_coord(result);
    return result;
}

/**************************************************
 *           Generating Surface Normals           *
 *                                                *
 *  Generate the surface normals of the objects   *
 *  using the cross product between two vectors   *
 *  that lie on the surface (plane) of interest.  *
 *  Recall that the direction of the normal to a  *
 *  surface follows the Right Hand Rule.          *
 *                                                *
 *************************************************/

vector<GLfloat> join_vectors (vector<GLfloat> A, vector<GLfloat> B) {
    vector<GLfloat> AB;
    AB.insert( AB.end(), A.begin(), A.end() );
    AB.insert( AB.end(), B.begin(), B.end() );
    return AB;
}

// Performs the cross product between two vectors
vector<GLfloat> cross_product(vector<GLfloat> A, vector<GLfloat> B) {
    vector<GLfloat> C;
    
    C = {
        A[1]*B[2] - A[2]*B[1],
        A[2]*B[0] - A[0]*B[2],
        A[0]*B[1] - A[1]*B[0]
      };

    return C;
}

// void generate_vector(vector<GLfloat> points, vector<GLfloat> &vector, int starting_position) {
//     for(int j = 0; j < 3; j++)  {
//         vector.push_back(points[i + starting_position + j] - points[j]);
//     }
// }

GLfloat find_magnitude(vector<GLfloat> A) { 
    return sqrt((A[0] * A[0]) + (A[1] * A[1]) + (A[2] * A[2]));
}

vector<GLfloat> make_unit(vector<GLfloat> A) {
    GLfloat magnitude_of_A = find_magnitude(A);
    A[0] = (A[0])/(magnitude_of_A);
    A[1] = (A[1])/(magnitude_of_A);
    A[2] = (A[2])/(magnitude_of_A);
    return A;
}

vector<GLfloat> vector_subtraction(vector<GLfloat> A, vector<GLfloat> B) {

    vector<GLfloat> C = {
        A[0] - B[0],
        A[1] - B[1],
        A[2] - B[2]
    };
    return C;
}

// Generates the normals to each surface (plane)
// *** May need to handle turning this into a unit normal vector by div by magnitude
vector<GLfloat> generate_normals(vector<GLfloat> points) {
    vector<GLfloat> normals;
    
    int three_components_of_vertex = 3;
    int points_in_plane = 4;
    vector<GLfloat> vector_1;
    vector<GLfloat> vector_2;

    for(int i = 0; i < points.size(); i = i + 12) {
        // generate_vector(points, vector_1, 3);
        // generate_vector(points, vector_2, 9);
        // *** Make a helper function possibly that makes vectors
        vector<GLfloat> q_0 = { points[i], points[i + 1], points[i + 2]};
        vector<GLfloat> q_1 = { points[i + 3], points[i + 4], points[i + 5]};
        vector<GLfloat> q_3 = { points[i + 9], points[i + 10], points[i + 11]};
        
        vector_1 = vector_subtraction(q_1, q_0);
        vector_2 = vector_subtraction(q_3, q_0);
        normals = join_vectors(normals, make_unit(cross_product(vector_1, vector_2)));
        normals = join_vectors(normals, make_unit(cross_product(vector_1, vector_2)));
        normals = join_vectors(normals, make_unit(cross_product(vector_1, vector_2)));
        normals = join_vectors(normals, make_unit(cross_product(vector_1, vector_2)));
    }
    return normals;
}

/**************************************************
 *         Shading via Lighting and Color         *
 *                                                *
 *  Generate the set of colors for each face of   *
 *  the planes that compose the objects in the    *
 *  scene. Based on the light source and surface  *
 *  normals, render the colors of the objects by  *
 *  applying the shading equation.                *
 *                                                *
 *************************************************/

// Initializes the base color of a plane to a single color
vector<GLfloat> init_base_color(GLfloat r, GLfloat g, GLfloat b) {
    vector<GLfloat> base_color = {
        r,   g,   b,
        r,   g,   b,
        r,   g,   b,
        r,   g,   b
    };
    return base_color;
}

// Initializes the base color of a plane by specifying the color of each point
vector<GLfloat> init_base_color(GLfloat r0, GLfloat g0, GLfloat b0, GLfloat r1, GLfloat g1, GLfloat b1,
                                GLfloat r2, GLfloat g2, GLfloat b2, GLfloat r3, GLfloat g3, GLfloat b3) {
    // This enables OpenGL to use interpolation for (Gouraud) shading the plane
    vector<GLfloat> base_color = {
        r0,   g0,   b0,
        r1,   g1,   b1,
        r2,   g2,   b2,
        r3,   g3,   b3
    };
    return base_color;
}


vector<GLfloat> generate_h(vector<GLfloat> &light_source, vector<GLfloat> &camera, vector<GLfloat> points) {
    GLfloat light_source_magnitude = find_magnitude(light_source);
    vector<GLfloat> v = vector_subtraction(points, camera);
    GLfloat v_magnitude = find_magnitude(v);

    GLfloat total_magnitude = v_magnitude + light_source_magnitude;
    vector<GLfloat> h = {
        (light_source[0] + v[0])/total_magnitude, 
        (light_source[1] + v[1])/total_magnitude, 
        (light_source[2] + v[2])/total_magnitude
    };
    return h; 
        
}

// Performs the dot product between two vectors
GLfloat dot_product(vector<GLfloat> A, vector<GLfloat> B) {
    GLfloat dot_product = ((A[0] * B[0]) + (A[1] * B[1]) + (A[2] * B[2]));
    return dot_product;
}

// Generates the colors of a set of surfaces based on the light source,
// surface normals, and base color of the surface
ObjectModel apply_shading(ObjectModel object_model, vector<GLfloat> light_source, vector<GLfloat> camera) {
    GLfloat amb = 0.65; 
    GLfloat diff = 0.35; 
    GLfloat spec = 0.11;

    vector<GLfloat> model_points = object_model.get_points();

    vector<GLfloat> base_colors = object_model.get_base_colors();

    vector<GLfloat> normals = object_model.get_normals();

    vector<GLfloat> colors;


    for(int i = 0; i < normals.size(); i+=3) {
        vector<GLfloat> current_points = { model_points[i], model_points[i + 1], model_points[i + 2] };
        vector<GLfloat> h = generate_h(light_source, camera, current_points);
        vector<GLfloat> normal = {normals[i], normals[i + 1], normals[i + 2]};
        GLfloat light_dot_product = dot_product(normal, light_source);
        GLfloat h_dot_product = dot_product(normals, h);
        
        colors.push_back(base_colors[i] * (amb + diff * (light_dot_product)) + (spec * base_colors[i] * h_dot_product));
        colors.push_back(base_colors[i + 1] * (amb + diff * (light_dot_product)) + (spec * base_colors[i + 1] * h_dot_product));
        colors.push_back(base_colors[i + 2] * (amb + diff * (light_dot_product)) + (spec * base_colors[i + 2] * h_dot_product));
    }

    object_model.set_colors(colors);

    return object_model;
}

// Allows for ambience (a vector of 3 values), diffusion (vector of 3 x n values) and specular (vector of 3 x n values)
// as input to the shading equation
ObjectModel apply_shading(ObjectModel object_model, vector<GLfloat> light_source, vector<GLfloat> camera,
                          vector<GLfloat> amb, vector<GLfloat> diff, vector<GLfloat> spec) {
    vector<GLfloat> colors;
    
    object_model.set_colors(colors);
    return object_model;    
}

/**************************************************
 *            Camera and World Modeling           *
 *                                                *
 *  Create a scene by applying transformations to *
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
    gluPerspective(50.0, 1.0, 2.0, 70.0);
    gluLookAt(35.0, 15.0, -30.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0); // 2, 6, 5
}

vector<GLfloat> color_cube(vector<GLfloat> &colors, float r, float g, float b) {
    colors = join_vectors(colors, init_base_color(r, g, b));
    colors = join_vectors(colors, init_base_color(r, g, b));
    colors = join_vectors(colors, init_base_color(r, g, b));
    colors = join_vectors(colors, init_base_color(r, g, b));
    colors = join_vectors(colors, init_base_color(r, g, b));
    colors = join_vectors(colors, init_base_color(r, g, b));

    return colors;

}

ObjectModel make_chair() {

    vector<GLfloat> unit_cube = to_homogenous_coord(build_cube());

    vector<GLfloat> collection_of_chair_pieces;

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

    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_leg_1);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_leg_2);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_leg_3);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_leg_4);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_seat);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_back_1);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_back_2);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_back_3);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_back_4);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_back_5);
    collection_of_chair_pieces = join_vectors(collection_of_chair_pieces, chair_back_6);

    chair_1.set_points(to_cartesian_coord(collection_of_chair_pieces));

    vector<GLfloat> colors;
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);
    color_cube(colors, .746, .5, .25);

    chair_1.set_base_colors(colors);

    chair_1.set_normals(generate_normals(chair_1.get_points()));

    chair_1 = apply_shading(chair_1, light_source, camera);
    return chair_1;
}

ObjectModel make_table() {

    vector<GLfloat> unit_cube = to_homogenous_coord(build_cube());
    
    vector<GLfloat> collection_of_table_pieces;

    vector<GLfloat> table_leg_1 = mat_mult(translation_matrix(-3.0f, 2.5f, 4.25f), mat_mult(scaling_matrix(0.5, 9.0f, 0.5f), unit_cube));
    vector<GLfloat> table_leg_2 = mat_mult(translation_matrix(3.0f, 2.5f, 4.25f), mat_mult(scaling_matrix(0.5f, 9.0f, 0.5f), unit_cube));
    vector<GLfloat> table_leg_3 = mat_mult(translation_matrix(-3.0f, 2.5f, -4.25f), mat_mult(scaling_matrix(0.5, 9.0f, 0.5f), unit_cube));
    vector<GLfloat> table_leg_4 = mat_mult(translation_matrix(3.0f, 2.5f, -4.25f), mat_mult(scaling_matrix(0.5f, 9.0f, 0.5f), unit_cube));
    vector<GLfloat> table_top = mat_mult(translation_matrix(0.0f, 7.0f, 0.0f), mat_mult(scaling_matrix(8.0f, 0.5f, 10.0f), unit_cube));

    collection_of_table_pieces = join_vectors(collection_of_table_pieces, table_leg_1);
    collection_of_table_pieces = join_vectors(collection_of_table_pieces, table_leg_2);
    collection_of_table_pieces = join_vectors(collection_of_table_pieces, table_leg_3);
    collection_of_table_pieces = join_vectors(collection_of_table_pieces, table_leg_4);
    collection_of_table_pieces = join_vectors(collection_of_table_pieces, table_top);

    vector<GLfloat> colors;
    table.set_points(to_cartesian_coord(collection_of_table_pieces));
    colors = color_cube(colors, .398, .199, 0);
    colors = color_cube(colors, .398, .199, 0);
    colors = color_cube(colors, .398, .199, 0);
    colors = color_cube(colors, .398, .199, 0);
    colors = color_cube(colors, .5, 0, 0);
    
    table.set_base_colors(colors);
    table.set_normals(generate_normals(table.get_points()));
    table = apply_shading(table, light_source, camera);
    return table;
}

ObjectModel make_carpet() {
    vector<GLfloat> colors;
    vector<GLfloat> unit_cube = to_homogenous_coord(build_cube());
    unit_cube =  to_cartesian_coord(mat_mult(translation_matrix(0, -2.5, 0), mat_mult(scaling_matrix(20, .2, 40), unit_cube))),
    carpet.set_points((unit_cube));
    carpet.set_normals(generate_normals(carpet.get_points()));
    color_cube(colors, .199, .398, 1);
    carpet.set_base_colors(colors);
    carpet = apply_shading(carpet, light_source, camera);

    return carpet;
}

ObjectModel make_lamp() {
    vector<GLfloat> unit_cube = to_homogenous_coord(build_cube());
    vector<GLfloat> collection_of_lamp_pieces;

    vector<GLfloat> leg = mat_mult(translation_matrix(0.0f, 6.0f, 0.0f), mat_mult(scaling_matrix(0.5, 16.0f, 0.5f), unit_cube));
    vector<GLfloat> base = mat_mult(translation_matrix(0.0f, -2.0f, 0.0f), mat_mult(scaling_matrix(4.0f, 0.5f, 4.0f), unit_cube));
    vector<GLfloat> shade_1 = mat_mult(translation_matrix(0.0f, 14.0f, -2.0f), mat_mult(scaling_matrix(4.2f, 5.0f, 0.3f), unit_cube));
    vector<GLfloat> shade_2 = mat_mult(translation_matrix(0.0f, 14.0f, 2.0f), mat_mult(scaling_matrix(4.2f, 5.0f, 0.3f), unit_cube));
    vector<GLfloat> shade_3 = mat_mult(translation_matrix(2.0f, 14.0f, 0.0f), mat_mult(scaling_matrix(0.3f, 5.0f, 4.2f), unit_cube));
    vector<GLfloat> shade_4 = mat_mult(translation_matrix(-2.0f, 14.0f, 0.0f), mat_mult(scaling_matrix(0.3f, 5.0f, 4.2f), unit_cube));
    
    collection_of_lamp_pieces = join_vectors(collection_of_lamp_pieces, leg);
    collection_of_lamp_pieces = join_vectors(collection_of_lamp_pieces, base);
    collection_of_lamp_pieces = join_vectors(collection_of_lamp_pieces, shade_1);
    collection_of_lamp_pieces = join_vectors(collection_of_lamp_pieces, shade_2);
    collection_of_lamp_pieces = join_vectors(collection_of_lamp_pieces, shade_3);
    collection_of_lamp_pieces = join_vectors(collection_of_lamp_pieces, shade_4);

    vector<GLfloat> colors;
    lamp.set_points(to_cartesian_coord(collection_of_lamp_pieces));
    colors = color_cube(colors, .5, .5, 0.5);
    colors = color_cube(colors, 0, 0, 0);
    colors = color_cube(colors, .5, 0, 0);
    colors = color_cube(colors, .5, 0, 0);
    colors = color_cube(colors, .5, 0, 0);
    colors = color_cube(colors, .5, 0, 0);
    
    lamp.set_base_colors(colors);
    lamp.set_normals(generate_normals(lamp.get_points()));
    lamp = apply_shading(lamp, light_source, camera);
    return lamp;
}

ObjectModel make_cup() {
    vector<GLfloat> unit_cube = to_homogenous_coord(build_cube());
    vector<GLfloat> collection_of_cup_pieces;

    vector<GLfloat> left = mat_mult(translation_matrix(0.4f, 8.5f, 0.0f), mat_mult(scaling_matrix(0.2, 2.5f, 1.0f), unit_cube));
    vector<GLfloat> right = mat_mult(translation_matrix(-0.4f, 8.5f, 0.0f), mat_mult(scaling_matrix(0.2, 2.5f, 1.0f), unit_cube));
    vector<GLfloat> front = mat_mult(translation_matrix(0.0f, 8.5f, 0.4f), mat_mult(scaling_matrix(1.0, 2.5f, 0.2f), unit_cube));
    vector<GLfloat> back = mat_mult(translation_matrix(0.0f, 8.5f, -0.4f), mat_mult(scaling_matrix(1.0, 2.5f, 0.2f), unit_cube));
    vector<GLfloat> bottom = mat_mult(translation_matrix(0.0f, 6.5f, 0.0f), mat_mult(scaling_matrix(1.0, 0.2f, 1.0f), unit_cube));
    vector<GLfloat> handle_1 = mat_mult(translation_matrix(0.0f, 9.0f, -0.8f), mat_mult(scaling_matrix(0.4f, 0.2f, 0.6f), unit_cube));
    vector<GLfloat> handle_2 = mat_mult(translation_matrix(0.0f, 8.5f, -1.0f), mat_mult(scaling_matrix(0.4f, 1.0f, 0.2f), unit_cube));
    vector<GLfloat> handle_3 = mat_mult(translation_matrix(0.0f, 8.0f, -0.8f), mat_mult(scaling_matrix(0.4f, 0.2f, 0.6f), unit_cube));
    
    collection_of_cup_pieces = join_vectors(collection_of_cup_pieces, left);
    collection_of_cup_pieces = join_vectors(collection_of_cup_pieces, right);
    collection_of_cup_pieces = join_vectors(collection_of_cup_pieces, front);
    collection_of_cup_pieces = join_vectors(collection_of_cup_pieces, back);
    collection_of_cup_pieces = join_vectors(collection_of_cup_pieces, bottom);
    collection_of_cup_pieces = join_vectors(collection_of_cup_pieces, handle_1);
    collection_of_cup_pieces = join_vectors(collection_of_cup_pieces, handle_2);
    collection_of_cup_pieces = join_vectors(collection_of_cup_pieces, handle_3);

    vector<GLfloat> colors;
    cup.set_points(to_cartesian_coord(collection_of_cup_pieces));
    colors = color_cube(colors, 1, 1, 1);
    colors = color_cube(colors, 1, 1, 1);
    colors = color_cube(colors, 1, 1, 1);
    colors = color_cube(colors, 1, 1, 1);
    colors = color_cube(colors, 1,1,1);
    colors = color_cube(colors, 1,1,1);
    colors = color_cube(colors, 1,1,1);
    colors = color_cube(colors, 1,1,1);
    
    cup.set_base_colors(colors);
    cup.set_normals(generate_normals(cup.get_points()));
    cup = apply_shading(cup, light_source, camera);
    return cup;
}

// Construct the scene using objects built from cubes/prisms
GLfloat* init_scene() {

    vector<GLfloat> objects_vector;
    vector<GLfloat> chair_points = to_homogenous_coord(chair_1.get_points());
    chair_1.set_points(to_cartesian_coord(mat_mult(translation_matrix(0.0f, 0.0f, -8.0f), mat_mult(rotation_matrix_y(-20), chair_points))));
    chair_2.set_points(to_cartesian_coord(mat_mult(translation_matrix(0.0f, 0.0f, 6.0f), mat_mult(rotation_matrix_y(70), chair_points))));
    
    vector<GLfloat> lamp_points = to_homogenous_coord(lamp.get_points());
    lamp.set_points(to_cartesian_coord(mat_mult(translation_matrix(5.0f, 0.0f, -15.0f), lamp_points)));

    for(int i = 0; i < chair_1.get_points().size(); i++){
        objects_vector.push_back(chair_1.get_points()[i]);
    } 

    for(int i = 0; i < chair_2.get_points().size(); i++){
        objects_vector.push_back(chair_2.get_points()[i]);
    }

    for(int i = 0; i < table.get_points().size(); i++){
        objects_vector.push_back(table.get_points()[i]);
    }

    for(int i = 0; i < carpet.get_points().size(); i++) {
        objects_vector.push_back(carpet.get_points()[i]);
    }

    for(int i = 0; i < lamp.get_points().size(); i++) {
        objects_vector.push_back(lamp.get_points()[i]);
    }

    for(int i = 0; i < cup.get_points().size(); i++) {
        objects_vector.push_back(cup.get_points()[i]);
    }

    GLfloat* objects = vector2array(objects_vector);
    return objects;
}

// Construct the color mapping of the scene
GLfloat* init_color() {
    vector<GLfloat> chair_1_colors = chair_1.get_colors();
    vector<GLfloat> chair_2_colors = chair_2.get_colors();
    vector<GLfloat> table_colors = table.get_colors();
    vector<GLfloat> carpet_colors = carpet.get_colors();
    vector<GLfloat> lamp_colors = lamp.get_colors();
    vector<GLfloat> cup_colors = cup.get_colors();
    vector<GLfloat> all_colors = join_vectors(chair_1_colors, chair_2_colors);
    all_colors = join_vectors(all_colors, table_colors);
    all_colors = join_vectors(all_colors, carpet_colors);
    all_colors = join_vectors(all_colors, lamp_colors);
    all_colors = join_vectors(all_colors, cup_colors);
    GLfloat* array_of_colors = vector2array(all_colors);
    return array_of_colors;
    
}

float theta = 0.0;

void display_func() {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glRotatef(theta, 0.0, 1.0, 0.0);
    // glRotatef(theta, 1.0, 0.0, 0.0);

    
    glVertexPointer(3,          
                    GL_FLOAT,  
                    0,         
                    objects);

    glColorPointer(3,
                   GL_FLOAT,
                   0,
                   colors);

    glDrawArrays(GL_QUADS, 0, 6*4*42);
    
    
    glFlush();         
    glutSwapBuffers();
}

void idle_func() {
    theta = theta+0.3;
    display_func();
}

int main (int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(800, 600);

    glutCreateWindow("Assignment 3");
    setup();
    init_camera();
    chair_1 = make_chair();
    chair_2 = make_chair();
    make_carpet();
    make_table();
    make_lamp();
    make_cup();
    objects = init_scene();
    colors = init_color();
    
    glutDisplayFunc(display_func);

    glutIdleFunc(idle_func);
    glutMainLoop();
    

    delete objects;
    delete colors;

    return 0;
}

