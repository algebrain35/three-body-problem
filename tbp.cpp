#include <cmath>
#include <iostream>
#include "raylib.h"
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include<GLFW/glfw3.h>
using namespace std;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Vector4d;


#define MIDPOINT(p1, p2, p3) (p1 + p2 + p3) / 3
#define SCALE 50
#define Vector3Normalize(v1) v1 / v1.norm()
double MAXDISTANCE(Vector3d p1, Vector3d p2, Vector3d p3, Vector3d C) {
  return max(
      max(sqrt(pow((p1 - C).x(), 2) + pow((p1 - C).z(), 2)),
      sqrt(pow((p2 - C).x(), 2) + pow((p2 - C).z(),2))),
      sqrt(pow((p3 - C).x(), 2) + pow((p3 - C).z(),2))
      );
}

Vector3 raylib_vec3(Vector3d v) {
  return Vector3{(float) (v.x()), (float) (v.y()), (float) (v.z())};
}


void get_min_max(Vector3d* v_min, Vector3d* v_max, string fname) {
  ifstream infile (fname);
  assert(infile.is_open());
  double x_max = -INFINITY, y_max = -INFINITY, z_max = -INFINITY;
  double x_min = INFINITY, y_min = INFINITY, z_min = INFINITY;

  while(!infile.eof()) {
    double x, y, z;
    infile >> x >> y >> z;
    if(x > x_max) {
      x_max = x;
    }

    if(y > y_max) {
      y_max = y;
    }

    if(z > z_max) {
      z_max = z;
    }

    if(x < x_min) {
      x_min = x;
    }

    if(y < y_min) {
      y_min = y;
    }

    if(z < z_min) {
      z_min = z;
    }

  }
  infile.close();
  *v_min = Vector3d{x_min, y_min, z_min};
  *v_max = Vector3d{x_max, y_max, z_max};
};
void next_position_state(Matrix3d& points, ifstream& infile) {
  assert(infile.is_open());  
  for(int i = 0; i < 3 && !infile.eof(); i++) {
      double x, y, z;
      infile >> x >> y >> z;
      points.row(i) << x, y, z;
  }
}

void load_positions(Matrix3d m[], ifstream& infile) {
  assert(infile.is_open());
  Matrix3d tmp;
  int idx = 0;
  while(!infile.eof()) {
      for(int i = 0; i < 3 && !infile.eof(); i++) {
        double x, y, z;
        infile >> x >> y >> z;
        tmp.row(i) << x, y, z;
      }
      m[idx++] = tmp;
      
  }

}
int main() {
    const string filename = "tbp_state_5.txt";
    Vector3d v_min, v_max;
    
    get_min_max(&v_min, &v_max, filename);
    Vector3d midpoint = SCALE * (v_max + v_min) / 2.0;
    float diag = (float)SCALE * (v_max - v_min).norm();

    
    ifstream tbp_file(filename);
    assert(tbp_file.is_open());

    InitWindow(1200, 800, "Three Body Problem Simulation");
    SetTargetFPS(60);

    
    Camera3D camera = { 0 };
    camera.position = {(float)midpoint.x(), (float)midpoint.y(), (float)(midpoint.z() + diag * 1.5f)};
    camera.target   = raylib_vec3(midpoint);
    camera.up       = { 0.0f, 1.0f, 0.0f };
    camera.fovy     = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    Matrix3d points;

    
    while (!WindowShouldClose() && !tbp_file.eof()) {
        next_position_state(points, tbp_file);
        UpdateCamera(&camera, CAMERA_FREE);

        BeginDrawing();
            ClearBackground(BLACK);
            BeginMode3D(camera);
                DrawSphere(raylib_vec3(SCALE * points.col(0)), 10.0f, YELLOW);
                DrawSphere(raylib_vec3(SCALE * points.col(1)), 2.5f, BLUE);
                
                DrawSphere(raylib_vec3(SCALE * points.col(2)), 1.0f, GRAY);

                
                DrawGrid(20, diag / 10.0f);

            EndMode3D();

            
            DrawFPS(10, 10);
            DrawText("Use Mouse + WASD to fly", 10, 30, 20, DARKGRAY);
        EndDrawing();
    }

    tbp_file.close();
    CloseWindow();
    return 0;
}

