#include <cstdlib>
#include <cstdio>
#if defined(_WIN32) // windows necessary setup for Delfem2 Library
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>
// include utility functions from dlm2 library
#include "delfem2/dtri2_v2dtri.h"  // dlfm2's triangle mesh generation
#include "delfem2/msh_topology_uniform.h"  // dlfm2's line mesh from triangle mesh
#include "delfem2/glfw/viewer2.h"


/**
distance of 2 points in 2D
**/
template<typename T>
double distance2(const T *p0, const T *p1) {
  return sqrt((p0[0] - p1[0]) * (p0[0] - p1[0]) + (p0[1] - p1[1]) * (p0[1] - p1[1]));
}

/**
 *
 * @param[out] W - constraint function C
 * @param[out] dW - differentiation of function W w.r.t. positions
 * @param[in] ap - positions of the end points of an edge
 * @param[in] Len - initial length of the spring
 * @param[in] stiffness - the stiffness of the spring
 */

void W_Wd_Spring(
    double &W,
    double dW[2][2],
    const double ap[2][2],
    const double Len,
    double stiffness) {
  constexpr unsigned nnode = 2;
  constexpr unsigned ndim = 2;
  const double len = distance2(ap[0], ap[1]);
  const double u01[ndim] = {
      (ap[1][0] - ap[0][0]) / len,
      (ap[1][1] - ap[0][1]) / len};
  const double C = (len - Len);       // the length difference
  const double dC[nnode][ndim] = {    // gradient of C
      {-u01[0], -u01[1]},
      {+u01[0], +u01[1]}};

  W = C;                               // distance constraint function, same with original paper formula 
  // W = 0.5 * stiffness * C * C;      // we can also use potential energy instead of constraint function 
  // W = C * stiffness                 // also we can add stiffness on the constraint function

  for (unsigned int ino = 0; ino < nnode; ++ino) {
    for (unsigned int idim = 0; idim < ndim; ++idim) {
      //dW[ino][idim] = stiffness * dC[ino][idim] * C;
      //dW[ino][idim] = dC[ino][idim] * stiffness;
      dW[ino][idim] = dC[ino][idim];
    }
  }
}

void StepTimePbdSpring2(
    std::vector<double> &aXY,
    std::vector<double> &aUV,
    const std::vector<double> &aXY0, // initial position 
    const std::vector<unsigned int> &aLine,
    double dt,
    double stiffness,
    const double *gravity,
    const std::vector<double> &aMassPointInv) {  // simulation
  const unsigned int np = static_cast<unsigned int>(aXY.size() / 2);  // number of points
  const std::vector<double> aXYt = aXY;
  // equation of motion
  for (unsigned int ip = 0; ip < np; ++ip) {
    if (aMassPointInv[ip] < 1.0e-5) { continue; }
    aXY[ip * 2 + 0] += dt * aUV[ip * 2 + 0] + dt * dt * gravity[0];
    aXY[ip * 2 + 1] += dt * aUV[ip * 2 + 1] + dt * dt * gravity[1];
  }
  // update position by gravity (one step of explicit Euler)
  for (unsigned int il = 0; il < aLine.size() / 2; ++il) {
    unsigned int ip0 = aLine[il * 2 + 0];
    unsigned int ip1 = aLine[il * 2 + 1];
    const double Len = distance2(aXY0.data() + ip0 * 2, aXY0.data() + ip1 * 2);
    const double ap[2][2] = {
        {aXY[ip0 * 2 + 0], aXY[ip0 * 2 + 1]},
        {aXY[ip1 * 2 + 0], aXY[ip1 * 2 + 1]}};
    double w, dw[2][2];

    W_Wd_Spring(
        w, dw,
        ap, Len, stiffness);

    double mpi0 = aMassPointInv[ip0];
    double mpi1 = aMassPointInv[ip1];
    double deno = 0.0;
    deno += mpi0 * (dw[0][0] * dw[0][0] + dw[0][1] * dw[0][1]);
    deno += mpi1 * (dw[1][0] * dw[1][0] + dw[1][1] * dw[1][1]);
    if (deno < 1.0e-5) { continue; }

    aXY[ip0 * 2 + 0] -= w * dw[0][0] / deno;
    aXY[ip0 * 2 + 1] -= w * dw[0][1] / deno;
    aXY[ip1 * 2 + 0] -= w * dw[1][0] / deno;
    aXY[ip1 * 2 + 1] -= w * dw[1][1] / deno;
    // update position to meet the constraint 
  }
  for (unsigned int ip = 0; ip < np; ++ip) {
    if (aMassPointInv[ip] < 1.0e-5) {
      aUV[ip * 2 + 0] = 0.0;
      aUV[ip * 2 + 1] = 0.0;
      aXY[ip * 2 + 0] = aXYt[ip * 2 + 0];
      aXY[ip * 2 + 1] = aXYt[ip * 2 + 1];
      continue;
    }
    // update velocity using new position 
    aUV[ip * 2 + 0] = (aXY[ip * 2 + 0] - aXYt[ip * 2 + 0]) / dt;
    aUV[ip * 2 + 1] = (aXY[ip * 2 + 1] - aXYt[ip * 2 + 1]) / dt;
  }
}


// using Delfem2 Lib to generate mesh 
void SettingUpSimulation(
    std::vector<double> &aXY,
    std::vector<unsigned int> &aLine,
    std::vector<double> &aMassPointInv) {
  const std::vector<std::vector<double> > aLoop = {{-1, 0, +1, 0, +1, +0.3, -1, +0.3}};
  std::vector<delfem2::CDynPntSur> aPo2D;
  std::vector<delfem2::CDynTri> aETri;
  std::vector<delfem2::CVec2d> aVec2;
  delfem2::GenMesh(
      aPo2D, aETri, aVec2,
      aLoop, 0.07, 0.07);
  std::vector<unsigned int> aTri;
  delfem2::CMeshTri2D(
      aXY, aTri,
      aVec2, aETri);

  delfem2::MeshLine_MeshElem(
      aLine,
      aTri.data(), aTri.size() / 3,
	  delfem2::MESHELEM_TRI, aXY.size() / 2);
  aMassPointInv.resize(aXY.size() / 2, 1.0);
  for (unsigned int ip = 0; ip < aXY.size() / 2; ++ip) {
    if (aXY[ip * 2 + 0] < -0.99 || aXY[ip * 2 + 0] > +0.99) {
      aMassPointInv[ip] = 0;
    }
  }
}

// =============================================

// print out error
static void error_callback(
	[[maybe_unused]] int error, 
	const char *description) {
  fputs(description, stderr);
}

void DrawMeshLine2(
    const std::vector<double> &aXY,
    const std::vector<unsigned int> &aLine) {
  glBegin(GL_LINES);
  for (unsigned int it = 0; it < aLine.size() / 2; ++it) {
    ::glVertex2dv(aXY.data() + aLine[it * 2 + 0] * 2);
    ::glVertex2dv(aXY.data() + aLine[it * 2 + 1] * 2);
  }
  glEnd();
}

int main() {
  std::vector<double> aXY;
  std::vector<unsigned int> aLine;
  std::vector<double> aMassPointInv;
  SettingUpSimulation(aXY, aLine, aMassPointInv);
  const std::vector<double> aXY0 = aXY;  // initial position
  std::vector<double> aUV(aXY.size(), 0.0);
  const double gravity[2] = {0., -10.0};
  const double stiffness = 1.0e+5;

  double dt = 1.0 / 60.0;
  //
  delfem2::glfw::CViewer2 viewer;
  {
    viewer.view_height = 1.0;
    viewer.title = "PBD spring - distance constraint";
    viewer.trans[1] = 0.3f;
  }
  glfwSetErrorCallback(error_callback);
  if (!glfwInit()) { exit(EXIT_FAILURE); }
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  // ------
  viewer.OpenWindow();

  while (!glfwWindowShouldClose(viewer.window)) {
    StepTimePbdSpring2(
        aXY, aUV,
        aXY0, aLine, dt, stiffness, gravity, aMassPointInv);

    viewer.DrawBegin_oldGL();
    // initial configuration (red color)
    glColor3f(1.f, 0.f, 0.f);
    DrawMeshLine2(aXY0, aLine);
    // PBD deformed configuration (black color)
    glColor3f(0.f, 0.f, 0.f);
    DrawMeshLine2(aXY, aLine);
    viewer.SwapBuffers();
    glfwPollEvents();
  }
  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}