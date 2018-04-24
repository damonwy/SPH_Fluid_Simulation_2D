//------------------------------------------------
//  Project 4:
//  SPH simulation
//  Ying Wang
//
//-------------------------------------------------

#include <cmath>
#include "CmdLineFind.h"
#include <OpenGL/gl.h>   // OpenGL itself.
#include <OpenGL/glu.h>  // GLU support library.
#include <GLUT/glut.h>
#include <Eigen/Dense>

#include <iostream>
#include <OpenImageIO/imageio.h>

#include "sph.h"
#include "Particle.h"

using namespace Eigen;
using namespace std;
using namespace lux;
OIIO_NAMESPACE_USING

SPH *sph_demo;
bool capture_screen;
bool is_on;
int frame;

string captured_file_basename;
int xmouse_prev, ymouse_prev;

//--------------------------------------------------------
//
//  Initialization routines
//
//  
// Initialize all of the fields to zero
void Initialize( float *data, int size, float value )
{
#pragma omp parallel for
   for(int i=0;i<size;i++ ) { data[i] = value; }
}

void setNbCores( int nb )
{
  // omp_set_num_threads( nb );
}

//----------------------------------------------------
//
//  GL and GLUT callbacks
//
//----------------------------------------------------

void setInitialProjectionPosition(){
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0f, 500, 500, 0.0f, 0.0f, 1.0f);
}

void draw(){
    glPointSize(4.0f);
    glBegin(GL_POINTS);
    for (int i = 0; i < (int)sph_demo->particles.size(); i++){
      Vector2f x = sph_demo->particles[i]->x;
      Vector3f color = sph_demo->particles[i]->color;
      glColor3f(color(0), color(1), color(2));
      glVertex3f(x(0), x(1), 0.0f);
    }

    glEnd();
    glFlush();
}


void cbDisplay( void )
{
   glClear(GL_COLOR_BUFFER_BIT );
   draw();
   glutSwapBuffers();
   if( capture_screen )
   {
      std::stringstream os; os<<frame;
      string dispframe = os.str();
      if( frame < 1000 ){ dispframe = "0" + dispframe; }
      if( frame < 100 ){ dispframe = "0" + dispframe; }
      if( frame < 10 ){ dispframe = "0" + dispframe; }
      string fname = captured_file_basename + "." + dispframe + ".exr";
      //writeOIIOImage( fname.c_str(), display_map ); 
      cout << "Frame written to file " << fname << endl;
   }
}

// animate and display new result
void cbIdle()
{
   if(is_on){
    sph_demo->updateFluid();
   } 

   glutPostRedisplay();
   frame++;
}

void cbOnKeyboard( unsigned char key, int x, int y )
{
   switch (key) 
   {
      case 'c':
      capture_screen = !capture_screen;
      break;

      case 'g':
      sph_demo->grav *= 0.9f;
      cout << "decrease gravity constant to "<<sph_demo->grav<<endl;
      break;

      case 'G':
      sph_demo->grav *= 1.0/0.9f;
      cout << "increase gravity constant to "<<sph_demo->grav<<endl;
      break;

      case 'q':
      sph_demo->dt *= 0.9f;
      cout << "decreaese time step to "<<sph_demo->dt<<endl;
      break;

      case 'Q':
      sph_demo->dt *= 1.0/0.9f;
      cout << "increase time step to "<<sph_demo->dt<<endl;
      break;

      case 'v':
      sph_demo->viscosity_coeff *= 0.9f;
      cout << "decrease viscosity_coeff to " << sph_demo->viscosity_coeff << endl;
      break;

      case 'V':
      sph_demo->viscosity_coeff *= 1.0/0.9f;
      cout << "increase viscosity_coeff to " << sph_demo->viscosity_coeff << endl;
      break;

      case 'e':
      sph_demo->epsilon *= 0.9f;
      cout << "decrease epsilon to " << sph_demo->epsilon << endl;
      break;

      case 'E':
      sph_demo->epsilon *= 1.0/0.9f;
      cout << "increase epsilon to " << sph_demo->epsilon << endl;
      break;

      case 'h':
      sph_demo->h *= 0.9f;
      cout << "decrease h to " << sph_demo->h << endl;
      break;

      case 'H':
      sph_demo->h *= 1.0/0.9f;
      cout << "increase h to " << sph_demo->h << endl;
      break;

      case 'u':
      sph_demo-> upsilon*= 0.9f;
      cout << "decrease upsilon to " << sph_demo->upsilon << endl;
      break;

      case 'U':
      sph_demo->upsilon *= 1.0/0.9f;
      cout << "increase upsilon to " << sph_demo->upsilon << endl;
      break;

      case 'd':
      sph_demo->den_bar *= 0.9f;
      cout << "decrease den_bar to " << sph_demo->den_bar << endl;
      break;

      case 'D':
      sph_demo->den_bar *= 1.0/0.9f;
      cout << "increase den_bar to " << sph_demo->den_bar << endl;
      break;

      case 'p':
      sph_demo->pres_bar *= 0.9f;
      cout << "decrease pres_bar to " << sph_demo->pres_bar << endl;
      break;

      case 'P':
      sph_demo->pres_bar *= 1.0/0.9f;
      cout << "increase pres_bar to " << sph_demo->pres_bar << endl;
      break;

      case 'w':
      sph_demo->wallsticky *= 0.9f;
      cout << "decrease wallsticky to " << sph_demo->wallsticky << endl;
      break;

      case 'W':
      sph_demo->wallsticky *= 1.0/0.9f;
      cout << "increase wallsticky to " << sph_demo->wallsticky << endl;
      break;

      case ' ':
      if(is_on){
         is_on = false;
         cout << "stop simulation" << endl;
      }else{
         is_on = true;
         cout << "start simulation" << endl;
      }
      break;

      default:
      break;
   }
}

void cbMouseDown( int button, int state, int x, int y )
{
   if( button != GLUT_LEFT_BUTTON ) { return; }
   if( state != GLUT_DOWN ) { return; }
   xmouse_prev = x;
   ymouse_prev = y;
}

void cbMouseMove( int x, int y )
{
   xmouse_prev = x;
   ymouse_prev = y;
}

void PrintUsage()
{
   cout << "sph_2d keyboard choices\n";
   cout << "space   start/stop simulation\n";

   cout << "c       toggles screen capture on/off\n";
   cout << "d       decrease density_bar\n";
   cout << "D       increase density_bar\n";

   cout << "h       decrease h\n";
   cout << "H       increase h\n";

   cout << "q       decrease size of the time step\n";
   cout << "Q       increase size of the time step\n";

   cout << "g       decrease gravity constant\n";
   cout << "G       increase gravity constant\n";

   cout << "v       decrease viscosity coeff\n";
   cout << "V       increase viscosity coeff\n";

   cout << "p       decrease pressure_bar\n";
   cout << "P       increase pressure_bar\n";

   cout << "w       decrease wallsticky\n";
   cout << "W       increase wallsticky\n";

   cout << "u       decrease upsilon\n";
   cout << "U       increase upsilon\n";

   cout << "e       decrease epsilon\n";
   cout << "E       increase epsilon\n";
}

//---------------------------------------------------

int main(int argc, char** argv)
{
    frame = 1;
    CmdLineFind clf( argc, argv );
        
    capture_screen = clf.findFlag("-capture");
    captured_file_basename = clf.find("-fname", "sphbeginning" );

    //setNbCores(4);
    clf.usage("-h");
    clf.printFinds();
    PrintUsage();
     
    is_on = false;
    double num_particles = 1500;
    float init_x = 300.0f;
    float init_y = 400.0f;
    float h = 20;
    float eps = 0.0008f;
    float vis = 100.0f;
    float base_density = 0.0017f;
    float power = 3.0f;
    float base_pressure = 1.0f;
    float dt = 0.01f;

    Vector2f window_llc, window_urc;
    window_llc << 0.0f, 0.0f;
    window_urc << 500.0f, 500.0f;

    sph_demo = new SPH(num_particles, init_x, init_y, h, eps, vis, base_density, base_pressure, power, dt);
    sph_demo->LLC = window_llc;
    sph_demo->URC = window_urc;

   // GLUT routines
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
   glutInitWindowSize( 500, 500 );

   // Open a window 
   char title[] = "sph2d Demo";
   glutCreateWindow( title );
   
   glClearColor( 1,1,1,1 );
   setInitialProjectionPosition();
   glutDisplayFunc(&cbDisplay);
   glutIdleFunc(&cbIdle);
   
   glutKeyboardFunc(&cbOnKeyboard);
   glutMouseFunc( &cbMouseDown );
   glutMotionFunc( &cbMouseMove );
  
   glutMainLoop();
   return 1;
};
