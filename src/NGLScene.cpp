#include <QMouseEvent>
#include <QGuiApplication>

#include "NGLScene.h"
#include <ngl/NGLInit.h>
#include <ngl/VAOPrimitives.h>
#include <ngl/ShaderLib.h>
#include <ngl/Util.h>
#include <iostream>


const auto ColourShader = "ColourShader";
const auto GridViz = "GridViz";
const auto PressureViz = "PressureViz";
const auto SolidShader = "SolidShader";

NGLScene::NGLScene(QWidget *_parent) : QOpenGLWidget(_parent)
{

}

////////// Functions called on a button click //////////

// initialise the simulation with the setting values from GUI and start timer.
void NGLScene::initialise()
{
  m_fluid = std::make_unique<Fluid>();
  m_fluid->initialise(m_resolutionX, m_resolutionY, m_gridSize, m_timeStep, ngl::Vec3(m_forceX, m_forceY, 0), m_initRight, m_initLeft, m_initTop, m_initBottom, m_type);
  startTimer(10);
  update();
}

void NGLScene::step()
{
  m_fluid->simulate();
  update();
}

void NGLScene::start()
{
  sim=true;
}

void NGLScene::stop()
{
  sim=false;
}

////////// Functions called on a value change //////////

void NGLScene::setResolutionX(int i)
{
  m_resolutionX = i;
  update();
}

void NGLScene::setResolutionY(int i)
{
  m_resolutionY = i;
  update();
}

void NGLScene::setGridSize(double d)
{
  m_gridSize = d;
  update();
}

void NGLScene::setTimeStep(double d)
{
  m_timeStep = d;
  update();
}

void NGLScene::setForceX(double d)
{
  m_forceX = d;
  update();
}

void NGLScene::setForceY(double d)
{
  m_forceY = d;
  update();
}

void NGLScene::setInitRight(int i)
{
  m_initRight = i;
  update();
}

void NGLScene::setInitLeft(int i)
{
  m_initLeft = i;
  update();
}

void NGLScene::setInitTop(int i)
{
  m_initTop = i;
  update();
}

void NGLScene::setInitBottom(int i)
{
  m_initBottom = i;
  update();
}

void NGLScene::setType(int i)
{
  m_type = i;
  update();
}

void NGLScene::setFlu(bool b)
{
  m_flu = b;
  update();
}

void NGLScene::setVel(bool b)
{
  m_vel = b;
  update();
}

void NGLScene::setPre(bool b)
{
  m_pre = b;
  update();
}

NGLScene::~NGLScene()
{
  std::cout<<"Shutting down NGL, removing VAO's and Shaders\n";
}

void NGLScene::resizeGL(int _w , int _h)
{
  m_win.width  = static_cast<int>( _w * devicePixelRatio() );
  m_win.height = static_cast<int>( _h * devicePixelRatio() );
}

void NGLScene::initializeGL()
{
  ngl::NGLInit::initialize();
  // Set the background colour to white.
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  // Enable depth testing for drawing.
  glEnable(GL_DEPTH_TEST);
  // Load shaders
  ngl::ShaderLib::loadShader(ColourShader, "shaders/ColourVertex.glsl", "shaders/ColourFragment.glsl");
  ngl::ShaderLib::loadShader(PressureViz, "shaders/PressureVizVertex.glsl", "shaders/PressureVizFragment.glsl");
  ngl::ShaderLib::loadShader(GridViz,"shaders/VectorVizVertex.glsl","shaders/VectorVizFragment.glsl","shaders/VectorVizGeometry.glsl");
  ngl::ShaderLib::loadShader(SolidShader,"shaders/SolidVertex.glsl","shaders/SolidFragment.glsl","shaders/SolidGeometry.glsl");
  // initialise the simulation.
  initialise();
}

void NGLScene::timerEvent ( QTimerEvent *_event)
{
  // Simulate one step on timer event when the simulation has been started.
  if(sim)
  {
    m_fluid->simulate();
    update();
  }
}

void NGLScene::paintGL()
{
  // Clear the screen and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glViewport(0,0,m_win.width,m_win.height);
  m_fluid->render(m_win.width,m_win.height,m_flu,m_vel,m_pre);
}

//----------------------------------------------------------------------------------------------------------------------

void NGLScene::keyPressEvent(QKeyEvent *_event)
{
  // this method is called every time the main window recives a key event.
  // we then switch on the key value and set the camera in the GLWindow
  switch (_event->key())
  {
  // escape key to quite
  case Qt::Key_Escape : QGuiApplication::exit(EXIT_SUCCESS); break;
  case Qt::Key_Space :
      m_win.spinXFace=0;
      m_win.spinYFace=0;
      m_modelPos.set(ngl::Vec3::zero());

  break;
  default : break;
  }
  // finally update the GLWindow and re-draw

    update();
}
