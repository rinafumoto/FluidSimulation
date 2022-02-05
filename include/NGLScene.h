#ifndef NGLSCENE_H_
#define NGLSCENE_H_
#include <ngl/Vec3.h>
#include <memory>
#include "Fluid.h"
#include "WindowParams.h"
// this must be included after NGL includes else we get a clash with gl libs
#include <QOpenGLWidget>
//----------------------------------------------------------------------------------------------------------------------
/// @file NGLScene.h
/// @brief this class inherits from the Qt OpenGLWindow and allows us to use NGL to draw OpenGL
/// @author Jonathan Macey
/// @version 1.0
/// @date 10/9/13
/// Revision History :
/// This is an initial version used for the new NGL6 / Qt 5 demos
/// @class NGLScene
/// @brief our main glwindow widget for NGL applications all drawing elements are
/// put in this file
//----------------------------------------------------------------------------------------------------------------------

class NGLScene : public QOpenGLWidget
{
    Q_OBJECT
  public:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor for our NGL drawing class
    /// @param [in] parent the parent window to the class
    //----------------------------------------------------------------------------------------------------------------------
    NGLScene(QWidget *_parent);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief dtor must close down ngl and release OpenGL resources
    //----------------------------------------------------------------------------------------------------------------------
    ~NGLScene() override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief the initialize class is called once when the window is created and we have a valid GL context
    /// use this to setup any default GL stuff
    //----------------------------------------------------------------------------------------------------------------------
    void initializeGL() override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this is called everytime we want to draw the scene
    //----------------------------------------------------------------------------------------------------------------------
    void paintGL() override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this is called everytime we resize the window
    //----------------------------------------------------------------------------------------------------------------------
    void resizeGL(int _w, int _h) override;

public slots :

    //----------------------------------------------------------------------------------------------------------------------
	/// @brief a slot to initialise the simulation with the setting values from GUI
    //----------------------------------------------------------------------------------------------------------------------
    void initialise();

    //----------------------------------------------------------------------------------------------------------------------
	/// @brief a slot to simulate one step
    //----------------------------------------------------------------------------------------------------------------------
    void step();

    //----------------------------------------------------------------------------------------------------------------------
	/// @brief a slot to start the simulation
    //----------------------------------------------------------------------------------------------------------------------
    void start();

    //----------------------------------------------------------------------------------------------------------------------
	/// @brief a slot to stop the simulation
    //----------------------------------------------------------------------------------------------------------------------
    void stop();
    
    //----------------------------------------------------------------------------------------------------------------------
	/// @brief a slot to set the values from GUI
	/// @param i the integer value to set
	/// @param f the float value to set
	/// @param b the boolean value to set
    //----------------------------------------------------------------------------------------------------------------------
    void setResolutionX(int i);
    void setResolutionY(int i);
    void setGridSize(double f);
    void setTimeStep(double f);
    void setForceX(double f);
    void setForceY(double f);
    void setInitRight(int i);
    void setInitLeft(int i);
    void setInitTop(int i);
    void setInitBottom(int i);
    void setType(int i);
    void setFlu(bool b);
    void setVel(bool b);
    void setPre(bool b);

private:

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief Qt Event called when a key is pressed
    /// @param [in] _event the Qt event to query for size etc
    //----------------------------------------------------------------------------------------------------------------------
    void keyPressEvent(QKeyEvent *_event) override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called every time a mouse is moved
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mouseMoveEvent (QMouseEvent * _event ) override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse button is pressed
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mousePressEvent ( QMouseEvent *_event) override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse button is released
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void mouseReleaseEvent ( QMouseEvent *_event ) override;
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the timer event occurs
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void timerEvent ( QTimerEvent *_event) override;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief this method is called everytime the mouse wheel is moved
    /// inherited from QObject and overridden here.
    /// @param _event the Qt Event structure
    //----------------------------------------------------------------------------------------------------------------------
    void wheelEvent( QWheelEvent *_event) override;


    /// @brief windows parameters for mouse control etc.
    WinParams m_win;
    /// position for our model
    ngl::Vec3 m_modelPos;


    std::unique_ptr<Fluid> m_fluid;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief variables to store the simulation settings from GUI
    //----------------------------------------------------------------------------------------------------------------------
    int m_type = 0;
    bool m_flu = 1;
    bool m_vel = 0;
    bool m_pre = 0;
    int m_resolutionX = 10;
    int m_resolutionY = 10;
    double m_gridSize = 0.1f;
    double m_timeStep = 0.01f;
    double m_forceX = 0.0f;
    double m_forceY = -9.81f;
    int m_initRight = 6;
    int m_initLeft = 4;
    int m_initTop = 6;
    int m_initBottom = 4;

    //----------------------------------------------------------------------------------------------------------------------
    /// @brief The flag used to start and stop the simulation
    //----------------------------------------------------------------------------------------------------------------------
    bool sim = false;
};



#endif
