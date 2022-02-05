#include "MainWindow.h"
#include "ui_MainWindow.h"

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent)
  , ui(new Ui::MainWindow)
{
  ui->setupUi(this);
  m_gl = new NGLScene(this);
  ui->m_mainWindowGridLayout->addWidget(m_gl);

  // connect the GUI buttons to slots
  connect(ui->m_initialise,&QPushButton::clicked,m_gl,&NGLScene::initialise);
  connect(ui->m_start,&QPushButton::clicked,m_gl,&NGLScene::start);
  connect(ui->m_stop,&QPushButton::clicked,m_gl,&NGLScene::stop);
  connect(ui->m_step,&QPushButton::clicked,m_gl,&NGLScene::step);
  // connect the GUI dropdown value to a slot
  connect(ui->m_type,SIGNAL(currentIndexChanged(int)),m_gl,SLOT(setType(int)));
  // connect the GUI checkboxs value to slots
  connect(ui->m_flu,&QCheckBox::stateChanged,m_gl,&NGLScene::setFlu);
  connect(ui->m_vel,&QCheckBox::stateChanged,m_gl,&NGLScene::setVel);
  connect(ui->m_pre,&QCheckBox::stateChanged,m_gl,&NGLScene::setPre);
  // connect the other GUI value to slots
  connect(ui->m_resolutionX,SIGNAL(valueChanged(int)),m_gl,SLOT(setResolutionX(int)));
  connect(ui->m_resolutionY,SIGNAL(valueChanged(int)),m_gl,SLOT(setResolutionY(int)));
  connect(ui->m_gridsize,SIGNAL(valueChanged(double)),m_gl,SLOT(setGridSize(double)));
  connect(ui->m_timestep,SIGNAL(valueChanged(double)),m_gl,SLOT(setTimeStep(double)));
  connect(ui->m_forceX,SIGNAL(valueChanged(double)),m_gl,SLOT(setForceX(double)));
  connect(ui->m_forceY,SIGNAL(valueChanged(double)),m_gl,SLOT(setForceY(double)));
  connect(ui->m_right,SIGNAL(valueChanged(int)),m_gl,SLOT(setInitRight(int)));
  connect(ui->m_left,SIGNAL(valueChanged(int)),m_gl,SLOT(setInitLeft(int)));
  connect(ui->m_top,SIGNAL(valueChanged(int)),m_gl,SLOT(setInitTop(int)));
  connect(ui->m_bottom,SIGNAL(valueChanged(int)),m_gl,SLOT(setInitBottom(int)));
}

MainWindow::~MainWindow()
{
  delete ui;
  delete m_gl;
}

