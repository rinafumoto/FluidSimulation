#include "MainWindow.h"

#include <QApplication>
#include <QSurfaceFormat>

int main(int argc, char *argv[])
{
  QSurfaceFormat format;
  format.setMajorVersion(4);
  format.setMinorVersion(6);
  format.setProfile(QSurfaceFormat::CoreProfile);
  format.setSamples(4);
  QSurfaceFormat::setDefaultFormat(format);
  QApplication a(argc, argv);
  MainWindow w;
  w.show();
  return a.exec();
}
