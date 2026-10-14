///@file
///@note If this fails to run with a segfault or "Invalid Device", it's likely due to a bug in Qt<->Wayland. Set the QT_QPA_PLATFORM=xcb environment variable to work around it.

#include <QApplication>
#include "mainwindow.h"

using namespace Waves;

/// Entry point to program
int main(int argc, char *argv[])
{
	QCoreApplication::setAttribute(Qt::AA_ShareOpenGLContexts);
	QApplication app(argc, argv);
	app.setQuitOnLastWindowClosed(true);

	mainwindow window;
	window.show();

	return app.exec();
}
