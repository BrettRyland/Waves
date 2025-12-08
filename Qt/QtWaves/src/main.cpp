///@file

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
