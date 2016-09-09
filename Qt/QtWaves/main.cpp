#include <QApplication>

#include "mainwindow.h"

using namespace Waves;

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);
	app.setQuitOnLastWindowClosed(true);
	app.setAttribute(Qt::AA_ShareOpenGLContexts);

	mainwindow window;
	window.show();

	return app.exec();
}
