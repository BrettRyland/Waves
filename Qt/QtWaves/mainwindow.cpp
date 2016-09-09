#include <QDebug>

#include "mainwindow.h"

namespace Waves {
	mainwindow::mainwindow(QWidget *parent) :
		QMainWindow(parent),
		ui(new Ui::mainwindowClass)
	{
		ui->setupUi(this);
		// Set up various signal/slot connections
		QObject::connect(this, &mainwindow::on_ICButton_clicked, ui->openGLWidget, &OGLWidget::change_initial_conditions);
		QObject::connect(this, &mainwindow::on_BCButton_clicked, ui->openGLWidget, &OGLWidget::change_boundary_conditions);
		QObject::connect(this, &mainwindow::on_resetviewButton_clicked, ui->openGLWidget, &OGLWidget::reset_view);
		QObject::connect(this, &mainwindow::on_quitButton_clicked, ui->openGLWidget, &OGLWidget::quit);
		QObject::connect(ui->openGLWidget, &OGLWidget::toggle_fullscreen, this, &mainwindow::on_fullscreenButton_clicked);
		QObject::connect(ui->openGLWidget, &OGLWidget::toggle_paused_integrator, this, &mainwindow::update_pauseButton);
	}

	mainwindow::~mainwindow()
	{
		delete ui;
	}

	// Fullscreen button
	void mainwindow::on_fullscreenButton_clicked()
	{
		// TODO Handle fullscreening of openGL widget.
		qDebug() << "toggle fullscreen";
	}

	// Pause/UnPause button
	void mainwindow::on_pauseButton_clicked()
	{
		emit ui->openGLWidget->toggle_paused_integrator();
		update_pauseButton();
	}

	// Update the text on the pause button to reflect the paused state.
	void mainwindow::update_pauseButton()
	{
		if (ui->openGLWidget->is_paused())
			ui->pauseButton->setText(QApplication::translate("Waves::mainwindowClass", "UnPause", 0));
		else
			ui->pauseButton->setText(QApplication::translate("Waves::mainwindowClass", "Pause", 0));
	}

}
