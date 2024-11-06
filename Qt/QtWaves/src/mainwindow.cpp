///@file

#include <QDebug>
#include <QTimer>
#include "mainwindow.h"

namespace Waves
{
	mainwindow::mainwindow(QWidget *parent) : QMainWindow(parent),
																						ui(new Ui::mainwindowClass)
	{
		ui->setupUi(this);
		// Set up various signal/slot connections
		QObject::connect(this, &mainwindow::on_reset_view_Button_clicked, ui->openGLWidget, &OGLWidget::reset_view);
		QObject::connect(this, &mainwindow::on_quit_Button_clicked, ui->openGLWidget, &OGLWidget::quit);
		QObject::connect(ui->openGLWidget, &OGLWidget::toggle_fullscreen, this, &mainwindow::on_fullscreen_Button_clicked);
		QObject::connect(ui->openGLWidget, &OGLWidget::notify_paused_state, this, &mainwindow::update_pause_Button);
		QObject::connect(ui->openGLWidget, &OGLWidget::notify_IC_changed, [=](int ic)
										 { ui->IC_comboBox->setCurrentIndex(ic); });
		QObject::connect(ui->openGLWidget, &OGLWidget::notify_BC_changed, [=](int bc)
										 { ui->BC_comboBox->setCurrentIndex(bc); });
		QObject::connect(ui->openGLWidget, &OGLWidget::notify_time_step, [=](double dt)
										 { ui->timestep_doubleSpinBox->setValue(dt); });
		QObject::connect(ui->dissipation_doubleSpinBox, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [=](double d)
										 { emit ui->openGLWidget->change_dissipation(d); });
		QObject::connect(ui->timestep_doubleSpinBox, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [=](double d)
										 { emit ui->openGLWidget->change_timestep(d); });
		QObject::connect(ui->IC_comboBox, static_cast<void (QComboBox::*)(int)>(&QComboBox::activated), [=](int index)
										 { emit ui->openGLWidget->change_initial_conditions(index); });
		QObject::connect(ui->BC_comboBox, static_cast<void (QComboBox::*)(int)>(&QComboBox::activated), [=](int index)
										 { emit ui->openGLWidget->change_boundary_conditions(index); });
		QObject::connect(ui->openGLWidget, &OGLWidget::update_time, [=]()
										 {
											ui->time_lineEdit->setText(QString("Time: %1s").arg(QString::number(ui->openGLWidget->get_time())));
											ui->stats_lineEdit->setText(QString("FPS: %1, Steps/Frame: %2").arg(QString::number(ui->openGLWidget->get_fps(), 'f', 1), QString::number(ui->openGLWidget->get_spf(), 'f', 1))); });
		QObject::connect(ui->waveSpeed_doubleSpinBox, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [=](double d)
										 { emit ui->openGLWidget->change_wave_speed(d); });
		QObject::connect(ui->heightScale_doubleSpinBox, static_cast<void (QDoubleSpinBox::*)(double)>(&QDoubleSpinBox::valueChanged), [=](double d)
										 { emit ui->openGLWidget->change_height_scale(d); });
		// Set the initial state of the pause button
		update_pause_Button();
		ui->fullscreen_Button->setDisabled(true);
	}

	mainwindow::~mainwindow()
	{
		delete ui;
	}

	// Fullscreen button
	void mainwindow::on_fullscreen_Button_clicked()
	{
		// TODO Handle fullscreening of openGL widget.
		// qDebug() << "toggle fullscreen";
	}

	// Pause/UnPause button
	void mainwindow::on_pause_Button_clicked()
	{
		emit ui->openGLWidget->toggle_paused_integrator();
	}

	void mainwindow::on_reset_Button_clicked()
	{
		emit ui->openGLWidget->reset_integrator();
	}

	// Update the text on the pause button to reflect the paused state.
	void mainwindow::update_pause_Button()
	{
		if (ui->openGLWidget->is_paused())
			ui->pause_Button->setText(QApplication::translate("Waves::mainwindowClass", "UnPause", 0));
		else
			ui->pause_Button->setText(QApplication::translate("Waves::mainwindowClass", "Pause", 0));
	}

}
