#pragma once
///@file

#include <QMainWindow>
#include "ui_mainwindow.h"

namespace Waves {
	/** The mainwindow class for UI elements.
	*/
	class mainwindow : public QMainWindow
	{
		Q_OBJECT

	public:
		explicit mainwindow(QWidget *parent = Q_NULLPTR); ///< Constructor
		~mainwindow(); ///< Destructor

	public slots:
		void on_reset_Button_clicked(); ///< Reset the integrator to the current initial conditions
		void on_pause_Button_clicked(); ///< Handle the pause button being clicked
		void update_pause_Button(); ///< Update the text on the pause button to reflect the current state
		void on_fullscreen_Button_clicked(); ///< Handle the fullscreen button being clicked

	signals:
		// We pass various events directly to the underlying OpenGLWidget
		void on_reset_view_Button_clicked(); ///< Reset view signal
		void on_quit_Button_clicked(); ///< Quit signal

	private:
		Ui::mainwindowClass *ui; ///< Main ui handle

	};
}
