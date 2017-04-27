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
		void on_pauseButton_clicked(); ///< Handle the pause button being clicked
		void update_pauseButton(); ///< Update the text on the pause button to reflect the current state
		void on_fullscreenButton_clicked(); ///< Handle the fullscreen button being clicked

	signals:
		// We pass various events directly to the underlying OpenGLWidget
		void on_ICButton_clicked(); ///< Change initial conditions signal
		void on_BCButton_clicked(); ///< Change boundary conditions signal
		void on_resetviewButton_clicked(); ///< Reset view signal
		void on_quitButton_clicked(); ///< Quit signal

	private:
		Ui::mainwindowClass *ui; ///< Main ui handle

	};
}
