#pragma once
#include <QMainWindow>

#include "ui_mainwindow.h"

namespace Waves {
	/** The mainwindow class for UI elements.
	*/
	class mainwindow : public QMainWindow
	{
		Q_OBJECT

	public:
		explicit mainwindow(QWidget *parent = Q_NULLPTR);
		~mainwindow();

	public slots:
		void on_pauseButton_clicked();
		void update_pauseButton();
		void on_fullscreenButton_clicked();

	signals:
		// We pass various events directly to the underlying OpenGLWidget
		void on_ICButton_clicked();
		void on_BCButton_clicked();
		void on_resetviewButton_clicked();
		void on_quitButton_clicked();

	private:
		Ui::mainwindowClass *ui;

	};
}
