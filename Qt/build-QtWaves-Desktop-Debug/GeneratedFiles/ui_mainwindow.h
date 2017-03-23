/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.7.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "oglwidget.h"

namespace Waves {

class Ui_mainwindowClass
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    OGLWidget *openGLWidget;
    QVBoxLayout *verticalLayout;
    QPushButton *ICButton;
    QPushButton *BCButton;
    QPushButton *pauseButton;
    QPushButton *fullscreenButton;
    QPushButton *resetviewButton;
    QSpacerItem *verticalSpacer;
    QPushButton *quitButton;
    QToolBar *mainToolBar;

    void setupUi(QMainWindow *Waves__mainwindowClass)
    {
        if (Waves__mainwindowClass->objectName().isEmpty())
            Waves__mainwindowClass->setObjectName(QStringLiteral("Waves__mainwindowClass"));
        Waves__mainwindowClass->setWindowModality(Qt::ApplicationModal);
        Waves__mainwindowClass->resize(930, 630);
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(1);
        sizePolicy.setVerticalStretch(1);
        sizePolicy.setHeightForWidth(Waves__mainwindowClass->sizePolicy().hasHeightForWidth());
        Waves__mainwindowClass->setSizePolicy(sizePolicy);
        centralWidget = new QWidget(Waves__mainwindowClass);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        openGLWidget = new OGLWidget(centralWidget);
        openGLWidget->setObjectName(QStringLiteral("openGLWidget"));
        QSizePolicy sizePolicy1(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy1.setHorizontalStretch(1);
        sizePolicy1.setVerticalStretch(1);
        sizePolicy1.setHeightForWidth(openGLWidget->sizePolicy().hasHeightForWidth());
        openGLWidget->setSizePolicy(sizePolicy1);
        openGLWidget->setMinimumSize(QSize(1, 1));
        openGLWidget->setSizeIncrement(QSize(1, 1));
        openGLWidget->setBaseSize(QSize(1, 1));

        horizontalLayout->addWidget(openGLWidget);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        ICButton = new QPushButton(centralWidget);
        ICButton->setObjectName(QStringLiteral("ICButton"));

        verticalLayout->addWidget(ICButton);

        BCButton = new QPushButton(centralWidget);
        BCButton->setObjectName(QStringLiteral("BCButton"));

        verticalLayout->addWidget(BCButton);

        pauseButton = new QPushButton(centralWidget);
        pauseButton->setObjectName(QStringLiteral("pauseButton"));

        verticalLayout->addWidget(pauseButton);

        fullscreenButton = new QPushButton(centralWidget);
        fullscreenButton->setObjectName(QStringLiteral("fullscreenButton"));

        verticalLayout->addWidget(fullscreenButton);

        resetviewButton = new QPushButton(centralWidget);
        resetviewButton->setObjectName(QStringLiteral("resetviewButton"));

        verticalLayout->addWidget(resetviewButton);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        quitButton = new QPushButton(centralWidget);
        quitButton->setObjectName(QStringLiteral("quitButton"));

        verticalLayout->addWidget(quitButton);


        horizontalLayout->addLayout(verticalLayout);

        Waves__mainwindowClass->setCentralWidget(centralWidget);
        mainToolBar = new QToolBar(Waves__mainwindowClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        Waves__mainwindowClass->addToolBar(Qt::TopToolBarArea, mainToolBar);

        retranslateUi(Waves__mainwindowClass);

        QMetaObject::connectSlotsByName(Waves__mainwindowClass);
    } // setupUi

    void retranslateUi(QMainWindow *Waves__mainwindowClass)
    {
        Waves__mainwindowClass->setWindowTitle(QApplication::translate("Waves::mainwindowClass", "QtGL", Q_NULLPTR));
        ICButton->setText(QApplication::translate("Waves::mainwindowClass", "Initial Condition", Q_NULLPTR));
        BCButton->setText(QApplication::translate("Waves::mainwindowClass", "Boundary Condition", Q_NULLPTR));
        pauseButton->setText(QApplication::translate("Waves::mainwindowClass", "Pause", Q_NULLPTR));
        fullscreenButton->setText(QApplication::translate("Waves::mainwindowClass", "Fullscreen", Q_NULLPTR));
        resetviewButton->setText(QApplication::translate("Waves::mainwindowClass", "Reset view", Q_NULLPTR));
        quitButton->setText(QApplication::translate("Waves::mainwindowClass", "Quit", Q_NULLPTR));
    } // retranslateUi

};

} // namespace Waves

namespace Waves {
namespace Ui {
    class mainwindowClass: public Ui_mainwindowClass {};
} // namespace Ui
} // namespace Waves

#endif // UI_MAINWINDOW_H
