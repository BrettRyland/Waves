/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.15.13
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLineEdit>
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
    QComboBox *IC_comboBox;
    QComboBox *BC_comboBox;
    QPushButton *pause_Button;
    QDoubleSpinBox *timestep_doubleSpinBox;
    QLineEdit *time_lineEdit;
    QDoubleSpinBox *dissipation_doubleSpinBox;
    QPushButton *fullscreen_Button;
    QPushButton *reset_view_Button;
    QSpacerItem *verticalSpacer;
    QPushButton *quit_Button;
    QToolBar *mainToolBar;

    void setupUi(QMainWindow *Waves__mainwindowClass)
    {
        if (Waves__mainwindowClass->objectName().isEmpty())
            Waves__mainwindowClass->setObjectName(QString::fromUtf8("Waves__mainwindowClass"));
        Waves__mainwindowClass->setWindowModality(Qt::ApplicationModal);
        Waves__mainwindowClass->resize(1024, 600);
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(1);
        sizePolicy.setVerticalStretch(1);
        sizePolicy.setHeightForWidth(Waves__mainwindowClass->sizePolicy().hasHeightForWidth());
        Waves__mainwindowClass->setSizePolicy(sizePolicy);
        centralWidget = new QWidget(Waves__mainwindowClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        openGLWidget = new OGLWidget(centralWidget);
        openGLWidget->setObjectName(QString::fromUtf8("openGLWidget"));
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
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        IC_comboBox = new QComboBox(centralWidget);
        IC_comboBox->addItem(QString());
        IC_comboBox->addItem(QString());
        IC_comboBox->addItem(QString());
        IC_comboBox->addItem(QString());
        IC_comboBox->addItem(QString());
        IC_comboBox->addItem(QString());
        IC_comboBox->setObjectName(QString::fromUtf8("IC_comboBox"));

        verticalLayout->addWidget(IC_comboBox);

        BC_comboBox = new QComboBox(centralWidget);
        BC_comboBox->addItem(QString());
        BC_comboBox->addItem(QString());
        BC_comboBox->addItem(QString());
        BC_comboBox->addItem(QString());
        BC_comboBox->addItem(QString());
        BC_comboBox->addItem(QString());
        BC_comboBox->setObjectName(QString::fromUtf8("BC_comboBox"));

        verticalLayout->addWidget(BC_comboBox);

        pause_Button = new QPushButton(centralWidget);
        pause_Button->setObjectName(QString::fromUtf8("pause_Button"));

        verticalLayout->addWidget(pause_Button);

        timestep_doubleSpinBox = new QDoubleSpinBox(centralWidget);
        timestep_doubleSpinBox->setObjectName(QString::fromUtf8("timestep_doubleSpinBox"));
        timestep_doubleSpinBox->setAlignment(Qt::AlignCenter);
        timestep_doubleSpinBox->setDecimals(4);
        timestep_doubleSpinBox->setMinimum(0.000000000000000);
        timestep_doubleSpinBox->setMaximum(0.100000000000000);
        timestep_doubleSpinBox->setSingleStep(0.000100000000000);
        timestep_doubleSpinBox->setValue(0.016700000000000);

        verticalLayout->addWidget(timestep_doubleSpinBox);

        time_lineEdit = new QLineEdit(centralWidget);
        time_lineEdit->setObjectName(QString::fromUtf8("time_lineEdit"));
        time_lineEdit->setAlignment(Qt::AlignCenter);
        time_lineEdit->setReadOnly(true);

        verticalLayout->addWidget(time_lineEdit);

        dissipation_doubleSpinBox = new QDoubleSpinBox(centralWidget);
        dissipation_doubleSpinBox->setObjectName(QString::fromUtf8("dissipation_doubleSpinBox"));
        dissipation_doubleSpinBox->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        dissipation_doubleSpinBox->setDecimals(3);
        dissipation_doubleSpinBox->setMaximum(1.000000000000000);
        dissipation_doubleSpinBox->setSingleStep(0.001000000000000);

        verticalLayout->addWidget(dissipation_doubleSpinBox);

        fullscreen_Button = new QPushButton(centralWidget);
        fullscreen_Button->setObjectName(QString::fromUtf8("fullscreen_Button"));

        verticalLayout->addWidget(fullscreen_Button);

        reset_view_Button = new QPushButton(centralWidget);
        reset_view_Button->setObjectName(QString::fromUtf8("reset_view_Button"));

        verticalLayout->addWidget(reset_view_Button);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        quit_Button = new QPushButton(centralWidget);
        quit_Button->setObjectName(QString::fromUtf8("quit_Button"));

        verticalLayout->addWidget(quit_Button);


        horizontalLayout->addLayout(verticalLayout);

        Waves__mainwindowClass->setCentralWidget(centralWidget);
        mainToolBar = new QToolBar(Waves__mainwindowClass);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        Waves__mainwindowClass->addToolBar(Qt::TopToolBarArea, mainToolBar);

        retranslateUi(Waves__mainwindowClass);

        IC_comboBox->setCurrentIndex(0);
        BC_comboBox->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(Waves__mainwindowClass);
    } // setupUi

    void retranslateUi(QMainWindow *Waves__mainwindowClass)
    {
        Waves__mainwindowClass->setWindowTitle(QCoreApplication::translate("Waves::mainwindowClass", "QtGL", nullptr));
        IC_comboBox->setItemText(0, QCoreApplication::translate("Waves::mainwindowClass", "Single frequency", nullptr));
        IC_comboBox->setItemText(1, QCoreApplication::translate("Waves::mainwindowClass", "Continuous spectrum", nullptr));
        IC_comboBox->setItemText(2, QCoreApplication::translate("Waves::mainwindowClass", "Localised hump", nullptr));
        IC_comboBox->setItemText(3, QCoreApplication::translate("Waves::mainwindowClass", "Two localised humps", nullptr));
        IC_comboBox->setItemText(4, QCoreApplication::translate("Waves::mainwindowClass", "Off-centre localised hump", nullptr));
        IC_comboBox->setItemText(5, QCoreApplication::translate("Waves::mainwindowClass", "A localised wave", nullptr));

#if QT_CONFIG(tooltip)
        IC_comboBox->setToolTip(QCoreApplication::translate("Waves::mainwindowClass", "Initial conditions", nullptr));
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        IC_comboBox->setStatusTip(QString());
#endif // QT_CONFIG(statustip)
#if QT_CONFIG(whatsthis)
        IC_comboBox->setWhatsThis(QCoreApplication::translate("Waves::mainwindowClass", "Initial conditions", nullptr));
#endif // QT_CONFIG(whatsthis)
        BC_comboBox->setItemText(0, QCoreApplication::translate("Waves::mainwindowClass", "Periodic square", nullptr));
        BC_comboBox->setItemText(1, QCoreApplication::translate("Waves::mainwindowClass", "Dirichlet, square", nullptr));
        BC_comboBox->setItemText(2, QCoreApplication::translate("Waves::mainwindowClass", "Dirichlet, circle", nullptr));
        BC_comboBox->setItemText(3, QCoreApplication::translate("Waves::mainwindowClass", "Dirichlet, circle with a cusp", nullptr));
        BC_comboBox->setItemText(4, QCoreApplication::translate("Waves::mainwindowClass", "Dirichlet, intersecting circles", nullptr));
        BC_comboBox->setItemText(5, QCoreApplication::translate("Waves::mainwindowClass", "Double slit (mixed BCs)", nullptr));

#if QT_CONFIG(tooltip)
        BC_comboBox->setToolTip(QCoreApplication::translate("Waves::mainwindowClass", "Boundary Conditions", nullptr));
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(whatsthis)
        BC_comboBox->setWhatsThis(QCoreApplication::translate("Waves::mainwindowClass", "Boundary Conditions", nullptr));
#endif // QT_CONFIG(whatsthis)
        pause_Button->setText(QCoreApplication::translate("Waves::mainwindowClass", "Pause", nullptr));
        timestep_doubleSpinBox->setPrefix(QCoreApplication::translate("Waves::mainwindowClass", "Timestep: ", nullptr));
        timestep_doubleSpinBox->setSuffix(QCoreApplication::translate("Waves::mainwindowClass", "s", nullptr));
        time_lineEdit->setPlaceholderText(QString());
#if QT_CONFIG(tooltip)
        dissipation_doubleSpinBox->setToolTip(QCoreApplication::translate("Waves::mainwindowClass", "Artificial dissipation", nullptr));
#endif // QT_CONFIG(tooltip)
        dissipation_doubleSpinBox->setPrefix(QCoreApplication::translate("Waves::mainwindowClass", "Dissipation factor: ", nullptr));
        fullscreen_Button->setText(QCoreApplication::translate("Waves::mainwindowClass", "Fullscreen", nullptr));
        reset_view_Button->setText(QCoreApplication::translate("Waves::mainwindowClass", "Reset view", nullptr));
        quit_Button->setText(QCoreApplication::translate("Waves::mainwindowClass", "Quit", nullptr));
    } // retranslateUi

};

} // namespace Waves

namespace Waves {
namespace Ui {
    class mainwindowClass: public Ui_mainwindowClass {};
} // namespace Ui
} // namespace Waves

#endif // UI_MAINWINDOW_H
