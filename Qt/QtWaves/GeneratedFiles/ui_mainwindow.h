/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.9.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
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
            Waves__mainwindowClass->setObjectName(QStringLiteral("Waves__mainwindowClass"));
        Waves__mainwindowClass->setWindowModality(Qt::ApplicationModal);
        Waves__mainwindowClass->resize(1024, 600);
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
        IC_comboBox = new QComboBox(centralWidget);
        IC_comboBox->setObjectName(QStringLiteral("IC_comboBox"));

        verticalLayout->addWidget(IC_comboBox);

        BC_comboBox = new QComboBox(centralWidget);
        BC_comboBox->setObjectName(QStringLiteral("BC_comboBox"));

        verticalLayout->addWidget(BC_comboBox);

        pause_Button = new QPushButton(centralWidget);
        pause_Button->setObjectName(QStringLiteral("pause_Button"));

        verticalLayout->addWidget(pause_Button);

        timestep_doubleSpinBox = new QDoubleSpinBox(centralWidget);
        timestep_doubleSpinBox->setObjectName(QStringLiteral("timestep_doubleSpinBox"));
        timestep_doubleSpinBox->setAlignment(Qt::AlignCenter);
        timestep_doubleSpinBox->setDecimals(4);
        timestep_doubleSpinBox->setMinimum(0);
        timestep_doubleSpinBox->setMaximum(0.1);
        timestep_doubleSpinBox->setSingleStep(0.0001);
        timestep_doubleSpinBox->setValue(0.0167);

        verticalLayout->addWidget(timestep_doubleSpinBox);

        time_lineEdit = new QLineEdit(centralWidget);
        time_lineEdit->setObjectName(QStringLiteral("time_lineEdit"));
        time_lineEdit->setAlignment(Qt::AlignCenter);
        time_lineEdit->setReadOnly(true);

        verticalLayout->addWidget(time_lineEdit);

        dissipation_doubleSpinBox = new QDoubleSpinBox(centralWidget);
        dissipation_doubleSpinBox->setObjectName(QStringLiteral("dissipation_doubleSpinBox"));
        dissipation_doubleSpinBox->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);
        dissipation_doubleSpinBox->setDecimals(3);
        dissipation_doubleSpinBox->setMaximum(1);
        dissipation_doubleSpinBox->setSingleStep(0.001);

        verticalLayout->addWidget(dissipation_doubleSpinBox);

        fullscreen_Button = new QPushButton(centralWidget);
        fullscreen_Button->setObjectName(QStringLiteral("fullscreen_Button"));

        verticalLayout->addWidget(fullscreen_Button);

        reset_view_Button = new QPushButton(centralWidget);
        reset_view_Button->setObjectName(QStringLiteral("reset_view_Button"));

        verticalLayout->addWidget(reset_view_Button);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        quit_Button = new QPushButton(centralWidget);
        quit_Button->setObjectName(QStringLiteral("quit_Button"));

        verticalLayout->addWidget(quit_Button);


        horizontalLayout->addLayout(verticalLayout);

        Waves__mainwindowClass->setCentralWidget(centralWidget);
        mainToolBar = new QToolBar(Waves__mainwindowClass);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        Waves__mainwindowClass->addToolBar(Qt::TopToolBarArea, mainToolBar);

        retranslateUi(Waves__mainwindowClass);

        IC_comboBox->setCurrentIndex(0);
        BC_comboBox->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(Waves__mainwindowClass);
    } // setupUi

    void retranslateUi(QMainWindow *Waves__mainwindowClass)
    {
        Waves__mainwindowClass->setWindowTitle(QApplication::translate("Waves::mainwindowClass", "QtGL", Q_NULLPTR));
        IC_comboBox->clear();
        IC_comboBox->insertItems(0, QStringList()
         << QApplication::translate("Waves::mainwindowClass", "Single frequency", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Continuous spectrum", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Localised hump", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Two localised humps", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Off-centre localised hump", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "A localised wave", Q_NULLPTR)
        );
#ifndef QT_NO_TOOLTIP
        IC_comboBox->setToolTip(QApplication::translate("Waves::mainwindowClass", "Initial conditions", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        IC_comboBox->setStatusTip(QString());
#endif // QT_NO_STATUSTIP
#ifndef QT_NO_WHATSTHIS
        IC_comboBox->setWhatsThis(QApplication::translate("Waves::mainwindowClass", "Initial conditions", Q_NULLPTR));
#endif // QT_NO_WHATSTHIS
        BC_comboBox->clear();
        BC_comboBox->insertItems(0, QStringList()
         << QApplication::translate("Waves::mainwindowClass", "Periodic square", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Dirichlet, square", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Dirichlet, circle", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Dirichlet, circle with a cusp", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Dirichlet, intersecting circles", Q_NULLPTR)
         << QApplication::translate("Waves::mainwindowClass", "Double slit (mixed BCs)", Q_NULLPTR)
        );
#ifndef QT_NO_TOOLTIP
        BC_comboBox->setToolTip(QApplication::translate("Waves::mainwindowClass", "Boundary Conditions", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_WHATSTHIS
        BC_comboBox->setWhatsThis(QApplication::translate("Waves::mainwindowClass", "Boundary Conditions", Q_NULLPTR));
#endif // QT_NO_WHATSTHIS
        pause_Button->setText(QApplication::translate("Waves::mainwindowClass", "Pause", Q_NULLPTR));
        timestep_doubleSpinBox->setPrefix(QApplication::translate("Waves::mainwindowClass", "Timestep: ", Q_NULLPTR));
        timestep_doubleSpinBox->setSuffix(QApplication::translate("Waves::mainwindowClass", "s", Q_NULLPTR));
        time_lineEdit->setPlaceholderText(QString());
#ifndef QT_NO_TOOLTIP
        dissipation_doubleSpinBox->setToolTip(QApplication::translate("Waves::mainwindowClass", "Artificial dissipation", Q_NULLPTR));
#endif // QT_NO_TOOLTIP
        dissipation_doubleSpinBox->setPrefix(QApplication::translate("Waves::mainwindowClass", "Dissipation factor: ", Q_NULLPTR));
        fullscreen_Button->setText(QApplication::translate("Waves::mainwindowClass", "Fullscreen", Q_NULLPTR));
        reset_view_Button->setText(QApplication::translate("Waves::mainwindowClass", "Reset view", Q_NULLPTR));
        quit_Button->setText(QApplication::translate("Waves::mainwindowClass", "Quit", Q_NULLPTR));
    } // retranslateUi

};

} // namespace Waves

namespace Waves {
namespace Ui {
    class mainwindowClass: public Ui_mainwindowClass {};
} // namespace Ui
} // namespace Waves

#endif // UI_MAINWINDOW_H
