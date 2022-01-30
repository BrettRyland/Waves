/********************************************************************************
** Form generated from reading UI file 'oglwidget.ui'
**
** Created by: Qt User Interface Compiler version 5.12.8
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_OGLWIDGET_H
#define UI_OGLWIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QOpenGLWidget>

QT_BEGIN_NAMESPACE

class Ui_OGLWidget
{
public:

    void setupUi(QOpenGLWidget *OGLWidget)
    {
        if (OGLWidget->objectName().isEmpty())
            OGLWidget->setObjectName(QString::fromUtf8("OGLWidget"));
        OGLWidget->setWindowModality(Qt::ApplicationModal);
        OGLWidget->resize(800, 600);

        retranslateUi(OGLWidget);

        QMetaObject::connectSlotsByName(OGLWidget);
    } // setupUi

    void retranslateUi(QOpenGLWidget *OGLWidget)
    {
        OGLWidget->setWindowTitle(QApplication::translate("OGLWidget", "OGLWidget", nullptr));
    } // retranslateUi

};

namespace Ui {
    class OGLWidget: public Ui_OGLWidget {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_OGLWIDGET_H
