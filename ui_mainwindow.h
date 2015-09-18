/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
<<<<<<< HEAD
** Created: Wed Dec 3 20:42:57 2014
=======
** Created: Wed Dec 3 20:43:58 2014
>>>>>>> 8bbf8275bf7ac810912bf2fb7eb5da11d2b9e72b
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "glpanel.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionExit;
    QAction *actionReset;
    QAction *actionReset_Everything;
    QWidget *centralWidget;
    GLPanel *GLWidget;
    QFrame *parameterFrame;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QGroupBox *groupBox;
    QGroupBox *groupBox_2;
    QRadioButton *addVelocityRadio;
    QRadioButton *addDensityRadio;
    QLabel *label_14;
    QLineEdit *densityRadiusEdit;
    QLabel *label_15;
    QLineEdit *velocityRadiusEdit;
    QLabel *label;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QLabel *label_6;
    QLabel *label_7;
    QLabel *label_8;
    QLabel *label_9;
    QLabel *label_10;
    QLineEdit *rodDensityEdit;
    QLineEdit *rodStretchEdit;
    QLineEdit *rodBendEdit;
    QLineEdit *rodSegmentsEdit;
    QLineEdit *ropeDensityEdit;
    QLineEdit *ropeBendEdit;
    QLineEdit *ropeSegmentsEdit;
    QLabel *label_11;
    QLabel *label_12;
    QLabel *label_13;
    QLineEdit *massEdit;
    QLineEdit *maxSpringDistEdit;
    QLabel *label_16;
    QLineEdit *springStiffnessEdit;
    QLabel *label_17;
    QLineEdit *maxStrainEdit;
    QLineEdit *dampingStiffnessEdit;
    QLabel *label_18;
    QCheckBox *isFixedCheckBox;
    QPushButton *startSimulationButton;
    QLabel *label_19;
    QLineEdit *velocityMagnitudeEdit;
    QLabel *label_20;
    QLineEdit *densityMagnitudeEdit;
    QLabel *timeStepLabel;
    QLineEdit *timeStepEdit;
    QLabel *timeStepLabel_2;
    QLineEdit *newtonTolEdit;
    QLabel *timeStepLabel_3;
    QLineEdit *newtonMaxItersEdit;
    QGroupBox *groupBox_3;
    QRadioButton *ropeButton;
    QRadioButton *springButton;
    QRadioButton *flexibleRodButton;
    QRadioButton *rigidRodButton;
    QLabel *label_21;
    QLineEdit *viscosityConstantEdit;
    QLabel *label_22;
    QLineEdit *diffusionConstantEdit;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QMenu *menuScene;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1200, 800);
        actionExit = new QAction(MainWindow);
        actionExit->setObjectName(QString::fromUtf8("actionExit"));
        actionReset = new QAction(MainWindow);
        actionReset->setObjectName(QString::fromUtf8("actionReset"));
        actionReset_Everything = new QAction(MainWindow);
        actionReset_Everything->setObjectName(QString::fromUtf8("actionReset_Everything"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        GLWidget = new GLPanel(centralWidget);
        GLWidget->setObjectName(QString::fromUtf8("GLWidget"));
        GLWidget->setGeometry(QRect(10, 0, 731, 731));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(GLWidget->sizePolicy().hasHeightForWidth());
        GLWidget->setSizePolicy(sizePolicy);
        parameterFrame = new QFrame(centralWidget);
        parameterFrame->setObjectName(QString::fromUtf8("parameterFrame"));
        parameterFrame->setGeometry(QRect(749, -1, 441, 731));
        parameterFrame->setFrameShape(QFrame::StyledPanel);
        parameterFrame->setFrameShadow(QFrame::Raised);
        verticalLayoutWidget = new QWidget(parameterFrame);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(10, 0, 431, 731));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        groupBox = new QGroupBox(verticalLayoutWidget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox_2 = new QGroupBox(groupBox);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        groupBox_2->setGeometry(QRect(9, 29, 411, 121));
        addVelocityRadio = new QRadioButton(groupBox_2);
        addVelocityRadio->setObjectName(QString::fromUtf8("addVelocityRadio"));
        addVelocityRadio->setGeometry(QRect(20, 30, 116, 22));
        addDensityRadio = new QRadioButton(groupBox_2);
        addDensityRadio->setObjectName(QString::fromUtf8("addDensityRadio"));
        addDensityRadio->setGeometry(QRect(20, 50, 116, 22));
        label_14 = new QLabel(groupBox_2);
        label_14->setObjectName(QString::fromUtf8("label_14"));
        label_14->setGeometry(QRect(173, 11, 111, 17));
        densityRadiusEdit = new QLineEdit(groupBox_2);
        densityRadiusEdit->setObjectName(QString::fromUtf8("densityRadiusEdit"));
        densityRadiusEdit->setGeometry(QRect(288, 7, 81, 27));
        label_15 = new QLabel(groupBox_2);
        label_15->setObjectName(QString::fromUtf8("label_15"));
        label_15->setGeometry(QRect(173, 39, 111, 17));
        velocityRadiusEdit = new QLineEdit(groupBox_2);
        velocityRadiusEdit->setObjectName(QString::fromUtf8("velocityRadiusEdit"));
        velocityRadiusEdit->setGeometry(QRect(288, 35, 81, 27));
        label = new QLabel(groupBox);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(90, 360, 151, 17));
        QFont font;
        font.setBold(true);
        font.setWeight(75);
        label->setFont(font);
        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(190, 580, 111, 17));
        label_3->setFont(font);
        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(190, 430, 111, 17));
        label_4->setFont(font);
        label_5 = new QLabel(groupBox);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(10, 400, 66, 17));
        label_6 = new QLabel(groupBox);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(10, 430, 141, 17));
        label_7 = new QLabel(groupBox);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(190, 460, 66, 17));
        label_8 = new QLabel(groupBox);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(190, 520, 66, 17));
        label_9 = new QLabel(groupBox);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setGeometry(QRect(190, 490, 71, 17));
        label_10 = new QLabel(groupBox);
        label_10->setObjectName(QString::fromUtf8("label_10"));
        label_10->setGeometry(QRect(190, 550, 81, 17));
        rodDensityEdit = new QLineEdit(groupBox);
        rodDensityEdit->setObjectName(QString::fromUtf8("rodDensityEdit"));
        rodDensityEdit->setGeometry(QRect(270, 460, 81, 27));
        rodStretchEdit = new QLineEdit(groupBox);
        rodStretchEdit->setObjectName(QString::fromUtf8("rodStretchEdit"));
        rodStretchEdit->setGeometry(QRect(270, 490, 81, 27));
        rodBendEdit = new QLineEdit(groupBox);
        rodBendEdit->setObjectName(QString::fromUtf8("rodBendEdit"));
        rodBendEdit->setGeometry(QRect(270, 520, 81, 27));
        rodSegmentsEdit = new QLineEdit(groupBox);
        rodSegmentsEdit->setObjectName(QString::fromUtf8("rodSegmentsEdit"));
        rodSegmentsEdit->setGeometry(QRect(270, 550, 81, 27));
        ropeDensityEdit = new QLineEdit(groupBox);
        ropeDensityEdit->setObjectName(QString::fromUtf8("ropeDensityEdit"));
        ropeDensityEdit->setGeometry(QRect(271, 606, 81, 27));
        ropeBendEdit = new QLineEdit(groupBox);
        ropeBendEdit->setObjectName(QString::fromUtf8("ropeBendEdit"));
        ropeBendEdit->setGeometry(QRect(271, 636, 81, 27));
        ropeSegmentsEdit = new QLineEdit(groupBox);
        ropeSegmentsEdit->setObjectName(QString::fromUtf8("ropeSegmentsEdit"));
        ropeSegmentsEdit->setGeometry(QRect(271, 666, 81, 27));
        label_11 = new QLabel(groupBox);
        label_11->setObjectName(QString::fromUtf8("label_11"));
        label_11->setGeometry(QRect(190, 640, 71, 17));
        label_12 = new QLabel(groupBox);
        label_12->setObjectName(QString::fromUtf8("label_12"));
        label_12->setGeometry(QRect(190, 610, 66, 17));
        label_13 = new QLabel(groupBox);
        label_13->setObjectName(QString::fromUtf8("label_13"));
        label_13->setGeometry(QRect(190, 670, 81, 17));
        massEdit = new QLineEdit(groupBox);
        massEdit->setObjectName(QString::fromUtf8("massEdit"));
        massEdit->setGeometry(QRect(60, 400, 81, 27));
        maxSpringDistEdit = new QLineEdit(groupBox);
        maxSpringDistEdit->setObjectName(QString::fromUtf8("maxSpringDistEdit"));
        maxSpringDistEdit->setGeometry(QRect(60, 450, 81, 27));
        label_16 = new QLabel(groupBox);
        label_16->setObjectName(QString::fromUtf8("label_16"));
        label_16->setGeometry(QRect(5, 526, 81, 17));
        springStiffnessEdit = new QLineEdit(groupBox);
        springStiffnessEdit->setObjectName(QString::fromUtf8("springStiffnessEdit"));
        springStiffnessEdit->setGeometry(QRect(120, 492, 50, 27));
        label_17 = new QLabel(groupBox);
        label_17->setObjectName(QString::fromUtf8("label_17"));
        label_17->setGeometry(QRect(5, 496, 101, 17));
        maxStrainEdit = new QLineEdit(groupBox);
        maxStrainEdit->setObjectName(QString::fromUtf8("maxStrainEdit"));
        maxStrainEdit->setGeometry(QRect(120, 522, 50, 27));
        dampingStiffnessEdit = new QLineEdit(groupBox);
        dampingStiffnessEdit->setObjectName(QString::fromUtf8("dampingStiffnessEdit"));
        dampingStiffnessEdit->setGeometry(QRect(120, 552, 51, 27));
        label_18 = new QLabel(groupBox);
        label_18->setObjectName(QString::fromUtf8("label_18"));
        label_18->setGeometry(QRect(7, 555, 111, 17));
        isFixedCheckBox = new QCheckBox(groupBox);
        isFixedCheckBox->setObjectName(QString::fromUtf8("isFixedCheckBox"));
        isFixedCheckBox->setGeometry(QRect(190, 390, 97, 22));
        startSimulationButton = new QPushButton(groupBox);
        startSimulationButton->setObjectName(QString::fromUtf8("startSimulationButton"));
        startSimulationButton->setGeometry(QRect(10, 210, 181, 27));
        label_19 = new QLabel(groupBox);
        label_19->setObjectName(QString::fromUtf8("label_19"));
        label_19->setGeometry(QRect(149, 122, 141, 17));
        velocityMagnitudeEdit = new QLineEdit(groupBox);
        velocityMagnitudeEdit->setObjectName(QString::fromUtf8("velocityMagnitudeEdit"));
        velocityMagnitudeEdit->setGeometry(QRect(297, 119, 81, 27));
        label_20 = new QLabel(groupBox);
        label_20->setObjectName(QString::fromUtf8("label_20"));
        label_20->setGeometry(QRect(150, 95, 141, 17));
        densityMagnitudeEdit = new QLineEdit(groupBox);
        densityMagnitudeEdit->setObjectName(QString::fromUtf8("densityMagnitudeEdit"));
        densityMagnitudeEdit->setGeometry(QRect(297, 92, 81, 27));
        timeStepLabel = new QLabel(groupBox);
        timeStepLabel->setObjectName(QString::fromUtf8("timeStepLabel"));
        timeStepLabel->setGeometry(QRect(10, 250, 81, 21));
        timeStepEdit = new QLineEdit(groupBox);
        timeStepEdit->setObjectName(QString::fromUtf8("timeStepEdit"));
        timeStepEdit->setGeometry(QRect(170, 250, 51, 21));
        timeStepLabel_2 = new QLabel(groupBox);
        timeStepLabel_2->setObjectName(QString::fromUtf8("timeStepLabel_2"));
        timeStepLabel_2->setGeometry(QRect(10, 280, 141, 21));
        newtonTolEdit = new QLineEdit(groupBox);
        newtonTolEdit->setObjectName(QString::fromUtf8("newtonTolEdit"));
        newtonTolEdit->setGeometry(QRect(170, 280, 51, 21));
        timeStepLabel_3 = new QLabel(groupBox);
        timeStepLabel_3->setObjectName(QString::fromUtf8("timeStepLabel_3"));
        timeStepLabel_3->setGeometry(QRect(10, 310, 131, 21));
        newtonMaxItersEdit = new QLineEdit(groupBox);
        newtonMaxItersEdit->setObjectName(QString::fromUtf8("newtonMaxItersEdit"));
        newtonMaxItersEdit->setGeometry(QRect(170, 310, 51, 21));
        groupBox_3 = new QGroupBox(groupBox);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        groupBox_3->setGeometry(QRect(10, 580, 151, 121));
        ropeButton = new QRadioButton(groupBox_3);
        ropeButton->setObjectName(QString::fromUtf8("ropeButton"));
        ropeButton->setGeometry(QRect(0, 90, 116, 22));
        springButton = new QRadioButton(groupBox_3);
        springButton->setObjectName(QString::fromUtf8("springButton"));
        springButton->setGeometry(QRect(0, 30, 116, 22));
        flexibleRodButton = new QRadioButton(groupBox_3);
        flexibleRodButton->setObjectName(QString::fromUtf8("flexibleRodButton"));
        flexibleRodButton->setGeometry(QRect(0, 70, 116, 22));
        rigidRodButton = new QRadioButton(groupBox_3);
        rigidRodButton->setObjectName(QString::fromUtf8("rigidRodButton"));
        rigidRodButton->setGeometry(QRect(0, 50, 116, 22));
        label_21 = new QLabel(groupBox);
        label_21->setObjectName(QString::fromUtf8("label_21"));
        label_21->setGeometry(QRect(149, 149, 141, 17));
        viscosityConstantEdit = new QLineEdit(groupBox);
        viscosityConstantEdit->setObjectName(QString::fromUtf8("viscosityConstantEdit"));
        viscosityConstantEdit->setGeometry(QRect(298, 146, 81, 27));
        label_22 = new QLabel(groupBox);
        label_22->setObjectName(QString::fromUtf8("label_22"));
        label_22->setGeometry(QRect(150, 177, 141, 17));
        diffusionConstantEdit = new QLineEdit(groupBox);
        diffusionConstantEdit->setObjectName(QString::fromUtf8("diffusionConstantEdit"));
        diffusionConstantEdit->setGeometry(QRect(298, 174, 81, 27));

        verticalLayout->addWidget(groupBox);

        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1200, 25));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuScene = new QMenu(menuBar);
        menuScene->setObjectName(QString::fromUtf8("menuScene"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuBar->addAction(menuScene->menuAction());
        menuFile->addAction(actionExit);
        menuScene->addAction(actionReset);
        menuScene->addAction(actionReset_Everything);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Final Project - Fluids", 0, QApplication::UnicodeUTF8));
        actionExit->setText(QApplication::translate("MainWindow", "Exit", 0, QApplication::UnicodeUTF8));
        actionReset->setText(QApplication::translate("MainWindow", "Clear Scene", 0, QApplication::UnicodeUTF8));
        actionReset_Everything->setText(QApplication::translate("MainWindow", "Reset Everything", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("MainWindow", "Simulation Variables", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Addition Selector", 0, QApplication::UnicodeUTF8));
        addVelocityRadio->setText(QApplication::translate("MainWindow", "Velocity", 0, QApplication::UnicodeUTF8));
        addDensityRadio->setText(QApplication::translate("MainWindow", "Density", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("MainWindow", "Density Radius :", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("MainWindow", "Velocity Radius :", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Mass Spring System", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "Rope Settings", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("MainWindow", "Rod Settings", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("MainWindow", "Mass :", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("MainWindow", "Max Connector Dist :", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("MainWindow", "Density :", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("MainWindow", "Bend K :", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("MainWindow", "Stretch K :", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("MainWindow", "Segments :", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("MainWindow", "Bend K :", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("MainWindow", "Density :", 0, QApplication::UnicodeUTF8));
        label_13->setText(QApplication::translate("MainWindow", "Segments :", 0, QApplication::UnicodeUTF8));
        label_16->setText(QApplication::translate("MainWindow", "Max Strain :", 0, QApplication::UnicodeUTF8));
        label_17->setText(QApplication::translate("MainWindow", "Base Stiffness:", 0, QApplication::UnicodeUTF8));
        label_18->setText(QApplication::translate("MainWindow", "Damp Stiffness :", 0, QApplication::UnicodeUTF8));
        isFixedCheckBox->setText(QApplication::translate("MainWindow", "Is Fixed", 0, QApplication::UnicodeUTF8));
        startSimulationButton->setText(QApplication::translate("MainWindow", "Start Simulation", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("MainWindow", "Velocity Magnitude :", 0, QApplication::UnicodeUTF8));
        label_20->setText(QApplication::translate("MainWindow", "Density Magnitude :", 0, QApplication::UnicodeUTF8));
        densityMagnitudeEdit->setText(QString());
        timeStepLabel->setText(QApplication::translate("MainWindow", "Time Step:", 0, QApplication::UnicodeUTF8));
        timeStepLabel_2->setText(QApplication::translate("MainWindow", "Newton Tolerance :", 0, QApplication::UnicodeUTF8));
        timeStepLabel_3->setText(QApplication::translate("MainWindow", "Newton Max Iters :", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("MainWindow", "Connector Type", 0, QApplication::UnicodeUTF8));
        ropeButton->setText(QApplication::translate("MainWindow", "Rope", 0, QApplication::UnicodeUTF8));
        springButton->setText(QApplication::translate("MainWindow", "Spring", 0, QApplication::UnicodeUTF8));
        flexibleRodButton->setText(QApplication::translate("MainWindow", "Flexible Rod", 0, QApplication::UnicodeUTF8));
        rigidRodButton->setText(QApplication::translate("MainWindow", "Rigid Rod", 0, QApplication::UnicodeUTF8));
        label_21->setText(QApplication::translate("MainWindow", "Viscosity K:", 0, QApplication::UnicodeUTF8));
        label_22->setText(QApplication::translate("MainWindow", "Diffusion K:", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        menuScene->setTitle(QApplication::translate("MainWindow", "Scene", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
