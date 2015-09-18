#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>

class Controller;
struct SimParameters;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(Controller &cont, int fps, QWidget *parent = 0);
    ~MainWindow();   

public slots:
    void setUIFromParameters(const SimParameters &params);

private slots:
    void updateGL();

    void on_actionExit_triggered();

    void on_actionReset_Everything_triggered();

    void on_actionReset_triggered();

    void on_addDensityRadio_clicked();

    void on_addVelocityRadio_clicked();

    void on_viscosityConstantEdit_editingFinished();

    void on_diffusionConstantEdit_editingFinished();

    void on_startSimulationButton_clicked();

    void on_timeStepEdit_editingFinished();

    void on_newtonTolEdit_editingFinished();

    void on_newtonMaxItersEdit_editingFinished();

    void on_springStiffnessEdit_editingFinished();

    void on_maxStrainEdit_editingFinished();

    void on_dampingStiffnessEdit_editingFinished();

    void on_massEdit_editingFinished();

    void on_maxSpringDistEdit_editingFinished();

    void on_isFixedCheckBox_clicked();

    void on_densityRadiusEdit_editingFinished();
    void on_densityMagnitudeEdit_editingFinished();
    void on_velocityRadiusEdit_editingFinished();
    void on_velocityMagnitudeEdit_editingFinished();

    //-------------------------------

    void on_springButton_clicked();
    void on_rigidRodButton_clicked();
    void on_flexibleRodButton_clicked();
    void on_ropeButton_clicked();

    void on_rodDensityEdit_editingFinished();
    void on_rodStretchEdit_editingFinished();
    void on_rodBendEdit_editingFinished();
    void on_rodSegmentsEdit_editingFinished();
    void on_ropeDensityEdit_editingFinished();
    void on_ropeBendEdit_editingFinished();
    void on_ropeSegmentsEdit_editingFinished();

private:
    Controller &cont_;
    Ui::MainWindow *ui;
    bool simRunning_;
    QTimer renderTimer_;

    void setParametersFromUI();
};

#endif // MAINWINDOW_H
