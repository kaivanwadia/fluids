#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "simparameters.h"
#include "controller.h"

MainWindow::MainWindow(Controller &cont, int fps, QWidget *parent) :
    QMainWindow(parent),
    cont_(cont),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->GLWidget->setController(&cont);
    simRunning_ = true;
    connect(&renderTimer_, SIGNAL(timeout()), this, SLOT(updateGL()));
    renderTimer_.start(1000/fps);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionExit_triggered()
{
    close();
}

void MainWindow::setParametersFromUI()
{
    SimParameters params;

    params.simRunning = simRunning_;
    if(ui->addDensityRadio->isChecked())
        params.clickMode = SimParameters::CM_ADDDENSITY;
    else if(ui->addVelocityRadio->isChecked())
        params.clickMode = SimParameters::CM_ADDVELOCITY;

    params.densityRadius = ui->densityRadiusEdit->text().toInt();
    params.densityMagnitude = ui->densityMagnitudeEdit->text().toDouble();
    params.velocityRadius = ui->velocityRadiusEdit->text().toInt();
    params.velocityMagnitude = ui->velocityMagnitudeEdit->text().toDouble();
    params.viscosityFluid = ui->viscosityConstantEdit->text().toDouble();
    params.diffusionConstant = ui->diffusionConstantEdit->text().toDouble();

    params.timeStep = ui->timeStepEdit->text().toDouble();
    params.NewtonTolerance = ui->newtonTolEdit->text().toDouble();
    params.NewtonMaxIters = ui->newtonMaxItersEdit->text().toInt();

    params.springStiffness = ui->springStiffnessEdit->text().toDouble();
    params.maxSpringStrain = ui->maxStrainEdit->text().toDouble();
    params.dampingStiffness = ui->dampingStiffnessEdit->text().toDouble();

    params.particleMass = ui->massEdit->text().toDouble();
    params.particleFixed = ui->isFixedCheckBox->isChecked();
    params.maxSpringDist = ui->maxSpringDistEdit->text().toDouble();

    if(ui->springButton->isChecked())
        params.connector = SimParameters::CT_SPRING;
    else if(ui->rigidRodButton->isChecked())
        params.connector = SimParameters:: CT_RIGID_ROD;
    else if(ui->flexibleRodButton->isChecked())
        params.connector = SimParameters:: CT_FLEXIBLE_ROD;
    else if(ui->ropeButton->isChecked())
        params.connector = SimParameters:: CT_ROPE;

    params.rodDensity = ui->rodDensityEdit->text().toDouble();
    params.rodStretchStiffness = ui->rodStretchEdit->text().toDouble();
    params.rodBendingStiffness = ui->rodBendEdit->text().toDouble();
    params.rodSegments = ui->rodSegmentsEdit->text().toInt();

    params.ropeDensity = ui->ropeDensityEdit->text().toDouble();
    params.ropeBend = ui->ropeBendEdit->text().toDouble();
    params.ropeSegments = ui->ropeSegmentsEdit->text().toInt();

    setUIFromParameters(params);
    QMetaObject::invokeMethod(&cont_, "updateParameters", Q_ARG(SimParameters, params));
}

void MainWindow::setUIFromParameters(const SimParameters &params)
{
    if(params.simRunning)
    {
        ui->startSimulationButton->setText(QString("Pause Simulation"));
        simRunning_ = true;
    }
    else
    {
        ui->startSimulationButton->setText(QString("Start Simulation"));
        simRunning_ = false;
    }

    switch(params.clickMode)
    {
        case SimParameters::CM_ADDDENSITY:
            ui->addDensityRadio->setChecked(true);
            break;
        case SimParameters::CM_ADDVELOCITY:
            ui->addVelocityRadio->setChecked(true);
            break;
    }
    ui->densityRadiusEdit->setText(QString::number(params.densityRadius));
    ui->densityMagnitudeEdit->setText(QString::number(params.densityMagnitude));
    ui->velocityRadiusEdit->setText(QString::number(params.velocityRadius));
    ui->velocityMagnitudeEdit->setText(QString::number(params.velocityMagnitude));
    ui->viscosityConstantEdit->setText(QString::number(params.viscosityFluid));
    ui->diffusionConstantEdit->setText(QString::number(params.diffusionConstant));

    ui->timeStepEdit->setText(QString::number(params.timeStep));
    ui->newtonTolEdit->setText(QString::number(params.NewtonTolerance));
    ui->newtonMaxItersEdit->setText(QString::number(params.NewtonMaxIters));

    ui->springStiffnessEdit->setText(QString::number(params.springStiffness));
    ui->maxStrainEdit->setText(QString::number(params.maxSpringStrain));
    ui->dampingStiffnessEdit->setText(QString::number(params.dampingStiffness));

    ui->massEdit->setText(QString::number(params.particleMass));
    ui->isFixedCheckBox->setChecked(params.particleFixed);
    ui->maxSpringDistEdit->setText(QString::number(params.maxSpringDist));

    if(params.connector == SimParameters:: CT_SPRING)
        ui->springButton->setChecked(true);
    else if(params.connector == SimParameters:: CT_RIGID_ROD)
        ui->rigidRodButton->setChecked(true);
    else if(params.connector == SimParameters:: CT_FLEXIBLE_ROD)
        ui->flexibleRodButton->setChecked(true);
    else if(params.connector == SimParameters:: CT_ROPE)
        ui->ropeButton->setChecked(true);

    ui->rodDensityEdit->setText(QString::number(params.rodDensity));
    ui->rodStretchEdit->setText(QString::number(params.rodStretchStiffness));
    ui->rodBendEdit->setText(QString::number(params.rodBendingStiffness));
    ui->rodSegmentsEdit->setText(QString::number(params.rodSegments));

    ui->ropeDensityEdit->setText(QString::number(params.ropeDensity));
    ui->ropeBendEdit->setText(QString::number(params.ropeBend));
    ui->ropeSegmentsEdit->setText(QString::number(params.ropeSegments));
}

void MainWindow::updateGL()
{
    ui->GLWidget->update();
}

void MainWindow::on_actionReset_Everything_triggered()
{
    QMetaObject::invokeMethod(&cont_, "reset");
}

void MainWindow::on_actionReset_triggered()
{
    QMetaObject::invokeMethod(&cont_, "clearScene");
}

void MainWindow::on_addDensityRadio_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_addVelocityRadio_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_viscosityConstantEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_diffusionConstantEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_densityRadiusEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_densityMagnitudeEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_velocityRadiusEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_velocityMagnitudeEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_startSimulationButton_clicked()
{
    simRunning_ = !simRunning_;
    setParametersFromUI();
}

void MainWindow::on_timeStepEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_newtonTolEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_newtonMaxItersEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_springStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_maxStrainEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_dampingStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_massEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_maxSpringDistEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_isFixedCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_springButton_clicked()
{
    setParametersFromUI();
}
void MainWindow::on_rigidRodButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_flexibleRodButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_ropeButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_rodDensityEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_rodStretchEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_rodBendEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_rodSegmentsEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_ropeDensityEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_ropeBendEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_ropeSegmentsEdit_editingFinished()
{
    setParametersFromUI();
}
