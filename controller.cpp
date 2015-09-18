#include "controller.h"
#include "mainwindow.h"
#include "simulation.h"
#include <QDebug>
#include <iostream>

using namespace std;

Controller::Controller(int fps) : QThread(), mw_(NULL), fps_(fps)
{
}

Controller::~Controller()
{
    delete sim_;
}

void Controller::initialize(MainWindow *mw)
{
    mw_ = mw;
    sim_ = new Simulation(params_);

    this->dragXOld = 0;
    this->dragYOld = 0;
    this->dragOn = false;
}

void Controller::run()
{
    reset();
    connect(&simtimer_, SIGNAL(timeout()), this, SLOT(simTick()));
    simtimer_.start(1000/fps_);
    exec();
}

void Controller::reset()
{
    params_ = SimParameters();
    QMetaObject::invokeMethod(mw_, "setUIFromParameters", Q_ARG(SimParameters, params_));
    clearScene();
}

void Controller::clearScene()
{
    sim_->clearScene();
}

void Controller::updateParameters(SimParameters params)
{
    params_ = params;
}

void Controller::render()
{
    sim_->render();
}

void Controller::mouseDrag(double x, double y)
{
    if (!this->dragOn)
    {
        this->dragXOld = x;
        this->dragYOld = y;
        return;
    }
    double velToAdd = params_.velocityMagnitude;
    double velX,velY;

    velX = (x - this->dragXOld) * velToAdd;
    velY = -(y - this->dragYOld) * velToAdd;

    this->dragXOld = x;
    this->dragYOld = y;

    if(params_.simRunning)
    {
        switch(params_.clickMode)
        {

            case SimParameters::CM_ADDVELOCITY:
                sim_->addVelocity(x, y, velX, velY);
//                cout<<"Veloctiy"<<endl;
                break;
            case SimParameters::CM_ADDDENSITY:
                sim_->addDensity(x,y);
//                cout<<"Density"<<endl;
                break;
        }
    }
//    cout<<"Mouse Moved :"<<this->dragOn<<endl;
}

void Controller::resetDrag()
{
    this->dragXOld = 0;
    this->dragYOld = 0;
    this->dragOn = false;
//    cout<<"Mouse Released :"<<this->dragOn<<endl;
}

void Controller::leftMouseClicked(double x, double y)
{
    if (!this->dragOn)
    {
        this->dragOn = true;
    }
//    cout<<"Left Mouse Clicked :"<<this->dragOn<<endl;
}

void Controller::rightMouseClicked(double x, double y)
{
    sim_->addParticle(x, y);
//    cout<<"Right Mouse Clicked :"<<this->dragOn<<endl;
}

void Controller::simTick()
{
    if(params_.simRunning)
        sim_->takeSimulationStep();
}
