#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "simparameters.h"
#include <QThread>
#include <QTimer>

class MainWindow;
class Simulation;

class Controller : public QThread
{
    Q_OBJECT

public:
    Controller(int fps);
    virtual ~Controller();
    void initialize(MainWindow *mw);
    void render();

public slots:
    void reset();
    void clearScene();
    void updateParameters(SimParameters params);
    void leftMouseClicked(double x, double y);
    void rightMouseClicked(double x, double y);
    void mouseDrag(double x, double y);
    void resetDrag();

    void simTick();

protected:
    virtual void run();

private:
    MainWindow *mw_;
    Simulation *sim_;
    SimParameters params_;
    double dragXOld;
    double dragYOld;
    bool dragOn;

    int fps_;
    QTimer simtimer_;
};

#endif // CONTROLLER_H
