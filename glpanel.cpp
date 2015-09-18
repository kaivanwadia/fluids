#include "glpanel.h"
#include "controller.h"
#include <QMouseEvent>

GLPanel::GLPanel(QWidget *parent) :
    QGLWidget(parent)
{
    cont_ = NULL;
}

void GLPanel::setController(Controller *cont)
{
    cont_ = cont;
}

void GLPanel::resizeGL(int , int )
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-1, 1, -1, 1, 0, 1);

    glMatrixMode(GL_MODELVIEW);

    glDisable(GL_DEPTH_TEST);

    glClearColor(1.0, 1.0, 1.0, 0.0);
}

void GLPanel::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT);
    cont_->render();
}

void GLPanel::mouseMoveEvent(QMouseEvent *me)
{
    double x = -1.0 + 2.0*double(me->x())/double(width());
    double y = 1.0 - 2.0*double(me->y())/double(height());

    QMetaObject::invokeMethod(cont_, "mouseDrag", Q_ARG(double, x), Q_ARG(double, y));
}

void GLPanel::mousePressEvent(QMouseEvent *me)
{
    if(me->button() == Qt::LeftButton)
    {
        double x = -1.0 + 2.0*double(me->x())/double(width());
        double y = 1.0 - 2.0*double(me->y())/double(height());
        QMetaObject::invokeMethod(cont_, "leftMouseClicked", Q_ARG(double, x), Q_ARG(double, y));
    }
    else if (me->button() == Qt::RightButton)
    {
        double x = -1.0 + 2.0*double(me->x())/double(width());
        double y = 1.0 - 2.0*double(me->y())/double(height());
        QMetaObject::invokeMethod(cont_, "rightMouseClicked", Q_ARG(double, x), Q_ARG(double, y));
    }
}

void GLPanel::mouseReleaseEvent(QMouseEvent *me)
{
    if(me->button() == Qt::LeftButton)
    {
        QMetaObject::invokeMethod(cont_, "resetDrag");
    }
}
