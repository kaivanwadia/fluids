/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Wed Dec 3 20:18:30 2014
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      34,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      19,   12,   11,   11, 0x0a,
      54,   11,   11,   11, 0x08,
      65,   11,   11,   11, 0x08,
      91,   11,   11,   11, 0x08,
     129,   11,   11,   11, 0x08,
     156,   11,   11,   11, 0x08,
     185,   11,   11,   11, 0x08,
     215,   11,   11,   11, 0x08,
     258,   11,   11,   11, 0x08,
     301,   11,   11,   11, 0x08,
     336,   11,   11,   11, 0x08,
     370,   11,   11,   11, 0x08,
     405,   11,   11,   11, 0x08,
     445,   11,   11,   11, 0x08,
     486,   11,   11,   11, 0x08,
     521,   11,   11,   11, 0x08,
     563,   11,   11,   11, 0x08,
     593,   11,   11,   11, 0x08,
     632,   11,   11,   11, 0x08,
     661,   11,   11,   11, 0x08,
     700,   11,   11,   11, 0x08,
     742,   11,   11,   11, 0x08,
     782,   11,   11,   11, 0x08,
     825,   11,   11,   11, 0x08,
     851,   11,   11,   11, 0x08,
     879,   11,   11,   11, 0x08,
     910,   11,   11,   11, 0x08,
     934,   11,   11,   11, 0x08,
     970,   11,   11,   11, 0x08,
    1006,   11,   11,   11, 0x08,
    1039,   11,   11,   11, 0x08,
    1076,   11,   11,   11, 0x08,
    1113,   11,   11,   11, 0x08,
    1147,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0params\0"
    "setUIFromParameters(SimParameters)\0"
    "updateGL()\0on_actionExit_triggered()\0"
    "on_actionReset_Everything_triggered()\0"
    "on_actionReset_triggered()\0"
    "on_addDensityRadio_clicked()\0"
    "on_addVelocityRadio_clicked()\0"
    "on_viscosityConstantEdit_editingFinished()\0"
    "on_diffusionConstantEdit_editingFinished()\0"
    "on_startSimulationButton_clicked()\0"
    "on_timeStepEdit_editingFinished()\0"
    "on_newtonTolEdit_editingFinished()\0"
    "on_newtonMaxItersEdit_editingFinished()\0"
    "on_springStiffnessEdit_editingFinished()\0"
    "on_maxStrainEdit_editingFinished()\0"
    "on_dampingStiffnessEdit_editingFinished()\0"
    "on_massEdit_editingFinished()\0"
    "on_maxSpringDistEdit_editingFinished()\0"
    "on_isFixedCheckBox_clicked()\0"
    "on_densityRadiusEdit_editingFinished()\0"
    "on_densityMagnitudeEdit_editingFinished()\0"
    "on_velocityRadiusEdit_editingFinished()\0"
    "on_velocityMagnitudeEdit_editingFinished()\0"
    "on_springButton_clicked()\0"
    "on_rigidRodButton_clicked()\0"
    "on_flexibleRodButton_clicked()\0"
    "on_ropeButton_clicked()\0"
    "on_rodDensityEdit_editingFinished()\0"
    "on_rodStretchEdit_editingFinished()\0"
    "on_rodBendEdit_editingFinished()\0"
    "on_rodSegmentsEdit_editingFinished()\0"
    "on_ropeDensityEdit_editingFinished()\0"
    "on_ropeBendEdit_editingFinished()\0"
    "on_ropeSegmentsEdit_editingFinished()\0"
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainWindow *_t = static_cast<MainWindow *>(_o);
        switch (_id) {
        case 0: _t->setUIFromParameters((*reinterpret_cast< const SimParameters(*)>(_a[1]))); break;
        case 1: _t->updateGL(); break;
        case 2: _t->on_actionExit_triggered(); break;
        case 3: _t->on_actionReset_Everything_triggered(); break;
        case 4: _t->on_actionReset_triggered(); break;
        case 5: _t->on_addDensityRadio_clicked(); break;
        case 6: _t->on_addVelocityRadio_clicked(); break;
        case 7: _t->on_viscosityConstantEdit_editingFinished(); break;
        case 8: _t->on_diffusionConstantEdit_editingFinished(); break;
        case 9: _t->on_startSimulationButton_clicked(); break;
        case 10: _t->on_timeStepEdit_editingFinished(); break;
        case 11: _t->on_newtonTolEdit_editingFinished(); break;
        case 12: _t->on_newtonMaxItersEdit_editingFinished(); break;
        case 13: _t->on_springStiffnessEdit_editingFinished(); break;
        case 14: _t->on_maxStrainEdit_editingFinished(); break;
        case 15: _t->on_dampingStiffnessEdit_editingFinished(); break;
        case 16: _t->on_massEdit_editingFinished(); break;
        case 17: _t->on_maxSpringDistEdit_editingFinished(); break;
        case 18: _t->on_isFixedCheckBox_clicked(); break;
        case 19: _t->on_densityRadiusEdit_editingFinished(); break;
        case 20: _t->on_densityMagnitudeEdit_editingFinished(); break;
        case 21: _t->on_velocityRadiusEdit_editingFinished(); break;
        case 22: _t->on_velocityMagnitudeEdit_editingFinished(); break;
        case 23: _t->on_springButton_clicked(); break;
        case 24: _t->on_rigidRodButton_clicked(); break;
        case 25: _t->on_flexibleRodButton_clicked(); break;
        case 26: _t->on_ropeButton_clicked(); break;
        case 27: _t->on_rodDensityEdit_editingFinished(); break;
        case 28: _t->on_rodStretchEdit_editingFinished(); break;
        case 29: _t->on_rodBendEdit_editingFinished(); break;
        case 30: _t->on_rodSegmentsEdit_editingFinished(); break;
        case 31: _t->on_ropeDensityEdit_editingFinished(); break;
        case 32: _t->on_ropeBendEdit_editingFinished(); break;
        case 33: _t->on_ropeSegmentsEdit_editingFinished(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 34)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 34;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
