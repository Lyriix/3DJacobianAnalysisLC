#pragma once

#ifndef MY_WINDOW_HPP
#define MY_WINDOW_HPP



#include <QMainWindow>


//forward declaration
namespace Ui
{
class MainWindow;
}
class myWidgetGL;

/** Declaration of the Window class */
class myWindow: public QMainWindow
{
    Q_OBJECT

public:

    myWindow(QWidget *parent=NULL);
    ~myWindow();



private slots:

    /** Quit the application */
    void action_quit();
    /** Enable the drawing of the meshes */
    void action_draw();
    /** Set the Wireframe mode for the meshes */
    void action_wireframe();
    /** Enable the deformation to the deformation state */
    void action_deform_to_baseline();
    void action_deform_to_iop1();
    void action_deform_to_iop2();
    void action_deform_to_recovery();

private:

    /** Layout for the Window */
    Ui::MainWindow *ui;
    /** The OpenGL Widget */
    myWidgetGL *glWidget;





};

#endif
