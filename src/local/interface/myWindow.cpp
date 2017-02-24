#include "myWindow.hpp"

#include "myWidgetGL.hpp"
#include "../../lib/common/error_handling.hpp"
#include "ui_mainwindow.h"

#include <iostream>

myWindow::myWindow(QWidget *parent)
    :QMainWindow(parent),ui(new Ui::MainWindow)
{
    try
    {
        //Setup window layout
        ui->setupUi(this);

        //Create openGL context
        QGLFormat qglFormat;
        qglFormat.setVersion(1,2);

        //Create OpenGL Widget renderer
        glWidget=new myWidgetGL(qglFormat);

        //Add the OpenGL Widget into the layout
        ui->layout_scene->addWidget(glWidget);
    }
    catch(cpe::exception_cpe const& e)
    {
        std::cout<<std::endl<<e.report_exception()<<std::endl;
    }

    //Connect slot and signals
    connect(ui->quit,SIGNAL(clicked()),this,SLOT(action_quit()));
    connect(ui->draw,SIGNAL(clicked()),this,SLOT(action_draw()));
    connect(ui->wireframe,SIGNAL(clicked()),this,SLOT(action_wireframe()));
    connect(ui->baseline,SIGNAL(clicked()),this,SLOT(action_deform_to_baseline()));
    connect(ui->iop1,SIGNAL(clicked()),this,SLOT(action_deform_to_iop1()));
    connect(ui->iop2,SIGNAL(clicked()), this, SLOT(action_deform_to_iop2()));
    connect(ui->recovery,SIGNAL(clicked()), this, SLOT(action_deform_to_recovery()));
}

myWindow::~myWindow()
{}

void myWindow::action_quit()
{
    close();
}

void myWindow::action_draw()
{
    glWidget->change_draw_state();
}

void myWindow::action_wireframe()
{
    bool const state_wireframe=ui->wireframe->isChecked();
    glWidget->wireframe(state_wireframe);
}

void myWindow::action_deform_to_baseline()
{
    glWidget->deform_to_baseline();
}
void myWindow::action_deform_to_iop1()
{
    glWidget->deform_to_iop1();
}
void myWindow::action_deform_to_iop2()
{
    glWidget->deform_to_iop2();
}

void myWindow::action_deform_to_recovery()
{
    glWidget->deform_to_recovery();
}


