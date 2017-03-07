
/** TP 4ETI - CPE Lyon - 2013/2014 */

#pragma once

#ifndef SCENE_HPP
#define SCENE_HPP

#include <GL/gl.h>
#include <GL/glew.h>

#include "../../lib/3d/mat3.hpp"
#include "../../lib/3d/vec3.hpp"
#include "../../lib/mesh/mesh.hpp"
#include "../../lib/opengl/mesh_opengl.hpp"
#include "../../lib/interface/camera_matrices.hpp"

#include <vector>


class myWidgetGL;

class scene
{
public:

    scene();



    /** \brief Method called only once at the beginning (load off files ...) */
    void load_scene();

    /** \brief Method called at every frame */
    void draw_scene();

    /** Method called to draw or not the tube */
    void change_draw_tube_state ();
    /** Method called to draw or not the grid */
    void change_draw_grid_state();

    /** Set the pointer to the parent Widget */
    void set_widget(myWidgetGL* widget_param);

    /** Generate  a surface */
    void generate_scene(float xmin, float xmax, float ymin, float ymax, int Nu, int Nv, float r, float x_center, float y_center);

    /** Get / Set for the Nu Nv Tube */
    int& set_Nu_tube() ;
    int& set_Nv_tube() ;
    /** Get / Set for Nu,Nv,Nw grid */
    int& set_Nu_grid() ;
    int& set_Nv_grid() ;
    int& set_Nw_grid() ;
    /** Generate a tube */
    void generate_tube(float r, float x_center, float y_center, float Nu, float Nv);
    cpe::vec3 getArcEnCielColor(int j, int inter);

    /** generate a grid */
    void generate_grid(int Nu, int Nv, int Nw);
    /** Generate a grid */





    /** Structure to analyze the csv "per slice" */
    struct deformationArrays {
        std::vector<float> baseline;  //mean of baseline a and baseline b
        std::vector<float> iop1;      //mean of iop1 and iop2
        std::vector<float> iop2;      //...
        std::vector<float> recovery;  //...
    };

    /** Analyze the deformation field from the Csv "per slice" */
    void analyzeCsv(const std::vector<std::vector<std::string>>& csvfile, deformationArrays &deform);

    /** Method called at every frame to perform animations on a draw object*/
    void animation();


    /** Method called by animation at every frame to perform the chosen animation on a TUBE */
    void applyDeformation(std::vector<float>& deformation, int axis);
    /** Method called by animation_grid to apply a deformation on a given GRID */
    void applyGridDeformation(std::vector<float> &deformation, int axis);


    /** Methods called to start animations */
    void setBaselineAnimation(){baselineAnimation = true;}
    void setIop1Animation(){iop1animation = true;}
    void setIop2Animation(){ iop2animation = true;}
    void setRecoveryAnimation() { recoveryAnimation = true; }


private:

    /** Load a texture from a given file and returns its id */
    GLuint load_texture_file(std::string const& filename);

    /** Access to the parent object */
    myWidgetGL* pwidget;

    /** Default id for the texture (white texture) */
    GLuint texture_default;

    /** The id of the shader do draw meshes */
    GLuint shader_program_id;

    /** Size of the Tube */
    int Nu_tube = 20;
    int Nv_tube = 600;
    /** Size of the grid */
    int Nu_grid ;
    int Nv_grid ;
    int Nw_grid ;

    // Data of the scene
    /** Boolean to chose which meshes display */
    bool draw_tube_state = false;
    bool draw_grid_state = false;

    cpe::mesh mesh_test;
    cpe::mesh_opengl mesh_test_opengl;

    cpe::mesh mesh_tube;
    cpe::mesh_opengl mesh_tube_opengl;

    cpe::mesh mesh_tube_displayed;
    cpe::mesh_opengl mesh_tube_displayed_opengl;

    cpe::mesh grid_d;
    cpe::mesh_opengl grid_d_opengl;

    cpe::mesh grid;
    cpe::mesh_opengl grid_opengl;

    /** Structure to analyze the csv */

    deformationArrays deformationYX;
    deformationArrays deformationYY;
    deformationArrays deformationYZ;

    float dmax = 0.0f;
    float dmin = 0.0f;

    bool baselineAnimation = false;
    bool iop1animation = false;
    bool iop2animation = false;
    bool recoveryAnimation = false;
    /** Variable for Dragon */
    float T = 1.0f;
    float tps=0.0f;
    float op = 1.0f;

};

#endif
