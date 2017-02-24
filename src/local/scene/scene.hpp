
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

    /** Set the pointer to the parent Widget */
    void set_widget(myWidgetGL* widget_param);

    /** Generate  a surface */

    void generate_scene(float xmin, float xmax, float ymin, float ymax, int Nu, int Nv, float r, float x_center, float y_center);

    /** Generate a tube */
    void generate_tube(float r, float x_center, float y_center, float Nu, float Nv);

    /** Analyze the deformation field from the Csv */

    void analyzeCsv(const std::vector<std::vector<std::string>>& csvfile);

    /** Structure to analyze the csv */
    struct deformationArrays {
        std::vector<float> baseline;  //mean of baseline a and baseline b
        std::vector<float> iop1;      //mean of iop1 and iop2
        std::vector<float> iop2;      //...
        std::vector<float> recovery;  //...
    };
    void animation(deformationArrays &deformation);
    void applyDeformation( std::vector<float>& deformation);

private:

    /** Load a texture from a given file and returns its id */
    GLuint load_texture_file(std::string const& filename);

    /** Access to the parent object */
    myWidgetGL* pwidget;

    /** Default id for the texture (white texture) */
    GLuint texture_default;

    /** The id of the shader do draw meshes */
    GLuint shader_program_id;


    // Data of the scene


    cpe::mesh mesh_test;
    cpe::mesh_opengl mesh_test_opengl;

    cpe::mesh mesh_tube;
    cpe::mesh_opengl mesh_tube_opengl;


    /** Structure to analyze the csv */
    //We are gonna compute the means
    deformationArrays deformation;

    bool baselineAnimation = false;
    /** Variable for Dragon */
    float tps=0.0f;
    float op = 1.0f;

};

#endif
