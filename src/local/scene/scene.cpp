#include <typeinfo>

#include <GL/glew.h>

#include "scene.hpp"

#include "../../lib/opengl/glutils.hpp"
#include "../../lib/perlin/perlin.hpp"
#include "../../lib/interface/camera_matrices.hpp"
#include "../interface/myWidgetGL.hpp"
#include "../../lib/mesh/mesh_io.hpp"

#include "../../lib/common/error_handling.hpp"
#include "../../lib/common/exception_cpe.hpp"
#include "../../readCSV.hpp"

#include <cmath>
#include <string>
#include <sstream>


using namespace cpe;


void scene::load_scene()
{

    //*****************************************//
    // Preload default structure               //
    //*****************************************//
    texture_default = load_texture_file("data/white.jpg");
    shader_program_id = read_shader("shaders/shader_mesh.vert",
                                    "shaders/shader_mesh.frag"); PRINT_OPENGL_ERROR();

    std::string csvName("/home/charly/workspace/OphtalmologyProject/JacobianAnalysis/3D/Modelisation/data/yvolution.csv");

    std::vector<std::vector<std::string>> bcd = readCsvFile(csvName);
    analyzeCsv(bcd);
    //For debugging only
    /*for (const auto & bc : bcd)
    {
        for(auto const & a : bc )
            std::cout << a ;
        std::cout <<std::endl;
    }*/


    //*****************************************//
    // Generate user defined mesh              //
    //*****************************************//

    float r=0.2f;
    float h=1.0f;
    int Nu=20; int Nv= 600;
    //scene::generate_scene(0.0f, 3.0f, 0.0f, 3.0f, Nu, Nv,r,0.8f,1.5f);
    scene::generate_tube(r,0.0f,0.0f,Nu,Nv);
    //baselineAnimation = true;
    //iop1animation = false;

}

void scene::draw_scene()
{


    //Setup uniform parameters
    glUseProgram(shader_program_id);                                                                           PRINT_OPENGL_ERROR();

    //Get cameras parameters (modelview,projection,normal).
    camera_matrices const& cam=pwidget->camera();

    //Set Uniform data to GPU
    glUniformMatrix4fv(get_uni_loc(shader_program_id,"camera_modelview"),1,false,cam.modelview.pointer());     PRINT_OPENGL_ERROR();
    glUniformMatrix4fv(get_uni_loc(shader_program_id,"camera_projection"),1,false,cam.projection.pointer());   PRINT_OPENGL_ERROR();
    glUniformMatrix4fv(get_uni_loc(shader_program_id,"normal_matrix"),1,false,cam.normal.pointer());           PRINT_OPENGL_ERROR();


    //Draw the meshes

   // mesh_tube_opengl.draw();
    mesh_tube_displayed_opengl.draw();

    animation(deformation);



}


scene::scene()
    :shader_program_id(0)
{}


GLuint scene::load_texture_file(std::string const& filename)
{
    return pwidget->load_texture_file(filename);
}

void scene::set_widget(myWidgetGL* widget_param)
{
    pwidget=widget_param;
}


void scene::generate_tube(float r, float x_center, float y_center, float Nu, float Nv)
{
    for ( int j = 0 ; j< Nv ; j++)
    {
        for( int i = 0 ; i < Nu ; i++)
        {

            float theta = i/(Nu)*2*M_PI ;
            float x = x_center + r*cos(theta);
            float y = y_center + r*sin(theta);
            float z = j/Nv;
            mesh_tube.add_vertex(vec3(x,y,z));
            mesh_tube.add_color(vec3(j/Nv,1-j/Nv,1-j/Nv));

        }
    }
    std::cout << "size mesh tube " << mesh_tube.size_vertex() << std::endl;

    for ( int j = 0 ; j< Nv-1 ; j++)
    {
        for( int i = 0 ; i < Nu-1 ; i++)
        {
            mesh_tube.add_triangle_index({i+Nu*j, (i+1)+Nu*(j+1), i+Nu*(j+1)});
            mesh_tube.add_triangle_index({i+Nu*j, (i+1)+Nu*(j+1), (i+1)+Nu*j});
        }
    }


    mesh_tube.fill_empty_field_by_default();
    mesh_tube_displayed = mesh_tube;
    mesh_tube_opengl.fill_vbo(mesh_tube);
    mesh_tube_displayed_opengl.fill_vbo(mesh_tube_displayed);
}

void scene::analyzeCsv(const std::vector<std::vector<std::string>>& csvfile)
{
    for(auto const & row : csvfile)
    {
        for(int i = 0; i < row.size() ; i+=2)
        {

            if( i == 0 )
                deformation.baseline.push_back(
                            (std::stof(row[1]) + std::stof(row[2])) / 2
                        - (std::stof(row[i]) + std::stof(row[1+i])) / 2 );
            if( i == 2 )
                deformation.iop1.push_back( (std::stof(row[1]) + std::stof(row[2])) / 2
                        - (std::stof(row[i]) + std::stof(row[1+i])) / 2);
            if( i == 4 )
                deformation.iop2.push_back( (std::stof(row[1]) + std::stof(row[2])) / 2
                        - (std::stof(row[i]) + std::stof(row[1+i])) / 2);
            if( i == 6 )
                deformation.recovery.push_back( (std::stof(row[1]) + std::stof(row[2])) / 2
                        - (std::stof(row[i]) + std::stof(row[1+i])) / 2);

            //std::cout << (std::stof(row[1]) + std::stof(row[2])) / 2
              //      - (std::stof(row[i]) + std::stof(row[1+i])) / 2 << std::endl;
            //std::cout << row.size() << std::endl;
            //std::cout << typeid(row[i]).name() << std::endl;
        }
    }

}

void scene::animation(deformationArrays &deformation)
{
    if(baselineAnimation){
        //std::cout << " to baseline" << std::endl;
        applyDeformation(deformation.baseline, baselineAnimation);
    }

    if(iop1animation)
    {
        //std::cout << "to IOP1" <<std::endl;
        applyDeformation(deformation.iop1, iop1animation);
    }
    if(iop2animation)
        applyDeformation(deformation.iop2, iop2animation);
    if(recoveryAnimation)
        applyDeformation(deformation.recovery, recoveryAnimation);
    //std::cout << baselineAnimation;
}

void scene::applyDeformation( std::vector<float> &deformation, bool &animationb)
{

    float timeStep = 0.02f;
    float T = 1.0f;

    //Parcourir les deformations
    //appliquer les deformations Y selon l'axe 0 : Nv
    int j=0;

    //temps is represented by tps in scene::scene.hpp

    //mesh_tube.vertex(100).x() += timeStep*deformation[j]/100;


     //std::cout << deformation.at(8) << std::endl;

    for( int j = 0 ; j < 600 ; j ++){
        for( int i=0; i< 20 ; i++){
            int ind = i + 20*(j) ;
            float theta = static_cast<float>(i) / 20* 2* M_PI;
            //mesh_tube_displayed.vertex(ind).x() +=  (deformation.at(j)*cos(theta)/1000);
            //mesh_tube_displayed.vertex(ind).y() += (deformation.at(j)*sin(theta)/1000);

           /* mesh_tube_displayed.vertex(ind).x() =
                    mesh_tube_displayed.vertex(ind).x()
                    - ( mesh_tube_displayed.vertex(ind).x() - mesh_tube.vertex(ind).x())
                    + (deformation.at(j)*cos(theta)/10000)*tps;
            mesh_tube_displayed.vertex(ind).y() =
                    mesh_tube_displayed.vertex(ind).y()
                    - ( mesh_tube_displayed.vertex(ind).y() - mesh_tube.vertex(ind).y())
                    + (deformation.at(j)*sin(theta)/10000)*tps;
                    */
            mesh_tube_displayed.vertex(ind).z() =
                                mesh_tube_displayed.vertex(ind).z()
                                - ( mesh_tube_displayed.vertex(ind).z() - mesh_tube.vertex(ind).z())
                                + (deformation.at(j)*cos(theta)/10000)*tps;

            //std::cout << cos(theta) << " " << sin(theta) << std::endl;
            //std::cout << ind << std::endl;
        }
    }
    tps+=timeStep;




    mesh_tube_displayed_opengl.fill_vbo(mesh_tube_displayed);
    //std::cout << T-tps << std::endl;
    if(T - tps < 10e-8){
        animationb = false;
        tps = 0.0f;
    }
}
