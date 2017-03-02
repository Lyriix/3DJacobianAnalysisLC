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


#include <cmath>
#include <string>
#include <sstream>


using namespace cpe;
int& scene::set_Nu_tube() {return Nu_tube;}
int& scene::set_Nv_tube() {return Nv_tube;}
int& scene::set_Nu_grid() {return Nu_grid;}
int& scene::set_Nv_grid() {return Nv_grid;}
int& scene::set_Nw_grid() {return Nw_grid;}
scene::deformationArrays& scene::set_deformationArrays(){return deformation;}

void scene::load_scene()
{

    //*****************************************//
    // Preload default structure               //
    //*****************************************//
    texture_default = load_texture_file("data/white.jpg");
    shader_program_id = read_shader("shaders/shader_mesh.vert",
                                    "shaders/shader_mesh.frag"); PRINT_OPENGL_ERROR();
    //*****************************************//
    // Generate user defined mesh              //
    //*****************************************//
    float r=0.6f;
    float h=1.0f;

    Nu_tube = 20;
    Nv_tube = 600;
    Nu_grid = 100;
    Nv_grid = 600;
    Nw_grid = 100;

    scene::generate_tube(r,0.0f,0.0f,Nu_tube,Nv_tube);
    scene::generate_grid(Nu_grid,Nv_grid,Nw_grid);
    std::string loc = "/home/charly/workspace/OphtalmologyProject/DeformationFieldAnalysis/Modelisation(PerVoxelwithItk)/data/InverseWarp/";
    std::cout << "ccccc" << std::endl;
    itkObject.read_file(loc + "den500101InverseWarp.nii.gz", "baselinea");
    itkObject.read_file(loc + "den500211InverseWarp.nii.gz", "baselineb");
    itkObject.read_file(loc + "den530121InverseWarp.nii.gz", "iop1a");
    itkObject.read_file(loc + "den530231InverseWarp.nii.gz", "iop1b");
    itkObject.read_file(loc + "den550141InverseWarp.nii.gz", "iop2a");
    itkObject.read_file(loc + "den550251InverseWarp.nii.gz", "iop2b");
    itkObject.read_file(loc + "den590161InverseWarp.nii.gz", "recoverya");
    itkObject.read_file(loc + "den590271InverseWarp.nii.gz", "recoveryb");

    itkObject.analyze_files();

    deformation.baseline = itkObject.get_deformations_baseline();
    deformation.iop1 = itkObject.get_deformations_iop1();
    deformation.iop2 = itkObject.get_deformations_iop2();
    deformation.recovery = itkObject.get_deformations_recovery();


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
    if(draw_tube_state)
        mesh_tube_displayed_opengl.draw();
    if(draw_grid_state)
        grid_d_opengl.draw();

    animation(deformation);
}

void scene::change_draw_tube_state()
{
    draw_tube_state=!draw_tube_state;
}
void scene::change_draw_grid_state()
{
    draw_grid_state=!draw_grid_state;
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

void scene::generate_grid(int Nu, int Nv, int Nw)
{
    float u=0, v=0, w=0;
    for(int k=0; k<Nw; k++)
    {
        w = (float)k/Nw;
        for( int j=0; j<Nv; j++)
        {
            v = (float)j/Nv;
            for( int i=0; i<Nu; i++)
            {
                u = (float)i/Nu;
                grid_d.add_vertex(vec3(u,v,w));
                grid_d.add_color(getArcEnCielColor(k,Nw));
            }
        }
    }
    std::cout << "Size mesh grid : "<<grid_d.size_vertex() << std::endl;
    for(int k=0; k<Nw-1; k++)
    {
        for( int j=0; j<Nv-1; j++)
        {
            for( int i=0; i<Nu-1; i++)
            {
                //faces
                grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+1+Nu*j+Nv*Nu*k, i+1+Nu*(1+j)+Nv*Nu*k});
                grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+Nu*(1+j)+Nv*Nu*k, i+1+Nu*(1+j)+Nv*Nu*k});
                //bottom
                grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+1+Nu*j+Nv*Nu*k, i+1+Nu*(j)+Nv*Nu*(1+k)});
                grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+Nu*(j)+Nv*Nu*(1+k), i+1+Nu*(j)+Nv*Nu*(1+k)});
                //Left Sides
                grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+Nu*(1+j)+Nv*Nu*k, i+Nu*(1+j)+Nv*Nu*(1+k)});
                grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+Nu*(j)+Nv*Nu*(1+k), i+Nu*(1+j)+Nv*Nu*(1+k)});
            }
        }
    }
    //Last left faces OR right faces // x = Nu
    for( int j=0; j<Nv-1; j++)
    {
        for( int k=0; k<Nw-1; k++)
        {
            int i = Nu -1;
            grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+Nu*(1+j)+Nv*Nu*k, i+Nu*(1+j)+Nv*Nu*(1+k)});
            grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+Nu*(j)+Nv*Nu*(1+k), i+Nu*(1+j)+Nv*Nu*(1+k)});
        }
    }
    //Last bottom OR upper // y = Nv
    for( int i=0; i<Nu-1; i++)
    {
        for( int k=0; k<Nw-1; k++)
        {
            int j = Nv -1;
            grid_d.add_triangle_index({i+Nu*j+Nv*Nv*k, i+1+Nu*j+Nv*Nu*k, i+1+Nu*(j)+Nv*Nu*(1+k)});
            grid_d.add_triangle_index({i+Nu*j+Nv*Nv*k, i+Nu*(j)+Nv*Nu*(1+k), i+1+Nu*(j)+Nv*Nu*(1+k)});
        }
    }
    //last faces OR back // z = Nw
    for( int i=0; i<Nu-1; i++)
    {
        for( int j=0; j<Nv-1; j++)
        {
            int k = Nw -1;
            grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+1+Nu*j+Nv*Nu*k, i+1+Nu*(1+j)+Nv*Nu*k});
            grid_d.add_triangle_index({i+Nu*j+Nv*Nu*k, i+Nu*(j+1)+Nv*Nu*k, i+1+Nu*(1+j)+Nv*Nu*k});
        }
    }

    std::cout << "grid connectivity " << grid_d.size_connectivity() << std::endl;
    grid_d.fill_empty_field_by_default();
    grid = grid_d;
    grid_opengl.fill_vbo(grid);
    grid_d_opengl.fill_vbo(grid_d);
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
            mesh_tube.add_color(getArcEnCielColor(j,Nv));
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

vec3 scene::getArcEnCielColor(int j, int inter) //intervalle = Nv
{
    //vec3 color;
    int x = (1530/inter)*j; //1530 = 6* 255 = max de colors pour cet algo
    float r, g, b;
    if ( x >= 0 && x<255 )
    {
        r=255;
        g=x;
        b=0;
    }
    if ( x >= 255 && x < 510)
    {
        r = 510 - x;
        g = 255;
        b = 0;
    }
    if( x >= 510 && x < 765)
    {
        r = 0;
        g = 255;
        b = x - 510;
    }
    if( x >= 765 && x < 1020)
    {
        r = 0;
        g = 1020 - x;
        b = 255;
    }
    if( x >= 1020 && x < 1275)
    {
        r = x - 1020;
        g = 0;
        b = 255;
    }
    if( x >= 1275 && x <= 1530)
    {
        r = 255;
        g = 0;
        b = 1530 - x ;
    }
    return vec3(r/255,g/255,b/255);
}


void scene::animation(deformationArrays &deformation)
{
    if(baselineAnimation)
    {
        if(draw_tube_state)
            applyDeformation(deformation.baseline, baselineAnimation);
        if(draw_grid_state)
            applyGridDeformation(deformation.baseline, baselineAnimation);
    }

    if(iop1animation)
    {
        if(draw_tube_state)
            applyDeformation(deformation.iop1, iop1animation);
        if(draw_grid_state)
            applyGridDeformation(deformation.iop1, iop1animation);
    }
    if(iop2animation)
    {
        if(draw_tube_state)
            applyDeformation(deformation.iop2, iop2animation);
        if(draw_grid_state)
            applyGridDeformation(deformation.iop2, iop2animation);
    }
    if(recoveryAnimation)
    {
        if(draw_tube_state)
            applyDeformation(deformation.recovery, recoveryAnimation);
        if(draw_grid_state)
            applyGridDeformation(deformation.recovery, recoveryAnimation);
    }
}

void scene::applyDeformation( std::vector<itk::Vector<float,3>> &deformation, bool &animationb)
{

    float timeStep = 0.02f;
    float T = 1.0f;

    //Parcourir les deformations
    //appliquer les deformations Y selon l'axe 0 : Nv

    //temps is represented by tps in scene::scene.hpp
    for( int j = 0 ; j < Nv_tube ; j ++){
        for( int i=0; i< Nu_tube ; i++){
            int ind = i + Nu_tube*(j) ;
            float theta = static_cast<float>(i) / Nu_tube* 2* M_PI;
            //mesh_tube_displayed.vertex(ind).x() +=  (deformation.at(j)*cos(theta)/1000);
            //mesh_tube_displayed.vertex(ind).y() += (deformation.at(j)*sin(theta)/1000);
            /*
            mesh_tube_displayed.vertex(ind).x() =
                    mesh_tube_displayed.vertex(ind).x()
                    - (( mesh_tube_displayed.vertex(ind).x() - mesh_tube.vertex(ind).x())
                    + deformation.at(j)*cos(theta)/100000)*tps;
            mesh_tube_displayed.vertex(ind).y() =
                    mesh_tube_displayed.vertex(ind).y()
                    - (( mesh_tube_displayed.vertex(ind).y() - mesh_tube.vertex(ind).y())
                    + deformation.at(j)*sin(theta)/100000)*tps;
                    */
            /*
            mesh_tube_displayed.vertex(ind).z() =
                                mesh_tube_displayed.vertex(ind).z()
                                -( ( mesh_tube_displayed.vertex(ind).z() - mesh_tube.vertex(ind).z())
                                + deformation.at(j)*cos(theta)/10000)*tps;*/
        }
    }
    tps+=timeStep;

    mesh_tube_displayed_opengl.fill_vbo(mesh_tube_displayed);
    if(T - tps < 10e-8){
        animationb = false;
        tps = 0.0f;
    }
}

void scene::applyGridDeformation(std::vector<itk::Vector<float,3>> &deformation, bool &animationb)
{
    /** deformation is a 3-D vector with X,Y,Z component of the deformation field at the voxel (i,j,k) */

    float timeStep = 0.2f;
    float T = 1.0f;


    std::vector<int> size_of_deformation_field = itkObject.get_deformation_size();

    float facteurx = size_of_deformation_field[0] / Nu_grid;
    float facteury = size_of_deformation_field[1] / Nv_grid;
    float facteurz = size_of_deformation_field[2] / Nw_grid;

    for(int i=0 ; i < Nu_grid; i++) //chq slice sur y on change tout les voxels
    {
        float defy = 0.0f;
        for(int j=0; j<Nv_grid; j++ )
        {
            for(int k=0; k<Nw_grid; k++)
            {

                for(int x=0 ; x < facteurx; x++)
                {
                    for(int y=0; y<facteury; y++)
                    {
                        for(int z=0; z<facteurz; z++)
                        {
                            //We need to ressample the deformations to fit the grid
                            int ind = (i+x) + (j+y)*size_of_deformation_field[0] + (k+z)*size_of_deformation_field[0]*size_of_deformation_field[1];
                            defy += deformation.at(ind)[1];
                        }
                    }
                   // std::cout << defy << std::endl;
                    int ind_grid = i + j*Nu_grid + k*Nu_grid*Nv_grid;
                    grid_d.vertex(ind_grid).y() +=
                            -( (grid_d.vertex(ind_grid).y()-grid.vertex(ind_grid).y()) +defy )*tps/10000000;
                }
                //std::cout<< i <<" " << j <<" "<< k <<std::endl;
                /* if(i==10 && j==4 && k==25 ){
                std::cout << deformation.at(i + j*size_of_deformation_field[0] + k*size_of_deformation_field[0]*size_of_deformation_field[1])[0] << std::endl;
                std::cout << deformation.at(i + j*size_of_deformation_field[0] + k*size_of_deformation_field[0]*size_of_deformation_field[1])[1] << std::endl;
                std::cout << deformation.at(i + j*size_of_deformation_field[0] + k*size_of_deformation_field[0]*size_of_deformation_field[1])[2] << std::endl;
                std::cout << deformation.size() << std::endl ;
                */


            }
        }
    }

    tps+=timeStep;
    grid_d_opengl.fill_vbo(grid_d);
    if(T - tps < 10e-2){
        animationb = false;
        tps = 0.0f;
    }

}
