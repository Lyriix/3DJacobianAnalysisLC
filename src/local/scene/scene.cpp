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
int& scene::set_Nu_tube() {return Nu_tube;}
int& scene::set_Nv_tube() {return Nv_tube;}
int& scene::set_Nu_grid() {return Nu_grid;}
int& scene::set_Nv_grid() {return Nv_grid;}
int& scene::set_Nw_grid() {return Nw_grid;}

void scene::load_scene()
{

    //*****************************************//
    // Preload default structure               //
    //*****************************************//
    texture_default = load_texture_file("data/white.jpg");
    shader_program_id = read_shader("shaders/shader_mesh.vert",
                                    "shaders/shader_mesh.frag"); PRINT_OPENGL_ERROR();

    std::string csvYXName("/home/charly/workspace/OphtalmologyProject/DeformationFieldAnalysis/Modelisation(PerSlice)/data/YX.csv");
    std::string csvYYName("/home/charly/workspace/OphtalmologyProject/DeformationFieldAnalysis/Modelisation(PerSlice)/data/YY.csv");
    std::string csvYZName("/home/charly/workspace/OphtalmologyProject/DeformationFieldAnalysis/Modelisation(PerSlice)/data/YZ.csv");

    std::vector<std::vector<std::string>> csvYX = readCsvFile(csvYXName);
    std::vector<std::vector<std::string>> csvYY = readCsvFile(csvYYName);
    std::vector<std::vector<std::string>> csvYZ = readCsvFile(csvYZName);
    analyzeCsv(csvYX, deformationYX);
    analyzeCsv(csvYY, deformationYY);
    analyzeCsv(csvYZ, deformationYZ);

    //*****************************************//
    // Generate user defined mesh              //
    //*****************************************//

    float r=0.6f;
    float h=1.0f;
    //int Nu=20; int Nv= 600;
    Nu_tube = 20;
    Nv_tube = 600;
    Nu_grid = 40;
    Nv_grid = 40;
    Nw_grid = 40;
    scene::generate_tube(r,0.0f,0.0f,Nu_tube,Nv_tube);
    scene::generate_grid(Nu_grid,Nv_grid,Nw_grid);
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

    //mesh_tube_opengl.draw();
    if(draw_tube_state)
    {
        mesh_tube_displayed_opengl.draw();

    }
    if(draw_grid_state)
        grid_d_opengl.draw();
    animation();
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
/* //TEST
void scene::generate_grid(int Nu, int Nv, int Nw)
{
    grid_d.add_vertex(vec3(0.0f,0.0f,0.0f));
    grid_d.add_vertex(vec3(0.0f,1.0f,0.0f));
    grid_d.add_vertex(vec3(1.0f,1.0f,0.0f));

    grid_d.add_vertex(vec3(1.0f,0.0f,0.0f));

    grid_d.add_vertex(vec3(0.0f,0.0f,1.0f));

    grid_d.add_triangle_index({0,1,2});
    grid_d.add_triangle_index({0,2,3});
    grid_d.add_triangle_index({0,1,4});
    grid_d.fill_empty_field_by_default();
    grid_d_opengl.fill_vbo(grid_d);
}
*/

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
                //grid_d.add_color(getArcEnCielColor(k,Nw));
                //grid_d.add_color(vec3(0.7f,0.7f,0.7f));
                // grid_d.add_color(vec3(1.0f,0.0f,0.0f));
                // std::cout << " vec3 " <<u << " " << v << " " << w<< std::endl;
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
                //std::cout<<"ajout de "<< i+Nv*j+Nv*Nu*k <<" "<< i+1+Nv*j+Nv*Nu*k <<" "<< i+1+Nv*(1+j)+Nv*Nu*k<<std::endl;
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
    std::cout << "Connectivity of grid_d : "<<grid_d.size_connectivity() << std::endl;

    grid_d.fill_empty_field_by_default();
    grid_d_opengl.fill_vbo(grid_d);
    grid = grid_d;
    grid.fill_empty_field_by_default();
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
            //mesh_tube.add_color(getArcEnCielColor(j,Nv));
            mesh_tube.add_color(vec3(0.7f,0.7f,0.7f));
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

void scene::analyzeCsv(const std::vector<std::vector<std::string>>& csvfile, deformationArrays& deform)
{

    for(auto const & row : csvfile)
    {
        for(int i = 0; i < 4 ; i++)
        {

            if( i==0 )
                deform.iop1.push_back(std::stof(row[i]));
            if( i==1 )
                deform.iop2.push_back(std::stof(row[i]));
            if( i==2 )
                deform.recovery.push_back(std::stof(row[i]));
            if( i==3 )
                deform.baseline.push_back(0.0f);

            if(std::stof(row[i]) > dmax )
                dmax = std::stof(row[i]);
            if(std::stof(row[i]) < dmin )
                dmin = std::stof(row[i]);

            /*if( i == 0 )
                deformation.baseline.push_back(
                            (std::stof(row[0]) + std::stof(row[1])) / 2
                        - (std::stof(row[i]) + std::stof(row[1+i])) / 2 );
            if( i == 2 )
                deformation.iop1.push_back( (std::stof(row[0]) + std::stof(row[1])) / 2
                        - (std::stof(row[i]) + std::stof(row[1+i])) / 2);
            if( i == 4 )
                deformation.iop2.push_back( (std::stof(row[0]) + std::stof(row[1])) / 2
                        - (std::stof(row[i]) + std::stof(row[1+i])) / 2);
            if( i == 6 )
                deformation.recovery.push_back( (std::stof(row[0]) + std::stof(row[1])) / 2
                        - (std::stof(row[i]) + std::stof(row[1+i])) / 2);
                        */

            //std::cout << (std::stof(row[1]) + std::stof(row[2])) / 2
            //      - (std::stof(row[i]) + std::stof(row[1+i])) / 2 << std::endl;
            //std::cout << row.size() << std::endl;
            //std::cout << typeid(row[i]).name() << std::endl;
        }
    }
}

void scene::animation()
{
    /*if(baselineAnimation)
    {
        if(draw_tube_state)
            applyDeformation(deformation.baseline, baselineAnimation);
        if(draw_grid_state)
            applyGridDeformation(deformation.baseline, baselineAnimation);
    }*/
    // std::cout << iop1animation << std::endl;
    if(baselineAnimation)
    {
        if(draw_tube_state)
        {
            applyDeformation(deformationYX.baseline, 0);
            applyDeformation(deformationYY.baseline, 1);
            applyDeformation(deformationYZ.baseline, 2);
        }
        if(draw_grid_state)
        {
            applyGridDeformation(deformationYX.baseline,0);
            applyGridDeformation(deformationYY.baseline,1);
            applyGridDeformation(deformationYZ.baseline,2);
        }
        if(T-tps < 10e-6)
        {
            baselineAnimation = false;
            tps = 0.0f;
        }
    }
    if(iop1animation)
    {
        if(draw_tube_state)
        {
            applyDeformation(deformationYX.iop1, 0);
            applyDeformation(deformationYY.iop1, 1);
            applyDeformation(deformationYZ.iop1, 2);

        }
        if(draw_grid_state)
        {
            applyGridDeformation(deformationYX.iop1,0);
            applyGridDeformation(deformationYY.iop1,1);
            applyGridDeformation(deformationYZ.iop1,2);
            iop1animation = false;
        }
        if(T-tps < 10e-6)
        {
            iop1animation = false;
            tps = 0.0f;
        }
    }
    if(iop2animation)
    {
        if(draw_tube_state)
        {
            applyDeformation(deformationYX.iop2, 0);
            applyDeformation(deformationYY.iop2, 1);
            applyDeformation(deformationYZ.iop2, 2);
        }
        if(draw_grid_state)
        {
            applyGridDeformation(deformationYX.iop2,0);
            applyGridDeformation(deformationYY.iop2,1);
            applyGridDeformation(deformationYZ.iop2,2);
        }
        if(T-tps < 10e-6)
        {
            iop2animation = false;
            tps = 0.0f;
        }
    }
    if(recoveryAnimation)
    {
        if(draw_tube_state)
        {
            applyDeformation(deformationYX.recovery, 0);
            applyDeformation(deformationYY.recovery, 1);
            applyDeformation(deformationYZ.recovery, 2);
        }
        if(draw_grid_state)
        {
            applyGridDeformation(deformationYX.recovery,0);
            applyGridDeformation(deformationYY.recovery,1);
            applyGridDeformation(deformationYZ.recovery,2);
        }
        if(T-tps < 10e-6)
        {
            recoveryAnimation = false;
            tps = 0.0f;
        }
    }
    /*if(iop2animation)
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
    }*/
}

void scene::applyDeformation(std::vector<float> &deformation, int axis)
{

    float timeStep = 0.02f;
    //float T = 1.0f;
    //std::cout << "dmax " << dmax << std::endl;
    //std::cout << "dmin " << dmin << std::endl;
    float facteurCouleur = (dmax - dmin) / (dmax + dmin);
    //Parcourir les deformations
    //appliquer les deformations Y selon l'axe 0 : Nv

    //temps is represented by tps in scene::scene.hpp

    //std::cout << deformation.size() << std::endl;
    for( int j = 0 ; j < 600 ; j ++){
        for( int i=0; i< 20 ; i++){
            int ind = i + 20*(j) ;
            float theta = static_cast<float>(i) / 20* 2* M_PI;

            if(axis == 0)
                mesh_tube_displayed.vertex(ind).x() = mesh_tube_displayed.vertex(ind).x()
                        - (( mesh_tube_displayed.vertex(ind).x() - mesh_tube.vertex(ind).x())
                           + deformation.at(j)/100000)*tps;
            if( axis == 1 )
                mesh_tube_displayed.vertex(ind).y() = mesh_tube_displayed.vertex(ind).y()
                        - (( mesh_tube_displayed.vertex(ind).y() - mesh_tube.vertex(ind).y())
                           + deformation.at(j)/100000)*tps;
            if( axis == 2){
                mesh_tube_displayed.vertex(ind).z() = mesh_tube_displayed.vertex(ind).z()
                        - (( mesh_tube_displayed.vertex(ind).z() - mesh_tube.vertex(ind).z())
                           + deformation.at(j)/100000)*tps;

            float c = std::sqrt(std::pow(mesh_tube_displayed.vertex(ind).x(),2.0f)
                                + std::pow(mesh_tube_displayed.vertex(ind).y(),2.0f)
                                + std::pow(mesh_tube_displayed.vertex(ind).z(),2.0f));
            // std::cout << c << std::endl;
            c = (c-dmin)/(dmax -dmin);
            mesh_tube_displayed.color(ind)[0] = c;
             mesh_tube_displayed.color(ind)[2] = 1-c;
}
            /*float c = (std::abs(deformation.at(j)) - dmin) / (dmax - dmin);
            mesh_tube_displayed.color(ind) = vec3(c,0.0f, c);*/
            //std::cout<<deformation.at(j) <<std::endl;
            //mesh_tube_displayed.vertex(ind).x() +=  (deformation.at(j)*cos(theta)/1000);
            //mesh_tube_displayed.vertex(ind).y() += (deformation.at(j)*sin(theta)/1000);

            /*mesh_tube_displayed.vertex(ind).x() =
                    mesh_tube_displayed.vertex(ind).x()
                    - (( mesh_tube_displayed.vertex(ind).x() - mesh_tube.vertex(ind).x())
                    + deformation.at(j)*cos(theta)/100000)*tps;
            mesh_tube_displayed.vertex(ind).y() =
                    mesh_tube_displayed.vertex(ind).y()
                    - (( mesh_tube_displayed.vertex(ind).y() - mesh_tube.vertex(ind).y())
                    + deformation.at(j)*sin(theta)/100000)*tps;*/
            /*
            mesh_tube_displayed.vertex(ind).z() =
                                mesh_tube_displayed.vertex(ind).z()
                                -( ( mesh_tube_displayed.vertex(ind).z() - mesh_tube.vertex(ind).z())
                                + deformation.at(j)*cos(theta)/10000)*tps;*/
        }
    }
    tps+=timeStep;

    std::cout << mesh_tube_displayed.color(2000)[0] << " ";

    mesh_tube_displayed_opengl.fill_vbo(mesh_tube_displayed);
    /*if(T - tps < 10e-8){
        animationb = false;
        tps = 0.0f;
    }*/
}

void scene::applyGridDeformation(std::vector<float> &deformation, int axis)
{
    float timeStep = 1.0f;
    //float T = 1.0f;
    float facteur = 600 / Nv_grid; //600 is the number of row in the csv file (hard coded I know)

    for(int j=0 ; j < Nv_grid; j++) //chq slice sur y on change tout les voxels
    {
        float def=0.0f;
        for(int n = 0 ; n < floor(facteur) ; n++)
            def += deformation.at(n+j);
        std::cout << j << " " << def << std::endl;
        for(int i=0; i<Nu_grid; i++ )
        {
            for(int k=0; k<Nw_grid; k++)
            {
                int ind = i + j*Nu_grid + k*Nu_grid*Nv_grid;

                if(axis == 0)
                    grid_d.vertex(ind).x() = grid_d.vertex(ind).x()
                            - (( grid_d.vertex(ind).x() - grid.vertex(ind).x())
                               + def)*tps;

                if(axis == 1)
                    grid_d.vertex(ind).y() = grid_d.vertex(ind).y()
                            - (( grid_d.vertex(ind).y() - grid.vertex(ind).y())
                               + def)*tps;
                if(axis == 2)
                    grid_d.vertex(ind).z() = grid_d.vertex(ind).z()
                            - (( grid_d.vertex(ind).z() - grid.vertex(ind).z())
                               + def)*tps;
            }
        }
    }
    tps+=timeStep;
    grid_d_opengl.fill_vbo(grid_d);

}
