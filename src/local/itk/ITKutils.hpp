/** Modelisation of LC deformation - NYU Tandon School - 2017 */

#pragma once

#ifndef ITKUTILS_HPP
#define ITKUTILS_HPP

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkVector.h>
#include <string>

class ITKutils
{
    struct deformationArrays {
        std::vector<itk::Vector<float,3>> baseline;  //mean of baseline a and baseline b
        std::vector<itk::Vector<float,3>> iop1;      //mean of iop1 and iop2
        std::vector<itk::Vector<float,3>> iop2;      //...
        std::vector<itk::Vector<float,3>> recovery;  //...
    };

public:

    /** \brief constructor */
    ITKutils();
    typedef float PixelT;
    typedef itk::Image<itk::Vector<PixelT,3>,3> ImageT;
    typedef itk::ImageFileReader<ImageT> ReaderT;
    /** \brief Method called to read a file with itk and store it in an image */
    void read_file(std::string filename, std::__cxx11::string nameImage);

    /** \brief method to get the size of the images x y z */
    std::vector<int> get_deformation_size();

    /** \brief Method to get the deformation field */
    std::vector<itk::Vector<float,3>> get_deformations_baseline(void);
    std::vector<itk::Vector<float,3>> get_deformations_iop1(void);
    std::vector<itk::Vector<float,3>> get_deformations_iop2(void);
    std::vector<itk::Vector<float,3>> get_deformations_recovery(void);

    /** \brief Method called to analyze the readed fileS for a specific use here */
    void analyze_files();




private:
    /** itk type definitions */
    //int dim;
    std::vector<int> deformation_size;


    /** \brief structure with each readed deformation field */
    struct deformationFieldImages{
        ImageT::Pointer baselinea;
        ImageT::Pointer baselineb;
        ImageT::Pointer iop1a;
        ImageT::Pointer iop1b;
        ImageT::Pointer iop2a;
        ImageT::Pointer iop2b;
        ImageT::Pointer recoverya;
        ImageT::Pointer recoveryb;
    };

    deformationFieldImages deformationImages;



    deformationArrays deformationMeans;
};


#endif
