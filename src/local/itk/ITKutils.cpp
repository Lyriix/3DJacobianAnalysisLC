
#include "ITKutils.hpp"

#include "itkImageSliceConstIteratorWithIndex.h"


ITKutils::ITKutils()
    :deformationImages()
{}

std::vector<int> ITKutils::get_deformation_size() {return deformation_size;}
std::vector<itk::Vector<float,3>> ITKutils::get_deformations_baseline(void){ return deformationMeans.baseline;}
std::vector<itk::Vector<float,3>> ITKutils::get_deformations_iop1(void){ return deformationMeans.iop1;}
std::vector<itk::Vector<float,3>> ITKutils::get_deformations_iop2(void){ return deformationMeans.iop2;}
std::vector<itk::Vector<float,3>> ITKutils::get_deformations_recovery(void){ return deformationMeans.recovery;}


//void ITKutils::read_file(std::string filename, ITKutils::ImageT::Pointer& image)
void ITKutils::read_file(std::string filename, std::string nameImage)
{
    std::cout << "reading " << filename <<std::endl;
    ITKutils::ReaderT::Pointer reader = ITKutils::ReaderT::New();
    reader->SetFileName(filename);
    try{
        reader->Update();
    }
    catch(itk::ExceptionObject &ex){
        std::cerr << "Error while reading file : " << filename << std::endl;
        std::cerr << ex << std::endl;
    }

    if(std::strcmp(nameImage.c_str(),"baselinea")==0)
    {
        deformationImages.baselinea = reader->GetOutput();
        ImageT::RegionType region = reader->GetOutput()->GetLargestPossibleRegion(); //We just need to get the size once ; did it here for an efficiency problem
        ImageT::SizeType size = region.GetSize();
        std::cout << "size : " << size << std::endl;
        deformation_size.push_back(size[0]);
        deformation_size.push_back(size[1]);
        deformation_size.push_back(size[2]);
        ImageT::IndexType idx;
        idx[0] = 174; idx[1] = 599; idx[2] = 249;
        std::cout << deformationImages.baselinea->GetPixel(idx) << std::endl;
    }
    if(std::strcmp(nameImage.c_str(),"baselineb")==0)
        deformationImages.baselineb = reader->GetOutput();
    if(std::strcmp(nameImage.c_str(),"iop1a")==0)
        deformationImages.iop1a = reader->GetOutput();
    if(std::strcmp(nameImage.c_str(),"iop1b")==0)
        deformationImages.iop1b = reader->GetOutput();
    if(std::strcmp(nameImage.c_str(),"iop2a")==0)
        deformationImages.iop2a = reader->GetOutput();
    if(std::strcmp(nameImage.c_str(),"iop2b")==0)
        deformationImages.iop2b = reader->GetOutput();
    if(std::strcmp(nameImage.c_str(),"recoverya")==0)
        deformationImages.recoverya = reader->GetOutput();
    if(std::strcmp(nameImage.c_str(),"recoveryb")==0)
        deformationImages.recoveryb = reader->GetOutput();





}


void ITKutils::analyze_files()
{
    //On parcours chaque pixel de l'imae et on store ces pixels dans des vecteurs
    // de sorte que ind = i + Nu + NvNuk
    std::cout << "computing means" << std::endl;
    typedef itk::ImageSliceConstIteratorWithIndex< ImageT
    > SliceIteratorType;

    ImageT::IndexType Index = {{10, 4, 25}};
    std::cout<< deformationImages.baselinea->GetPixel(Index) << std::endl;
    std::cout<< deformationImages.baselineb->GetPixel(Index) << std::endl;
    int dir0 = 0;
    int dir1 = 1;
    SliceIteratorType baseait(deformationImages.baselinea, deformationImages.baselinea->GetRequestedRegion());
    baseait.SetFirstDirection(dir0); // Slice throught the Y axis
    baseait.SetSecondDirection(dir1); // Linear against X axis
    SliceIteratorType basebit(deformationImages.baselineb, deformationImages.baselineb->GetRequestedRegion());
    basebit.SetFirstDirection(dir0); // Slice throught the Y axis
    basebit.SetSecondDirection(dir1); // Linear against X axis
    SliceIteratorType iop1ait(deformationImages.iop1a, deformationImages.iop1a->GetRequestedRegion());
    iop1ait.SetFirstDirection(dir0); // Slice throught the Y axis
    iop1ait.SetSecondDirection(dir1); // Linear against X axis
    SliceIteratorType iop1bit(deformationImages.iop1b, deformationImages.iop1b->GetRequestedRegion());
    iop1bit.SetFirstDirection(dir0); // Slice throught the Y axis
    iop1bit.SetSecondDirection(dir1); // Linear against X axis
    SliceIteratorType iop2ait(deformationImages.iop2a, deformationImages.iop2a->GetRequestedRegion());
    iop2ait.SetFirstDirection(dir0); // Slice throught the Y axis
    iop2ait.SetSecondDirection(dir1); // Linear against X axis
    SliceIteratorType iop2bit(deformationImages.iop2b, deformationImages.iop2b->GetRequestedRegion());
    iop2bit.SetFirstDirection(dir0); // Slice throught the Y axis
    iop2bit.SetSecondDirection(dir1); // Linear against X axis
    SliceIteratorType recoveryait(deformationImages.recoverya, deformationImages.recoverya->GetRequestedRegion());
    recoveryait.SetFirstDirection(dir0); // Slice throught the Y axis
    recoveryait.SetSecondDirection(dir1); // Linear against X axis
    SliceIteratorType recoverybit(deformationImages.recoveryb, deformationImages.recoveryb->GetRequestedRegion());
    recoverybit.SetFirstDirection(dir0); // Slice throught the Y axis
    recoverybit.SetSecondDirection(dir1); // Linear against X axis


    baseait.GoToBegin();
    basebit.GoToBegin();
    iop1ait.GoToBegin();
    iop1bit.GoToBegin();
    iop2ait.GoToBegin();
    iop2bit.GoToBegin();
    recoveryait.GoToBegin();
    recoverybit.GoToBegin();

    while(!baseait.IsAtEnd())
    {
        while( !baseait.IsAtEndOfSlice())
        {
            while( !baseait.IsAtEndOfLine())
            {
                deformationMeans.baseline.push_back( (baseait.Get() + basebit.Get()) / 2 );
                deformationMeans.iop1.push_back( (iop1ait.Get() + iop1bit.Get()) / 2 );
                deformationMeans.iop2.push_back( (iop2ait.Get() + iop2bit.Get()) / 2 );
                deformationMeans.recovery.push_back( (recoveryait.Get() + recoverybit.Get()) / 2 );
                ++baseait;++basebit;
                ++iop1ait; ++iop1bit;
                ++iop2ait; ++iop2bit;
                ++recoveryait; ++recoverybit;
            }
            baseait.NextLine(); basebit.NextLine();iop1ait.NextLine();iop1bit.NextLine();
            iop2ait.NextLine();iop2bit.NextLine(); recoveryait.NextLine();recoverybit.NextLine();
        }
        baseait.NextSlice(); basebit.NextSlice();iop1ait.NextSlice();iop1bit.NextSlice();
        iop2ait.NextSlice();iop2bit.NextSlice(); recoveryait.NextSlice();recoverybit.NextSlice();
    }
}
