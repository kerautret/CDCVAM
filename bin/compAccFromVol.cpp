#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>


#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/readers/PointListReader.h>

#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/io/writers/VolWriter.h>

#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/writers/MeshWriter.h>
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageContainerBySTLMap.h>
#include <DGtal/kernel/sets/DigitalSetFromMap.h>


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

// #include <pcl/features/integral_image_normal.h>
// #include <pcl/features/normal_3d.h>
// #include <pcl/point_types.h>


// ToDO : mettre à jour avec la nouvelle version et la classe NormalAccumulato

// #include "TubeAnalyseHelper.h"


using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

typedef ImageContainerBySTLVector<Z3i::Domain, unsigned int> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef ImageContainerBySTLVector<Z3i::Domain,  Z3i::RealPoint> ImageVector;


typedef ImageContainerBySTLMap<Z3i::Domain,double> DistanceImage;
typedef DigitalSetFromMap<DistanceImage> AcceptedPointSet;
typedef Z3i::Domain::Predicate DomainPredicate;






/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("input,i", po::value<std::string>(), "input vol file.")
  ("output,o", po::value<std::string>()->default_value("accumulation.vol"), "output accumulation vol file.")
  ("outputR,R", po::value<std::string>()->default_value("radius.vol"), "output radius vol file.")
  ("reScaleOutRadiusImg,s", "auto scale out image values between 0 255.")
  ("exportNormal,e", po::value<std::string>()->default_value("normal.sdp"), "export normals estimated from gradients.")
  ("radius,r", po::value<double>()->default_value(5), "maximal radius of accumulation analysis.")
  ("expNormalMinG,m", po::value<double>()->default_value(10.0), "min value of gradient to export the normals.")
  ("minGrad,g", po::value<double>()->default_value(10.0), "min radient to consider accumulation.");;


  
  bool parseOK = true;
  po::variables_map vm;
  
  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  } catch (const std::exception &ex) {
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
    parseOK = false;
  }
  
  po::notify(vm);
  
  if (vm.count("help") || argc <= 1 || !parseOK || !vm.count("input")) {
    trace.info() << "Compute accumulation from vol and compute the radius (median value) "
    << "of all normal participating to the accumulation" << std::endl
    << "Options: " << std::endl
    << general_opt << "\n";
    return 0;
  }
  
  double radius = vm["radius"].as<double>();
  string inputFile = vm["input"].as<std::string>();
  string outputFile = vm["output"].as<std::string>();
  string outputFileRadius = vm["outputR"].as<std::string>();
  string normalFileName = vm["exportNormal"].as<std::string>();
  double expNormalMinG = vm["expNormalMinG"].as<double>();
  double minGrad = vm["minGrad"].as<double>();

  // 1) Reading input image
  Image3DChar inputImage = GenericReader<Image3DChar>::import(vm["input"].as<std::string>());
  ImageVector imageVF(inputImage.domain());

  // 2) compute gradient
  AccuGradientFieldHelper::compGradientVects(inputImage, imageVF);
  

  // 3) compute accumulation
  Image3D imageAcc(inputImage.domain());
  ImageVector mainDir (inputImage.domain());


  AccuGradientFieldHelper::computeNormalIntersectImage(imageAcc, mainDir, radius, imageVF, minGrad, "", false);
  
   unsigned int maxVal = 0;
  // 4)  export normals 
  std::ofstream outNormals;
  outNormals.open(normalFileName);
  double scaleV = radius;
   for(typename Image3DChar::Domain::ConstIterator it = imageAcc.domain().begin();
      it != imageAcc.domain().end(); it++){
    if(imageAcc(*it)>maxVal){
      maxVal = imageAcc(*it);
    }
    Z3i::RealPoint pt = *it;
    Z3i::RealPoint pt2 = pt+(imageVF(*it).getNormalized());

    if(imageAcc(*it)> expNormalMinG){
      outNormals << pt[0] <<  " "  << pt[1] << " " << pt[2] << std::endl;
      outNormals << pt2[0] <<  " "  << pt2[1] << " " << pt2[2] << std::endl;
    }
   }
   outNormals.close();
  

  
  
   
   trace.info()<< "export accumulation image .... on domain: " << imageAcc.domain();
   
   
    trace.info() << "Saving in " << outputFile;
    typedef functors::Rescaling<unsigned int, unsigned char> ScaleFct;
    ScaleFct  scaleFct (0 ,maxVal, 0, 255);
    VolWriter<Image3D, ScaleFct >::exportVol(outputFile, imageAcc, scaleFct);
    trace.info() << "[done]" << std::endl;
  
    trace.info() << "Saving in " << outputFileRadius;
    // if(vm.count("reScaleOutRadiusImg")){
    //   typedef functors::Rescaling<unsigned char, unsigned char> ScaleFct;
    //   ScaleFct  scaleFctR (0 ,radius, 0, 255);
    //   VolWriter<Image3DChar, ScaleFct >::exportVol(outputFileRadius, imageRadius, scaleFctR);
    //   trace.info() << "[done]" << std::endl;
    // }else{
    //   imageRadius >> outputFileRadius;
    // }
  
  
  
  
  return 0;
}

