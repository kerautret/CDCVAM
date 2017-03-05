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
#include <DGtal/io/writers/LongvolWriter.h>


#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/writers/MeshWriter.h>
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageContainerBySTLMap.h>
#include <DGtal/kernel/sets/DigitalSetFromMap.h>


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <pcl/features/integral_image_normal.h>
#include <pcl/features/normal_3d.h>
#include <pcl/point_types.h>



#include "AccumulatorHelper.h"
#include "NormalAccumulator.h"


using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain,  double> ImageDouble;







/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
  ("help,h", "display this message")
  ("input,i", po::value<std::string>(), "input sdp points (x y z format).")
  ("outputAcc,o", po::value<std::string>()->default_value("accumulation.vol"), "output accumulation vol file.")
  ("outputConf,c", po::value<std::string>()->default_value("confidence.longvol"), "output confidence longvol file.")
  ("outputRad,R", po::value<std::string>()->default_value("radius.vol"), "output radius vol file.")  
  ("radiusEstimator,e", po::value<std::string>()->default_value("min"),  "use: {min (default), max, mean, median} to estimate the radius") 
      
  ("neighborhoodDistance,d",po::value<double>(), "use a neighborhood distance to estimate the normals." )
  ("neighborhoodNum,n",po::value<unsigned int>(), "set a number of neighbors to define the neighborhood to estimate the normals." )

  ("exportNormalSegments", po::value<std::string>()->default_value("normalPCL.sdp"), "export normal segments estimated from PCL (extremity of normals vectors).")
  ("exportNormals", po::value<std::string>()->default_value("normalSegPCL.normals"), "export normals estimated from PCL (unit vectors).")

    ("exportPointAssociations", po::value<unsigned int>(), "export in sdp (voxelAssos.sdp) the associated points contributing to the accumulation score (for voxels with score < threshold). Each pair point of the association are exported sequentially.")
  ("radius,r", po::value<double>()->default_value(5), "radius of accumulation analysis.")
  ("maxValOutConf", po::value<DGtal::uint64_t>()->default_value(255), "set MAX scale of confidence out image: 0 .. 1 -> 0 ..  MAX.")
  ("maxValOutRad", po::value<DGtal::uint64_t>()->default_value(255), "set MAX scale of radius out image: 0 .. 1 -> 0 ..  MAX.")
  ("autoScaleRad", po::value<std::string>(), "auto scale out radius image values between 0 and 255.")
  ("autoScaleConf", po::value<std::string>(), "auto scale out confidence image values between 0 and 255.")
  ("autoScaleAcc", po::value<std::string>(), "auto scale out accumulation image values between 0 and 255.");
  
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
    trace.info() << "Compute accumulation from point of cloud and compute the radius (median value) "
    << "of all normal participating to the accumulation" << std::endl
    << "Options: " << std::endl
    << general_opt << "\n";
    return 0;
  }
  
  double radius = vm["radius"].as<double>();
  string inputFile = vm["input"].as<std::string>();
  string outputFileAcc = vm["outputAcc"].as<std::string>();
  string outputFileRad = vm["outputRad"].as<std::string>();
  string outputFileConf = vm["outputConf"].as<std::string>();

  string typeStat = vm["radiusEstimator"].as<string>();


  DGtal::uint64_t outConfidenceMax = vm["maxValOutConf"].as<DGtal::uint64_t>();
  DGtal::uint64_t outRadiusMax = vm["maxValOutRad"].as<DGtal::uint64_t>();
  

  // 1) Reading input point of cloud and constructoing the cloud

  DGtal::trace.info() << "Reading imput ... " << std::endl;
  std::vector<Z3i::RealPoint> setOfPt = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFile);
  DGtal::trace.info() << "[done]" << std::endl;

  
  // 2) create and init an NormalAccumulator

  NormalAccumulator normAcc(radius, typeStat);
  if(vm.count("neighborhoodDistance")){
    double d = vm["neighborhoodDistance"].as<double>();
    normAcc.initFromPointCloud(setOfPt, d);
  }else{
    //ne.setKSearch(vm["neighborhoodNum"].as<unsigned int>());
    trace.error() << "you need to use either --neighborhoodDistance or --neighborhoodNum options" << std::endl; 
  }
  
  
  
  // 3) check if export
  if(vm.count("exportNormalSegments")){
    std::ofstream outNormals;
    string normalSegFileName = vm["exportNormalSegments"].as<std::string>();
    outNormals.open(normalSegFileName);
    double scaleV = 2.0;
    std::vector<Z3i::RealPoint> vPts = normAcc.getNormalOrigins();
    std::vector<Z3i::RealPoint> vNormals = normAcc.getNormalField();
    for (unsigned int i = 0; i < vPts.size() ; i++) {
      outNormals << vPts[i][0] << " " << vPts[i][1] << " "<< vPts[i][2] << std::endl;
      outNormals <<vPts[i][0]+vNormals[i][0]*scaleV << " "<<vPts[i][1]+vNormals[i][1]*scaleV  
                 << " "<< vPts[i][2]+vNormals[i][2]*scaleV << std::endl;
    }  
    outNormals.close();
  }

  if(vm.count("exportNormals")){
    std::ofstream outNormals;
    string normalSegFileName = vm["exportNormals"].as<std::string>();
    outNormals.open(normalSegFileName);
    std::vector<Z3i::RealPoint> vNormals = normAcc.getNormalField();
    for (unsigned int i = 0; i < vNormals.size() ; i++) {
      Z3i::RealPoint n = vNormals[i].getNormalized();
      outNormals << n[0] << " " << n[1] << " "<< n[2] << std::endl;
    }  
    outNormals.close();
  }

  

  
  // 3) Generate accumulation image
  normAcc.computeAccumulation();
  NormalAccumulator::Image3D &imageAccumulation = normAcc.getAccumulationImage();

  trace.info() << "Saving accumulation image in " << outputFileAcc;
  LongvolWriter<Image3D>::exportLongvol(outputFileAcc, imageAccumulation);
  trace.info() << "[done]" << std::endl;

  if (vm.count("autoScaleAcc")){
    trace.info() << "Saving accumulation (auto scale 0 " << normAcc.getMaxAccumulation() 
                 <<"  -> [0, 255]) image in " << vm["autoScaleAcc"].as<std::string>() << " ... ";
    typedef functors::Rescaling<DGtal::uint64_t, unsigned char> ScaleFct;
    ScaleFct  scaleFct (0.0 ,normAcc.getMaxAccumulation(), 0, 255);
    VolWriter<Image3D,ScaleFct>::exportVol(vm["autoScaleAcc"].as<std::string>(), imageAccumulation, true,  scaleFct);
    trace.info() << "[done]" << std::endl;
  }
  
  
  // 4) Compute confidence image 
  normAcc.computeConfidence();
  ImageDouble imageConfidance = normAcc.getConfidenceImage();
  trace.info() << "Saving confidence image in " << outputFileConf << " ... ";
  typedef functors::Rescaling<double, DGtal::uint64_t> ScaleFctD;
  ScaleFctD scaleFct(0.0, 1.0, 0, outConfidenceMax);
  LongvolWriter<ImageDouble, ScaleFctD>::exportLongvol(outputFileConf, imageConfidance, true, scaleFct);
  trace.info() << "[done]" << std::endl;
  
  if(vm.count("autoScaleConfidence"))
    {
      std::string outNameAutoConf = vm["autoScaleConf"].as<std::string>();
      trace.info() << "Saving confidence (auto scale 0 1 -> [0, 255]) image in " << outNameAutoConf << " ... ";
      typedef functors::Rescaling<double, unsigned char> ScaleFctD;
      ScaleFctD  scaleFct (0.0 ,1.0, 0, 255);
      VolWriter<ImageDouble,ScaleFctD>::exportVol(outNameAutoConf, imageConfidance, true, scaleFct);
      trace.info() << "[done]" << std::endl;    
    }
  

  // 5) Compute Radius image
  normAcc.computeRadiusFromConfidence();
  ImageDouble imageRadius = normAcc.getRadiusImage();  
  trace.info() << "Saving radius image in " << outputFileRad << " ... ";
  ScaleFctD scaleFct2(0.0, normAcc.getMaxRadius(), 0, outRadiusMax);
  LongvolWriter<ImageDouble, ScaleFctD>::exportLongvol(outputFileRad, imageRadius, true, scaleFct2);
  trace.info() << "[done]" << std::endl;
 if(vm.count("autoScaleRadius"))
    {
      std::string outNameAutoRad = vm["autoScaleRad"].as<std::string>();
      trace.info() << "Saving confidence (auto scale 0 " << normAcc.getMaxRadius()
                   << " -> 0 255) image in " << outNameAutoRad << " ... ";
      typedef functors::Rescaling<double, unsigned char> ScaleFctD;
      ScaleFctD  scaleFct (0.0 ,normAcc.getMaxRadius(), 0, 255);
      VolWriter<ImageDouble,ScaleFctD>::exportVol(outNameAutoRad, imageRadius, true, scaleFct);
      trace.info() << "[done]" << std::endl;    
    }      

 if(vm.count("exportPointAssociations"))
   {
     std::ofstream outAssociations;
     outAssociations.open("voxelAssos.sdp");     
     unsigned int th = vm["exportPointAssociations"].as<unsigned int>();
     for( const auto &p: imageAccumulation.domain())
       {
         if (imageAccumulation(p)>th)
           {
             outAssociations << "# associations of point: (" << p[0] << " " << p[1] << " "<< p[2]<< std::endl;
             NormalAccumulator::PointContainer neighbors = normAcc.getAssociatedPoints(p);
             for(const auto & pa: neighbors)
               {
                 outAssociations << p[0] << " " << p[1] << " " << p[2] << std::endl;
                 outAssociations << pa[0] << " " << pa[1] << " " << pa[2] << std::endl;                 
               }
           }         
       }
   }
 return 0;
}

