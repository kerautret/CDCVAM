#include <iostream>
#include <stdio.h>
#include <stdlib.h>


#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/io/writers/LongvolWriter.h>

#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/writers/MeshWriter.h>

#include <DGtal/images/ImageContainerBySTLVector.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "NormalAccumulator.h"


using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain,  double> ImageDouble;

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
    ("input,i", po::value<std::string>(), "input mesh.")
    ("outputAcc,o", po::value<std::string>()->default_value("accumulation.longvol"), "output accumulation longvol file.")
    ("outputConf,c", po::value<std::string>()->default_value("confidence.longvol"), "output confidence longvol file.")
    ("outputRad,R", po::value<std::string>()->default_value("radius.longvol"), "output radius vol file.")
    ("invertNormal", "invert input normal.")
    ("radiusEstimator,e", po::value<std::string>()->default_value("min"), 
     "use: {min (default), max, mean, median} to estimate the radius") 
    ("maxValOutConf", po::value<DGtal::uint64_t>()->default_value(255),
     "set MAX scale of confidence out image: 0 .. 1 -> 0 ..  MAX.")
    ("maxValOutRad", po::value<DGtal::uint64_t>()->default_value(255),
     "set MAX scale of radius out image: 0 .. 1 -> 0 ..  MAX.")
    ("autoScaleAcc", po::value<std::string>(), "auto scale out accumulation image values between 0 and 255.")
    ("autoScaleConf", po::value<std::string>(), "auto scale out condidence image values between 0 and 255.")
    ("radius,r", po::value<double>()->default_value(5), "radius of accumulation analysis.");
  
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
    trace.info() << "Compute mesh accumulation from a mesh and compute the radius (median value) "
                 << "of all faces participating to the accumulation" << std::endl
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
  bool invertNormal = vm.count("invertNormal");
  unsigned int outConfidenceMax = vm["maxValOutConf"].as<DGtal::uint64_t>();
  
  
  // 1) Reading input mesh
  Mesh<Z3i::RealPoint> aMesh(true);
  aMesh << inputFile;
    
  // 2) Init an NormalAccumulator
  NormalAccumulator normAcc(radius, typeStat);
  normAcc.initFromMesh(aMesh, invertNormal);

  // 3) Generate accumulation image
  normAcc.computeAccumulation();
  NormalAccumulator::Image3D &imageAccumulation = normAcc.getAccumulationImage();

  trace.info() << "Saving accumulation image in " << outputFileAcc;
  LongvolWriter<Image3D>::exportLongvol(outputFileAcc, imageAccumulation);
  trace.info() << "[done]" << std::endl;

  if (vm.count("autoScaleAcc")){
    trace.info() << "Saving accumulation (auto scale 0 " << normAcc.getMaxAccumulation() 
                 <<"  -> [0, 255]) image in " << outputFileAcc << " ... ";
    typedef functors::Rescaling<DGtal::uint64_t, unsigned char> ScaleFct;
    ScaleFct  scaleFct (0.0 ,normAcc.getMaxAccumulation(), 0, 255);
    VolWriter<Image3D,ScaleFct>::exportVol(vm["autoScaleAcc"].as<std::string>(), imageAccumulation, true, scaleFct);
    trace.info() << "[done]" << std::endl;
  }
  
  
  // 4) Compute confidence image 
  normAcc.computeConfidence();
  ImageDouble imageConfidence = normAcc.getConfidenceImage();
  trace.info() << "Saving confidence image in " << outputFileConf << " ... ";
  typedef functors::Rescaling<double, DGtal::uint64_t> ScaleFctD;
  ScaleFctD scaleFct(0.0, 1.0, 0, vm["maxValOutConf"].as<DGtal::uint64_t>());
  LongvolWriter<ImageDouble, ScaleFctD>::exportLongvol(outputFileConf, imageConfidence, true, scaleFct);
  trace.info() << "[done]" << std::endl;
  
  if(vm.count("autoScaleConf"))
    {
      std::string outNameAutoConfidence = vm["autoScaleConf"].as<std::string>();
      trace.info() << "Saving confidence (auto scale 0 1 -> 0 255) image in " << outNameAutoConfidence << " ... ";
      typedef functors::Rescaling<double, unsigned char> ScaleFctD;
      ScaleFctD  scaleFct (0.0 ,1.0, 0, 255);
      VolWriter<ImageDouble,ScaleFctD>::exportVol(outNameAutoConfidence, imageConfidence, true, scaleFct);
      trace.info() << "[done]" << std::endl;    
    }
  

  // 5) Compute Radius image
  normAcc.computeRadiusFromConfidence();
  ImageDouble imageRadius = normAcc.getRadiusImage();  
  trace.info() << "Saving radius image in " << outputFileRad << " ... ";
  ScaleFctD scaleFct2(0.0, normAcc.getMaxRadius(), 0, vm["maxValOutRad"].as<DGtal::uint64_t>());
  LongvolWriter<ImageDouble, ScaleFctD>::exportLongvol(outputFileRad, imageRadius, true, scaleFct2);
  trace.info() << "[done]" << std::endl;
 if(vm.count("autoScaleRad"))
    {
      std::string outNameAutoRadius = vm["autoScaleRad"].as<std::string>();
      trace.info() << "Saving confidence (auto scale 0 " << normAcc.getMaxRadius()
                   << " ->[0,255]) image in " << outNameAutoRadius << " ... ";
      typedef functors::Rescaling<double, unsigned char> ScaleFctD;
      ScaleFctD  scaleFct (0.0 ,normAcc.getMaxRadius(), 0, 255);
      VolWriter<ImageDouble,ScaleFctD>::exportVol(outNameAutoRadius, imageRadius, true, scaleFct);
      trace.info() << "[done]" << std::endl;    
    }    
  return 0;
}

