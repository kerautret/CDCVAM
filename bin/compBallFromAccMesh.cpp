#include <iostream>
#include <fstream>
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
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageContainerBySTLMap.h>
#include <DGtal/kernel/sets/DigitalSetFromMap.h>


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "NormalAccumulator.h"


using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef ImageContainerBySTLVector<Z3i::Domain,  Z3i::RealPoint> ImageVector;
typedef ImageContainerBySTLVector<Z3i::Domain,  double> ImageDouble;

typedef ImageContainerBySTLMap<Z3i::Domain,double> DistanceImage;





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
    ("output,o", po::value<std::string>()->default_value("out.sdp"), "output file representing sdp with radius.")
    ("invertNormal", "invert input normal.")
    ("radiusEstimDef,e", po::value<std::string>()->default_value("min"), 
     "use: {min (default), max, mean, median} to estimate the radius")
    ("addColor", po::value<std::vector<unsigned int> >()->multitoken(), "insert Color in the sdp export.")
    ("minConfidence,m",po::value<double>()->default_value(0.5), "set the minimal confidance rate of ball extraction.")
    ("minAccumulation,a",po::value<unsigned int>()->default_value(0), "set the minimal accumulation of ball extraction.")
    ("minRadius",po::value<double>()->default_value(0), "min threshold on radius image.")
    ("maxRadius",po::value<double>()->default_value(std::numeric_limits<double>::max()), "max threshold on radius image.")
    ("radius,r", po::value<double>()->default_value(5), "maximal radius of accumulation analysis.");

  
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
    trace.info() << "Compute ball centers with radius defined from accumulation image. " << std::endl
                 << "Options: " << std::endl
                 << general_opt << "\n";
    return 0;
  }
  
  double radius = vm["radius"].as<double>();
  string inputFile = vm["input"].as<std::string>();
  string outputFile = vm["output"].as<std::string>();
  bool invertNormal = vm.count("invertNormal");
  double minConfidence = vm["minConfidence"].as<double>();
  unsigned int minAccumulation = vm["minAccumulation"].as<unsigned int>();
  string typeStat = vm["radiusEstimDef"].as<string>();
  double minRadius = vm["minRadius"].as<double>();
  double maxRadius = vm["maxRadius"].as<double>();  
  
  Color colorAdded(250, 250, 250);
  if(parseOK && vm.count("addColor")){
    std::vector<unsigned int> vcol= vm["addColor"].as<std::vector<unsigned int > >();
    if(vcol.size()<4){
      trace.error() << " Not enough parameter: color specification should contains four elements: red, green, blue and alpha values "
                    << "(Option --addColor ignored). "  << std::endl;
    }
    colorAdded.setRGBi(vcol[0], vcol[1], vcol[2], vcol[3]);
  }
  
  // 1) Reading input mesh
  Mesh<Z3i::RealPoint> aMesh(true);
  aMesh << inputFile;

  // 2) Init an NormalAccumulator
  NormalAccumulator normAcc(radius, typeStat);
  normAcc.initFromMesh(aMesh, invertNormal);
  
  // 3) compute accumulation with confidence rate 
  normAcc.computeRadiusFromConfidence();
  
  
  // 4) Extract centers from resulting acc and extract radius
  std::ofstream out;
  out.open(outputFile);
  NormalAccumulator::Image3DDouble &imageConf = normAcc.getConfidenceImage();
  NormalAccumulator::Image3DDouble &imageRadius = normAcc.getRadiusImage();
  NormalAccumulator::Image3D &imageAccumulation = normAcc.getAccumulationImage();
  
  unsigned int nb = 0;
  for(auto &v: imageConf.domain()){
    if(imageRadius(v)>0  && imageConf(v)>minConfidence && 
       imageAccumulation(v) > minAccumulation && imageRadius(v) < maxRadius
       && imageRadius(v) > minRadius){
      out << v[0] << " " << v[1] <<  " " << v[2] << " " << imageRadius(v);
      if (vm.count("addColor")){
        out << " " << (unsigned int) colorAdded.red() << " " 
            << (unsigned int) colorAdded.green() << " " << (unsigned int) colorAdded.blue();
      }
      out << std::endl;
      nb++;
    }
  } 
  DGtal::trace.info() << "total exported:" << nb << std::endl;
  out.close();

    
  
  
  return 0;
}
