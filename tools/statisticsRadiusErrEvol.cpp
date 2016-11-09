#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageContainerBySTLMap.h>
#include <DGtal/kernel/sets/DigitalSetFromMap.h>
#include <DGtal/kernel/sets/DigitalSetBySTLSet.h>
#include "DGtal/geometry/volumes/distance/FMM.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include <DGtal/math/Statistic.h>

#include "NormalAccumulator.h"



using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


typedef ImageContainerBySTLVector<Z3i::Domain, bool> Image3DMaker;
typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, double> Image3DDouble;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef ImageContainerBySTLVector<Z3i::Domain,double> DistanceImage;

typedef DigitalSetBySTLSet<Z3i::Domain> AcceptedPointSet;
typedef Z3i::Domain::Predicate DomainPredicate;
typedef L2FirstOrderLocalDistance<DistanceImage, AcceptedPointSet> DistanceMeasure;
typedef FMM<DistanceImage, AcceptedPointSet, DomainPredicate, DistanceMeasure> TFMM;

typedef FMM<DistanceImage, AcceptedPointSet, Z3i::DigitalSet, DistanceMeasure> TFMM2;







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
    ("output,o", po::value<std::string>()->default_value("result"), "the output base name.")
    ("invertNormal,n", "invert normal to apply accumulation.")
    ("estimRadiusType", po::value<std::string>()->default_value("mean"), "set the type of the"
                                               "radius estimation (mean, min, median or max).")
    ("radiusCompare,R", po::value<double>()->default_value(10.0), "radius used to compute the accumulation.")
    ("radiusComputeAcc,r", po::value<double>()->default_value(10.0), "radius used to compute the accumulation.")
    ("maxSampledOutValues,m", po::value<unsigned int>()->default_value(100), "set the maximal out values resulting of the stats.");
  
  
  
  bool parseOK = true;
  po::variables_map vm;
  try
  {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  }
  catch (const std::exception &ex)
  {
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
    parseOK = false;
  }
  po::notify(vm);
  if ( !parseOK || vm.count("help") || argc <= 1 || !vm.count("input") )
  {
    trace.info() << "Compute radius statistics obtained on confidence and accumulation estimation" << std::endl
                 << "Options: " << std::endl
                 << general_opt << std::endl;
    return 0;
  }


  // Reading parameters:
  std::string inputMeshName = vm["input"].as<std::string>();
  std::string outputName = vm["output"].as<std::string>();
  double radiusCompare = vm["radiusCompare"].as<double>();
  double radiusComputeAcc = vm["radiusComputeAcc"].as<double>();
  bool invertNormal = vm.count("invertNormal");
  std::string estimRadiusType = vm["estimRadiusType"].as<std::string>();
  unsigned int maxSampledOutValues = vm["maxSampledOutValues"].as<unsigned int>();
  
  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;
  DGtal::trace.info() << "Step 1: Reading input mesh ... ";

  // 1) Reading input mesh:
  DGtal::Mesh<DGtal::Z3i::RealPoint> aMesh;
  aMesh << inputMeshName;
  DGtal::trace.info() << " [done] " << std::endl;  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;  

  // 2) Init accumulator:
  DGtal::trace.info() << "Step 2: Init accumulator ... ";
  NormalAccumulator acc(radiusComputeAcc, estimRadiusType);
  acc.initFromMesh(aMesh, invertNormal);
  DGtal::trace.info() << " [done] " << std::endl;  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;  
  

  

  // 3) Compute accumulation and confidence
  DGtal::trace.info() << "Step 3: Compute accumulation/confidence ... ";
  acc.computeAccumulation();
  acc.computeRadiusFromOrigins();
  Image3DDouble imageRadiusAcc = acc.getRadiusImage();
  acc.computeRadiusFromConfidence();
  Image3DDouble imageRadiusConf = acc.getRadiusImage();

  Image3DDouble &imageConfidence = acc.getConfidenceImage();
  Image3D &imageAccumulation = acc.getAccumulationImage();
  
  double maxAccumulation = (double) (acc.getMaxAccumulation());
  // init empty vectors
  std::vector<DGtal::Statistic<double>> statsErrorAcc, statsErrorConf; 


  // represent error as pair with first: accumulation, second error (used to apply sort) 
  std::vector<std::pair<double, double > > vectAccWithError;
  std::vector<std::pair<double, double > > vectConfWithError;
  // compute stats
  for(auto const &p: imageRadiusAcc.domain())
    {
      double valConf = imageConfidence(p);
      double valAcc = imageAccumulation(p)/maxAccumulation;
      if(valConf>0)
        {
          double absErr = std::abs(imageRadiusConf(p)-radiusCompare);
          std::pair<double, double > confAndErr;
          confAndErr.first = imageConfidence(p);
          confAndErr.second = absErr;
          vectConfWithError.push_back(confAndErr);
        }
      if(valAcc>0){
        double absErr = std::abs(imageRadiusAcc(p)-radiusCompare);
        std::pair<double, double > accAndErr;
        accAndErr.first = imageAccumulation(p);
        accAndErr.second = absErr;
        vectAccWithError.push_back(accAndErr);
      }
    }

  
   // sorting values:
   DGtal::trace.info() << "Sorting acc errors from acc values:" ;
   std::sort(vectAccWithError.begin(),  vectAccWithError.end(), [](std::pair<double,double> a,std::pair<double,double> b ) -> bool { return a.first>b.first;});
   DGtal::trace.info() << " [done] \n Sorting acc errors from acc values:" ;
   std::sort(vectConfWithError.begin(),  vectConfWithError.end(), [](std::pair<double,double> a,std::pair<double,double> b ) -> bool { return a.first>b.first;});
   DGtal::trace.info() << " [done] " << std::endl;
  

  
   stringstream ssConfResName;
   ssConfResName<< outputName << "Confidence.dat";
   stringstream ssAccResName;
   ssAccResName<< outputName << "Accumulation.dat";
  
  
   std::ofstream outConf, outAcc;
   outConf.open(ssConfResName.str().c_str());
   outAcc.open(ssAccResName.str().c_str());

   outConf << "# statistics on repartion of radius errors from confidence: i Nun_of_error" << std::endl;
   outConf << "# radius used to compute acc: " << radiusComputeAcc << "radius to compare val: " << radiusCompare << std::endl;
 
   outAcc << "# statistics on repartition of radius errors from accumulation: i Nun_of_error"  << std::endl;
   outAcc << "# radius used to compute acc: " << radiusComputeAcc << "radius to compare val: " << radiusCompare << std::endl;
 

  for(unsigned int i = 0; i < maxSampledOutValues && i< vectConfWithError.size() ; i++)
  {
    outConf << i << " " <<  vectConfWithError[i].second << std::endl;
  }
  for(unsigned int i = 0; i < maxSampledOutValues && i< vectAccWithError.size(); i++)
  {
    outAcc << i << " " <<  vectAccWithError[i].second << std::endl;
  }

  
   outAcc.close();
   outConf.close();
    
  
  DGtal::trace.info() << "[Done]" <<std::endl;

  return 0;
}



