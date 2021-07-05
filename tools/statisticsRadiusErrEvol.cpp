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


#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include <DGtal/math/Statistic.h>

#include "NormalAccumulator.h"

#include "CLI11.hpp"


using namespace std;
using namespace DGtal;


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
  // Reading parameters:
  std::string inputMeshName;
  std::string outputName {"result"};
  bool invertNormal {false};
  std::string estimRadiusType {"mean"};

  double radiusCompare {10.0};
  double radiusComputeAcc {10.0};
  unsigned int maxSampledOutValues {100};
  
  CLI::App app;
  app.description("Compute radius statistics obtained on confidence and accumulation estimation");
  app.add_option("-i,--input,1", inputMeshName, "input mesh" )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputName, "the output base name.", true );
  app.add_option("--invertNormal,-n", invertNormal, "invert normal to apply accumulation.");
  app.add_option("--estimRadiusType,-e", estimRadiusType,  "set the type of theradius estimation (mean, min, median or max). ", true)
  -> check(CLI::IsMember({"max", "min", "mean", "median"}));
  app.add_option("--radiusCompare,-R", radiusCompare, "radius used to compute the accumulation.", true);
  app.add_option("--radiusComputeAcc,-r", radiusComputeAcc, "radius used to compute the accumulation.", true);
  app.add_option("--maxSampledOutValues,-m", radiusComputeAcc, "set the maximal out values resulting of the stats.", true);

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  
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



