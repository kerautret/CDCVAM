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
  std::string inputMeshName ;
  std::string outputName {"result"};
  bool invertNormal {false};

  double radiusCompare {10.0};
  double radiusComputeAcc {10.0};

  double confidenceTh {0.5};
  double accumulationTh {0.5};
  
  CLI::App app;
  app.description("Computes centricity statistics obtained on confidence and accumulation estimation.");
  app.add_option("-i,--input,1", inputMeshName, "input mesh" )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputName, "the output base name.", true );
  app.add_option("--invertNormal,-n", invertNormal, "invert normal to apply accumulation.");
  app.add_option("--radiusCompare,-r", radiusCompare, "reference Radius to apply comparisons.", true);
  app.add_option("--radiusComputeAcc,-d", radiusComputeAcc, "distance (d_{acc}) used to compute the accumulation.", true);
  app.add_option("--confidenceTh,-c", confidenceTh, "threshold confidence image to apply stat.", true);
  app.add_option("--accumulationTh,-a", confidenceTh, "threshold accumulation image to apply stat.", true);
    
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
  // IMPORTANT we use min in order to have the image radius equivalent to an approximation of the hausdorff distance.
  NormalAccumulator acc(radiusComputeAcc, "min");
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
  double stepRes = 50.0;
  for (unsigned int i = 0; i<= stepRes; i++)
  {
    statsErrorConf.push_back(DGtal::Statistic<double>());
    statsErrorAcc.push_back(DGtal::Statistic<double>());
  }
  
  // compute stats
  DGtal::trace.progressBar(0, imageRadiusAcc.size());
  for(auto const &p: imageRadiusAcc.domain())
    {
      double valConf = imageConfidence(p);
      double valAcc = imageAccumulation(p)/maxAccumulation;
      double minDistance = imageRadiusAcc(p);
      
      if(valConf > confidenceTh )
      {
        double errorConf = (radiusCompare-minDistance)/radiusCompare;
        unsigned int index = (int)(valConf*stepRes);
        for(unsigned int i = 0; i <=index; i++){
          statsErrorConf[i].addValue(errorConf);
        }
      }
      if(valAcc > accumulationTh)
      {
        double errorAcc = (radiusCompare-minDistance)/radiusCompare;
        unsigned int index = (int)(valAcc*stepRes);
        for(unsigned int i = 0; i <=index; i++){
          statsErrorAcc[i].addValue(errorAcc);
        }
      }
    }
  
  stringstream ssConfResName;
  ssConfResName<< outputName << "Confidence.dat";
  stringstream ssAccResName;
  ssAccResName<< outputName << "Accumulation.dat";
  
  
  std::ofstream outConf, outAcc;
  outConf.open(ssConfResName.str().c_str());
  outAcc.open(ssAccResName.str().c_str());
  outConf << "#Stats generated from statisticsCentricity (from confidence) prog: MEAN MAX MIN VAR UVAR NB" << std::endl;
  outAcc << "#Stats generated from statisticsCentricity (from accumulation) prog: MEAN MAX MIN VAR UVAR NB" << std::endl;
  
  
  for(unsigned int i = 0; i < stepRes; i++)
  {
    statsErrorConf[i].terminate();
    statsErrorAcc[i].terminate();
    outConf << i << " " << statsErrorConf[i].mean() << " " <<  statsErrorConf[i].max() << " " << statsErrorConf[i].min()
    << " "  <<  statsErrorConf[i].variance() << " " << statsErrorConf[i].unbiasedVariance() << " "
    << statsErrorConf[i].samples() << std::endl;
    
    outAcc  << i << " " << statsErrorAcc[i].mean() << " " <<  statsErrorAcc[i].max() << " " << statsErrorAcc[i].min()
    << " "  <<  statsErrorAcc[i].variance() << " " << statsErrorAcc[i].unbiasedVariance() << " "
    << statsErrorAcc[i].samples() << std::endl;
  }
  outAcc.close();
  outConf.close();
  
  DGtal::trace.info() << "[Done]" <<std::endl;
  return 0;
}



