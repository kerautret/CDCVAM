#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageContainerBySTLMap.h"
#include "DGtal/kernel/sets/DigitalSetFromMap.h"
#include "DGtal/kernel/sets/DigitalSetBySTLSet.h"
#include "DGtal/geometry/volumes/distance/FMM.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/math/Statistic.h"

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
  std::string estimRadiusType {"mean"};
  double radius {10.0};
  
  
  CLI::App app;
  app.description("Compute radius statistics obtained on confidence and accumulation estimation");
  app.add_option("-i,--input,1", inputMeshName, "input mesh" )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputName, "the output base name.", true );
  app.add_option("--invertNormal,-n", invertNormal, "invert normal to apply accumulation.");
  app.add_option("--estimRadiusType", estimRadiusType,  "set the type of the"
                 "radius estimation (mean (default), min, median or max).", true)
    -> check(CLI::IsMember({"max", "min", "mean", "median"}));
  app.add_option("--radius,-R",radius, "radius used to compute the accumulation.", true);
  

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
  NormalAccumulator acc(radius, estimRadiusType);
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

  double stepRes = 100.0;
  for (unsigned int i = 0; i<= stepRes; i++)
    {
      statsErrorConf.push_back(DGtal::Statistic<double>());
      statsErrorAcc.push_back(DGtal::Statistic<double>());
    }
  
  // compute stats
  for(auto const &p: imageRadiusAcc.domain())
    {
      double valConf = imageConfidence(p);
      double valAcc = imageAccumulation(p)/maxAccumulation;
      
      if(valConf>0)
        {
          double errorConf = (imageRadiusConf(p)-radius)*(imageRadiusConf(p)-radius);
          statsErrorConf[(int)(valConf*stepRes)].addValue(errorConf);

        }
      
      if(valAcc>0){
        double errorAcc = (imageRadiusAcc(p)-radius)*(imageRadiusAcc(p)-radius);
        statsErrorAcc[(int)(valAcc*stepRes)].addValue(errorAcc);
      }

    }
  
  stringstream ssConfResName;
  ssConfResName<< outputName << "Confidence.dat";
  stringstream ssAccResName;
  ssAccResName<< outputName << "Accumulation.dat";
  
  
  std::ofstream outConf, outAcc;
  outConf.open(ssConfResName.str().c_str());
  outAcc.open(ssAccResName.str().c_str());

  outConf << "# statistics on Radius estimation on Confidence: thresholdIntervalle mean max min variance unbiased variance " << std::endl;
  outAcc << "# statistics on Radius estimation on Accumulation: thresholdIntervalle  mean max min variance unbiased variance " << std::endl;

  
  for(unsigned int i = 0; i < stepRes; i++)
    {
      statsErrorConf[i].terminate();
      statsErrorAcc[i].terminate();
      outConf << i << " " << statsErrorConf[i].mean() << " " <<  statsErrorConf[i].max() << " " << statsErrorConf[i].min()
              << " "  <<  statsErrorConf[i].variance() << " " << statsErrorConf[i].unbiasedVariance()
              << " " << statsErrorConf[i].samples() << std::endl;
      
      outAcc  << i << " " << statsErrorAcc[i].mean() << " " <<  statsErrorAcc[i].max() << " " << statsErrorAcc[i].min()
              << " "  <<  statsErrorAcc[i].variance() << " " << statsErrorAcc[i].unbiasedVariance()
              << " " << statsErrorAcc[i].samples() << std::endl;
    }
  outAcc.close();
  outConf.close();
    
  
  DGtal::trace.info() << "[Done]" <<std::endl;

  return 0;
}

