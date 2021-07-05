#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/MeshReader.h"
#include <DGtal/io/readers/PointListReader.h>
#include <DGtal/math/Statistic.h>

#include "CLI11.hpp"
#include "NormalAccumulator.h"

using namespace std;
using namespace DGtal;

typedef ImageContainerBySTLVector<Z3i::Domain, bool> Image3DMaker;
typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, double> Image3DDouble;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef ImageContainerBySTLVector<Z3i::Domain,double> DistanceImage;


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
  std::string outputGlobalStat {"resultGlobal"};
  double confidenceTh {0.5};


  double radiusCompare {10.0};
  double radiusComputeAcc {10.0};
  bool commonSet {false};
  std::string manualFixRefName {""};
  

  // parse command line using CLI ----------------------------------------------
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
  app.add_option("--outputGlobalStat", outputGlobalStat, "basename of the global result (mean max, min) (result appened to the end of file.", true);
 
  app.add_option("--manualFixRefName", manualFixRefName,  "use a manual reference set instead using confidence threshold", true);
  
  app.add_option("--confidenceTh", confidenceTh, "with this options errors are calculated on a same set defined from confidence threshold", true);
  app.add_option("--radiusCompare,-R", radiusCompare, "reference Radius to apply comparisons.", true);
  app.add_option("--radiusComputeAcc,-r", radiusComputeAcc, "radius used to compute the accumulation.", true);
  app.add_flag("--commonSet,-c", commonSet, "with this options errors are computed on a same set defined form confidence threshold 0.5.");

  
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
  std::vector<unsigned int> histoErrorAcc, histoErrorConf;
  DGtal::Statistic<double> statsErrorConf;
  DGtal::Statistic<double> statsErrorAcc;
  
  double stepRes = 200.0;
  for (unsigned int i = 0; i<= stepRes; i++)
  {
    histoErrorConf.push_back(0);
    histoErrorAcc.push_back(0);
  }
  
  // compute stats by thresholding common confidence.
  if(manualFixRefName == ""){
    for(auto const &p: imageRadiusAcc.domain())
    {
      double valConf = imageConfidence(p);
      double valAcc = imageAccumulation(p)/maxAccumulation;
      if((valConf>0 && !commonSet) || (valConf > confidenceTh  ))
      {
        double absErr = std::abs(imageRadiusConf(p)-radiusCompare);
        statsErrorConf.addValue(absErr);
        unsigned int index = (absErr/radiusCompare)*stepRes;
        histoErrorConf[index]++;
      }
      if((valAcc>0 && !commonSet) || (valConf > confidenceTh ))
      {
        double absErr = std::abs(imageRadiusAcc(p)-radiusCompare);
        statsErrorAcc.addValue(absErr);
        unsigned int index = (absErr/radiusCompare)*stepRes;
        histoErrorAcc[index]++;
      }
    }
  }else{
    // computes stat on manual sdp fix:
    std::vector<Z3i::Point> inputPoints = PointListReader<Z3i::Point>::getPointsFromFile(manualFixRefName);
    for( auto const &p: inputPoints)
    {
      double absErr = std::abs(imageRadiusConf(p)-radiusCompare);
      statsErrorConf.addValue(absErr);
      unsigned int index = (absErr/radiusCompare)*stepRes;
      histoErrorConf[index]++;

      double absErrA = std::abs(imageRadiusAcc(p)-radiusCompare);
      statsErrorAcc.addValue(absErrA);
      unsigned int indexA = (absErr/radiusCompare)*stepRes;
      histoErrorAcc[indexA]++;
      
    }
  }
  
  
  stringstream ssConfResName;
  ssConfResName<< outputName << "Confidence.dat";
  stringstream ssAccResName;
  ssAccResName<< outputName << "Accumulation.dat";
  
  std::ofstream outConf, outAcc;
  outConf.open(ssConfResName.str().c_str());
  outAcc.open(ssAccResName.str().c_str());
  
  outConf << "# statistics on repartion of radius errors from confidence: i Nun_of_error"
  << std::endl;
  outAcc << "# statistics on repartition of radius errors from accumulation: i Nun_of_error"
  << std::endl;
  
  for(unsigned int i = 0; i <= stepRes; i++)
  {
    outConf << i << " " <<  histoErrorConf[i] << std::endl;
    outAcc  << i << " " << histoErrorAcc[i] << std::endl;
  }
  outAcc.close();
  outConf.close();
  
  stringstream ssConf, ssAcc;
  ssConf <<  outputGlobalStat << "Confidence.dat";
  ssAcc << outputGlobalStat << "Accumulation.dat";
  std::ofstream outConfGlob, outAccGlob;
  std::ifstream inConfGlob; inConfGlob.open(ssConf.str().c_str());
  bool isNew = inConfGlob.peek() == std::ifstream::traits_type::eof();  inConfGlob.close();
  outConfGlob.open(ssConf.str().c_str(), std::ofstream::out | std::ofstream::app);
  outAccGlob.open(ssAcc.str().c_str(), std::ofstream::out | std::ofstream::app);
  if (isNew){
    outConfGlob << "# stats on confidence : Radius_ComputeAcc Radius_Compare Mean Min Max VAR UVAR " << std::endl;
    outAccGlob << "# stats on accumulation : Radius_ComputeAcc Radius_Compare Mean Min Max VAR UVAR " << std::endl;
    
  }
  
  statsErrorConf.terminate();
  statsErrorAcc.terminate();
  outConfGlob << radiusComputeAcc << " " << radiusCompare << " "<<  statsErrorConf.mean() << "  " << statsErrorConf.min()
  <<  "  " << statsErrorConf.max() << " " << statsErrorConf.variance() << " " <<statsErrorConf.unbiasedVariance() <<  std::endl;
  
  outAccGlob  << radiusComputeAcc << " " << radiusCompare << " "<<statsErrorAcc.mean() << " " << statsErrorAcc.min()
  <<  "  " << statsErrorAcc.max() <<  " " << statsErrorConf.variance() << " " << statsErrorConf.unbiasedVariance() << std::endl;
  
  outConfGlob.close();
  outAccGlob.close();
  
  DGtal::trace.info() << "[Done]" <<std::endl;
  
  return 0;
}
