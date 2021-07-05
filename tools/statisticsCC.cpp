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


typedef DGtal::DigitalSetBySTLSet<DGtal::Z3i::Domain> TSet;

using namespace std;
using namespace DGtal;


typedef ImageContainerBySTLVector<Z3i::Domain, bool> Image3DMaker;
typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, double> Image3DDouble;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef ImageContainerBySTLVector<Z3i::Domain,double> DistanceImage;
typedef ImageContainerBySTLVector<Z3i::Domain, Z3i::Point> Image3DPoint;


typedef DigitalSetBySTLSet<Z3i::Domain> AcceptedPointSet;
typedef Z3i::Domain::Predicate DomainPredicate;
typedef L2FirstOrderLocalDistance<DistanceImage, AcceptedPointSet> DistanceMeasure;
typedef FMM<DistanceImage, AcceptedPointSet, DomainPredicate, DistanceMeasure> TFMM;
typedef FMM<DistanceImage, AcceptedPointSet, Z3i::DigitalSet, DistanceMeasure> TFMM2;


Z3i::Point find(Image3DPoint &imageParent, const Z3i::Point &p)
{
  if (imageParent(p) != p)
    return find(imageParent, imageParent(p));
  return p;
}

void
unionF(Image3DPoint &imageParent, const Z3i::Point &p1,
                                  const Z3i::Point &p2)
{
  Z3i::Point xset = find(imageParent, p1);
  Z3i::Point yset = find(imageParent, p2);
  imageParent.setValue(xset, yset);
}


template<typename TImage>
unsigned int
computeCC_UF(const TImage &image, typename TImage::Value th )
{
  Z3i::Domain dom = image.domain();
  Image3DPoint imageParent (dom);
  for( const auto &p: dom ) imageParent.setValue(p, p);
  for( const auto &p: dom ){
    if (image(p) <= th){
      continue;
    }
    for(int i=-1; i<2; i++ )
    {
      for(int j=-1; j<2; j++ )
      {
        for(int k=-1; k<2; k++ )
        {
          Z3i::Point pv = p + Z3i::Point(i,j,k);
          if( p != pv && dom.isInside(pv) && image(pv) >= th)
          {
            unionF(imageParent, p, pv);
          }
        }
      }
    }
  }
  // counting elements insite the th and for which image(p)==p
  int nbComp = 0;
  for( const auto &p: dom ){
    if( p ==imageParent(p) && image(p) >= th ){
      nbComp++;
    }
  }
  return nbComp;
}







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
  int indexExportSDP {50};
  double radiusComputeAcc {10.0};
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Computes centricity on the nb of CC obtained on confidence and accumulation estimation.");
  app.add_option("-i,--input,1", inputMeshName, "input mesh." )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputName, "the output base name.", true);
  app.add_option("--invertNormal,-n", invertNormal, "invert normal to apply accumulation.");
  app.add_option("--indexExportSDP,-e",indexExportSDP, "the threhold to export the sdp file of elements" );
  app.add_option("--radiusComputeAcc,-d",radiusComputeAcc, "distance (d_{acc}) used to compute the accumulation.", true );
  
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
  
  ofstream outSDPAcc, outSDPConf;
  outSDPAcc.open("setAcc.sdp");
  outSDPConf.open("setConf.sdp");
  
  double maxAccumulation = (double) (acc.getMaxAccumulation());
  // init empty vectors
  
  std::vector<unsigned int > vectCCAcc, vectCCConf;
  double stepRes = 100.0;
  trace.progressBar(0, stepRes);
  // computing the set of confidence and accumulation
  TSet setAcc (imageConfidence.domain());
  TSet setConf (imageConfidence.domain());

  
  for (unsigned int i = 0; i<= stepRes; i++)
  {
    trace.progressBar(i, stepRes);
    double threshold = (stepRes-i)/stepRes;
    if(i==indexExportSDP){
      for(auto const &p: imageRadiusAcc.domain()){
        double valConf = imageConfidence(p);
        double valAcc = imageAccumulation(p)/maxAccumulation;
        if(valConf > threshold){
           outSDPConf << p[0] << " " << p[1] << "  " << p[2] << std::endl;
        }
        if(valAcc > threshold){
          outSDPAcc << p[0] << " " << p[1] << "  " << p[2] << std::endl;
        }
        
      }
    }
  
    unsigned int nbCCConf  = computeCC_UF(imageConfidence, threshold);
    unsigned int nbCCAcc = computeCC_UF(imageAccumulation, threshold*maxAccumulation);
    vectCCConf.push_back(nbCCConf);
    vectCCAcc.push_back(nbCCAcc);
    trace.info() << "nb cc conf: " << nbCCConf << std::endl;
    trace.info() << "nb cc acc: " << nbCCAcc << std::endl;

  }
  
  trace.progressBar(stepRes, stepRes);
  outSDPAcc.close();
  outSDPConf.close();
  
  
  stringstream ssConfResName;
  ssConfResName<< outputName << "Confidence.dat";
  stringstream ssAccResName;
  ssAccResName<< outputName << "Accumulation.dat";
  
  
  std::ofstream outConf, outAcc;
  outConf.open(ssConfResName.str().c_str());
  outAcc.open(ssAccResName.str().c_str());
  outConf << "#Stats generated from statisticsCC (from confidence) prog: Index_Th NB_CC" << std::endl;
  outAcc << "#Stats generated from statisticsCC (from accumulation) prog: Index_Th NB_CC" << std::endl;
  
  
  for(unsigned int i = 0; i < stepRes; i++)
  {
       outConf << i << " " << vectCCConf[stepRes-i-1]<< std::endl;
       outAcc << i << " " << vectCCAcc[stepRes-i-1]<< std::endl;
  }
  outAcc.close();
  outConf.close();
  
  
  
  DGtal::trace.info() << "[Done]" <<std::endl;

  return 0;
}



