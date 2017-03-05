#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageContainerBySTLMap.h>
#include <DGtal/kernel/sets/DigitalSetFromMap.h>
#include <DGtal/kernel/sets/DigitalSetBySTLSet.h>
#include "DGtal/geometry/volumes/distance/FMM.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "NormalAccumulator.h"
#include "CenterLineHelper.h"
#include "GeodesicGraphComputer.h"


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


  std::string inputMeshName = "test.off";

  double dilateDist = 2.0;
  double th = 0.5;
  double radius = 7.0;
  double radiusOpti = 5.5;
  bool invertNormal = false;
  std::string estimRadiusType = "mean";
  double deltaG = 4.0;
  unsigned int maxComp = 10;
  unsigned int minSizeCC = 0;
  
  
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
  acc.computeConfidence(false);
  //acc.computeRadiusFromConfidence();
  Image3DDouble &imageConfidence = acc.getConfidenceImage();
  // Image3DDouble &imageRadius = acc.getRadiusImage();

  
  DGtal::trace.info() << " [done] " << std::endl;  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;  
  

  
  // 4) Apply graph reconstruction from confidence image.
  // 4a) get initial sef of points:
  GeodesicGraphComputer::TSet aConfidenceSet(imageConfidence.domain());
  for (auto const &p: imageConfidence.domain())
    {
      if(imageConfidence(p)>=th)
        {
          aConfidenceSet.insert(p);
        }      
    }
  // 4b) set initial point: 
  DGtal::Z3i::Point startPoint; 

  // init from max radius
  //  startPoint = acc.getMaxRadiusPoint(); 

  // init from max accumulation:
  startPoint = acc.getMaxAccumulationPoint();   
  
  trace.info() << "Starting point: " << startPoint << std::endl;
  GeodesicGraphComputer gg(deltaG, aConfidenceSet, dilateDist,  acc.getDomain(), startPoint);
  //gg.setConfidenceImage(&(acc.getConfidenceImage()));
  gg.computeGraphFromGeodesic(maxComp, minSizeCC);
  

  
  // 4c) graph export no opti:
  std::ofstream resBrut;
  resBrut.open("resTestCenterLineHelperBrut.sdp");
  std::vector<DGtal::Z3i::Point> vertices = gg.getGraphVerticesRepresentant();
  for (auto const &p: vertices)
    {
      resBrut << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
  resBrut.close();
  


  //-------------------
  // test optimisation
  //-------------------
  std::vector<DGtal::Z3i::RealPoint> centerLineOpti = CenterLineHelper::optimizeCenterLineElasticForces(vertices,
                                                                                                        acc,  radiusOpti, 0.01);
  
  std::ofstream resOpti;
  resOpti.open("resTestCenterLineHelperOpti.sdp");
  for (auto const &p: centerLineOpti )
    {
      resOpti << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
  resOpti.close();
                                                                                                       
                                                                                                       

  return 0;
}
