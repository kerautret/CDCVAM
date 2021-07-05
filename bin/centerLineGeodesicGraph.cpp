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

#include "CLI11.hpp"

#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "NormalAccumulator.h"
#include "GeodesicGraphComputer.h"


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
  std::string outputName {"result"};
  std::string inputMeshName;
  double dilateDist {2.0};
  bool invertNormal {false};
  bool initFromMaxRadius {false};
  std::string exportConfident {""};
  std::string importPointsNormals {""};
  std::string importNormals {""};

  double th {0.5};
  double thAcc {0.0};
  double deltaG {3.0};
  unsigned int maxComp {100};
  unsigned int minSizeCC {0} ;
  std::string estimRadiusType {"mean"};
  double radius {10.0};
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Center line extraction from accumulation  and geodesic graph");
  
  app.add_option("-i,--input,1", inputMeshName, "input mesh." )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputName, "the output base name.", true);
  app.add_option("--dilateDist", dilateDist, "dilate distance of the confidence voxels.", true);
  app.add_flag("--invertNormal", invertNormal, "invert normal to apply accumulation.");
  app.add_option("--exportConfident", exportConfident, "to export the confident voxel as a sequence of points", true);
  app.add_option("--th,-t",th, "threshold in the confidence estimation (included).", true);
  app.add_option("--thAcc,-a",thAcc,  "threshold on the accumultion (not resampled).", true);
  app.add_option("--importPointsNormals", importPointsNormals, "import normals and source points.", true);
  app.add_option("--importNormals", importNormals, "Use imported normals instead the one computed from the mesh faces.", true);
  app.add_option("--deltaG,-g", deltaG, "the param to consider interval of distances"
                 "to reconstruct the graph from the geodesic.", true);
  
  app.add_flag("--initFromMaxRadius", initFromMaxRadius, "init from max radius instead max confidence.");
  app.add_option("--maxComp,-m", maxComp, "the maximal number of cc of the graph.", true );
  app.add_option("--minSizeCC,-C", minSizeCC, "the min size of the graph Connected Component.", true );
  app.add_option("--radiusEstimator,-e", estimRadiusType,  "set the type of the radius estimation (mean (default), min, median or max).", true)
  -> check(CLI::IsMember({"max", "min", "mean", "median"}));
  app.add_option("--radius,-R", radius, "radius used to compute the accumulation.", true );
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

  
  NormalAccumulator acc(radius, estimRadiusType);
  DGtal::trace.info() << "------------------------------------ "<< std::endl;
  DGtal::trace.info() << "Step 1: Reading input mesh ... ";
  
  if(importPointsNormals != "")
    {
      // 1) Reading input mesh:
      // format point A and B from normals
      std::vector<Z3i::RealPoint> vPtA, vPtB;
      std::string namePt = importPointsNormals;
      vPtA = PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(namePt);
      std::vector<unsigned int> index= {3,4,5};
      vPtB = PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(namePt, index);
      std::vector<Z3i::RealPoint> vectNormals;
      for(unsigned int i=0; i<vPtB.size(); i++)
        {
          vectNormals.push_back((vPtB[i]-vPtA[i]).getNormalized());
        }
      // 2) Init accumulator:
      DGtal::trace.info() << "Step 2: Init accumulator ... ";
      acc.initFromNormals(vPtA, vectNormals);
    }
  else
    {
      // 1) Reading input mesh:
      DGtal::Mesh<DGtal::Z3i::RealPoint> aMesh;
      aMesh << inputMeshName;
      DGtal::trace.info() << " [done] " << std::endl;  
      DGtal::trace.info() << "------------------------------------ "<< std::endl;  
      if(importNormals != ""){
        std::vector<DGtal::Z3i::RealPoint> vectorNorm =
          PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(importNormals);
        acc.initFromMeshAndNormals(aMesh, vectorNorm, invertNormal);
    
      }else{
        acc.initFromMesh(aMesh, invertNormal);
      }
    }
  DGtal::trace.info() << " [done] " << std::endl;  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;  


  // 3) Compute accumulation and confidence
  DGtal::trace.info() << "Step 3: Compute accumulation/confidence ... ";
  acc.computeAccumulation();
  acc.computeConfidence(false, 30);
  acc.computeRadiusFromConfidence();
  Image3DDouble imageConfidence = acc.getConfidenceImage();
  Image3D imageAccumulation = acc.getAccumulationImage();
  Image3DDouble imageRadius = acc.getRadiusImage();

  
  DGtal::trace.info() << "Exporting ... ";
  typedef DGtal::functors::Rescaling< double, unsigned char> ScaleFctD;
  ScaleFctD  confidencescale (0 , 1.0, 0, 255);
  DGtal::VolWriter<NormalAccumulator::Image3DDouble,ScaleFctD>::exportVol("confidence.vol",
                                                                          imageConfidence,
                                                                          true, confidencescale);
  if(exportConfident != ""){
    std::string confExportName = exportConfident;
    ofstream exstream;
    exstream.open(confExportName.c_str());
    for (auto &p: imageConfidence.domain()){
      if(imageConfidence(p)>=th){
        exstream<< p[0] << " " << p[1] << " " << p[2] << std::endl;
      }
    }
    exstream.close();
  }
  DGtal::trace.info() << " [done] " << std::endl;  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;  
  

  // 4) Apply graph reconstruction from confidence image.
  // 4a) get initial sef of points:
  GeodesicGraphComputer::TSet aConfidenceSet(imageConfidence.domain());
  for (auto const &p: imageConfidence.domain())
    {
      if(imageConfidence(p)>=th && imageAccumulation(p) > thAcc)
        {
          aConfidenceSet.insert(p);
        }      
    }
  // 4b) set initial point: 
  DGtal::Z3i::Point startPoint; 
  if(initFromMaxRadius){
    startPoint = acc.getMaxRadiusPoint(); 
  }else {
    startPoint = acc.getMaxAccumulationPoint();   
  }
  trace.info() << "Starting point: " << startPoint << std::endl;
  GeodesicGraphComputer gg(deltaG, aConfidenceSet, dilateDist,  acc.getDomain(), startPoint);
  gg.computeGraphFromGeodesic(maxComp, minSizeCC);
  
  
  // 4c) graph export:
  // export vertex
  stringstream ss;
  ss<<outputName << "Vertex.sdp";
  std::ofstream outVertex;
  outVertex.open(ss.str().c_str());
  std::vector<Z3i::Point> vertices = gg.getGraphVerticesRepresentant();
  for (auto const &p: vertices)
    {
      outVertex << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }
  outVertex.close();


  // export edges
  stringstream ss2;
  ss2<<outputName << "Edges.sdp";
  std::ofstream outEdges;
  outEdges.open(ss2.str().c_str());
  std::vector<std::pair<unsigned int, unsigned int> > edges;
  edges = gg.myGraph.getEdges();
  for (auto const &e: edges)
    {
      outEdges << e.first << " " << e.second << std::endl;
    }
  outEdges.close();
  
  trace.info() << "] [Done]."<< std::endl ;



  // export radius
  stringstream ss3;
  ss3<<outputName << "Radius.sdp";
  std::ofstream outRadius;
  outRadius.open(ss3.str().c_str());
  for (const auto& e: vertices)
   {
     outRadius << imageRadius(e) << std::endl;
   }
  outRadius.close();
  
  return 0;
}
