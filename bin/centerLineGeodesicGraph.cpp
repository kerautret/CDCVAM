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

#include "NormalAccumulator.h"
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
  po::options_description general_opt("Allowed options are: ");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string>(), "input mesh.")
    ("output,o", po::value<std::string>()->default_value("result"), "the output base name.")
    ("dilateDist", po::value<double>()->default_value(2.0), "dilate distance of the confidence voxels.")
    ("invertNormal,n", "invert normal to apply accumulation.")
    ("exportGeodesic,e",po::value<std::string>(), "to export the geodesic as a sequence of points")
    ("exportConfident,c",po::value<std::string>(), "to export the confident voxel as a sequence of points")
    ("importPointsNormals", po::value<std::string>(), "import normals and source points.")
    ("importNormals", po::value<std::string>(), "Use imported normals instead the one computed from the mesh faces.")
    
    ("deltaG,g",  po::value<double>()->default_value(3.0), "the param to consider interval of distances "
                                                            "to reconstruct the graph from the geodesic.")
    ("initFromMaxRadius", "init from max radius instead max confidence.")
    ("th,t", po::value<double>()->default_value(0.5), "threshold in the confidence estimation (included).")
    ("thAcc,a", po::value<double>()->default_value(0), "threshold on the accumultion (not resampled).")
    ("maxComp,m", po::value<unsigned int>()->default_value(100), "the maximal number of cc of the graph.")
    ("minSizeCC,C", po::value<unsigned int>()->default_value(0), "the min size of the graph Connected Component.")
    ("estimRadiusType", po::value<std::string>()->default_value("mean"), "set the type of the"
                                              "radius estimation (mean, min, median or max).")
    ("radius,R", po::value<double>()->default_value(10.0), "radius used to compute the accumulation.");
  

  
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
  if ( !parseOK || vm.count("help") || argc <= 1 || (!vm.count("input")&&!vm.count("importPointsNormals") ))
  {
    trace.info() << "Center line extraction from accumulation  and geodesic graph" << std::endl
                 << "Options: " << std::endl
                 << general_opt << std::endl;
    return 0;
  }


  // Reading parameters:
  
  std::string outputName = vm["output"].as<std::string>();
  double dilateDist = vm["dilateDist"].as<double>();
  double th = vm["th"].as<double>();
  double thAcc = vm["thAcc"].as<double>();
  double radius = vm["radius"].as<double>();
  bool invertNormal = vm.count("invertNormal");
  std::string estimRadiusType = vm["estimRadiusType"].as<std::string>();
  double deltaG = vm["deltaG"].as<double>();
  unsigned int maxComp = vm["maxComp"].as<unsigned int>();
  unsigned int minSizeCC = vm["minSizeCC"].as<unsigned int>();
  bool importNormals = vm.count("importNormals");
  
  NormalAccumulator acc(radius, estimRadiusType);
  DGtal::trace.info() << "------------------------------------ "<< std::endl;
  DGtal::trace.info() << "Step 1: Reading input mesh ... ";
  
  if(vm.count("importPointsNormals"))
    {
      // 1) Reading input mesh:
      // format point A and B from normals
      std::vector<Z3i::RealPoint> vPtA, vPtB;
      std::string namePt = vm["importPointsNormals"].as<std::string>();
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
      std::string inputMeshName = vm["input"].as<std::string>();
      DGtal::Mesh<DGtal::Z3i::RealPoint> aMesh;
      aMesh << inputMeshName;
      DGtal::trace.info() << " [done] " << std::endl;  
      DGtal::trace.info() << "------------------------------------ "<< std::endl;  
      if(importNormals){
        std::vector<DGtal::Z3i::RealPoint> vectorNorm =
          PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(vm["importNormals"].as<std::string>());
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
  if(vm.count("exportConfident")){
  
    std::string confExportName = vm["exportConfident"].as<std::string>();
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
  if(vm.count("initFromMaxRadius")){
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
