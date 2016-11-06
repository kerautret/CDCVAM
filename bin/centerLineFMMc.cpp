#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
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

#include "NormalAccumulator.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


typedef ImageContainerBySTLVector<Z3i::Domain, bool> Image3DMaker;
typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, double> Image3DDouble;
typedef ImageContainerBySTLVector<Z3i::Domain,double> DistanceImage;

typedef DigitalSetBySTLSet<Z3i::Domain> AcceptedPointSet;
typedef Z3i::Domain::Predicate DomainPredicate;
typedef L2SecondOrderLocalDistance<DistanceImage,AcceptedPointSet> DistanceMeasure;
typedef FMM<DistanceImage, AcceptedPointSet, DomainPredicate, DistanceMeasure> TFMM;




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
    ("output,o", po::value<std::string>(), "the output base name.")
    ("distanceSearch,d", po::value<double>(), "distance search for tracking next (to determine next ATP).")
    ("minRadius,m", po::value<double>(), "min radius to track center line.")
    ("maxRadius,M", po::value<double>(), "max radius to track center line.")
    ("distanceMerge,r", po::value<double>(), "distance to consider merged points in seg ATP.")
    ("invertNormal,n", "invert normal to apply accumulation.")
    ("reconsWithRadiusSort", "graph reconstruction by using decreasing radius size else it use confidence.")
    ("initFromMaxRadius", "init from max radius instead max confidence.")
    ("th,t", po::value<double>(), "threshold in the confidence estimation.")
    ("estimRadiusType", po::value<std::string>()->default_value("mean"), "set the type of the radius estimation (mean, min, median or max).")
    ("radius,R", po::value<double>(), "radius used to compute the accumulation.");
  

  
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
    trace.info() << "Center line extraction from accumulation volumetric image" << std::endl
                 << "Options: " << std::endl
                 << general_opt << std::endl;
    return 0;
  }


  // Reading parameters:
  std::string inputMeshName = vm["input"].as<std::string>();
  std::string outputName = vm["output"].as<std::string>();
  double distanceSearch = vm["distanceSearch"].as<double>();
  double distanceMerge = vm["distanceMerge"].as<double>();
  double th = vm["th"].as<double>();
  double radius = vm["radius"].as<double>();
  double minRadius = 0;
  double maxRadius = std::numeric_limits<double>::max();
  bool invertNormal = vm.count("invertNormal");
  std::string estimRadiusType = vm["estimRadiusType"].as<std::string>();
  if(vm.count("minRadius")){
    minRadius = vm["minRadius"].as<double>();
  }

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
  acc.computeConfidence();
  Image3DDouble imageConfidence = acc.getConfidenceImage();
  Image3DMaker imageVisited (imageConfidence.domain());
  for (auto const &p: imageVisited.domain()){ imageVisited.setValue(p, false);}

  DGtal::trace.info() << "Exporting ... ";
  typedef DGtal::functors::Rescaling< double, unsigned char> ScaleFctD;
  ScaleFctD  confidencescale (0 , 1.0, 0, 255);
  DGtal::VolWriter<NormalAccumulator::Image3DDouble,ScaleFctD>::exportVol("confidence.vol",
                                                                          imageConfidence,
                                                                          confidencescale);
  DGtal::trace.info() << " [done] " << std::endl;  
  DGtal::trace.info() << "------------------------------------ "<< std::endl;  
  

  acc.computeRadiusFromConfidence();
  Image3DDouble imageRadius = acc.getRadiusImage();

  // 4) Starting to track from maximal radius point 
  DGtal::trace.info() << "Step 4: Starting track ... " << std::endl;
  
  DGtal::Z3i::Point startPoint; 
  if(vm.count("initFromMaxRadius")){
    startPoint = acc.getMaxRadiusPoint(); 
  }else {
    startPoint = acc.getMaxAccumulationPoint();   
  }
  DGtal::trace.info() << "Start point:" << startPoint << std::endl;
  
  // Tracking Steps: init with maximal accumulation point:
  std::deque<DGtal::Z3i::Point> vectATP;
  std::vector<DGtal::Z3i::Point> centerLineVertex;
  std::vector< std::pair<int,int> > centerLineEdges;
  std::map<DGtal::Z3i::Point, unsigned int> vertexLabels;
  
 
  vectATP.push_back(startPoint);
  vertexLabels[startPoint] = 0;
  imageVisited.setValue(startPoint,true);
  centerLineVertex.push_back(startPoint);   
  DGtal::trace.info() << "tracking in progress: [ ";
  // Tracking Steps: main loop 
  // => while exits ATP
  unsigned int nb = 0;

  
  while(vectATP.size() != 0){
    nb++;
    DGtal::trace.info() << ".";
    DGtal::Z3i::Point p = vectATP.front(); vectATP.pop_front();
    imageVisited.setValue(p,true);
    // Tracking Steps: looking for extension of actual tracking point
    // => apply FFM
    DistanceImage distanceImage( imageConfidence.domain() );
    AcceptedPointSet set( imageConfidence.domain() );
    set.insert(p);
    TFMM fmm( distanceImage, set, imageConfidence.domain().predicate(), 
              imageConfidence.domain().size(), distanceSearch );
    double locMax = 0.0;
    double lastDist = 0.0;
    DGtal::Z3i::Point lastPt;
    std::vector<DGtal::Z3i::Point> listCandidates;
    while ( fmm.computeOneStep( lastPt, lastDist ) )
      {
        if ( imageConfidence(lastPt) > th &&
             imageRadius(lastPt) > minRadius &&
             imageRadius(lastPt) < maxRadius &&(lastPt-p).norm()>distanceMerge  && !imageVisited(lastPt)) {
          listCandidates.push_back(lastPt);
        }
         imageVisited.setValue(lastPt, true);
      }
    // sorting candidates from confidence
    if(vm.count("reconsWithRadiusSort")){
      //sorting candidates from Radius:
            std::sort(listCandidates.begin(), listCandidates.end(), 
                      [imageRadius](const DGtal::Z3i::Point &a,
                                    const DGtal::Z3i::Point &b) -> bool 
                      {
                        return imageRadius(a)<imageRadius(b);
                      });
    }else{
      std::sort(listCandidates.begin(), listCandidates.end(),
                [imageConfidence](const DGtal::Z3i::Point &a,
                                  const DGtal::Z3i::Point &b) -> bool 
                {
                return imageConfidence(a)<imageConfidence(b);
                });
      
      
    }
    
    // Take candidates which are not included in another newATP 
    std::vector<DGtal::Z3i::Point> newATP;
    while( ! listCandidates.empty()){
      DGtal::Z3i::Point anATP = listCandidates.back();
      listCandidates.pop_back();
      // check if new candidates  are not included in the neightboorhood of an existing ATP:
      bool isNotIncluded = true;
      for(auto const &pp: newATP){
        isNotIncluded = isNotIncluded && (anATP-pp).norm() > distanceMerge;
      }
      for(auto const &pp: centerLineVertex){
        isNotIncluded = isNotIncluded && (anATP-pp).norm() > distanceMerge;
      }
      if(isNotIncluded){
        if(!imageVisited(anATP)){
          trace.warning() << "Not visited..." << anATP;
        }
        newATP.push_back(anATP);
      }
    }
    unsigned int pos = 0;
    for (auto const & pNew: newATP){
      if(vertexLabels.count(pNew) == 0){
        vertexLabels[pNew]=centerLineVertex.size();
        centerLineVertex.push_back(pNew);
      }
      centerLineEdges.push_back( std::pair<int,int>(vertexLabels[p], 
                                                    vertexLabels[pNew] ) );
      vectATP.push_back(pNew);
      pos++;
    }
  }
  DGtal::trace.info() << "[done]" << std::endl;
  

  // export vertex
  stringstream ss;
  ss<<outputName << "Vertex.sdp";
  std::ofstream outVertex;
  outVertex.open(ss.str().c_str());
  std::sort(centerLineVertex.begin(), centerLineVertex.end(), 
            [vertexLabels](const DGtal::Z3i::Point &a, const DGtal::Z3i::Point &b) -> bool {
              auto valA = vertexLabels.find(a)->second;
              auto valB = vertexLabels.find(b)->second;
              return valA < valB;} ); 
  for (unsigned int i = 0; i< centerLineVertex.size(); i++)
    {
      auto v = centerLineVertex[i];
      outVertex << v[0] << "  " << v[1] << " " << v[2] << " " << vertexLabels[v] << std::endl;
    }
  outVertex.close();

  // export radius

  stringstream ss3;
  ss3<<outputName << "Radius.sdp";
  std::ofstream outRadius;
  outRadius.open(ss3.str().c_str());
  for (const auto& e: centerLineVertex)
   {
     outRadius << std::max(imageRadius(e), 2.0) << std::endl;
   }
  outRadius.close();
  
  

  // export edges
  stringstream ss2;
  ss2<<outputName << "Edges.sdp";
  std::ofstream outEdges;
  outEdges.open(ss2.str().c_str());
  
  for (const auto& e: centerLineEdges)
   {
     outEdges << e.first << "  " << e.second << std::endl;
   }
  outEdges.close();
  
 

  return 0;
}
