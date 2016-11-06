#include <iostream>
#include <stdio.h>
#include <stdlib.h>


#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/readers/PointListReader.h>
#include <DGtal/io/writers/MeshWriter.h>
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageContainerBySTLMap.h>
#include <DGtal/kernel/sets/DigitalSetFromMap.h>


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>


#include "AccumulatorHelper.h"
#include "NormalAccumulator.h"




using namespace std;
using namespace DGtal;
namespace po = boost::program_options;




typedef ImageContainerBySTLVector<Z3i::Domain, unsigned int> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef ImageContainerBySTLVector<Z3i::Domain,  Z3i::RealPoint> ImageVector;





typedef ImageContainerBySTLMap<Z3i::Domain,double> DistanceImage;
typedef DigitalSetFromMap<DistanceImage> AcceptedPointSet;
typedef Z3i::Domain::Predicate DomainPredicate;



template<typename TImage>
TImage
applyMedian(const TImage &anImage, unsigned int size){
  unsigned int imageSize = (anImage.domain().upperBound()[0]-anImage.domain().lowerBound()[0])*
  (anImage.domain().upperBound()[1]-anImage.domain().lowerBound()[1])*
  (anImage.domain().upperBound()[2]-anImage.domain().lowerBound()[2]);
  
  trace.progressBar(0, imageSize);
  typedef typename TImage::Domain::ConstIterator ImageDomIterator;
  typedef   typename Z3i::Point Point3D;
  TImage imageRes (anImage.domain());
  unsigned int cpt = 0;
  
  for(ImageDomIterator it = anImage.domain().begin(); it != anImage.domain().end(); ++it){
    trace.progressBar(cpt, imageSize);
    Point3D pt = *it;
    std::vector<unsigned char> vectVal;
    for (int k = -(int)size; k <= (int) size; k++) {
      for (int l = -(int)size; l <= (int) size; l++) {
        for (int m = -(int)size; m <= (int) size; m++) {
          Point3D p (((int)(pt[0]))+k, (int)(pt[1])+l, ((int)(pt[2]))+m);
          if(anImage.domain().isInside(p) && (p-pt).norm() < size ){
            vectVal.push_back(anImage(p));
          }
        }
      }
    }
    std::sort(vectVal.begin(), vectVal.end());
    imageRes.setValue(pt, vectVal[vectVal.size()/2.0]);
    cpt++;
  }
  return imageRes;
}






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
    ("output,o", po::value<std::string>(), "output mesh.")
    ("minAccumulation,m", po::value<double>()->default_value(0.0), "min accumulation or confidence to consider display.")
    ("cutoff", po::value<double>()->default_value(10000.0), "max cut off display.")
    ("estimRadiusType", po::value<std::string>()->default_value("mean"), "set the type of the radius estimation (mean, min, median or max).")
    ("displayFromConfidence", "flag to display value from confident image insted default accumulation.")
    ("invertNormal", "invert input normal.")
    ("filterFacesNoAcc", "filter faces with no accumulation.")
    ("importNormals", po::value<std::string>(), "Use imported normals instead the one computed from the mesh faces.")
    ("filterFacesMaxRadius,R", po::value<double>()->default_value(100000.0), "filter faces with radius greater than a threshold.")
    ("filterFacesMinRadius", po::value<double>()->default_value(0.0), "filter faces with radius less than a threshold.")
    ("segmentColor", "filter and display mesh faces in color (to be used with --segmentTh) .")
    ("segmentThreshold", po::value<double>()->default_value(0.0), "set the segmentation thershold.")
    ("radius,r", po::value<double>()->default_value(5), "maximal radius of accumulation analysis.");
  
  
  bool parseOK = true;
  po::variables_map vm;
  
  
  try {
    po::store(po::parse_command_line(argc, argv, general_opt), vm);
  } catch (const std::exception &ex) {
    trace.info() << "Error checking program options: " << ex.what() << std::endl;
    parseOK = false;
  }

  po::notify(vm);

  if (vm.count("help") || argc <= 1 || !parseOK || !vm.count("input")) {
    trace.info() << "Display on mesh accumulation values from normal vector" << std::endl <<
                 "Options: " << std::endl
                 << general_opt << "\n";
    return 0;
  }
  
  double radius = vm["radius"].as<double>();
  string inputFile = vm["input"].as<std::string>();
  string outputFile = vm["output"].as<std::string>();
  bool invertNormal = vm.count("invertNormal");
  std::string estimRadiusType = vm["estimRadiusType"].as<std::string>();
  bool useConfident = vm.count("displayFromConfidence");
  double minAccumulation = vm["minAccumulation"].as<double>();
  double maxRadiusFilter = vm["filterFacesMaxRadius"].as<double>();
  double minRadiusFilter = vm["filterFacesMinRadius"].as<double>();
  bool segmentColor = vm.count("segmentColor");
  double segmentThreshold = vm["segmentThreshold"].as<double>();
  double cutoff = vm["cutoff"].as<double>();
  bool importNormals = vm.count("importNormals");
  
  
  DGtal::Color colorInf = DGtal::Color::Blue;
  DGtal::Color colorSup = DGtal::Color::Yellow;
  DGtal::Color colorNoAcc = DGtal::Color(100, 100, 100, 30);
  
  // 1) Reading input mesh
  Mesh<Z3i::RealPoint> aMesh(true);
  aMesh << inputFile;

  // 2) Generate accumulation image
  std::pair<Z3i::RealPoint, Z3i::RealPoint> bb = aMesh.getBoundingBox();

  ImageVector imageVect (Z3i::Domain(bb.first, bb.second));
  trace.info()<< "Computing accumulation image .... on domain: " << bb;

  NormalAccumulator normAcc (radius, estimRadiusType);
  if(importNormals){
    std::vector<DGtal::Z3i::RealPoint> vectorNorm =
      PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(vm["importNormals"].as<std::string>());
    normAcc.initFromMeshAndNormals(aMesh, vectorNorm, invertNormal);
    
  }else{
    normAcc.initFromMesh(aMesh, invertNormal);
  }
  normAcc.computeAccumulation();
  normAcc.computeConfidence();
  trace.info() << "[done]" << std::endl;
  NormalAccumulator::Image3DDouble &imageConfidence = normAcc.getConfidenceImage();
  NormalAccumulator::Image3D &imageAccumulation = normAcc.getAccumulationImage();
  
  
  //imageAccu = applyMedian(imageAccu, 1);
  
  // 3) Compute for each face its maximal accumulation value and change the output mesh colors from its value.
  // Get the max value:
  
  double maxVal = 0;
  HueShadeColorMap<double > map(0.0, radius);
  std::vector<unsigned int> faceToBeRemoved;
  if (useConfident){
      DGtal::trace.info() << "compute face accumulation from Confident image...";
  }else{
    DGtal::trace.info() << "compute face accumulation from Accumumation image...";
  }
  for(unsigned int i = 0; i< aMesh.nbFaces(); i++){
    trace.progressBar(i, aMesh.nbFaces());
    double val = 0;
    if (useConfident){
      val = AccumulatorHelper::getFaceMaxAccDist(imageConfidence, aMesh, i, radius, minAccumulation, invertNormal );
    }else{
      val = AccumulatorHelper::getFaceMaxAccDist(imageAccumulation, aMesh, i, radius, minAccumulation, invertNormal );
    }
    if(val> maxVal){
            maxVal = val;
    }
    if(segmentColor){
      if (val == 0){
        aMesh.setFaceColor(i, colorNoAcc );
      }
      else if(val < segmentThreshold){
        aMesh.setFaceColor(i, colorInf );
      }
      else  {
        aMesh.setFaceColor(i, colorSup );
      }
    }else  if((val==0 && vm.count("filterFacesNoAcc"))||
              (vm.count("filterFacesMinRadius") && val  < minRadiusFilter)||
              (vm.count("filterFacesMaxRadius") && val > maxRadiusFilter)){
      faceToBeRemoved.push_back(i);
    }
    if(!segmentColor){      
      
      aMesh.setFaceColor(i, map(std::min(cutoff,val)));
    }
  }
  
  DGtal::trace.info() << "[done]" << std::endl;
  trace.info()<< "Val max dist accumulation: " << maxVal << std::endl;
  if(!segmentColor && (vm.count("filterFacesNoAcc") ||vm.count("filterFacesMinRadius")
                       || vm.count("filterFacesMaxRadius")) ){
    aMesh.removeFaces(faceToBeRemoved);
  }
  
  aMesh >> outputFile  ;
  return 0;
}











