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

#include "CLI11.hpp"


#include "AccumulatorHelper.h"
#include "NormalAccumulator.h"


using namespace std;
using namespace DGtal;




typedef ImageContainerBySTLVector<Z3i::Domain, unsigned int> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3DChar;
typedef ImageContainerBySTLVector<Z3i::Domain,  Z3i::RealPoint> ImageVector;


typedef ImageContainerBySTLMap<Z3i::Domain,double> DistanceImage;
typedef DigitalSetFromMap<DistanceImage> AcceptedPointSet;
typedef Z3i::Domain::Predicate DomainPredicate;

typedef Mesh<Z3i::RealPoint>::FaceStorage FaceStorage;

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

  string inputFile;
  string outputFile {"result.off"};
  double minAccumulation {0.0};
  double radius {5.0};
  double cutoff {10000.0};
  std::string estimRadiusType {"mean"};
  bool invertNormal {false};
  bool useConfident {false};
  bool filterFacesNoAcc {false};
  string importNormals {""};
  double maxRadiusFilter {100000.0};

  double minRadiusFilter {0.0};
  bool segmentColor {false};
  double segmentThreshold {0.0};
  
  
  DGtal::Color colorInf = DGtal::Color::Blue;
  DGtal::Color colorSup = DGtal::Color::Yellow;
  DGtal::Color colorNoAcc = DGtal::Color(100, 100, 100, 30);
  
  
  
  // parse command line using CLI ----------------------------------------------
  // parse command line using CLI ----------------------------------------------
  CLI::App app;

  app.description("Display on mesh accumulation values from normal vector ");
  app.add_option("-i,--input,1", inputFile, "input mesh." )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2", outputFile, "output mesh.", true);
  app.add_option("--minAccumulation,-m", minAccumulation, "min accumulation or confidence to consider display.");
  app.add_option("--radius,-r", radius, "maximal radius of accumulation analysis.", true);
  app.add_option("--cutoff", radius, "max cut off display.", true);

  app.add_option("--radiusEstimator,-e", estimRadiusType,  "set the type of the radius estimation (mean (default), min, median or max).", true)
  -> check(CLI::IsMember({"max", "min", "mean", "median"}));
  app.add_option("--importNormals", importNormals,"Use imported normals instead the one computed from the mesh faces." );
  auto maxRadiusFilterOpt = app.add_option("--filterFacesMaxRadius,-R",maxRadiusFilter, "filter faces with radius greater than a threshold.", true);
  auto minRadiusFilterOpt = app.add_option("--filterFacesMinRadius", minRadiusFilter,"filter faces with radius less than a threshold.", true);
  app.add_option("--segmentThreshold", segmentThreshold,"set the segmentation thershold.", true);
  app.add_flag("--invertNormal", invertNormal, "invert input normal.");
  app.add_flag("--displayFromConfidence", useConfident, "flag to display value from confident image insted default accumulation.");
  app.add_flag("--filterFacesNoAcc", filterFacesNoAcc, "filter faces with no accumulation." );
  app.add_flag("--segmentColor", segmentColor, "filter and display mesh faces in color (to be used with --segmentTh).");
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

    
   
  // 1) Reading input mesh
  Mesh<Z3i::RealPoint> aMesh(true);
  aMesh << inputFile;

  // 2) Generate accumulation image
  std::pair<Z3i::RealPoint, Z3i::RealPoint> bb = aMesh.getBoundingBox();

  ImageVector imageVect (Z3i::Domain(bb.first, bb.second));
  trace.info()<< "Computing accumulation image .... on domain: " << bb;

  NormalAccumulator normAcc (radius, estimRadiusType);
  if(importNormals != ""){
    std::vector<DGtal::Z3i::RealPoint> vectorNorm =
      PointListReader<DGtal::Z3i::RealPoint>::getPointsFromFile(importNormals);
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
  HueShadeColorMap<double > map(0.0, std::min(cutoff, radius));
  std::vector<Mesh<Z3i::RealPoint>::Index> faceToBeRemoved;
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
    }else  if((val==0 && filterFacesNoAcc)||
              (*minRadiusFilterOpt && val  < minRadiusFilter)||
              (*maxRadiusFilterOpt && val > maxRadiusFilter)){
      faceToBeRemoved.push_back(i);
    }
    if(!segmentColor){      
      
      aMesh.setFaceColor(i, map(std::min(cutoff,val)));
    }
  }
  
  DGtal::trace.info() << "[done]" << std::endl;
  trace.info()<< "Val max dist accumulation: " << maxVal << std::endl;
  if(!segmentColor && (filterFacesNoAcc || *minRadiusFilterOpt
                       || *maxRadiusFilterOpt) ){
    aMesh.removeFaces(faceToBeRemoved);
  }
  
  aMesh >> outputFile  ;
  return 0;
}
