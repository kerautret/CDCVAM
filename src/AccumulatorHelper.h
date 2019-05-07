#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/io/boards/Board2D.h"
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/kernel/BasicPointFunctors.h>
#include <DGtal/shapes/implicit/ImplicitBall.h>
#include "DGtal/geometry/curves/AlphaThickSegmentComputer.h"
#include <DGtal/shapes/EuclideanShapesDecorator.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/math/Statistic.h>


class AccumulatorHelper {


public:  


  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned int> Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> Image3DChar;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  DGtal::Z3i::RealPoint> ImageVector;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  std::vector<DGtal::Z3i::RealPoint> > ImagePointAssociation;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  std::vector<unsigned int> > ImageFaceAssociation;
  typedef DGtal::ConstImageAdapter<Image3D, DGtal::Z2i::Domain, DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>,
                                   unsigned int,  DGtal::functors::Identity >  ImageAdapterExtractor;



  /**
   *
   * Orient Normal vectors from max accumulation
   *
   **/
  template<typename TPoint>
  static DGtal::Z3i::Domain
  orientNormalsFromAcc(const std::vector<TPoint> &vPoints,
                       std::vector<TPoint> &vNormals, double maxRadius = 5.0)
  {
    // get bounding box of set of points
    TPoint pl, pu;
    for(auto &pt: vPoints)
      {
        if(pt[0] < pl[0]) pl[0] = pt[0];
        if(pt[1] < pl[1]) pl[1] = pt[1];
        if(pt[2] < pl[2]) pl[2] = pt[2];

        if(pt[0] > pu[0]) pu[0] = pt[0];
        if(pt[1] > pu[1]) pu[1] = pt[1];
        if(pt[2] > pu[2]) pu[2] = pt[2];
      }
    DGtal::Z3i::RealPoint p (maxRadius,maxRadius,maxRadius);
    Image3D imgAcc(Image3D::Domain(pl-p, pu+p));
    DGtal::trace.info() << imgAcc.domain() << std::endl;
    // compute accumulation in the two directions
    for(unsigned int i = 0; i < vPoints.size(); i++)
      {
        TPoint scanDir = vNormals[i].getNormalized();
        DGtal::Z3i::RealPoint centerPoint = vPoints[i] ;
        DGtal::Z3i::RealPoint currentPoint = centerPoint;
        DGtal::Z3i::RealPoint previousPoint;       
        DGtal::Z3i::RealPoint currentPointN = centerPoint;
        DGtal::Z3i::RealPoint previousPointN;       
        while((currentPoint - centerPoint).norm()<maxRadius){
          DGtal::Z3i::Point currentPointI = DGtal::Z3i::Point(currentPoint,DGtal::functors::Round<>());
          if(imgAcc.domain().isInside(currentPointI)
             && previousPoint != currentPoint){
            imgAcc.setValue(currentPointI, imgAcc(currentPointI)+1);
            previousPoint = currentPoint;
          }
          DGtal::Z3i::Point c = DGtal::Z3i::Point(currentPointN,DGtal::functors::Round<>());
          if(imgAcc.domain().isInside(c) &&
             previousPointN != currentPointN){
            imgAcc.setValue(c, imgAcc(c)+1);
            previousPointN = currentPointN;
          }
          previousPoint = currentPoint;
          currentPoint += scanDir;          
          previousPointN = currentPointN;
          currentPointN -= scanDir;
        }
      }  


    for(unsigned int i = 0; i < vPoints.size(); i++)
      {
        unsigned int somme = 0, sommeN = 0;
        
        TPoint scanDir = vNormals[i].getNormalized();
        DGtal::Z3i::RealPoint centerPoint = vPoints[i] ;
        DGtal::Z3i::RealPoint currentPoint = centerPoint;
        DGtal::Z3i::RealPoint previousPoint;       
        DGtal::Z3i::RealPoint currentPointN = centerPoint;
        DGtal::Z3i::RealPoint previousPointN;       
        while((currentPoint - centerPoint).norm()<maxRadius){
          DGtal::Z3i::Point currentPointI = DGtal::Z3i::Point(currentPoint,DGtal::functors::Round<>());
          if(imgAcc.domain().isInside(currentPointI) && previousPoint != currentPoint){
            somme+=imgAcc(currentPointI);
            previousPoint = currentPoint;
          }
          DGtal::Z3i::Point c = DGtal::Z3i::Point(currentPointN,DGtal::functors::Round<>());
          if(imgAcc.domain().isInside(c) && previousPointN != currentPointN){
            sommeN+=imgAcc(c);
            previousPointN = currentPointN;
          }
          previousPoint = currentPoint;
          currentPoint += scanDir;          
          previousPointN = currentPointN;
          currentPointN -= scanDir;
        }

        if(somme < sommeN){
          vNormals[i] = -vNormals[i];
        }
      }

    // Direct comparisons of acc to orient perhaps less robust (not sure)
    // for(unsigned int i = 0; i < vPoints.size(); i++)
    //     {
    //       TPoint p = vPoints[i];
    //       TPoint n = vNormals[i].getNormalized();
    //       TPoint pp = p+n*maxRadius;
    //       TPoint pm = p-n*maxRadius;
    //       if(imgAcc.domain().isInside(pm) && imgAcc.domain().isInside(pp) && 
    //          imgAcc(pp) < imgAcc(pm)){        
    //         vNormals[i] = -vNormals[i];
    //       }
    //     }

    return Image3D::Domain(pl-p, pu+p);
  }
  


  /**
   * From an a face id of a mesh and from its associated accumulation
   * image, return the maximal score found on the ray starting from
   * mesh center until le max radius.
   * 
   **/
  
  template<typename TImageAcc, typename TPoint>
  static double
  getFaceMaxAccDist(const TImageAcc &imageAcc, const DGtal::Mesh<TPoint> &aMesh, 
                    const unsigned int faceId, const double maxRadius, const double thresholdAcc,
                    bool invertNormal = false){

    typedef typename DGtal::Mesh<TPoint>::MeshFace Face;
    double resValue = 0.0;
    
    Face aFace = aMesh.getFace(faceId);
    
    TPoint p0 = aMesh.getVertex(aFace.at(0));
    TPoint p1 = aMesh.getVertex(aFace.at(2));
    TPoint p2 = aMesh.getVertex(aFace.at(1));
    TPoint scanDir = ((p1-p0).crossProduct(p2 - p0)).getNormalized();
    if (invertNormal){
      scanDir *= -1;
    }
    
    DGtal::Z3i::RealPoint centerPoint = (p0+p1+p2)/3.0;
    DGtal::Z3i::RealPoint currentPoint = centerPoint;
    typename TImageAcc::Value maxVal = 0;
    while((currentPoint - centerPoint).norm()<maxRadius){
      DGtal::Z3i::Point currentPointI = DGtal::Z3i::Point(currentPoint,DGtal::functors::Round<>());
              
      if(imageAcc.domain().isInside(currentPointI)){
        typename TImageAcc::Value v = imageAcc(currentPointI);
        if(v>maxVal && v> thresholdAcc){
          maxVal = v;
          resValue = (currentPoint-centerPoint).norm();
        }
      }
      currentPoint += scanDir;
    }
    return resValue;
  }




  
};





