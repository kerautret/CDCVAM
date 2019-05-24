#include <DGtal/math/Statistic.h>

#include "NormalAccumulator.h"
#include "AccumulatorHelper.h"

#if defined USE_PCL
#include <pcl/features/integral_image_normal.h>
#include <pcl/features/normal_3d.h>
#include <pcl/point_types.h>
#endif


using namespace std;






///////////////////////////////////////////////////////////////////////////////
// class NormalAccumulator
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////
//  Init from Mesh

void
NormalAccumulator::initFromMesh(const DGtal::Mesh<Point> &aMesh, bool invertNormal)
{  
  typedef typename DGtal::Mesh<Point>::MeshFace Face;
  // compute the normal field and vertex center
  if (myIsVerbose)
    { 
      DGtal::trace.info() << "starting init from mesh..." ;
    }
  for (unsigned int iFace = 0; iFace< aMesh.nbFaces(); iFace++){
     Face aFace = aMesh.getFace(iFace);
     if(aFace.size()>4){
       DGtal::trace.warning() << "ignoring face, not a triangular one: " << aFace.size() << " vertex" << std::endl;
       continue;
     }
     
     Point p0 = aMesh.getVertex(aFace.at(0));
     Point p1 = aMesh.getVertex(aFace.at(2));
     Point p2 = aMesh.getVertex(aFace.at(1));
     
     Vector n = ((p1-p0).crossProduct(p2 - p0)).getNormalized();
     if (invertNormal) {n.negate();}
     myNormalField.push_back(n);
     myNormalOrigins.push_back( (p0+p1+p2)/3.0 );
  }
  
  std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> bb = aMesh.getBoundingBox();
  myDomain = DGtal::Z3i::Domain(bb.first - DGtal::Z3i::RealPoint::diagonal(3*myRadius), 
                                bb.second + DGtal::Z3i::RealPoint::diagonal(3*myRadius));
  myAccumulationImage = NormalAccumulator::Image3D(myDomain);
  myConfidenceImage = NormalAccumulator::Image3DDouble(myDomain);
  myRadiusImage = NormalAccumulator::Image3DDouble(myDomain);
  myAssociationImage = NormalAccumulator::ImagePointAssociation(myDomain);
  if (myIsVerbose)
    { 
      DGtal::trace.info() << " [done]"<< std::endl 
                          << *this << std::endl;
    }
  myIsAccumulationComputed = false;
  myIsRadiusComputed = false;
  myIsConfidenceComputed = false;
  myIsAssociationCompFromConfidence = false;
  myMaxAccumulation = 0;
  myMaxRadius = 0;

  myIsInitialized = true;
}



void
NormalAccumulator::initFromMeshAndNormals(const DGtal::Mesh<Point> &aMesh,
                                          const std::vector<DGtal::Z3i::RealPoint> &vectNormals,
                                          bool invertNormal)
{  
  typedef typename DGtal::Mesh<Point>::MeshFace Face;
  assert(vectNormals.size()==aMesh.nbFaces());
  
  // compute the normal field and vertex center
  if (myIsVerbose)
    { 
      DGtal::trace.info() << "starting init from mesh..." ;
    }
  for (unsigned int iFace = 0; iFace< aMesh.nbFaces(); iFace++){
     Face aFace = aMesh.getFace(iFace);
     if(aFace.size()>4){
       DGtal::trace.warning() << "ignoring face, not a triangular one: " << aFace.size() << " vertex" << std::endl;
       continue;
     }
     
     Point p0 = aMesh.getVertex(aFace.at(0));
     Point p1 = aMesh.getVertex(aFace.at(2));
     Point p2 = aMesh.getVertex(aFace.at(1));
     
     Vector n = vectNormals[iFace];
     if (invertNormal) {n.negate();}
     myNormalField.push_back(n);
     myNormalOrigins.push_back( (p0+p1+p2)/3.0 );
  }
  
  std::pair<DGtal::Z3i::RealPoint, DGtal::Z3i::RealPoint> bb = aMesh.getBoundingBox();
  myDomain = DGtal::Z3i::Domain(bb.first - DGtal::Z3i::RealPoint::diagonal(3*myRadius), 
                                bb.second + DGtal::Z3i::RealPoint::diagonal(3*myRadius));
  myAccumulationImage = NormalAccumulator::Image3D(myDomain);
  myConfidenceImage = NormalAccumulator::Image3DDouble(myDomain);
  myRadiusImage = NormalAccumulator::Image3DDouble(myDomain);
  myAssociationImage = NormalAccumulator::ImagePointAssociation(myDomain);
  if (myIsVerbose)
    { 
      DGtal::trace.info() << " [done]"<< std::endl 
                          << *this << std::endl;
    }
  myIsAccumulationComputed = false;
  myIsRadiusComputed = false;
  myIsConfidenceComputed = false;
  myIsAssociationCompFromConfidence = false;
  myMaxAccumulation = 0;
  myMaxRadius = 0;

  myIsInitialized = true;
}






void
NormalAccumulator::initFromNormals(const std::vector<DGtal::Z3i::RealPoint> &aSetOfPoints,
                                   const std::vector<DGtal::Z3i::RealPoint> &aSetOfNormals,
                                   const bool invertNormal)
{
  // Simply insert input points and normals.
  if (myIsVerbose)
    { 
      DGtal::trace.info() << "starting init from points and normals..." ;
    }
  assert(aSetOfNormals.size()==aSetOfNormals.size());
  DGtal::Z3i::RealPoint lowPt, upPt;
  lowPt = aSetOfPoints[0];
  upPt = aSetOfPoints[1];
  for (unsigned int i = 0; i  < aSetOfPoints.size(); i++)
    {
      DGtal::Z3i::RealPoint n = aSetOfNormals[i];
       DGtal::Z3i::RealPoint p = aSetOfPoints[i];
      if (invertNormal) {n.negate();}
      myNormalField.push_back(n);
      myNormalOrigins.push_back( p );
      lowPt = lowPt.inf(p);
      upPt = upPt.sup(p);
    }
  
  
  myDomain = DGtal::Z3i::Domain(lowPt - DGtal::Z3i::RealPoint::diagonal(3*myRadius), 
                                upPt + DGtal::Z3i::RealPoint::diagonal(3*myRadius));
  myAccumulationImage = NormalAccumulator::Image3D(myDomain);
  myConfidenceImage = NormalAccumulator::Image3DDouble(myDomain);
  myRadiusImage = NormalAccumulator::Image3DDouble(myDomain);
  myAssociationImage = NormalAccumulator::ImagePointAssociation(myDomain);
  if (myIsVerbose)
    { 
      DGtal::trace.info() << " [done]"<< std::endl 
                          << *this << std::endl;
    }
  myIsAccumulationComputed = false;
  myIsRadiusComputed = false;
  myIsConfidenceComputed = false;
  myIsAssociationCompFromConfidence = false;
  myMaxAccumulation = 0;
  myMaxRadius = 0;
  myIsInitialized = true;

}





///////////////////////////////////////////
//  Init from Point cloud

#if defined USE_PCL

void 
NormalAccumulator::initFromPointCloud(const std::vector<Point> &aVectPoints,
                                      const double size, const bool autoOrient)
{
  myNormalOrigins = aVectPoints;
  // prepare cloud points
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);
  for (auto &pt: aVectPoints){
    cloud->push_back(pcl::PointXYZ(pt[0], pt[1], pt[2]));
  }
  // compute normal from PCL:
  pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
  ne.setInputCloud (cloud);
  pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ> ());
  ne.setSearchMethod (tree);
  pcl::PointCloud<pcl::Normal>::Ptr cloud_normals (new pcl::PointCloud<pcl::Normal>);
  ne.setRadiusSearch(size);
  ne.compute (*cloud_normals);
  
  // import normal:
  for(unsigned int i = 0; i < (*cloud_normals).size(); i++ ){
    myNormalField.push_back(DGtal::Z3i::RealPoint((*cloud_normals)[i].normal_x,
                                           (*cloud_normals)[i].normal_y,
                                           (*cloud_normals)[i].normal_z));
  }  
  if(autoOrient){
    myDomain = AccumulatorHelper::orientNormalsFromAcc(myNormalOrigins, myNormalField, myRadius);
  }
  myAccumulationImage = NormalAccumulator::Image3D(myDomain);
  myConfidenceImage = NormalAccumulator::Image3DDouble(myDomain);
  myRadiusImage = NormalAccumulator::Image3DDouble(myDomain);
  myAssociationImage = NormalAccumulator::ImagePointAssociation(myDomain);

  
  myIsAccumulationComputed = false;
  myIsRadiusComputed = false;
  myIsConfidenceComputed = false;
  myIsAssociationCompFromConfidence = false;
  myMaxAccumulation = 0;
  myMaxRadius = 0;

  myIsInitialized = true;
  
}


void 
NormalAccumulator::initFromPointCloud(const std::vector<Point> &aVectPoints, 
                                      const unsigned int nbNeighbohr)
{


}


#endif






void
NormalAccumulator::computeAccumulation(bool retainVertexAsso, bool verbose)
{
  assert(myIsInitialized);
  if(verbose)
    DGtal::trace.progressBar(0, myNormalField.size());

  for (unsigned int numNorm = 0; numNorm< myNormalField.size(); numNorm++){
    if(verbose)
      DGtal::trace.progressBar(numNorm, myNormalField.size());
    NormalAccumulator::Vector aNormal  = myNormalField[numNorm];    
    DGtal::Z3i::RealPoint centerPoint = myNormalOrigins[numNorm];
    DGtal::Z3i::RealPoint currentPoint = centerPoint;
    DGtal::Z3i::RealPoint previousPoint;
    
    while((currentPoint - centerPoint).norm()<myRadius){
      DGtal::Z3i::Point currentPointI = DGtal::Z3i::Point(currentPoint,DGtal::functors::Round<>());

      if(myDomain.isInside(currentPointI) && previousPoint != currentPoint){
        // retain for each voxel its associated faces contributing to the accumulation
        if(retainVertexAsso){
          std::vector<unsigned int> v = myAssociationImage(currentPointI);
          v.push_back(numNorm);        
          myAssociationImage.setValue(currentPointI, v);
        }
        myAccumulationImage.setValue(currentPointI, myAccumulationImage(currentPointI)+1);
        previousPoint = currentPoint;
        // update the max of accumulation
        if( myAccumulationImage(currentPointI)>myMaxAccumulation ) {
          myMaxAccumulation = myAccumulationImage(currentPointI);
          myMaxAccumulationPoint = currentPointI;
        }
      }
      previousPoint = currentPoint;
      currentPoint += aNormal;
    }
  }
  if(verbose)
    {
      DGtal::trace.progressBar(myNormalField.size(), myNormalField.size());
      DGtal::trace.info() << std::endl << "Max accumulation  value: "<< myMaxAccumulation << std::endl;
      DGtal::trace.info()  << "Max accumulation Point: "<< myMaxAccumulationPoint << std::endl;
    }
  myIsAccumulationComputed = true;
}





void
NormalAccumulator::computeRadiusFromOrigins()
{
  // Step 1: ensure that accumulation is computed.
  assert(myIsAccumulationComputed);

  // Step 2:  Recover for each voxel the list of point for which the normal contribute and compute radius from stats
  myMaxRadius = 0;
  for(auto &v: myDomain){
    myRadiusImage.setValue(v, 0);
  }
  double r;
  for(auto &v: myDomain){
    if(myAssociationImage(v).size() >= 2){
      DGtal::Statistic<double> stat(true);
      for(unsigned int ptId :myAssociationImage(v)){
        stat.addValue((myNormalOrigins[ptId] -v).norm());
      }
      stat.terminate();
      switch (myRadiusStatEstim)
        {
        case StatRadiusDef::min:
          r = stat.min();
          break;
        case StatRadiusDef::mean:
          r = stat.mean();  
          break;
        case  StatRadiusDef::max:
          r = stat.max();  
          break;
        case  StatRadiusDef::median:
          r = stat.median();  
          break;
        }
      
      myRadiusImage.setValue(v, r);
      if ( r > myMaxRadius ) 
        {
          myMaxRadiusPoint = v; 
          myMaxRadius = r;
        }
    }
  }
  myIsRadiusComputed = true;
}


void
NormalAccumulator::computeRadiusFromConfidence()
{
  if(!myIsAccumulationComputed){
    computeAccumulation();
  }
  if(!myIsAssociationCompFromConfidence){
    computeConfidence(true);
  }
  computeRadiusFromOrigins();
}




NormalAccumulator::PointContainer
NormalAccumulator::getAssociatedPoints(const DGtal::Z3i::Point &aVoxel)
{
  PointContainer result; 
  for(unsigned int ptId :myAssociationImage(aVoxel)){
    result.push_back(myNormalOrigins[ptId]);
  }
  return result;
}




void
NormalAccumulator::computeConfidence(bool updateVertexAsso, unsigned int minAcc )
{
  // Step 1: ensure that accumulation is computed and clean association image if needed.
  assert(myIsAccumulationComputed);
  if (updateVertexAsso){
    for (auto &p: myDomain){
      myAssociationImage.setValue(p, std::vector<unsigned int>());
    }
  }  
  
  // Step 2:  Apply second face scan to add 1 to the voxel which is the maximal value of the current scan.
  // stored in scoreConfidance.  
  Image3DDouble scoreConfidance (myDomain);
  for (unsigned int numNorm = 0; numNorm< myNormalField.size(); numNorm++){
    NormalAccumulator::Vector aNormal  = myNormalField[numNorm];    
    DGtal::Z3i::RealPoint centerPoint = myNormalOrigins[numNorm];
    DGtal::Z3i::RealPoint currentPoint = centerPoint;
    DGtal::Z3i::RealPoint previousPoint;
    
    DGtal::Z3i::Point aPtMaxAcc = DGtal::Z3i::Point(currentPoint,DGtal::functors::Round<>());
    unsigned int aMaxAcc = 0;
    
    while((currentPoint - centerPoint).norm()<myRadius){
      DGtal::Z3i::Point currentPointI = DGtal::Z3i::Point(currentPoint,DGtal::functors::Round<>());
      
      if(myDomain.isInside(currentPointI) && previousPoint != currentPoint){                
        unsigned int valAcc = myAccumulationImage(currentPointI);
        if( valAcc > aMaxAcc){
          aMaxAcc = valAcc;
          aPtMaxAcc = currentPointI;
        }        
        previousPoint = currentPoint;
      }
      previousPoint = currentPoint;
      currentPoint += aNormal;
    }
    if(aPtMaxAcc != centerPoint){
      scoreConfidance.setValue(aPtMaxAcc, scoreConfidance(aPtMaxAcc)+1);
      if (updateVertexAsso){
        std::vector<unsigned int> v = myAssociationImage(aPtMaxAcc);
        v.push_back(numNorm);        
        myAssociationImage.setValue(aPtMaxAcc, v);
      }
    }
  }
  
  // Step 3: Compute confidance image indicating the rate between accIsMax/acc
  for(auto &v: myDomain){
    if ( myAccumulationImage(v) > minAcc )
      myConfidenceImage.setValue(v, scoreConfidance(v)/(double)myAccumulationImage(v));
  }
  myIsAssociationCompFromConfidence = updateVertexAsso;
  myIsConfidenceComputed = true;
}




unsigned int 
NormalAccumulator::getMaxAccumulation() const {
  return myMaxAccumulation;
}

double
NormalAccumulator::getMaxRadius() const {
  return myMaxRadius;
}

DGtal::Z3i::Point
NormalAccumulator::getMaxRadiusPoint() const {
  return myMaxRadiusPoint;
}


DGtal::Z3i::Point
NormalAccumulator::getMaxAccumulationPoint() const{
  return myMaxAccumulationPoint;
}


NormalAccumulator::Image3D &
NormalAccumulator::getAccumulationImage() {
  return myAccumulationImage;
}


NormalAccumulator::Image3DDouble & 
NormalAccumulator::getConfidenceImage(){
  return myConfidenceImage;
}


NormalAccumulator::Image3DDouble & 
NormalAccumulator::getRadiusImage(){
  return myRadiusImage;
}


const 
NormalAccumulator::VectorContainer& 
NormalAccumulator::getNormalField () const{
  return myNormalField;
}

const 
NormalAccumulator::VectorContainer& 
NormalAccumulator::getNormalOrigins () const{
  return myNormalOrigins;
}

NormalAccumulator::Domain
NormalAccumulator::getDomain() const
{
  return myDomain;
}

/**
 * Overloads 'operator<<' for displaying objects of class 'NormalAccumulator'.
 * @param out the output stream where the object is written.
 * @param aColor the object of class 'NormalAccumulator' to write.
 * @return the output stream after the writing.
 */
void
NormalAccumulator::selfDisplay( std::ostream & out ) const {
  out << "----" << std::endl;
  DGtal::Z3i::Point dim = myDomain.upperBound()-myDomain.lowerBound();
  out << "NormalAccumulator: \ninitialized with:\n " << " \t# normals: " << myNormalField.size()
      << "\n\t# normal base: " << myNormalOrigins.size() << std::endl;    
  out << "Domain: \n\t" << "size:" << myDomain.size() << " (" << dim[0] << " " << dim[1] << " " << dim[2] << ")\n\t";
  out << "bounds: "<< myDomain.lowerBound() << " " << myDomain.upperBound() << std::endl;
  out << "\tmax accumulation: " << myMaxAccumulation  << " (" << myMaxAccumulationPoint << ")" << std::endl;
  out << "\t radius estim stat: " << myRadiusStatEstim   << std::endl;
  out << "----" << std::endl;
}


std::ostream&
operator<< ( std::ostream & out,
             const NormalAccumulator & aNormalAcc )
{
    aNormalAcc.selfDisplay ( out );
    return out;
}

