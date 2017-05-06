#if defined(NORMAL_ACCUMULATOR_RECURSES)
#error Recursive header files inclusion detected in Color.h
#else // defined(NORMAL_ACCUMULATOR_RECURSES)
/** Prevents recursive inclusion of headers. */
#define NORMAL_ACCUMULATOR_RECURSES

#if !defined NORMAL_ACCUMULATOR_h
/** Prevents repeated inclusion of headers. */
#define NORMAL_ACCUMULATOR_h




#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Mesh.h"
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/images/ImageContainerBySTLVector.h>


///////////////////////////////////////////////////////////////////////////////
// class NormalAccumulator
/**
 * Description of class 'NormalAccumulator' <p>
 * 
 * @brief Class to compute volumic accumulation from normal vector field.
 *
 *
 */




class NormalAccumulator{


  // types of main object
public:
  enum StatRadiusDef {mean, median, min, max};
  typedef typename DGtal::Z3i::Domain Domain; 
  typedef typename DGtal::Z3i::RealPoint Point; 
  typedef typename DGtal::Z3i::RealPoint Vector; 
  typedef std::vector<Point> PointContainer;
  typedef std::vector<Vector> VectorContainer;  
  
  // types of image containers:
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, DGtal::uint64_t> Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, double> Image3DDouble;
  // to recover the origin point which has contributed to a particular voxel.
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,  std::vector<unsigned int> > ImagePointAssociation;

  

  // ----------------------- Standard services ------------------------------
  
public:
  
  /**
   * Default constructor
   **/

  NormalAccumulator(const double aRadius,
                    const std::string &rEstimDef = "min"): myRadius(aRadius), 
                                                          myAccumulationImage(DGtal::Z3i::Domain()),
                                                          myConfidenceImage(DGtal::Z3i::Domain()),
                                                          myRadiusImage(DGtal::Z3i::Domain()),
                                                          myAssociationImage(DGtal::Z3i::Domain())
                                                                         

  {
    if(rEstimDef =="mean"){
      myRadiusStatEstim = NormalAccumulator::StatRadiusDef::mean; 
    }else if(rEstimDef =="median"){
      myRadiusStatEstim = NormalAccumulator::StatRadiusDef::median; 
    } else if(rEstimDef =="min"){
      myRadiusStatEstim = NormalAccumulator::StatRadiusDef::min; 
    } else if(rEstimDef =="max"){
      myRadiusStatEstim = NormalAccumulator::StatRadiusDef::max; 
    } 
  }

  // ----------------------- Interface --------------------------------------
  
  /**
   * Initialize a NormalAccumulator object from a mesh:
   * - compute the vertex center of each face
   * - construct the normal set from the each vertex
   * - compute the image domain
   **/
  void initFromMesh(const DGtal::Mesh<Point> &aMesh, bool invertNormal=false);



  /**
   * Initialize a NormalAccumulator object from a mesh and normals: (used for digital object)
   *
   * - compute the vertex center of each face
   * - construct the normal set from the each vertex
   * - compute the image domain
   * 
   **/
  void initFromMeshAndNormals(const DGtal::Mesh<Point> &aMesh,
                              const std::vector<Point> &vectNormals,
                              bool invertNormal=false);
  

  /**
   * Initialize a NormalAccumulator object from a set points and from
   * normals vectors
   * 
   * 
   **/

  void initFromNormals(const std::vector<DGtal::Z3i::RealPoint> &aSetOfPoints,
                       const std::vector<DGtal::Z3i::RealPoint> &aSetOfNormals,
                       const bool invertNormal=false);
  
  /**
   * @todo with testing if PCL is installed
   * 
   **/
#if defined USE_PCL



  void initFromPointCloud(const std::vector<Point> &aVectPoints, const double size=5.0,
                          bool autoOrient=true);



  void initFromPointCloud(const std::vector<Point> &aVectPoints, const unsigned int nbNeighbohr = 10);


#endif

  
  /**
   * @todo import normal computation from DGtal surface
   * 
   **/
  //void initFromDGtalSurface(const DGtal::Surface &aSurface);



  /**
   * Compute accumulation  from normal vectors and update the maximum accumulation value.
   *
   **/
  
  void computeAccumulation(bool retainVertexAsso = true, bool verbose = true);

  
  /**
   * Compute the confidence image
   * 
   * @param updateVertexAsso if true update the image association with only the confident vertex. 
   * 
   **/
  void computeConfidence(bool updateVertexAsso = true, unsigned int minAcc=1);
  
  
  /**
   * Compute the radius from all contributing normal origins
   *
   **/
  void computeRadiusFromOrigins();



  /**
   * Compute the radius from all contributing normal with confidance rate.
   *
   **/
  void computeRadiusFromConfidence();
  


  /**
   * Get the associated input points associated to a voxel. It can be
   * all the points with normal contributing to the accumulation or
   * the confidence (if the point association is updated whe when
   * calling computeConfidence and passing true on argument
   * updateVertexAsso).
   *
   * @param aVoxel the input voxel
   * @return the container with all point contributing to the condidence.
   **/

  PointContainer getAssociatedPoints(const DGtal::Z3i::Point &aVoxel); 
  
  
  /**
   * return the maximum accumulation value. 
   *
   **/  
  unsigned int getMaxAccumulation() const;


  /**
   * return the maximum radius value. 
   *
   **/  
  double getMaxRadius() const;


  /**
   * return the maximum radius point. 
   *
   **/  
  DGtal::Z3i::Point getMaxRadiusPoint() const;
  

  /**
   * Return the point of maximum accumulation value.
   *
   **/  
  DGtal::Z3i::Point getMaxAccumulationPoint() const;

  const VectorContainer& getNormalField ()const;
 
  const VectorContainer& getNormalOrigins ()const;
  
  
  Image3D & getAccumulationImage() ;
  

  Image3DDouble & getConfidenceImage() ;


  Image3DDouble & getRadiusImage() ;

  Domain getDomain() const;
  
  /**
   * self display method.
   **/
  void  selfDisplay( std::ostream & out) const;

  

  

  
  // protected attributes:
protected:

  VectorContainer myNormalField; // the objet normals
  PointContainer myNormalOrigins; // the origin of the object normals
  double myRadius; // the maximal radius of accumulation
  Domain myDomain; // the domain of the source objet.
  

  Image3D myAccumulationImage;
  Image3DDouble myConfidenceImage;
  Image3DDouble myRadiusImage;    
  ImagePointAssociation myAssociationImage; // to recover all normal origins contributing to an accumulation point.


private:
  bool myIsVerbose = true;
  bool myIsInitialized = false;
  bool myIsAccumulationComputed = false;
  bool myIsConfidenceComputed = false;
  bool myIsRadiusComputed = false;
  bool myIsAssociationCompFromConfidence = false;
  DGtal::Z3i::Point myMaxAccumulationPoint;
  unsigned int myMaxAccumulation = 0;
  double myMaxRadius = 0;
  DGtal::Z3i::Point myMaxRadiusPoint; 
  StatRadiusDef myRadiusStatEstim = StatRadiusDef::min;  
};



/**
 * Overloads 'operator<<' for displaying objects of class 'NormalAccumulator'.
 * @param out the output stream where the object is written.
 * @param aColor the object of class 'NormalAccumulator' to write.
 * @return the output stream after the writing.
 */
std::ostream&
operator<< ( std::ostream & out, const NormalAccumulator & aNormalAcc );




///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined NORMAL_ACCUMULATOR_h

#undef NORMAL_ACCUMULATOR_RECURSES
#endif // else defined(NORMAL_ACCUMULATOR_RECURSES)


