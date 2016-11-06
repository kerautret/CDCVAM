#if defined(GEODESIC_GRAPH_RECURSES)
#error Recursive header files inclusion detected in GeodesicGraphComputer.h
#else // defined(GEODESIC_GRAPH_RECURSES)
/** Prevents recursive inclusion of headers. */
#define GEODESIC_GRAPH_RECURSES

#if !defined GEODESIC_GRAPH_h
/** Prevents repeated inclusion of headers. */
#define GEODESIC_GRAPH_h




#include <iostream>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/Mesh.h>
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/kernel/sets/DigitalSetBySTLSet.h>
#include <DGtal/geometry/volumes/distance/FMM.h>
#include <DGtal/topology/SurfelAdjacency.h>


///////////////////////////////////////////////////////////////////////////////
// class NormalAccumulator
/**
 * Description of class 'GeodesicGraphComputer' <p>
 * 
 * @brief Class to construct a graph from geodesic propagation of input point set.
 *
 *
 */




class GeodesicGraphComputer{


public:  

  //////////////////////////////////////////////////////////////////////
  // Data related to graph representation and construction:

  class GGraph{
    
  public:

    /**
     * Add a new vertex to the graph and returns its associated id.
     **/

    unsigned int addVertex()
    {
      std::vector<unsigned int> adj;
      myAdjacencyList.push_back(adj);
      myAssociatedPoint.push_back(DGtal::Z3i::Point(0,0,0));
      myActiveVertex.push_back(true);
      return myAdjacencyList.size()-1;
    } 


    
    void addEdge(const unsigned int labelV1, const unsigned int labelV2)
    {
      myAdjacencyList.at(labelV1).push_back(labelV2);
      myAdjacencyList.at(labelV2).push_back(labelV1);
    } 

    
    unsigned int nbVertex()
    {
      return myAdjacencyList.size();
    } const 

    
    void associatePointToVertex(const unsigned int labelVertex,
                                const DGtal::Z3i::Point &p)
    {
      myAssociatedPoint[labelVertex] = p;
    }
    

    std::vector<std::pair<unsigned int, unsigned int> >
    getEdges()
    {
      std::vector<std::pair<unsigned int, unsigned int>> result;
      std::map<std::pair<unsigned int, unsigned int >, bool > mapEdges;
      for(unsigned int i=0 ; i<nbVertex(); i++)
        {
          for(unsigned int k = 0 ; k < myAdjacencyList[i].size(); k++)
            {
              std::pair<unsigned int, unsigned int> pair;
              pair.first = i<myAdjacencyList[i][k] ? i : myAdjacencyList[i][k];
              pair.second = i<myAdjacencyList[i][k] ? myAdjacencyList[i][k]: i;
              mapEdges[pair] = true;
            }
        }
      for (auto  e = mapEdges.begin(); e != mapEdges.end(); e++)
        {
          result.push_back(e->first);
        }
      return result;
    } const
    

    DGtal::Z3i::Point getVertexRepresentant(unsigned int label)
    {
      return myAssociatedPoint[label];
    } const
    
    
    std::vector<unsigned int> getAdjacentVertices(unsigned int vLabel)
    {
      assert(myAdjacencyList.size()>vLabel);
      return myAdjacencyList[vLabel];      
    } const 
    

    unsigned int getLastAdjacentVertex(unsigned int label)
    {
      assert(myAdjacencyList.size()>label);
      return myAdjacencyList[label][myAdjacencyList[label].size()-1];      
    } const
    
    
    void setAllVerticesInactive()
    {
      for(unsigned int i = 0; i < myActiveVertex.size(); i++)
        {
          myActiveVertex[i] = false;
        }
    }


    bool hasActiveVertices()
    {
      for(unsigned int i = 0; i < myActiveVertex.size(); i++)
        {
          if(myActiveVertex[i])
            {
              return true;
            }
        }
      return false;
    }

    
    std::vector<unsigned int> getActiveVertex()
    {
      std::vector<unsigned int> result;
      for(unsigned int i = 0; i < myActiveVertex.size(); i++)
        {
          if(myActiveVertex[i])
            {
              result.push_back(i);
            }
        }      
      return result;
    } const 

    
    bool isActiveVertex(const unsigned int label)
    {
      assert(label<myActiveVertex.size());
      return myActiveVertex[label];
    } const

    
    void unactiveVertex(const unsigned int label)
    {
      assert(label<myActiveVertex.size());
      myActiveVertex[label]=false;
    }

    
    void activeVertex(const unsigned int label)
    {
      assert(label<myActiveVertex.size());
      myActiveVertex[label]=true;
    }

    
    // basic but just need at the end of tracking of one graph component
    void updateAdjacency()
    {
      for (unsigned int i=0; i<myAdjacencyList.size(); i++)
        {
          std::vector<unsigned int> newAdj;
          for(auto const &e: myAdjacencyList[i])
            {
              if(e < nbVertex())
                {
                  newAdj.push_back(e);
                }
            }
          myAdjacencyList[i] = newAdj;
        }

    }
    
    std::vector<DGtal::Z3i::Point> myAssociatedPoint; // to have 3D points associated to a graph vertex
    std::vector<std::vector<unsigned int>> myAdjacencyList;  /// to represent the graph
    std::vector<bool> myActiveVertex; /// the vertices which contribute to extension
  };


  // End of class GGraph
  //////////////////////////////////////////////////////////////////////


  

    


  typedef DGtal::DigitalSetBySTLSet<DGtal::Z3i::Domain> TSet;   
  
  
  /////////////////////////////////////////
  // Data related to geodesic :

  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain,double> DistanceImage;
  typedef TSet AcceptedPointSet;
  typedef DGtal::Z3i::Domain::Predicate DomainPredicate;
  typedef DGtal::L2FirstOrderLocalDistance<DistanceImage, AcceptedPointSet> DistanceMeasure;
  typedef DGtal::FMM<DistanceImage, AcceptedPointSet, DomainPredicate, DistanceMeasure> TFMM_Dilate;
  typedef DGtal::FMM<DistanceImage, AcceptedPointSet, TSet, DistanceMeasure> TFMM_Geodes;

  
  
  // types of image containers:
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, DGtal::uint64_t> Image3D;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, double> Image3DDouble;
  typedef DGtal::ImageContainerBySTLVector<DGtal::Z3i::Domain, unsigned char> Image3DChar;  


  
private:


  ///////////////////////////////////////////
  // Internal data:
  
  TSet myInputPoints; /// the input points 
  TSet myDilatedPoints; 
  DGtal::Z3i::Domain myDomain; /// the domain associated to the input points
  double myGeodesicStep = 5.0;
  double myDilateDistance = 2.0;
  DGtal::Z3i::Point myStartPoint;
  Image3DDouble myDistanceImage;
  bool myMaintainImageDegree = false;
  Image3D myLabelImage; /// if the voxel ilabel associated to each CC comp on geodesic parts
  Image3DChar myDegreeImage; /// the degree of the graph being reconstructed. 
  double myMaxGeodesDistance = 0; /// the maximal distance obtained in the geodesic extension
public:
  Image3D myOrderGeodesImage; /// used only to illustrate geodesic progapation dans dist step



public:

  GGraph myGraph;
  
  // ----------------------- Standard services ------------------------------
  
public:
  
  /**
   * Default constructor
   **/

  GeodesicGraphComputer(const double aStep, const TSet &aPointSet, const double dilateDistance,
                const DGtal::Z3i::Domain aDomain, DGtal::Z3i::Point aStartPoint): 
                                                                                  myInputPoints(aPointSet),
                                                                                  myDilatedPoints(aDomain),
                                                                                  myDomain(aDomain),
                                                                                  myGeodesicStep(aStep),
                                                                                  myDilateDistance(dilateDistance),
                                                                                  myStartPoint(aStartPoint),
                                                                                  myDistanceImage(aDomain),
                                                                                  myLabelImage(aDomain),
                                                                                  myDegreeImage(aDomain),
                                                                                  myOrderGeodesImage(aDomain)
  {
    for(const auto &p: aDomain)
      {
        myLabelImage.setValue(p, -1);
      }
    for (const auto &p: myInputPoints)
      {
        myDilatedPoints.insert(p);
      }
  };
  
  // ----------------------- Interface --------------------------------------
  
  /**
   * Compute all the graph by iteratively apply a graph connected component construction 
   * @see computeGraphConnectedComponent
   *
   * @param[in] nbCompMax the maximal number of graph part 
   * @param[in] minSize the min size to consider the connected part of the graph.
   *
   **/
  
  void computeGraphFromGeodesic(const unsigned int nbCompMax = 1000,
                                const unsigned int minSize = 0);


  /**
   * Compute one connected component and return the number of source
   * point processed.
   *  
   * @param[in] initPoint the initial point to start the reconstruction.
   * @param[in] verbose to output processing information.
   **/
  
  unsigned int computeGraphConnectedComponent(const DGtal::Z3i::Point initPoint, bool verbose = true);

   
  /**
   *  Return the dilated set of points
   */

  TSet getDilatedSet();

  
  /**
   * Return the label image
   */
  
  const Image3D & getLabelImage();

    
  /**
   * Return the degree of the label image 
   **/
  const Image3DChar & getDegreeImage();

  const Image3DDouble & getDistanceImage();
  
  const double getMaxGeodesDistance();
  
  const std::vector<DGtal::Z3i::Point> getGraphVerticesRepresentant();


  /**
   * Set a mode to compute the image degree. Can be used to debug and
   * display the degree related information.
   * 
   * @param status if true it update the degree image. 
   * 
   **/
  
  void setModeComputeDegreeImage(bool status);
  
  
  /**
   * self display method.
   *
   * @param out the out stream.
   **/
  void  selfDisplay( std::ostream & out) const;

  

  // ----------------------- Internal use --------------------------------------

  /**
   * Use to transmit labels
   * @param[in] pt a point
   **/

 int getNeighborhoodLabel(const DGtal::Z3i::Point &pt);

  
  /**
   * Computes the number of connected component of a given set and
   * update data: (label image @ref myLabelImage and vertex degree
   * image @ref myDegreeImage).
   *
   * @param[in] aSet the input set of point.
   * @param[in] originLabel the label of the vertex which generate it in the graph construction 
   **/

  unsigned int computeCCnumberAndUpdateData(const TSet &aSet, const unsigned int originLabel) ;
  


  /**
   * Filters a set and retains only points which are labeled with a specific label.
   *  
   * @param[in] aSet an input set of voxels.
   * @param[in] label the label of the points to be filtered.
   *
   **/

  TSet filterSetFromLabel(const GeodesicGraphComputer::TSet &aSet, const unsigned int label) const;

    
  /**
   * Simple function to test if a set of surfel represent a hole or not.
   * @param[in] K a KhalimskySpaceND 
   * @param[in] SAdj an surfel adjacency
   * @param[in] vectSurfel the surfel to be tested if it correspond to an hole
   * @param[in] aSet the source voxel set.
   **/ 

  bool isHoleSurface(const DGtal::Z3i::KSpace &K,
                     const DGtal::SurfelAdjacency<3> &SAdj,
                     const std::vector<DGtal::Z3i::SCell> &vectSurfel,
                     const GeodesicGraphComputer::TSet &aSet) const ;



  /**
   *  To associate the vertex of the graph to a 3D representant.
   *  process all image labels and compute the mean representant.
   **/
  void computeVertex3DRepresentant();
  
  

  
  
};



/**
 * Overloads 'operator<<' for displaying objects of class 'NormalAccumulator'.
 * @param out the output stream where the object is written.
 * @param aColor the object of class 'NormalAccumulator' to write.
 * @return the output stream after the writing.
 */

std::ostream&
operator<< ( std::ostream & out, const GeodesicGraphComputer & aNormalAcc );




///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods 


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined GEODESIC_GRAPH_h

#undef GEODESIC_GRAPH_RECURSES
#endif // else defined(GEODESIC_GRAPH_RECURSES)


