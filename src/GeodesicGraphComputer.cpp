
#include "GeodesicGraphComputer.h"

#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/SurfelSetPredicate.h"


///////////////////////////////////////////////////////////////////////////////
// class GeodesicGraphComputer
///////////////////////////////////////////////////////////////////////////////




void
GeodesicGraphComputer::computeGraphFromGeodesic(const unsigned int nbCompMax,
                                                const unsigned int minSize )
{
  unsigned int nbComp = 0;
  // 1) First step: compute the dilatation of input set
  TFMM_Dilate fmm( myDistanceImage, myDilatedPoints, myDomain.predicate(), 
                   std::numeric_limits<int>::max(), myDilateDistance );
  fmm.compute();

  computeGraphConnectedComponent(myStartPoint);
  nbComp++;
  while(nbComp < nbCompMax &&  myInputPoints.size() != 0){
    unsigned int nbVertex = myGraph.nbVertex();
    
    int size = computeGraphConnectedComponent(*(myInputPoints.begin()));
    nbComp++;
    if (size < minSize){
      myGraph.myAdjacencyList.erase(myGraph.myAdjacencyList.begin()+nbVertex,
                                    myGraph.myAdjacencyList.end());
      myGraph.myActiveVertex.erase(myGraph.myActiveVertex.begin()+nbVertex,
                                    myGraph.myActiveVertex.end());
      myGraph.myAssociatedPoint.erase(myGraph.myAssociatedPoint.begin()+nbVertex,
                                    myGraph.myAssociatedPoint.end());
      myGraph.updateAdjacency();
    }
  }  
  DGtal::trace.info() << "Total detected " << nbComp << "components"  << std::endl;  
}



unsigned int
GeodesicGraphComputer::computeGraphConnectedComponent(const DGtal::Z3i::Point initPoint, bool verbose)
{
  
  
  // 2) Second step: apply geodesic:
  TSet geoSet(myDomain);
  geoSet.insert(initPoint);
  for(auto const &p: myDomain) myDistanceImage.setValue(p, 0.0);
  TFMM_Geodes geodesicFMM(myDistanceImage, geoSet, myDilatedPoints,
                          std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
  DGtal::Z3i::Point lastPt=initPoint;
  double lastDistance = 0.0;
  double distStartInterval = myGeodesicStep;
  unsigned int nbVertexIni = myGraph.nbVertex();
  unsigned int label = myGraph.addVertex();
  myLabelImage.setValue(initPoint, label);

  if (verbose) DGtal::trace.info() << "[." ;
  // apply initial geodesic
  unsigned int order = 1;
  while( lastDistance <= distStartInterval && geodesicFMM.computeOneStep(lastPt, lastDistance) )
    {
      myOrderGeodesImage.setValue(lastPt, order);
      myLabelImage.setValue(lastPt, label );
      myDegreeImage.setValue(lastPt, 1);
      myInputPoints.erase(lastPt);

    }
  
  distStartInterval += myGeodesicStep;
  // main extension algorithm
  while(myGraph.hasActiveVertices())
    {
      order++;
      myGraph.setAllVerticesInactive();
      GeodesicGraphComputer::TSet setNewGeo(myDilatedPoints.domain());
      
      // Graph Construction Step 1 : extension of connected points to the graph
      // - transmit labels
      
      while( lastDistance <= distStartInterval && geodesicFMM.computeOneStep(lastPt, lastDistance) )
        {
          if(lastDistance>myMaxGeodesDistance){
            myMaxGeodesDistance = lastDistance;
          }
          myOrderGeodesImage.setValue(lastPt, order);
          int labelNeigbor = GeodesicGraphComputer::getNeighborhoodLabel(lastPt);
          myInputPoints.erase(lastPt);
          if(labelNeigbor != -1){
            setNewGeo.insert(lastPt);
            myLabelImage.setValue(lastPt, labelNeigbor );
            myGraph.activeVertex(labelNeigbor);
          }else{
            int newLabel  = myGraph.addVertex();
            myLabelImage.setValue(lastPt, newLabel );
            myGraph.activeVertex(newLabel);
          }
        }
      
      
      // Graph Construction Step 2: add new vertices to the graph:
      for(auto  const & newPt: setNewGeo)
        {
          // - step 2a : each extended areas of setNewGeo produces a new
          //  graph vertex according their origins label which are
          //  transmitted on myLabelImage. (and mark the origin vertex as unactive).
          
          if(myGraph.isActiveVertex(myLabelImage(newPt)))
            { // First time found we create a new graph vertex
              if (verbose) DGtal::trace.info() << ".";
              myGraph.unactiveVertex(myLabelImage(newPt));
              unsigned int newLabel = myGraph.addVertex();
              myGraph.activeVertex(newLabel);
              myGraph.addEdge(newLabel, myLabelImage(newPt));
              myLabelImage.setValue(newPt, newLabel);
            }
          else
            {
              myLabelImage.setValue(newPt, myGraph.getLastAdjacentVertex(myLabelImage(newPt)));
            }
        }
      
      // Step 2b: we check if each node are connected or not: if not we update the graph.
      std::vector<unsigned int> vertexActive =  myGraph.getActiveVertex();
      for(unsigned int i =0; i < vertexActive.size(); i++){
        GeodesicGraphComputer::TSet aSet = filterSetFromLabel(setNewGeo, vertexActive[i]);
        computeCCnumberAndUpdateData(aSet, myGraph.getAdjacentVertices(vertexActive[i])[0]);
      }
      distStartInterval += myGeodesicStep;
    }
  myGraph.setAllVerticesInactive();
    
  if (verbose) DGtal::trace.info() << "]"<< std::endl;
  // Step 3: Final graph reconstruction: associate to each vertex its 3D associated point.
  myGraph.associatePointToVertex(nbVertexIni, initPoint);
  computeVertex3DRepresentant();
  myInputPoints.erase(initPoint);
  
  return myGraph.nbVertex() - nbVertexIni;
}




int
GeodesicGraphComputer::getNeighborhoodLabel(const DGtal::Z3i::Point &pt)
{
  for(int i = -1; i<=1; i++){
    for(int j = -1; j<=1; j++){
      for(int k = -1; k<=1; k++){
        DGtal::Z3i::Point neighbor = pt+ DGtal::Z3i::Point(i,j,k);
        if(myLabelImage.domain().isInside(neighbor) && myLabelImage(neighbor) != -1){
          return myLabelImage(neighbor);
        }          
      }
    }
  }
  return -1;
}




unsigned int
GeodesicGraphComputer::computeCCnumberAndUpdateData(const GeodesicGraphComputer::TSet &aSet,
                                                    const unsigned int originLabel)
{  
  DGtal::Z3i::KSpace K;
  K.init(aSet.domain().lowerBound(), aSet.domain().upperBound(), true);
  DGtal::SurfelAdjacency<3> SAdj( false );
  std::vector<std::vector<DGtal::Z3i::SCell> > vectConnectedSCell, vectConnectedSCellTmp;
  DGtal::Surfaces<DGtal::Z3i::KSpace>::extractAllConnectedSCell(vectConnectedSCellTmp, K, SAdj, aSet, false);
  
  unsigned int nbHoles = 0;
  for (unsigned int i = 0; i< vectConnectedSCellTmp.size(); i++ ){
    if(!isHoleSurface(K, SAdj, vectConnectedSCellTmp[i], aSet)){
      vectConnectedSCell.push_back(vectConnectedSCellTmp[i]);
    }else{
      nbHoles++;
    }
  }
  if (vectConnectedSCell.size() >= 2){
    if(myMaintainImageDegree){
      std::set<DGtal::Z3i::SCell> setSCell; for (auto const &s: vectConnectedSCell[0]) setSCell.insert(s);
      DGtal::functors::SurfelSetPredicate<std::set<DGtal::Z3i::SCell>,DGtal::Z3i::SCell> surfacePred (setSCell);
      DGtal::Surfaces<DGtal::Z3i::KSpace>::uFillInterior(K, surfacePred, myDegreeImage, vectConnectedSCell.size(),
                                                         false, false);
    }
    for (unsigned int i = 1; i< vectConnectedSCell.size(); i++ ){
      std::set<DGtal::Z3i::SCell> setSCell2; for (auto const &s: vectConnectedSCell[i]) setSCell2.insert(s);
      DGtal::functors::SurfelSetPredicate<std::set<DGtal::Z3i::SCell>,DGtal::Z3i::SCell> surfacePred2 (setSCell2);
      unsigned int newLabel = myGraph.addVertex();
      myGraph.addEdge(newLabel, originLabel);
      DGtal::Surfaces<DGtal::Z3i::KSpace>::uFillInterior(K, surfacePred2,
                                                         myLabelImage, newLabel,
                                                         false, false);
      if(myMaintainImageDegree){
        DGtal::Surfaces<DGtal::Z3i::KSpace>::uFillInterior(K, surfacePred2,
                                                           myDegreeImage,
                                                           vectConnectedSCell.size(),
                                                           false, false);
      }
    }
  }else {
    if(myMaintainImageDegree){
      std::set<DGtal::Z3i::SCell> setSCell; for (auto const &s: vectConnectedSCell[0]) setSCell.insert(s);
      DGtal::functors::SurfelSetPredicate<std::set<DGtal::Z3i::SCell>,DGtal::Z3i::SCell> surfacePred (setSCell);
      DGtal::Surfaces<DGtal::Z3i::KSpace>::uFillInterior(K, surfacePred, myDegreeImage, 1,
                                                         false, false);

    }
  }
  return vectConnectedSCell.size();
}




GeodesicGraphComputer::TSet
GeodesicGraphComputer::filterSetFromLabel(const GeodesicGraphComputer::TSet &aSet, unsigned int aLabel) const
{
  
  DGtal::Z3i::Point lowPt, upperPt;
  TSet preResult( aSet.domain() );
  lowPt = *(aSet.begin());
  upperPt = *(aSet.begin());

  for( auto const &p: aSet)
    {
      if(myLabelImage(p)==aLabel){
        preResult.insert(p);
        lowPt = lowPt.inf(p);
        upperPt = upperPt.sup(p);
      }
    }
  TSet result(  DGtal::Z3i::Domain(lowPt-DGtal::Z3i::Point::diagonal(1),
                                   upperPt+DGtal::Z3i::Point::diagonal(1)));
  for( auto const &p: preResult)
    {
       result.insert(p);
    }
  return result;
  
}
   





GeodesicGraphComputer::TSet
GeodesicGraphComputer::getDilatedSet()
{
  return myDilatedPoints;
}

const GeodesicGraphComputer::Image3D &
GeodesicGraphComputer::getLabelImage()
{
  return myLabelImage;
}

const GeodesicGraphComputer::Image3DChar &
GeodesicGraphComputer::getDegreeImage()
{
  return myDegreeImage;
}
const GeodesicGraphComputer::Image3DDouble &
GeodesicGraphComputer::getDistanceImage()
{
  return myDistanceImage;
}

const double
GeodesicGraphComputer::getMaxGeodesDistance()
{
  return myMaxGeodesDistance;
}

const std::vector<DGtal::Z3i::Point>
GeodesicGraphComputer::getGraphVerticesRepresentant()
{
  std::vector<DGtal::Z3i::Point> result;
  for (unsigned int i = 0; i < myGraph.nbVertex(); i++)
    {
      result.push_back(myGraph.getVertexRepresentant(i));
    }
  return result;
}


void
GeodesicGraphComputer::setModeComputeDegreeImage(bool status)
{
  myMaintainImageDegree = status;
}



bool
GeodesicGraphComputer::isHoleSurface(const DGtal::Z3i::KSpace &K,
                             const DGtal::SurfelAdjacency<3> &SAdj,
                             const std::vector<DGtal::Z3i::SCell> &vectSurfel,
                             const GeodesicGraphComputer::TSet &aSet) const
{
  DGtal::Z3i::SCell aSurfel = vectSurfel[0];
  // get orht Dir associated to the surfel
  auto dirVoxel = *(K.sOrthDirs(aSurfel));
  DGtal::Z3i::SCell aVoxel = K.sIncident(aSurfel, dirVoxel, !K.sSign(aSurfel));
  return aSet.find(K.sCoords(aVoxel)) == aSet.end();
}



void
GeodesicGraphComputer::computeVertex3DRepresentant()
{
  // for each graph vertex we construct their representant
  std::vector<DGtal::Z3i::Point> vertexRepresentant;
  /// used to compute the mean representant, for each indexed vertex, stores the number of voxel which belong to its label
  std::vector<unsigned int> vertexNb; 
  unsigned int nbVertex = myGraph.nbVertex();  
  for(unsigned int i=0; i<nbVertex; i++)
    {
        vertexRepresentant.push_back(DGtal::Z3i::Point(0,0,0));
        vertexNb.push_back(0);
    }
  
  // processing all label of the label image:
  for(auto const &p: myLabelImage.domain())
    {
      if(myLabelImage(p) != -1){
        vertexRepresentant[myLabelImage(p)] = vertexRepresentant[myLabelImage(p)]+p;
        vertexNb[myLabelImage(p)] = vertexNb[myLabelImage(p)] + 1;
        myLabelImage.setValue(p, -1);
      }
    }

  // compute the mean representant
  for(unsigned int i=0; i<nbVertex; i++)
    {
      if(vertexNb[i]!=0){
        vertexRepresentant[i] = vertexRepresentant[i]/(int)vertexNb[i];
      }
    }
  // Special case to handle (first vertex is 0 labeled as background voxel, but is a single point)
  vertexRepresentant[0] = myStartPoint;

  for(unsigned int i=0; i<nbVertex; i++)
    {
      if(vertexNb[i]!=0){
        myGraph.associatePointToVertex(i, vertexRepresentant[i]);
        myLabelImage.setValue(vertexRepresentant[i], -1);
      }
    }
  
}





/**
 * Overloads 'operator<<' for displaying objects of class 'GeodesicGraphComputer'.
 * @param out  the object of class 'GeodesicGraphComputer' to write.
 */
void
GeodesicGraphComputer::selfDisplay( std::ostream & out ) const {
  out << "----" << std::endl;
  DGtal::Z3i::Point dim = myDomain.upperBound()-myDomain.lowerBound();
  out << "GeodesicGraphComputer: \n Point set size: " << myInputPoints.size();
  out << "Domain: \n\t" << "size:" << myDomain.size() << " (" << dim[0] << " " << dim[1] << " " << dim[2] << ")\n\t";
  out << "----" << std::endl;
}




std::ostream&
operator<< ( std::ostream & out,
             const GeodesicGraphComputer & aGeodesicG )
{
    aGeodesicG.selfDisplay ( out );
    return out;
}

