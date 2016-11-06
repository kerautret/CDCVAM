#include <DGtal/io/readers/PointListReader.h>
#include <DGtal/io/colormaps/HueShadeColorMap.h>
#include <iostream>
#include <fstream>



#include "GeodesicGraphComputer.h"


using namespace DGtal;





int
main(int argc,char **argv)
{


  std::vector<Z3i::Point> input = PointListReader<Z3i::Point>::getPointsFromFile("/Users/kerautre/EnCours/TubeAnalyse/build/conf.sdp");
  Z3i::Domain d (Z3i::Point(496, 400, 288), Z3i::Point(602, 476, 392));

  // get an input set
  GeodesicGraphComputer::TSet aSet(d);
  for (auto &p: input){
    aSet.insert(p);
  }
  GeodesicGraphComputer gg(5.0, aSet, 2.0, d, Z3i::Point(531, 425, 302));
  gg.computeGraphFromGeodesic();

  
  std::ofstream outVertex, outEdges;
  outVertex.open("testGeodesicGraphVertex.sdp");
  outEdges.open("testGeodesicGraphEdges.sdp");
  
  std::vector<Z3i::Point> vertices = gg.getGraphVerticesRepresentant();
  for (auto const &p: vertices)
    {
      outVertex << p[0] << " " << p[1] << " " << p[2] << std::endl;
    }

  outVertex.close();
  std::vector<std::pair<unsigned int, unsigned int> > edges;
  edges = gg.myGraph.getEdges();
  for (auto const &e: edges)
    {
      outEdges << e.first << " " << e.second << std::endl;
    }
  outEdges.close();
  
  
}






