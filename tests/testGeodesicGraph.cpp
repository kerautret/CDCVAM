
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
  GeodesicGraphComputer gg(3.0, aSet, 2.0, d, Z3i::Point(531, 425, 302));
  gg.computeGraphFromGeodesic();
  DGtal::trace.info() << "max geodesic distance: " << gg.getMaxGeodesDistance() << std::endl;
  // Export dilatation
  GeodesicGraphComputer::TSet dilateSet = gg.getDilatedSet();  
  std::ofstream resDilate;
  resDilate.open("testDilate.sdp");
  for (auto const &p: dilateSet)
    {
      resDilate << p[0] << " " << p[1] << " " << p[2] << std::endl; 
    }
  resDilate.close();


  // Export geodesic labels
  const GeodesicGraphComputer::Image3D &imageLabel = gg.getLabelImage();
  std::ofstream geoLabel;
  geoLabel.open("testGeodes.sdp");
  unsigned int maxLabel = 0;
  for (auto const &p: imageLabel.domain()){
    if(imageLabel(p)>maxLabel){
      maxLabel = imageLabel(p);
    }
  }

  HueShadeColorMap<unsigned int > hueMap(0, maxLabel);

  for (auto const &p: imageLabel.domain())
    {
      if(imageLabel(p) != 0){
        auto color = hueMap(imageLabel(p));
        geoLabel << p[0] << " " << p[1] << " " << p[2];
        geoLabel << " " << (int) color.red() << " "
                 << (int) color.green() <<  " "
          << (int) color.blue() << std::endl;
      
      }
    }
  geoLabel.close();

  // Export graph degree
  const GeodesicGraphComputer::Image3DChar &imageDegree = gg.getDegreeImage();
  std::ofstream graphDegree;
  graphDegree.open("testDegree.sdp");
  for (auto const &p: imageDegree.domain())
    {
      if(imageDegree(p) != 0){
        graphDegree << p[0] << " " << p[1] << " " << p[2];
        if(imageDegree(p) == 1)
          {
            graphDegree << " 255 25 25 255" << std::endl;
          }
        else if(imageDegree(p) == 2)
          {
            graphDegree << " 25 25 255 255" << std::endl;
          }
        else if(imageDegree(p) == 3)
          {
            graphDegree << " 25 255 25 255" << std::endl;
          }
        else if(imageDegree(p) == 4)
          {
            graphDegree << " 25 255 255 255" << std::endl;
          }
        else if(imageDegree(p) == 5)
          {
            graphDegree << " 255 255 25 255" << std::endl;
          }
        else{
           graphDegree << " 25 25 25 255" << std::endl;
        }
      }
    }
  graphDegree.close();
  
}






