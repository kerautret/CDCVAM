#include <DGtal/io/viewers/Viewer3D.h>
#include <DGtal/io/readers/PointListReader.h>
#include <iostream>
#include <fstream>



#include "GeodesicGraphComputer.h"


using namespace DGtal;

typedef Viewer3D<> MyViewer;


typedef struct VoxelsAndGeodesicGraph
{ 
  std::vector<Z3i::Point> *vectVoxels;
  GeodesicGraphComputer *geoGraph;
  
}VoxelsAndGeodesicGraph;


int displayVoxelInfos( void* viewer, int name, void* data )
{
  VoxelsAndGeodesicGraph*  aData = (VoxelsAndGeodesicGraph*) data;
  std::vector<Z3i::Point>* vectVoxels = aData->vectVoxels;
  MyViewer *aViewer = (MyViewer *) viewer;
  GeodesicGraphComputer * gg = aData->geoGraph;
  std::stringstream ss ; ss << "Voxel label: " << (int) gg->getLabelImage()((*vectVoxels)[ name ]);
  ss << "degree "<< (int) gg->getDegreeImage()((*vectVoxels)[ name ])
     << " representant:" << gg->myGraph.getVertexRepresentant(gg->getLabelImage()((*vectVoxels)[ name ]));

  aViewer->setForegroundColor (QColor(255, 20, 20));
  aViewer->displayMessage(QString(ss.str().c_str()), 100000);
  aViewer->deleteCubeList(name);
  //(*aViewer).updateList(false);
  std::cout << "voxel: " << (*vectVoxels)[ name ] << " selected." << std::endl;
  return 0;
}


int
main(int argc,char **argv)
{

  QApplication app(argc, argv);
  MyViewer viewer;
  MyViewer viewer2;

  
  std::vector<Z3i::Point> input = PointListReader<Z3i::Point>::getPointsFromFile("/Users/kerautre/EnCours/TubeAnalyse/build/conf.sdp");
  Z3i::Domain d (Z3i::Point(496, 400, 288), Z3i::Point(602, 476, 392));

  // get an input set
  GeodesicGraphComputer::TSet aSet(d);
  for (auto &p: input){
    aSet.insert(p);
  }
  GeodesicGraphComputer gg(5.0, aSet, 2.0, d, Z3i::Point(531, 425, 302));
  gg.setModeComputeDegreeImage(true);
  gg.computeGraphFromGeodesic();

  // Export dilatation
  GeodesicGraphComputer::TSet dilateSet = gg.getDilatedSet();  
  std::ofstream resDilate;
  resDilate.open("testDilate.sdp");
  for (auto const &p: dilateSet)
    {
      resDilate << p[0] << " " << p[1] << " " << p[2] << std::endl; 
    }
  resDilate.close();



  // display graph degree
  std::vector<Z3i::Point> vectVoxelsDisplayed;
  const GeodesicGraphComputer::Image3DChar &imageDegree = gg.getDegreeImage();
  std::ofstream graphDegree;
  graphDegree.open("testDegree.sdp");
  int name = 0;
  for (auto const &p: imageDegree.domain())
    {
      if(imageDegree(p) != 0){
        
        if(imageDegree(p) == 1)
          {
            viewer<< DGtal::CustomColors3D(DGtal::Color(255, 25, 25, 255),DGtal::Color(255, 25, 25, 255)); 
          }
        else if(imageDegree(p) == 2)
          {
            viewer<< DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 255),DGtal::Color(25, 25, 255, 255)); 
          }
        else if(imageDegree(p) == 3)
          {
            viewer<< DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 255),DGtal::Color(25, 25, 255, 255)); 
          }
        else if(imageDegree(p) == 4)
          {
            viewer<< DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 255),DGtal::Color(25, 25, 255, 255)); 
          }
        else if(imageDegree(p) == 5)
          {
            viewer<< DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 255),DGtal::Color(25, 25, 255, 255)); 
          }
        else{
          viewer<< DGtal::CustomColors3D(DGtal::Color(25, 25, 25, 255),DGtal::Color(25, 25, 25, 255)); 
 
        }
        vectVoxelsDisplayed.push_back(p);
        viewer << DGtal::SetName3D(name++) << p; 
      }
    }

  graphDegree.close();
  VoxelsAndGeodesicGraph data;
  data.vectVoxels = &vectVoxelsDisplayed;
  data.geoGraph = &gg;

  viewer << DGtal::SetSelectCallback3D(displayVoxelInfos, &data, 0, vectVoxelsDisplayed.size());
  viewer << DGtal::CustomColors3D(DGtal::Color::Black, DGtal::Color::Black);
  viewer << Z3i::Point(531, 425, 302); 
  viewer<< MyViewer::updateDisplay;

  // test display graph vertex:
  
  std::vector<Z3i::Point> vertices = gg.getGraphVerticesRepresentant();
  for (auto const &p: vertices)
    {
      viewer2 << p;
    }

  viewer2<< MyViewer::updateDisplay;
  
  viewer2.show();
  viewer.show();
  return app.exec();

}






