#include <DGtal/io/viewers/Viewer3D.h>
#include <DGtal/io/readers/PointListReader.h>
#include <DGtal/io/readers/MeshReader.h>
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
  ss << "degree "<< (int) gg->getDegreeImage()((*vectVoxels)[ name ]);
//     << " representant:" << gg->myGraph.getVertexRepresentant(gg->getLabelImage()((*vectVoxels)[ name ]));

  aViewer->setForegroundColor (QColor(255, 20, 20));
  aViewer->displayMessage(QString(ss.str().c_str()), 100000);
  //aViewer->deleteCubeList(name);
  //(*aViewer).updateList(false);
  std::cout << "voxel: " << (*vectVoxels)[ name ] << " selected." << std::endl;
  return 0;
}


int
main(int argc,char **argv)
{

  QApplication app(argc, argv);
  
  MyViewer viewerStep0;
  MyViewer viewerStep1;
  MyViewer viewerStep2;
  MyViewer viewerStep3;
  MyViewer viewerStep4;
  MyViewer viewerStep5;
  
  std::vector<Z3i::Point> input = PointListReader<Z3i::Point>::getPointsFromFile("/Users/kerautre/EnCours/TubeAnalyse/build/conf2.sdp");
  DGtal::Mesh<Z3i::RealPoint> myMesh(false);
  myMesh << "/Users/kerautre/EnCours/TubeAnalyse/Samples2/bk1.off";
  
  Z3i::Domain d (Z3i::Point(496, 400, 288), Z3i::Point(602, 476, 392));

  // get an input set
  GeodesicGraphComputer::TSet aSet(d);
  for (auto &p: input){
    aSet.insert(p);
    if(p != Z3i::Point(531, 425, 302)){
      viewerStep0 << p;
    }
  }
  viewerStep0<< DGtal::CustomColors3D(DGtal::Color(255, 25, 25, 255),DGtal::Color(255, 25, 25, 255));
  viewerStep0 << Z3i::Point(531, 425, 302);
  GeodesicGraphComputer gg(3.0, aSet, 2, d, Z3i::Point(531, 425, 302));
  gg.setModeComputeDegreeImage(true);
  gg.computeGraphFromGeodesic();
  
  // test display graph vertex:
  std::vector<Z3i::Point> vertices = gg.getGraphVerticesRepresentant();
  for (auto const &p: vertices)
  {
    viewerStep5<< DGtal::CustomColors3D(DGtal::Color(255, 25, 25, 255),DGtal::Color(255, 25, 25, 255));
    viewerStep5 << p;
    viewerStep4<< DGtal::CustomColors3D(DGtal::Color(255, 25, 25, 255),DGtal::Color(255, 25, 25, 255));
    viewerStep4 << p;
    
  }
  // with edges:
  auto edges = gg.myGraph.getEdges();
  for (auto const &e: edges)
  {
    DGtal::Z3i::Point p1 = gg.myGraph.getVertexRepresentant(e.first);
    DGtal::Z3i::Point p2 = gg.myGraph.getVertexRepresentant(e.second);
    viewerStep5.setLineColor(DGtal::Color::Red);
    viewerStep5.addLine(p1, p2);
    
  }
  

  // Export dilatation
  GeodesicGraphComputer::TSet dilateSet = gg.getDilatedSet();  
  std::ofstream resDilate;
  resDilate.open("testDilate.sdp");
  DGtal::HueShadeColorMap<double> colorDistMap(0, 109 );

  for (auto const &p: dilateSet)
    {
      // to be working we need to remove the distance map init to 0... (else it is erased by tracking all cC)
      viewerStep1 <<  DGtal::CustomColors3D(colorDistMap(gg.getDistanceImage()(p)), colorDistMap(gg.getDistanceImage()(p)));
      
      
      resDilate << p[0] << " " << p[1] << " " << p[2] << std::endl;
      viewerStep1 << p;
      //viewerStep5<< DGtal::CustomColors3D(DGtal::Color(255, 225, 225, 15),DGtal::Color(225, 225, 225, 15));
      //viewerStep5<< p;
    }
  resDilate.close();

  // display labels:
  const GeodesicGraphComputer::Image3D &imageLabel = gg.getLabelImage();
  for(unsigned int i = 1; i<gg.myGraph.nbVertex(); i++){
    for (auto const &p: imageLabel.domain())
    {
      if(gg.myOrderGeodesImage(p)==i){
        if(i%2==0){
          viewerStep2<< DGtal::CustomColors3D(DGtal::Color(255, 255, 255, 255),DGtal::Color(255, 255, 255, 255));
          viewerStep4<< DGtal::CustomColors3D(DGtal::Color(255, 225, 225, 15),DGtal::Color(225, 225, 225, 25));
          viewerStep4<< p;
          
        }else {
          viewerStep2<< DGtal::CustomColors3D(DGtal::Color(155, 125, 125, 255),DGtal::Color(25, 25, 125, 255));
          viewerStep4<< DGtal::CustomColors3D(DGtal::Color(155, 125, 125, 15),DGtal::Color(25, 25, 125, 25));
          viewerStep4<< p;
          
        }
        viewerStep2 << p;
      }
    }
    
  }

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
            viewerStep3<< DGtal::CustomColors3D(DGtal::Color(255, 25, 25, 255),DGtal::Color(255, 25, 25, 255));
          }
        else if(imageDegree(p) == 2)
          {
            viewerStep3<< DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 255),DGtal::Color(25, 25, 255, 255));
          }
        else if(imageDegree(p) == 3)
          {
            viewerStep3<< DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 255),DGtal::Color(25, 25, 255, 255));
          }
        else if(imageDegree(p) == 4)
          {
            viewerStep3<< DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 255),DGtal::Color(25, 25, 255, 255));
          }
        else if(imageDegree(p) == 5)
          {
            viewerStep3<< DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 255),DGtal::Color(25, 25, 255, 255));
          }
        else{
          viewerStep3<< DGtal::CustomColors3D(DGtal::Color(25, 25, 25, 255),DGtal::Color(25, 25, 25, 255));
 
        }
        vectVoxelsDisplayed.push_back(p);
        viewerStep3 << DGtal::SetName3D(name++) << p;
      }
    }

  graphDegree.close();
  VoxelsAndGeodesicGraph data;
  data.vectVoxels = &vectVoxelsDisplayed;
  data.geoGraph = &gg;
  viewerStep0 << DGtal::CustomColors3D(DGtal::Color(25, 25, 255, 25),DGtal::Color(25, 25, 255, 25));
  viewerStep0 << myMesh;
  viewerStep3 << DGtal::SetSelectCallback3D(displayVoxelInfos, &data, 0, vectVoxelsDisplayed.size());
  viewerStep3 << DGtal::CustomColors3D(DGtal::Color::Black, DGtal::Color::Black);
  viewerStep3 << Z3i::Point(531, 425, 302);

  viewerStep0<< MyViewer::updateDisplay;
  viewerStep1<< MyViewer::updateDisplay;
  viewerStep2<< MyViewer::updateDisplay;
  viewerStep3 << MyViewer::updateDisplay;
  viewerStep4 << MyViewer::updateDisplay;
  viewerStep5 << MyViewer::updateDisplay;

  viewerStep0.setWindowTitle(QString("Step 0"));
  viewerStep0.show();
  viewerStep1.setWindowTitle(QString("Step 1"));
  viewerStep1.show();
  viewerStep2.setWindowTitle(QString("Step 2"));
  viewerStep2.show();
  viewerStep3.setWindowTitle(QString("Step 3"));
  viewerStep3.show();
  viewerStep4.setWindowTitle(QString("Step 4"));
  viewerStep4.show();
  viewerStep5.setWindowTitle(QString("Step 5"));
  viewerStep5.show();

  return app.exec();

}






