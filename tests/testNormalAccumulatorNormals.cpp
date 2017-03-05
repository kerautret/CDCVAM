#include <iostream>
#include <fstream>


#include <DGtal/io/readers/PointListReader.h>
#include "DGtal/io/writers/LongvolWriter.h"
#include <DGtal/io/readers/MeshReader.h>

#include "NormalAccumulator.h"

using namespace DGtal;
using namespace std;


int
main(int argc,char **argv)
{
  DGtal::Mesh<DGtal::Z3i::RealPoint>  aMesh;
  aMesh << "../Samples/tube1.off";

  std::vector<Z3i::RealPoint> normalsVect = PointListReader<Z3i::RealPoint>::getPointsFromFile("../Samples/tube1.normals");

  
  NormalAccumulator acc(8);
  acc.initFromMeshAndNormals(aMesh, normalsVect);
  acc.computeAccumulation();
  DGtal::LongvolWriter<NormalAccumulator::Image3D>::exportLongvol("testAccumulation.longvol",
                                                                  acc.getAccumulationImage());


  DGtal::trace.info() << acc;
  acc.computeConfidence();
  typedef DGtal::functors::Rescaling< double, DGtal::uint64_t> ScaleFctD;

  ScaleFctD  confidencescale (0 , 1.0, 0, 1000);
  DGtal::LongvolWriter<NormalAccumulator::Image3DDouble,ScaleFctD>::exportLongvol("testConfidence.longvol",
                                                                                  acc.getConfidenceImage(),
                                                                                  true, confidencescale);

  acc.computeRadiusFromOrigins();
  double maxRadius = acc.getMaxRadius();
  ScaleFctD  radiiScale (0 , maxRadius, 0, maxRadius * 1000);
  DGtal::LongvolWriter<NormalAccumulator::Image3DDouble,ScaleFctD>::exportLongvol("testRadius.longvol",
                                                                                  acc.getRadiusImage(),
                                                                                  true, radiiScale);
  
  DGtal::trace.info()<< "processing test accumulator 2" << std::endl;
  std::vector<DGtal::Z3i::RealPoint> listPts = PointListReader<Z3i::RealPoint>::getPointsFromFile("../Samples/sectionS0.1.sdp");
  std::vector<DGtal::Z3i::RealPoint> listNormals = PointListReader<Z3i::RealPoint>::getPointsFromFile("../Samples/sectionS0.1.normals");
  NormalAccumulator acc2(6);
  acc2.initFromNormals(listPts, listNormals);
  acc2.computeAccumulation();
  DGtal::trace.info() << "val max:" << acc2.getMaxAccumulation() << std::endl;
  DGtal::LongvolWriter<NormalAccumulator::Image3D>::exportLongvol("testAccumulation2.longvol",
                                                                  acc2.getAccumulationImage());

 

  
}






