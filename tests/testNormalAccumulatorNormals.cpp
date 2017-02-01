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
  aMesh << "/Users/kerautre/EnCours/TubeAnalyse/ExpeICPR/Discrete/tube1Vox300.off";

  std::vector<Z3i::RealPoint> normalsVect = PointListReader<Z3i::RealPoint>::getPointsFromFile("/Users/kerautre/EnCours/TubeAnalyse/ExpeICPR/Discrete/tube1Vox300.normals");
  
  NormalAccumulator acc(12);
  acc.initFromMeshAndNormals(aMesh, normalsVect);

  acc.computeAccumulation();
  DGtal::LongvolWriter<NormalAccumulator::Image3D>::exportLongvol("testAccumulation.lvol",
                                                                  acc.getAccumulationImage());

  DGtal::trace.info() << acc;
  acc.computeConfidence();
  typedef DGtal::functors::Rescaling< double, DGtal::uint64_t> ScaleFctD;

  ScaleFctD  confidencescale (0 , 1.0, 0, 1000);
  DGtal::LongvolWriter<NormalAccumulator::Image3DDouble,ScaleFctD>::exportLongvol("testConfidence.lvol",
                                                                                  acc.getConfidenceImage(),
                                                                                  true, confidencescale);

  acc.computeRadiusFromOrigins();
  double maxRadius = acc.getMaxRadius();
  ScaleFctD  radiiScale (0 , maxRadius, 0, maxRadius * 1000);
  DGtal::LongvolWriter<NormalAccumulator::Image3DDouble,ScaleFctD>::exportLongvol("testRadius.lvol",
                                                                                  acc.getRadiusImage(),
                                                                                  true, radiiScale);
  
  

  
}






