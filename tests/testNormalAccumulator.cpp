#include <iostream>
#include <fstream>

#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/writers/LongvolWriter.h"


#include "NormalAccumulator.h"

int
main(int argc,char **argv)
{

  DGtal::Mesh<DGtal::Z3i::RealPoint>  aMesh;
  aMesh << "../Samples2/bk7.off";
  NormalAccumulator acc(12);
  acc.initFromMesh(aMesh);

  acc.computeAccumulation();
  DGtal::LongvolWriter<NormalAccumulator::Image3D>::exportLongvol("testAccumulation.lvol",
                                                                  acc.getAccumulationImage());

  DGtal::trace.info() << acc;
  acc.computeConfidence();
  typedef DGtal::functors::Rescaling< double, DGtal::uint64_t> ScaleFctD;

  ScaleFctD  confidencescale (0 , 1.0, 0, 1000);
  DGtal::LongvolWriter<NormalAccumulator::Image3DDouble,ScaleFctD>::exportLongvol("testConfidence.lvol",
                                                                                  acc.getConfidenceImage(),
                                                                                  confidencescale);

  acc.computeRadiusFromOrigins();
  double maxRadius = acc.getMaxRadius();
  ScaleFctD  radiiScale (0 , maxRadius, 0, maxRadius * 1000);
  DGtal::LongvolWriter<NormalAccumulator::Image3DDouble,ScaleFctD>::exportLongvol("testRadius.lvol",
                                                                                  acc.getRadiusImage(),
                                                                                  radiiScale);
  
  

  
}






