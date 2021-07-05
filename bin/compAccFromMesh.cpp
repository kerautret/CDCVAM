#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>


#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/io/writers/LongvolWriter.h>

#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/writers/MeshWriter.h>

#include <DGtal/images/ImageContainerBySTLVector.h>

#include "CLI11.hpp"

#include "NormalAccumulator.h"


using namespace std;
using namespace DGtal;

typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain,  double> ImageDouble;

typedef Z3i::Domain::Predicate DomainPredicate;




/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  string inputFile;
  string outputFileAcc {"accumulation.longvol"};
  string outputFileConf {"confidence.longvol"};
  string outputFileRad {"radius.longvol"};
  string outputAccVectors {""};
  string outputConfVectors {""};
  string autoScaleAcc {""};
  string autoScaleConf {""};
  string autoScaleRad {""};

  unsigned int thAcc {50};
  double thConf {0.75} ;
  DGtal::uint64_t  maxValOutConf {255};
  DGtal::uint64_t  maxValOutRad {255};

  double radius {5.0};
  string typeStat {"min"};
  bool invertNormal {false};
  
  
  // parse command line using CLI ----------------------------------------------
     CLI::App app;
  app.description("Compute mesh accumulation from a mesh and compute the radius (median value) "
                  "of all faces participating to the accumulation" );
  app.add_option("-i,--input,1", inputFile, "input mesh." )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--outputAcc,-o,2",outputFileAcc,"output accumulation longvol file.", true);
  app.add_option("--outputRad,-R",outputFileConf,"output radius vol file.", true);
  app.add_option("--outputConf",outputFileConf,"output confidence longvol file.", true);

  app.add_option("--radius,-r", radius, "radius of accumulation analysis.", true);
  app.add_flag("--invertNormal", invertNormal, "export the accumulation vectors.");
  app.add_option("--radiusEstimator,-e",typeStat,  "use: {min (default), max, mean, median} to estimate the radius", true)
  -> check(CLI::IsMember({"max", "min", "mean", "median"}));
  app.add_option("--outputAccVectors", outputAccVectors, "export the accumulation vectors.", true);
  app.add_option("--outputConfVectors", outputConfVectors,"export the confidence vectors.", true );
  app.add_option("--maxThAccVectors", thAcc, "threshold the value of accumulations to export the accumulations vectors (used with outputAccVectors) .", true);
  app.add_option("--maxThConfVectors",thConf, "threshold the value of confidence to export the accumulations vectors (used with outputConfVectors) .", true);
  app.add_option("--maxValOutConf", maxValOutConf,"set MAX scale of confidence out image: 0 .. 1 -> 0 ..  MAX.", true);
  app.add_option("--maxValOutRad", maxValOutRad,"set MAX scale of radius out image: 0 .. 1 -> 0 ..  MAX.", true);
  app.add_option("--autoScaleAcc", autoScaleAcc,"auto scale out accumulation image values between 0 and 255.", true);
  app.add_option("--autoScaleConf", autoScaleConf,"auto scale out condidence image values between 0 and 255.", true);
  app.add_option("--autoScaleRad", autoScaleRad,"auto scale out condidence image values between 0 and 255.", true);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

 
  // 1) Reading input mesh
  Mesh<Z3i::RealPoint> aMesh(true);
  aMesh << inputFile;
    
  // 2) Init an NormalAccumulator
  NormalAccumulator normAcc(radius, typeStat);
  normAcc.initFromMesh(aMesh, invertNormal);

  // 3) Generate accumulation image
  normAcc.computeAccumulation();
  NormalAccumulator::Image3D &imageAccumulation = normAcc.getAccumulationImage();
  
  trace.info() << "Saving accumulation image in " << outputFileAcc;
  LongvolWriter<Image3D>::exportLongvol(outputFileAcc, imageAccumulation);
  trace.info() << "[done]" << std::endl;

  if (autoScaleAcc != ""){
    trace.info() << "Saving accumulation (auto scale 0 " << normAcc.getMaxAccumulation() 
                 <<"  -> [0, 255]) image in " << outputFileAcc << " ... ";
    typedef functors::Rescaling<DGtal::uint64_t, unsigned char> ScaleFct;
    ScaleFct  scaleFct (0.0 ,normAcc.getMaxAccumulation(), 0, 255);
    VolWriter<Image3D,ScaleFct>::exportVol(autoScaleAcc, imageAccumulation, true, scaleFct);
    trace.info() << "[done]" << std::endl;
  }

  // 3bis) Optionally export the vectors of the accumulation
 if(outputAccVectors != "")
   {
     std::string outName = outputAccVectors;
     ofstream fout; fout.open(outName.c_str(), ofstream::binary|ofstream::out);
     for(const auto &p: imageAccumulation.domain())
       {
         if(imageAccumulation(p)>thAcc)
           {
             NormalAccumulator::PointContainer setPt = normAcc.getAssociatedPoints(p);
             if(setPt.size()!=0){
               fout << p[0] << " " << p[1] << " " << p[2] << " " ;
               for(const auto &n: setPt)
                 {
                   Z3i::RealPoint v = p-n;
                   fout << v[0] << " " << v[1] << " " << v[2] << " "; 
                 }
               fout << std::endl;
             }
           }
       }
     fout.close();
   }
     
  
  // 4) Compute confidence image 
  normAcc.computeConfidence();
  ImageDouble imageConfidence = normAcc.getConfidenceImage();
  trace.info() << "Saving confidence image in " << outputFileConf << " ... ";
  typedef functors::Rescaling<double, DGtal::uint64_t> ScaleFctD;
  ScaleFctD scaleFct(0.0, 1.0, 0, maxValOutConf);
  LongvolWriter<ImageDouble, ScaleFctD>::exportLongvol(outputFileConf, imageConfidence, true, scaleFct);
  trace.info() << "[done]" << std::endl;
  
  if(autoScaleConf != "")
    {
      std::string outNameAutoConfidence = autoScaleConf;
      trace.info() << "Saving confidence (auto scale 0 1 -> 0 255) image in " << outNameAutoConfidence << " ... ";
      typedef functors::Rescaling<double, unsigned char> ScaleFctD;
      ScaleFctD  scaleFct (0.0 ,1.0, 0, 255);
      VolWriter<ImageDouble,ScaleFctD>::exportVol(outNameAutoConfidence, imageConfidence, true, scaleFct);
      trace.info() << "[done]" << std::endl;    
    }
  

  // 5) Compute Radius image
  normAcc.computeRadiusFromConfidence();
  ImageDouble imageRadius = normAcc.getRadiusImage();  
  trace.info() << "Saving radius image in " << outputFileRad << " ... ";
  ScaleFctD scaleFct2(0.0, normAcc.getMaxRadius(), 0, maxValOutRad);
  LongvolWriter<ImageDouble, ScaleFctD>::exportLongvol(outputFileRad, imageRadius, true, scaleFct2);
  trace.info() << "[done]" << std::endl;
 if(autoScaleRad != "")
    {
      std::string outNameAutoRadius = autoScaleRad;
      trace.info() << "Saving confidence (auto scale 0 " << normAcc.getMaxRadius()
                   << " ->[0,255]) image in " << outNameAutoRadius << " ... ";
      typedef functors::Rescaling<double, unsigned char> ScaleFctD;
      ScaleFctD  scaleFct (0.0 ,normAcc.getMaxRadius(), 0, 255);
      VolWriter<ImageDouble,ScaleFctD>::exportVol(outNameAutoRadius, imageRadius, true, scaleFct);
      trace.info() << "[done]" << std::endl;    
    }    


 // 5 bis) Optionally export the vectors of the confidence
 if(outputConfVectors != "")
   {
     std::string outName = outputConfVectors;
     ofstream fout; fout.open(outName.c_str(), ofstream::binary|ofstream::out);
     for(const auto &p: imageAccumulation.domain())
       {
         if(imageConfidence(p)>thConf)
           {
             NormalAccumulator::PointContainer setPt = normAcc.getAssociatedPoints(p);
             if(setPt.size()!=0){
               fout << p[0] << " " << p[1] << " " << p[2] << " " ;
               for(const auto &n: setPt)
                 {
                   Z3i::RealPoint v = p-n;
                   fout << v[0] << " " << v[1] << " " << v[2] << " "; 
                 }
               fout << std::endl;
             }

           }
       }
     fout.close();
   }
   
 
 return 0;
}

