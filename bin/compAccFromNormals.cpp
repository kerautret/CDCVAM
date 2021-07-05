#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>


#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include <DGtal/io/readers/GenericReader.h>
#include <DGtal/io/readers/PointListReader.h>


#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/io/writers/VolWriter.h>
#include <DGtal/io/writers/LongvolWriter.h>


#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/writers/MeshWriter.h>
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/images/ImageContainerBySTLMap.h>


#include "CLI11.hpp"


#include "AccumulatorHelper.h"
#include "NormalAccumulator.h"


using namespace std;
using namespace DGtal;


typedef ImageContainerBySTLVector<Z3i::Domain, DGtal::uint64_t> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain,  double> ImageDouble;


void checkFileFormat(const std::string &filename,
                     const std::string &reqExt,
                     const std::string &warningInstr)
{
   std::string extension = filename.substr( filename.find_last_of(".") + 1 );
  if (extension != reqExt)
  {
    trace.warning() << "You are specifying a wrong extension (" << extension << ") for you custom file name "<< filename << ", " << warningInstr << std::endl;  
  }
}



/**
 * @brief main function call
 *
 */
int main(int argc, char *const *argv)
{
  
  // parse command line using CLI ----------------------------------------------
     CLI::App app;
     std::string inputFileName;
    
  string inputFile;
  string outputFileAcc {"accumulation.longvol"} ;
  string outputFileConf {"confidence.longvol"};
  string outputFileRad {"radius.longvol"};
  string outputAccVectors {""};
  string outputConfVectors {""};
  string autoScaleRadius {""};
  string autoScaleConf {""};
  string autoScaleAcc {""};
  string exportNormalSegments {""};
  string exportNormals {""};
  bool outputPointAssociations {false};
  
  double radius {5.0};
  string typeStat {"min"};
  DGtal::uint64_t outConfidenceMax {255};
  DGtal::uint64_t outRadiusMax {255};
  DGtal::uint64_t minAccumulation {255};
  unsigned int thAcc {50};//= vm["maxThAccVectors"].as<unsigned int>();
  double thConf {0.75}; //vm["maxThConfVectors"].as<double>();

  
  // parse command line using CLI ----------------------------------------------
  app.description("Compute accumulation from point of cloud and compute the radius (median value) "
                  "of all normal participating to the accumulation");
  app.add_option("-i,--input,1", inputFile, "input sdp points (x y z format)." )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--outputAcc,-o,2", outputFileAcc, "output accumulation longvol file (to export in vol you can use the option --autoScaleAcc to autoscale).", true);
  app.add_option("--outputConf,-c", outputFileConf, "output confidence longvol file.", true);
  app.add_option("--outputRad,-R", outputFileRad, "output radius vol file.", true);
  app.add_option("--radius,-r", radius, "radius of accumulation analysis.", true);
  app.add_option("--radiusEstimator,-e",typeStat,  "use: {min (default), max, mean, median} to estimate the radius", true)
  -> check(CLI::IsMember({"max", "min", "mean", "median"}));
  app.add_option("--maxValOutConf",outConfidenceMax, "set MAX scale of confidence out image: 0 .. 1 -> 0 ..  MAX.", true);
  app.add_option("--maxValOutRad",outRadiusMax, "set MAX scale of radius out image: 0 .. 1 -> 0 ..  MAX.", true);
  
  app.add_option("--minAccumulation",outRadiusMax, "min value of the accumumation to consider the computation of the confidence.", true);
  app.add_option("--outputAccVectors",outputAccVectors, "export the accumulation vectors.");
  app.add_option("--outputConfVectors",outputAccVectors, "export the confidence vectors.");
  app.add_flag("--outputPointAssociations", outputPointAssociations, "output point associations instead vectors (used with outputAccVectors and outputConfVectors options). Each pair point of the association are exported sequentially.");
  
  app.add_option("--maxThAccVectors",thAcc, "threshold the value of accumulations to export the accumulations vectors (used with outputAccVectors).", true);
  app.add_option("--maxThConfVectors",thConf, "threshold the value of confidence to export the accumulations vectors (used with outputConfVectors) .", true);
    
  app.add_option("--autoScaleRadius",autoScaleRadius, "auto scale out radius image values between 0 and 255.");
  app.add_option("--autoScaleConf",autoScaleConf, "auto scale out confidence image values between 0 and 255.");
  app.add_option("--autoScaleAcc",autoScaleAcc, "auto scale out accumulation image values between 0 and 255.");
  app.add_option("--exportNormalSegments", exportNormalSegments, "export normal segments");
  app.add_option("--exportNormals", exportNormals, "export normals");

  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------


  
 
  

  // 1) Reading input point and normals

  DGtal::trace.info() << "Reading imput points and normals ... ";
  std::vector<unsigned int> vectIndPoint1 = {0, 1, 2};
  std::vector<Z3i::RealPoint> setOfPt1 = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFile, vectIndPoint1);
  std::vector<unsigned int> vectIndNorm = {6, 7, 8};
  std::vector<Z3i::RealPoint> setOfNormals = PointListReader<Z3i::RealPoint>::getPointsFromFile(inputFile, vectIndNorm);

  DGtal::trace.info() << "[done]" << std::endl;


  
  // 2) create and init an NormalAccumulator
  NormalAccumulator normAcc(radius, typeStat);
  normAcc.initFromNormals(setOfPt1, setOfNormals);  
  
  
  
  // 3) check if export
  if(exportNormalSegments != ""){
    std::ofstream outNormals;
    string normalSegFileName = exportNormalSegments;
    outNormals.open(normalSegFileName);
    double scaleV = 2.0;
    std::vector<Z3i::RealPoint> vPts = normAcc.getNormalOrigins();
    std::vector<Z3i::RealPoint> vNormals = normAcc.getNormalField();
    for (unsigned int i = 0; i < vPts.size() ; i++) {
      outNormals << vPts[i][0] << " " << vPts[i][1] << " "<< vPts[i][2] << std::endl;
      outNormals <<vPts[i][0]+vNormals[i][0]*scaleV << " "<<vPts[i][1]+vNormals[i][1]*scaleV  
                 << " "<< vPts[i][2]+vNormals[i][2]*scaleV << std::endl;
    }  
    outNormals.close();
  }

  if(exportNormals != ""){
    std::ofstream outNormals;
    string normalSegFileName = exportNormals;
    outNormals.open(normalSegFileName);
    std::vector<Z3i::RealPoint> vNormals = normAcc.getNormalField();
    for (unsigned int i = 0; i < vNormals.size() ; i++) {
      Z3i::RealPoint n = vNormals[i].getNormalized();
      outNormals << n[0] << " " << n[1] << " "<< n[2] << std::endl;
    }  
    outNormals.close();
  }

  

  
  // 3) Generate accumulation image
  normAcc.computeAccumulation();
  NormalAccumulator::Image3D &imageAccumulation = normAcc.getAccumulationImage();

  checkFileFormat(outputFileAcc,"longvol", "please change it to .longvol to avoid any problem (or use the option --autoScaleAcc to autoscale and to save in vol format) ");

  trace.info() << "Saving accumulation image in " << outputFileAcc;
  LongvolWriter<Image3D>::exportLongvol(outputFileAcc, imageAccumulation);
  trace.info() << " [done]" << std::endl;
  
  if (autoScaleAcc != ""){
    std::string outputFileAcc = autoScaleAcc;
    checkFileFormat(outputFileAcc,"vol", "please change it to .vol to avoid any problem.");
    trace.info() << "Saving accumulation (auto scale 0 " << normAcc.getMaxAccumulation() 
                 <<"  -> [0, 255]) image in " << outputFileAcc << " ... ";
    typedef functors::Rescaling<DGtal::uint64_t, unsigned char> ScaleFct;
    ScaleFct  scaleFct (0.0 ,normAcc.getMaxAccumulation(), 0, 255);
    VolWriter<Image3D,ScaleFct>::exportVol(outputFileAcc, imageAccumulation, true,  scaleFct);
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
               if(outputPointAssociations)
                 {
                   fout  << "# associations of point: (" << p[0] << " " << p[1] << " "<< p[2]<< std::endl;
                   NormalAccumulator::PointContainer neighbors = normAcc.getAssociatedPoints(p);
                   for(const auto & pa: neighbors)
                     {
                       fout << p[0] << " " << p[1] << " " << p[2] << std::endl;
                       fout << pa[0] << " " << pa[1] << " " << pa[2] << std::endl;                 
                     }
                 }
               else
                 {
                 
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
       }
     fout.close();
   }

  // 4) Compute confidence image 
  normAcc.computeConfidence(false, minAccumulation);
  ImageDouble imageConfidance = normAcc.getConfidenceImage();
  checkFileFormat(outputFileConf,"longvol", "please change it to .longvol to avoid any problem (or use the option --autoScaleConf to autoscale and to save in vol format)");


  trace.info() << "Saving confidence image in " << outputFileConf << " ... ";
  typedef functors::Rescaling<double, DGtal::uint64_t> ScaleFctD;
  ScaleFctD scaleFct(0.0, 1.0, 0, outConfidenceMax);
  LongvolWriter<ImageDouble, ScaleFctD>::exportLongvol(outputFileConf, imageConfidance, true, scaleFct);
  trace.info() << "[done]" << std::endl;
  
  if(autoScaleConf != "")
    {
      std::string outNameAutoConf = autoScaleConf;
      checkFileFormat(outNameAutoConf,"vol", "please change it to .vol to avoid any problem.");
      trace.info() << "Saving confidence (auto scale 0 1 -> [0, 255]) image in " << outNameAutoConf << " ... ";
      typedef functors::Rescaling<double, unsigned char> ScaleFctD;
      ScaleFctD  scaleFct (0.0 ,1.0, 0, 255);
      VolWriter<ImageDouble,ScaleFctD>::exportVol(outNameAutoConf, imageConfidance, true, scaleFct);
      trace.info() << "[done]" << std::endl;    
    }
  

  // 5) Compute Radius image
  normAcc.computeRadiusFromConfidence();
  ImageDouble imageRadius = normAcc.getRadiusImage();
  checkFileFormat(outputFileRad,"longvol", "please change it to .longvol to avoid any problem (or use the option --autoScaleRadius to autoscale and to save in vol format) ");
  trace.info() << "Saving radius image in " << outputFileRad << " ... ";
  ScaleFctD scaleFct2(0.0, normAcc.getMaxRadius(), 0, outRadiusMax);
  LongvolWriter<ImageDouble, ScaleFctD>::exportLongvol(outputFileRad, imageRadius, true, scaleFct2);
  trace.info() << "[done]" << std::endl;
  if(autoScaleRadius != "")
    {
      std::string outNameAutoRad = autoScaleRadius;
      checkFileFormat(outNameAutoRad,"vol", "please change it to .vol to avoid any problem.");
      trace.info() << "Saving confidence (auto scale 0 " << normAcc.getMaxRadius()
                   << " -> 0 255) image in " << outNameAutoRad << " ... ";
      typedef functors::Rescaling<double, unsigned char> ScaleFctD;
      ScaleFctD  scaleFct (0.0 ,normAcc.getMaxRadius(), 0, 255);
      VolWriter<ImageDouble,ScaleFctD>::exportVol(outNameAutoRad, imageRadius, true, scaleFct);
      trace.info() << "[done]" << std::endl;    
    }      


 // 5 bis) Optionally export the vectors of the confidence
 if(outputConfVectors != "")
   {
     std::string outName = outputConfVectors;
     ofstream fout; fout.open(outName.c_str(), ofstream::binary|ofstream::out);
     for(const auto &p: imageAccumulation.domain())
       {
         if(imageConfidance(p)>thConf)
           {
             NormalAccumulator::PointContainer setPt = normAcc.getAssociatedPoints(p);
             if(setPt.size()!=0){
               if(outputPointAssociations)
                 {
                   fout << "# associations of point: (" << p[0] << " " << p[1] << " "<< p[2]<< std::endl;
                   NormalAccumulator::PointContainer neighbors = normAcc.getAssociatedPoints(p);
                   for(const auto & pa: neighbors)
                     {
                       fout << p[0] << " " << p[1] << " " << p[2] << std::endl;
                       fout << pa[0] << " " << pa[1] << " " << pa[2] << std::endl;                 
                     }
                 }
               else
                 {
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
       }
     fout.close();
   }
 return 0;
}

