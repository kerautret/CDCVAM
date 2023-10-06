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
  string outputAccOrigins {"exportedAcc.dat"};
  unsigned int thAcc {50};
  double thConf {0.75};
  double radius {5.0};
  bool invertNormal {false};
  
  
  // parse command line using CLI ----------------------------------------------
  CLI::App app;
  app.description("Compute the accumulation origins from a mesh and export in dat format : X Y Z NbAcc NbAccM Id1 Id2 ... IdNbAccM ... NbAcc  "
                  "of all faces participating to the accumulation" );
  app.add_option("-i,--input,1", inputFile, "input mesh." )
      ->required()
      ->check(CLI::ExistingFile);
  app.add_option("--output,-o,2",outputAccOrigins,"output accumulation origins.", true);
 
  app.add_option("--radius,-r", radius, "radius of accumulation analysis.", true);
  app.add_flag("--invertNormal", invertNormal, "export the accumulation vectors.");
  app.add_option("--maxThAccVectors", thAcc, "threshold the value of accumulations used to decide to export the accumulations vectors (used with outputAccVectors) .", true);
  app.add_option("--maxThConfVectors",thConf, "threshold the value of confidence to export the accumulations vectors (used with outputConfVectors) .", true);
  
  app.get_formatter()->column_width(40);
  CLI11_PARSE(app, argc, argv);
  // END parse command line using CLI ----------------------------------------------

 
  // 1) Reading input mesh
  Mesh<Z3i::RealPoint> aMesh(true);
  aMesh << inputFile;
    
  // 2) Init an NormalAccumulator
  NormalAccumulator normAcc(radius, "mean");
  normAcc.initFromMesh(aMesh, invertNormal);

  // 3) Generate accumulation image
  normAcc.computeAccumulation();
  NormalAccumulator::Image3D &imageAccumulation = normAcc.getAccumulationImage();

 if(outputAccOrigins != "")
   {
     ofstream fout; fout.open(outputAccOrigins.c_str(), ofstream::binary|ofstream::out);
     for(const auto &p: imageAccumulation.domain())
       {
         if(imageAccumulation(p)>thAcc)
           {
             NormalAccumulator::PointIndexContainer setIndexPt = normAcc.getAssociatedIndexPoints(p);
             if(setIndexPt.size()!=0){
               fout << p[0] << " " << p[1] << " " << p[2] << " " << imageAccumulation(p) << " "  ;
               for(const auto &i: setIndexPt)
                 {
                   fout << i << " ";
                 }
               fout << std::endl;
             }
           }
       }
     fout.close();
   }
     
  
  
 return 0;
}

