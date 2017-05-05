#include <fstream>
#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/readers/PointListReader.h"
#include <DGtal/io/readers/TableReader.h>



#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <iostream>
#include <fstream>
#include <DGtal/math/linalg/SimpleMatrix.h>
#include <DGtal/math/linalg/EigenDecomposition.h>


using namespace DGtal;
using namespace std;
namespace po = boost::program_options;

typedef ImageContainerBySTLVector<Z3i::Domain, unsigned int> Image3D;
typedef ImageContainerBySTLVector<Z3i::Domain, Z3i::RealPoint> ImageVector;
typedef PointVector<9, double> Matrix3x3Point;
typedef SimpleMatrix<double, 3, 3 > CoVarianceMat;



 /**
   * Order natural (book reading from upper left to lower right)
   **/
  static CoVarianceMat
  getCoVarianceMatFrom(const Matrix3x3Point &aMatrixPt){
    CoVarianceMat res;
    for(unsigned int i = 0; i<3; i++){
      for(unsigned int j = 0; j<3; j++){
        res.setComponent(i, j, aMatrixPt[j+i*3] );
      }
    }
    return res;
  }



DGtal::Z3i::RealPoint getMainDirsCoVar(const std::vector<DGtal::Z3i::RealPoint> &accVect)
{
  Matrix3x3Point c ;
  for(const auto & v: accVect)
    {
      c[0] += (v[0]*v[0]);
      c[1] += (v[0]*v[1]);
      c[2] += (v[0]*v[2]);
      
      c[3] += (v[1]*v[0]);
      c[4] += (v[1]*v[1]);
      c[5] += (v[1]*v[2]);
         
      c[6] += (v[2]*v[0]);
      c[7] += (v[2]*v[1]);
      c[8] += (v[2]*v[2]);  
    }
  c = c/accVect.size();
  CoVarianceMat covar = getCoVarianceMatFrom(c);
  SimpleMatrix<double, 3, 3 > eVects;
  PointVector<3, double> eVals;
  DGtal::EigenDecomposition<3, double, CoVarianceMat>::getEigenDecomposition (covar, eVects, eVals);
  return eVects.column(0);
}




DGtal::Z3i::RealPoint getMainDirsAlgoLast(const std::vector<DGtal::Z3i::RealPoint> &accVect,
                                          double minSinAngle)
{
  DGtal::Z3i::RealPoint lastDir;
  DGtal::Z3i::RealPoint mainDir (0, 0, 0);
  bool first = true;
  unsigned int num = 0;
  for(const auto & v: accVect)
    {
      if(!first)
        {
          DGtal::Z3i::RealPoint d = lastDir.crossProduct(v.getNormalized());
          if(mainDir != DGtal::Z3i::RealPoint(0,0,0) && d.dot(mainDir)<0){
            d *=-1.0;
          }
          if(d.norm()>minSinAngle){
            mainDir += d;
            num++;
          }
        }
      else
        {
          first = false;
        }
      lastDir = v.getNormalized();
    }

  if (num !=0)
    {
      return mainDir.getNormalized();      
    }
  return mainDir;
}

  


int main(int argc, char **argv)

{

  // parse command line -------------------------------------------------------
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h", "display this message")
    ("input,i", po::value<std::string >(), "an input file representing the accumulating vectors which contribute in a voxel (format X Y Z V0x V0y V0z V1x V1y V1z ... . " )
    ("scale,s", po::value<double >()->default_value(2.0), "set the scale of the resultng vectors (ie. length of vector). " )
    ("useCovariance", "use covariance analysis")
    ("minAngle",po::value<double >()->default_value(0.01),  "minimal angle to consider two directions equals (used when useCovariance is not selected).")
    ("output,o", po::value<std::string>()->default_value("outMainDir.sdp"), "the resulting main dir vectors " );
  
  bool parseOK=true;
  po::variables_map vm;
  try
    {
      po::store(po::parse_command_line(argc, argv, general_opt), vm);
    }
  catch(const std::exception& ex)
    {
      parseOK=false;
      trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
  

  // check if min arguments are given and tools description ------------------
  po::notify(vm);
  if( !parseOK || vm.count("help")||argc<=1)
    {
      std::cout << "Usage: " << argv[0] << " [input]\n"
                << "It computes the main tube direction from the vectors contributing in the accumulation/confidence. \n"
                << "Two methods are proposed from the algorithm of accumulation and from the covariance matrix (with option useCovariance). \n"
                << general_opt << "\n"
                << "Typical use example:\n \t compMainDirs -i contribVect.sdp    -o vFieldCoVariance.sdp --useCovariance \n";
      return 0;
    }  
  if(! vm.count("input"))
    {
      trace.error() << " The file name was not defined" << endl;
      return 1;
    }


  //  recover the  args ----------------------------------------------------
  string inputFileName = vm["input"].as<string>();
  string outputFileName =  vm["output"].as<string>();
  double minSinAngle = sin(vm["minAngle"].as<double>());
                           
  // Reading input file
    trace.info() << "Reading input vectors...";
  std::vector<std::vector<double > > vectAccDirs =
    TableReader<double>::getLinesElementsFromFile(inputFileName);
  trace.info()<< "[Done]"<<std::endl;

  struct VectDir{Z3i::RealPoint pt; Z3i::RealPoint mainDir;};
  std::vector<VectDir> result;
  // construct reference point and its set of main dirs
  for(const auto &line: vectAccDirs)
    {
      if(line.size()>6)
        {
          VectDir vd; 
          Z3i::RealPoint p(line[0], line[1], line[2]);
          vd.pt = p;
          std::vector<Z3i::RealPoint> listV;
          for (unsigned int i = 3; i+2< line.size(); i+=3)
            {
              listV.push_back(Z3i::RealPoint(line[i], line[i+1], line[i+2]));
            }
          Z3i::RealPoint md;
          if(vm.count("useCovariance"))
            {
              md = getMainDirsCoVar(listV);
            }
          else{
            md = getMainDirsAlgoLast(listV, minSinAngle);
          }
          if(md != Z3i::RealPoint(0,0,0))
            {
              vd.mainDir = md;
              result.push_back(vd);
            }
        } 
    }

  
  ofstream fout;
  double scale = vm["scale"].as<double>();
  fout.open(outputFileName.c_str(), ofstream::out|ofstream::binary);
  for(const auto  &vd: result)
    {
      fout << vd.pt[0] << " " << vd.pt[1] << " " << vd.pt[2] << " "
           << vd.pt[0]+vd.mainDir[0]*scale <<  " " << vd.pt[1]+vd.mainDir[1]*scale << " " << vd.pt[2]+vd.mainDir[2]*scale << std::endl;       
    }
  fout.close();
  return 0;
}



