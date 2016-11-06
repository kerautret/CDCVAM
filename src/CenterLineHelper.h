#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Mesh.h"
#include "DGtal/io/boards/Board2D.h"
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/kernel/BasicPointFunctors.h>
#include <DGtal/shapes/implicit/ImplicitBall.h>
#include "DGtal/geometry/curves/AlphaThickSegmentComputer.h"
#include <DGtal/shapes/EuclideanShapesDecorator.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/math/Statistic.h>

#include "NormalAccumulator.h"


class CenterLineHelper {


public:  







  /**
   * Optimize center line  according attached  vertex
   *
   **/

  static
  std::vector<DGtal::Z3i::RealPoint>
  optimizeCenterLineElasticForces(const std::vector<DGtal::Z3i::Point> &aCenterLine, 
                                  NormalAccumulator &aNormalAcc,                 
                                  double aRadius, double epsilon=0.1){
    typedef DGtal::Z3i::RealPoint TPointOut;
    typedef DGtal::Z3i::Point TPointIn;
    std::vector<TPointOut> optiCenterLine;
  
  
    std::ofstream  outStat;
    outStat.open("optiStat.dat", std::ofstream::out);  
    std::vector<std::vector<NormalAccumulator::Point> > vectAssociations;
  
    for (unsigned int i = 0; i < aCenterLine.size(); i++){
      TPointOut ptCL (aCenterLine.at(i)[0], aCenterLine.at(i)[1], aCenterLine.at(i)[2]);
      vectAssociations.push_back(aNormalAcc.getAssociatedPoints(aCenterLine[i]));
      DGtal::trace.info() << "nb associated:" << aNormalAcc.getAssociatedPoints(aCenterLine[i]).size() << std::endl;
      optiCenterLine.push_back(ptCL);
    }
    
    double aStep = 0.0;
    double deltaE;
    unsigned int  num = 0;
    double previousTot = 0;
    bool  first = true;
    DGtal::trace.info() << "Starting optimisation with min precision diff:  " << epsilon <<  std::endl;
    while (first || deltaE > epsilon){
      num++;
      DGtal::trace.info() << "--------------------------------------" << std::endl;
      DGtal::trace.info() << "Iteration num " << num << "delta E: "<< deltaE <<  std::endl;
      double totalError = 0.0;
      for (unsigned int i = 0; i < aCenterLine.size(); i++){
        TPointOut ptCL (optiCenterLine.at(i)[0], optiCenterLine.at(i)[1], optiCenterLine.at(i)[2]);
        TPointOut sumForces; sumForces[0]=0; sumForces[1]=0; sumForces[2]=0;
        unsigned int nb = 0;
      
        for (unsigned int j = 0; j < vectAssociations[i].size(); j++){
          TPointOut ptMean = vectAssociations[i][j];
          TPointOut vecSM = ptMean - ptCL;
          double normPMoriente = vecSM.norm() - aRadius;  
          TPointOut vectPM = (vecSM/vecSM.norm())*normPMoriente;
          double error = normPMoriente*normPMoriente;
          totalError += error;
          sumForces += vectPM;
          nb++;
        }
      
        optiCenterLine[i] += sumForces/nb;
      }
      if(!first){
        outStat << std::fixed << std::setw( 20 ) << std::setprecision( 10 ) << num << " "<< totalError << std::endl;
      }

        
      DGtal::trace.info() << " Total error:" << totalError << std::endl;
      if (first) {
        deltaE = totalError;
        first = false;
      }else{
        deltaE = previousTot - totalError;
      }
      DGtal::trace.info() << " Total error:" << totalError/optiCenterLine.size() << std::endl;
      DGtal::trace.info() << " delta E:" << deltaE << std::endl;
      previousTot = totalError;    
    }
    return optiCenterLine;
  }





  
  
};


