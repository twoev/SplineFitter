
#include "Spline/BSpline.hh"

#include "TFile.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TClass.h"

#include "tclap/CmdLine.h"
#include "tclap/ArgException.h"
#include "tclap/CmdLineOutput.h"

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

int main(int argc, char **argv){
 
  TCLAP::CmdLine cmd("B spline fitting example.  Bug reports to James Monk <jmonk@cern.ch>", ' ', "0.1.0");
  TCLAP::ValueArg<std::string> inputArg("i", "inputfile", "name of ROOT file containing histogram", true, "", "string");
  TCLAP::ValueArg<std::string> histoArg("p", "plot", "name of plot to fit", true, "", "string");
  TCLAP::ValueArg<std::string> signalArg("s", "signal", "name of signal template plot to be included in fit", false, "", "string");
  TCLAP::ValueArg<std::string> covarianceFile("c", "covariance-file", "output name of file containing covariance matrix", false, "covariance.dat", "string");
  TCLAP::ValueArg<std::string> outputFile("o", "output-file", "Name of output file", false, "BSplineFit.root", "string");
  
  TCLAP::SwitchArg controlPtSwitch("d", "control-points", "Determine the de Boor control points of the curve", true);
  TCLAP::SwitchArg weightSwitch("w", "uniform-weights", "Use equal weights for each data point (instead of 1/er)", false);
  
  TCLAP::ValueArg<size_t> degree("n", "n-degrees", "Number of degrees per spline", false, 4, "int");
  TCLAP::ValueArg<int> nCoeffsArg("u", "uniform-coeffs", "Use a given number of uniform coefficients", true, 12, "int");
  TCLAP::ValueArg<std::string> knotsArg("k", "knots", "Use the given knots vector (comma separated or a file name)", true, "", "string");
  TCLAP::ValueArg<int> evalArg("e", "eval-pts", "The number of points to evaluate in producing the output fit (default 10,000)", false, 10000, "int");
  
  cmd.add(inputArg);
  cmd.add(histoArg);
  cmd.add(signalArg);
  cmd.add(covarianceFile);
  cmd.add(outputFile);
  cmd.add(controlPtSwitch);
  cmd.add(weightSwitch);
  cmd.add(degree);
  cmd.add(evalArg);

  cmd.xorAdd(nCoeffsArg, knotsArg);
  
  cmd.parse( argc, argv );
  
  Spline::DataSet dataset;
  Spline::DataSet signalSet;
  
  TFile *file = new TFile(inputArg.getValue().c_str(), "READ");
  TNamed *input = (TNamed*)file->Get(histoArg.getValue().c_str())->Clone();
  std::vector<std::string> path;
  boost::split(path, histoArg.getValue(), boost::is_any_of("/"));
  input->SetName(path.back().c_str());
  
  TH1 *signal = 0;
  if(signalArg.getValue() != ""){
    signal = (TH1*) file->Get(signalArg.getValue().c_str())->Clone();
  }
  
  // Any bin with contents below this threshold is not included - discards empty bins
  double threshold;
  
  if(input->IsA()->InheritsFrom("TH1")){
    
    std::cout<<"found TH1"<<std::endl;
    
    // The height of the highest bin
    TH1 *histo = (TH1*)input;
    threshold = histo->GetBinContent(histo->GetMaximumBin()) * 1.e-11;
    if(threshold < 0.){
      threshold = fabs(histo->GetBinContent(histo->GetMinimumBin())) * 1.e-11;
    }

    for(int bin = 1; bin != histo->GetNbinsX()+1; ++bin){
      double height = histo->GetBinContent(bin);
      if(fabs(height) < threshold) continue;
            
      double binEr = (weightSwitch.getValue())? 1.: histo->GetBinError(bin);
      dataset.push_back(Spline::DataPoint(histo->GetBinCenter(bin), height, binEr));
      if(signal){
        binEr = (weightSwitch.getValue())? 1.: signal->GetBinError(bin);
        signalSet.push_back(Spline::DataPoint(signal->GetBinCenter(bin), signal->GetBinContent(bin), binEr));
      }
    }
    
  }else if(input->IsA()->InheritsFrom("TGraph")){
    
    if(signal){
      throw std::runtime_error("signal fit not available when input is a TGraph");
    }
    
    TGraph *graph = (TGraph*)input;
    double *xs = graph->GetX();
    double *ys = graph->GetY();
    double *er = graph->GetEY();      
      
    for(int bin=0; bin != graph->GetN(); ++bin){
      double error = (er == 0 || weightSwitch.getValue())? 1.: er[bin];
      dataset.push_back(Spline::DataPoint(xs[bin], ys[bin], error));
    }
    
  }else{
    throw std::runtime_error("plot to fit is neither TH1 nor TGraph");
  }
  
  Spline::BSpline spline(degree.getValue());
  
  if(dataset.size() < 2) throw std::runtime_error("Error: Fewer than two data points to fit!");
  
  if(nCoeffsArg.isSet()){
  
    double rangeMin = 1.5 * dataset[0].x() - 0.5 * dataset[1].x();
    double rangeMax = 1.5 * dataset[dataset.size()-1].x() - 0.5 * dataset[dataset.size()-2].x();
    spline.setNCoefficientsUniform(nCoeffsArg.getValue(), rangeMin, rangeMax);
  }else{
    
    std::vector<double> knotVector;
    
    std::string delim;
    std::istream *buffer=0;
    if(knotsArg.getValue().find_first_of(",") == std::string::npos){
      delim='\n';
      buffer = new std::ifstream(knotsArg.getValue().c_str());
    }else{
      delim=',';
      buffer = new std::stringstream(knotsArg.getValue());
    }
    
    std::string knot;
    
    while(std::getline(*buffer, knot, *(delim.c_str()))){
      knotVector.push_back(boost::lexical_cast<double>(knot));
    }
    
    if(buffer) delete buffer;
    spline.setKnotVector(knotVector);
  }
  
  if(signalArg.getValue() == ""){
    spline.fitLeastSquares(dataset);
  }else{
    spline.fitLeastSquares(dataset, signalSet);
  }
  
  std::ofstream covOutput(covarianceFile.getValue().c_str());
    
  covOutput<<spline.covarianceMatrix();
  
  Spline::writeVector(std::cout<<std::endl<<"Knot vector: "<<std::endl, spline.knots());
  
  std::cout<<std::endl<<"chi squared of fit: "<< spline.chi2()<<std::endl;
  std::cout<<"chi squared per degree of freedom: "<< spline.chi2_dof()<<std::endl;
  
  std::cout<<std::endl<<"Coefficients of b spline fit:"<<std::endl;
  std::vector<double> coeffs = spline.coefficients();
  Spline::writeVector(std::cout,coeffs);
  
  std::cout<<std::endl<<"Covariance matrix of spline coefficients:"<<std::endl<<std::endl;
  
  std::cout<<spline.covarianceMatrix();  
    
  std::vector<double> xvec;
  std::vector<double> yvec;
  std::vector<double> erUp;
  std::vector<double> erDn;
    
  std::vector<std::vector<double> > basisFuncs(spline.nCoefficients(), std::vector<double>());
  
  double incr = (dataset.back().x() - dataset.front().x()) / (double)evalArg.getValue(); 
  
  for(double x = dataset.front().x(); x < dataset.back().x(); x+=incr){
    Spline::DataPoint pt = spline.evaluate(x);
    xvec.push_back(pt.x());
    yvec.push_back(pt.y());
    erUp.push_back(pt.y() + pt.error());
    erDn.push_back(pt.y() - pt.error());
    
    std::vector<double> basisTmp = spline.basisFunctions(x);
    for(size_t j = 0; j != basisTmp.size(); ++j){
      basisFuncs[j].push_back(basisTmp[j]);
    }
    
  }
  
  std::vector<double> chi2vec;
  std::vector<double> chi2x;
  chi2x.push_back(1.0);
  chi2vec.push_back(spline.chi2_dof());
  
  TGraph *chi2Graph = new TGraph(1, &(chi2x[0]), &(chi2vec[0]));
  
  TGraph *central = new TGraph(xvec.size(), &(xvec[0]),&(yvec[0]));
  TGraph *up = new TGraph(xvec.size(), &(xvec[0]),&(erUp[0]));
  TGraph *dn = new TGraph(xvec.size(), &(xvec[0]),&(erDn[0]));
  
  chi2Graph->SetName("chi2_dof");
  central->SetName("central");
  up->SetName("ErUp");
  dn->SetName("ErDn");
  
  TFile *outputRoot = new TFile(outputFile.getValue().c_str(), "RECREATE");
  
  chi2Graph->Write();
  central->Write();
  up->Write();
  dn->Write();
  input->Write();
  if(signal != 0){
    
    TH1 *histo = (TH1*)input;
    
    signal->Write();
    
    double signalContribution = spline.signalCoefficientSize();
    
    std::cout<<std::endl<<"Signal template coefficient = "<<signalContribution<<std::endl;
    
    std::vector<double> xFit;
    std::vector<double> yFit;
    std::vector<double> xerFit;
    std::vector<double> yerFit;
    
    for(int bin = 1; bin != histo->GetNbinsX(); ++bin){
      double x = histo->GetBinCenter(bin);
      if( x < dataset.front().x()) continue;
      if( x > dataset.back().x()) break;
      Spline::DataPoint pt = spline.evaluate(histo->GetBinCenter(bin));
      xFit.push_back(pt.x());
      yFit.push_back(pt.y() + signalContribution * signal->GetBinContent(bin));
      xerFit.push_back(0.);
      // Don't include the signal template error in the fit error because
      // 1) it is purely statistical from the MC 
      // 2) To include it properly you would have to vary the fit template so it propagates through the covariance matrix.
      //    The second is still possible, you just have to re-run with a different template :)
      yerFit.push_back(pt.error());
    }
    
    TGraphErrors *fullFit = new TGraphErrors(xFit.size(), &(xFit[0]), &(yFit[0]), &(xerFit[0]), &(yerFit[0]));
    fullFit->SetName("fullFit");
    fullFit->Write();
  }
  
  if(controlPtSwitch.getValue()){
  
    Spline::PointSet controlPoints = spline.controlPoints();
    
    std::vector<double> cx;
    std::vector<double> cy;
    
    for(Spline::PointSet::const_iterator pt = controlPoints.begin(); pt != controlPoints.end(); ++pt){
      cx.push_back(pt->x());
      cy.push_back(pt->y());    
    }
    
    TGraph *control = new TGraph(cx.size(), &(cx[0]), &(cy[0]));
    control->SetName("controlPts");
    control->Write();
  }
  
  std::vector<double> coefficientIndex;
  
  for(size_t ii=0; ii != spline.nCoefficients(); ++ii){
    coefficientIndex.push_back((double)ii);
  }
  
  TGraph *coeffGraph = new TGraph(spline.nCoefficients(), &(coefficientIndex[0]), &(spline.coefficients()[0]));
  coeffGraph->SetName("coefficients");
  coeffGraph->Write();
    
  std::vector<double> covX;
  std::vector<double> covY;
  std::vector<double> covVal;
  
  for(size_t ii=0; ii != spline.covarianceMatrix().nRows(); ++ii){
    for(size_t jj=0; jj != spline.covarianceMatrix().nColumns(); ++jj){
      covX.push_back(ii);
      covY.push_back(jj);
      covVal.push_back(spline.covarianceMatrix().element(ii, jj));
    }
  }
  
  TGraph2D *covarianceGraph = new TGraph2D(covVal.size(), &(covX[0]), &(covY[0]), &(covVal[0]));
  covarianceGraph->SetName("CovarianceMatrix");
  covarianceGraph->Write();
  
  std::vector<double> knotIndex;
  
  for(size_t ii=0; ii != spline.knots().size(); ++ii){
    knotIndex.push_back((double)ii);
  }
  
  TGraph *knotGraph = new TGraph(knotIndex.size(), &(knotIndex[0]), &(spline.knots()[0]));
  knotGraph->SetName("knots");
  knotGraph->Write();
  
  outputRoot->mkdir("bases")->cd();
  
  for(size_t coeff=0; coeff != spline.nCoefficients(); ++coeff){
    TGraph *base = new TGraph(xvec.size(), &(xvec[0]), &(basisFuncs[coeff][0]));
    std::string baseName = "B_" + boost::lexical_cast<std::string>(coeff);
    base->SetName(baseName.c_str());
    base->Write();
  }
  
  outputRoot->Close();
  
}
















