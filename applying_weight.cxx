/* Calculating weights using method from email with michael
 * Note that all of the datasets are D0Bar as this is where
 * the weight is applied to. 
 */

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include "RooDataSet.h" 
#include "RooRealVar.h"
#include "RooAbsArg.h"
#include "math.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooBifurGauss.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "RooPlot.h"

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
  os << "{ ";
  for (const auto& val: v) {
    os << val << (&val == &v.back() ? " }" : ", ");
  }
  return os;
}


template <typename T> // returns  contained within *.root files
std::unique_ptr<T>
get_rootfile(const char *filename, const char *obj_name)
{
  T* retVal = nullptr;
  try {
    TFile *file0 = TFile::Open(filename);
    if (!file0) throw std::logic_error("file_not_found");
    file0->GetObject(obj_name, retVal);
    return std::unique_ptr<T>(retVal);
  } catch (const std::exception& e) {
    delete retVal;
    throw e;
  }
}


void applying_weight(){

// Get ratio histogram
auto hratio = get_rootfile<TH1D>("/home/ppe/l/ldickson/MSci/detector_assymmetry_code/code_2_0_1_6/ProjectCode/asymm_hist_magup_2016.root","h3");

// Get original dataset
auto unweighteddata = get_rootfile<RooDataSet>("/home/ppe/l/ldickson/MSci/detector_assymmetry_code/code_2_0_1_6/ProjectCode/unweighted_magup16_D0Bar.root","data");

// Creating new dataset with added weights
  const RooArgSet* vars= unweighteddata->get(0);
  auto weight = RooRealVar("weight","weight",1);
  auto cp_vars = RooArgSet(*vars,"cp_vars"); 
  cp_vars.add(weight);
  RooDataSet weighteddata("weighteddata","weighteddata",cp_vars, "weight");

// Calculating weights for each data point
  for (unsigned i=0; i < unweighteddata->numEntries(); i++){
	const RooArgSet* args = unweighteddata->get(i);
	double bdtval = ((RooRealVar*)&(*args)["BDT"])->getVal();
	double weightval = hratio->GetBinContent(hratio->GetXaxis()->FindBin(bdtval)); 
	weighteddata.add(*args, weightval);
}

// Saving weighted dataset
TFile* fout = new TFile("D0Bar_weightapplied_magup16.root","RECREATE");
fout->WriteObject(&weighteddata,"data");
fout->Close();
delete fout;


}
