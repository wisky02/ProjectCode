
/*Takes the BDT dataset and calculates the weight for 
each BDT value from the ratio of the number of D0 and 
D0BAR candidates for a given BDT value. This weight 
is then outputed via a TTree which can be used to 
create a dataset with reweighted observables thus 
cancelling the detector asymmetry
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

template <typename T>
std::vector<T> get_bin_edges(unsigned nbins, std::vector<T>& data)
{
  std::vector<T> retVal;
  retVal.reserve(nbins + 1);
  const auto itBeg = data.begin();
  const auto itEnd = data.end();
  for (unsigned i = 0; i < nbins; ++i) {
    const auto it = itBeg + std::size_t(float(i) / float(nbins) * (itEnd
          - itBeg));
    std::nth_element(itBeg, it, itEnd);
    retVal.push_back(*it);
  }
  std::nth_element(itBeg, itEnd - 1, itEnd);
  retVal.push_back(std::nextafter(*(itEnd-1), std::numeric_limits<T>::infinity()));
  assert(retVal.size() == nbins + 1);
  return retVal;
}

template <typename T>
std::vector<T> bin_data(const std::vector<T>& bin_edges, const std::vector<T>& data)
{
  std::vector<T> retVal(bin_edges.size() - 1, 0);
  assert(bin_edges.size() == retVal.size() + 1);
  for (auto val: data) {
    const auto it = std::lower_bound(bin_edges.begin(),
				     bin_edges.end(), val);
    assert(bin_edges.begin() <= it && it < bin_edges.end());
    const auto idx = (val < *it) ? ((it - bin_edges.begin()) - 1) :
      (it - bin_edges.begin());
    assert(idx < retVal.size());
    ++retVal[idx];
  }
  return retVal;
}

template <typename T>
TH1D* make_TH1D(const char* name, const char* title,
		const std::vector<T>& bin_bounds,
		const std::vector<T>& bin_data)
{
  assert(bin_bounds.size() == bin_data.size() + 1);
  std::vector<double> bounds;
  bounds.reserve(bin_bounds.size());
  std::copy(bin_bounds.begin(), bin_bounds.end(), std::back_inserter(bounds));
  TH1D* retVal = new TH1D(name, title, bin_data.size(), bounds.data());
  retVal->Sumw2();
  double total = 0;
  for (unsigned i = 0; bin_data.size() != i; ++i) {
    total += bin_data[i];
    retVal->SetBinContent(i + 1, bin_data[i]);
  }
  retVal->SetEntries(total);
  return retVal;
}

void weight_calculation_magup_2016(){

  using namespace RooFit;
  gStyle->SetPalette(kBird);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  using namespace std; 

  ////Data Input\\\\

  // Getting combined data for D0
  TFile *_file1 = TFile::Open("/home/ppe/l/ldickson/MSci/detector_assymmetry_code/code_2_0_1_6/ProjectCode/BDT_ONLY_datasets/magup_2016/D0_Kpipi_magup_2016_BDT_ONLY.root");
  // Defining dataset and variables D0_M, deltam and ctau for D0
  RooDataSet* data_D0_comb = (RooDataSet*) _file1->Get("data");
  const RooArgSet* variables_D0 = data_D0_comb->get(0);
  RooRealVar* BDT_D0 = (RooRealVar*)&(*variables_D0)["BDT"];
  float ctau_float;

  // Getting combined data for D0_bar
  TFile *_file2 = TFile::Open("/home/ppe/l/ldickson/MSci/detector_assymmetry_code/code_2_0_1_6/ProjectCode/BDT_ONLY_datasets/magup_2016/D0BAR_Kpipi_magup_2016_BDT_ONLY.root");

  // Defining dataset and variables D0_M, deltam and ctau for D0_bar
  RooDataSet* data_D0BAR_comb = (RooDataSet*) _file2->Get("data");
  const RooArgSet* variables_D0BAR = data_D0BAR_comb->get(0);
  RooRealVar* BDT_D0BAR = (RooRealVar*)&(*variables_D0)["BDT"];
  float ctaubar_double;

  // Combining D0 and D0bar and sorting bins
  unsigned entries_D0 =  data_D0_comb->numEntries() ;
  unsigned entries_D0bar =  data_D0BAR_comb->numEntries();
  vector<float> BDTvalues_D0 ;
  vector<float> BDTvalues_D0bar ;
  vector<float> BDTvalues;
  float BDTvalue_D0 ;
  float BDTvalue_D0bar ;

  for(unsigned i = 0 ; i < entries_D0 ; ++i){
    BDTvalue_D0 = ((RooRealVar*)&(*data_D0_comb->get(i))["BDT"])->getVal() ;
    BDTvalues_D0.push_back(BDTvalue_D0) ;
  }

  for(unsigned q = 0 ; q < entries_D0bar ; ++q){
    BDTvalue_D0bar = ((RooRealVar*)&(*data_D0BAR_comb->get(q))["BDT"])->getVal() ;
    BDTvalues_D0bar.push_back(BDTvalue_D0bar) ;
  } 

  BDTvalues.reserve( BDTvalues_D0.size() + BDTvalues_D0bar.size() ); // preallocate memory
  BDTvalues.insert( BDTvalues.end(), BDTvalues_D0.begin(), BDTvalues_D0.end() );
  BDTvalues.insert( BDTvalues.end(), BDTvalues_D0bar.begin(), BDTvalues_D0bar.end() );

  unsigned total_entries = BDTvalues.size(); 
  const unsigned number_bins = 100 ;
  unsigned  number_per_bin = total_entries/number_bins ;

  vector<float> bin_edges = get_bin_edges(number_bins, BDTvalues);
  cout << "bin_edges=" << bin_edges << endl;

  ////COMBINING/SORTING VECTORS AND CREATING VARIABLE BINS\\\\

  //  Finding the signal strength between bins for BDTvalues_D0
  auto histdata_D0 = bin_data(bin_edges, BDTvalues_D0);
  cout << "histdata_D0=" << histdata_D0 << endl;
  // Finding the signal strength between bins for BDTvalues_D0bar
  auto histdata_D0bar = bin_data(bin_edges, BDTvalues_D0bar);
  cout << "histdata_D0bar=" << histdata_D0bar << endl;

  ////FITTING NEW TIME-DEPENDANT SIGNAL VALUES AND PUTTING INTO VECTOR////


  // Saving reweighting values for further analysis 
//  TFile* calculated_weights_file_2016_magup = new TFile("calculated_weights_file_2016_magup.root","RECREATE");
//  calculated_weights_file_2016_magup->Close();
//  delete calculated_weights_file_2016_magup;


  ////PLOTTING DATA ON HIST\\\\

  TH1D* D0_sig_hist = make_TH1D("D0_sig_hist","D0 signal hist", bin_edges, histdata_D0); 
  TH1D* D0bar_sig_hist = make_TH1D("D0bar_sig_hist","D0 BAR signal hist", bin_edges, histdata_D0bar);

//  auto c6 =  new TCanvas();
//  D0_sig_hist->Draw();

//  auto c7 = new TCanvas();
//  D0bar_sig_hist->Draw();

  TH1D *h3 = (TH1D*)D0bar_sig_hist->Clone("h3");
  TH1D *h4 = (TH1D*)D0_sig_hist->Clone("h4");
  h3->Divide(h4);

  auto c8 = new TCanvas();
  h3->Draw();

  TFile* hist_file = new TFile("asymm_hist_magup_2016.root","RECREATE");
  h3->Write();
  hist_file->Close();
  delete hist_file;

  ////EXTRACTING WEIGHTS\\\\

  vector<double> weights_vector1;
  Double_t bin_val_extract;
  for (unsigned i = 1; i <=number_bins; i++){
    bin_val_extract = h3->GetBinContent(i);
    weights_vector1.push_back(bin_val_extract);
  }
 
 ////FINAL OUTPUT\\\\

/* 
 TFile* fout = new TFile("weights_and_bin_edges_magup16.root","RECREATE");
 fout->WriteObject(&weights_vector1,"weights_vector");
 fout->WriteObject(&bin_edges,"bin_edges");
 fout->Close();
 delete fout;
*/

}
