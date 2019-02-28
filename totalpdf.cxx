/* Runs the totalpdf fit seeded with the MC data and saves it so that the other fits (experimental not final results 
 * which should be run in full) can be quicker
 */

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
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

// saving plots as pdf for checking on laptop
template <typename T>
void save_this_plot(const char* filename, T plot_to_save ){
TCanvas canv;
plot_to_save.Draw();
canv.SaveAs(filename);
}

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
T get_rootfile(const char *filename, const char *obj_name)
{
  try {
    TFile *file0 = TFile::Open(filename);
    if (!file0) throw std::logic_error("file_not_found");
    T retVal =(T)  file0->Get(obj_name);
    return(retVal);
  } catch (const std::exception& e) {
    throw e;
  }
}

const RooArgSet* get_RAS(RooDataSet* dataset){
  try{
    const RooArgSet* vars = dataset->get(0);
    if(!dataset) throw std::logic_error("dataset_does_not_exist");
    return(vars); } 
  catch (const std::exception& e){
  throw e;
  }
}

// extracts RooRealVar from RooDataSet
RooRealVar* get_RRV(RooDataSet* dataset,const char* variable_name){
  try{
    const RooArgSet* vars = dataset->get(0);
    if(!dataset){
      throw std::logic_error("dataset_does_not_exist");
    }
    RooRealVar* rrv_value = (RooRealVar*)&(*vars)[variable_name];
    return(rrv_value); } 
  catch (const std::exception& e){
  throw e;
  }
}

// returns value for RRV
float get_val(RooArgSet* vars, const char* val_name){
float val = ((RooRealVar*)&(*vars)[val_name])->getVal() ;
return(val);
}

template<typename T>
void save_file(const char* save_loc, T item_2_save){
  TFile* sfile = new TFile(save_loc,"RECREATE");
  item_2_save->Write();
  delete sfile;
}

RooPlot* RooPlot_make(RooRealVar* RRV, RooDataSet* data,RooAddPdf pdf){
  RooPlot* Frame = RRV->frame();
  data->plotOn(Frame);
  pdf.plotOn(Frame);
  return(Frame);
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

void set_const(RooArgSet* pdf_args,const char* const_var){
  ((RooRealVar*)&(*pdf_args)[const_var])->setConstant();
}

void totalpdf(){

  using namespace RooFit;
  gStyle->SetPalette(kBird);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  using namespace std; 
  
////Data Input\\\\ - Currently not working

  auto data_mass = get_rootfile<RooDataSet*>("/home/ppe/l/ldickson/MSci/2016_datasets/D0_M_1825_1905_magup_Kpipi_2016.root","data");
  RooRealVar* deltam = get_RRV(data_mass,"deltam");  
  cout << "RRV suc" <<endl;
  cout << "Finished loading data" << endl;
  
  // Loading MC data for seeding fits
  auto monte_vars = get_rootfile<RooArgSet*>("/home/ppe/l/ldickson/MSci/improving_fits_code/datasets/monte_vars_RAS.root","monte_var_values");
 
  // Defining MC variables for seeding fits 
  auto monte_cdel_1_float = get_val(monte_vars,"monte_cdel_1");
  auto monte_cdel_2_float  = get_val(monte_vars,"monte_cdel_2"); 
  auto monte_mean_bi_float  = get_val(monte_vars,"monte_mean_bi");
  auto monte_mean_del_float  = get_val(monte_vars,"monte_mean_del");
  auto monte_mean_del2_float = get_val(monte_vars,"monte_mean_del2");
  auto monte_sigma_L_bi_float = get_val(monte_vars,"monte_sigma_L_bi");
  auto monte_sigma_R_bi_float  = get_val(monte_vars,"monte_sigma_R_bi");
  auto monte_sigma_del_float  = get_val(monte_vars,"monte_sigma_del");
  auto monte_sigma_del2_float  = get_val(monte_vars,"monte_sigma_del2");
  cout << "Finished loading MC data" <<endl;
  
  ////SIGNAL: SHAPE PARAMETER\\\\

  // Creating guassians and bifurcated guassian signal pdfs for deltam 
  RooRealVar mean_del("mean_del","mean_del",monte_mean_del_float,140,150); 
  RooRealVar sigma_del("sigma_del","sigma_del",monte_sigma_del_float,0.,10.); 
  RooGaussian gauss_del("gauss_del","gauss_del",*deltam,mean_del,sigma_del); 

  RooRealVar mean_del2("mean_del2","mean_del2",monte_mean_del2_float,140,150); 
  RooRealVar sigma_del2("sigma_del2","sigma_del2",monte_sigma_del2_float,0.,10.); 
  RooGaussian gauss_del2("gauss_del2","gauss_del2",*deltam,mean_del2,sigma_del2);

  RooRealVar mean_bi("mean_bi","mean_bi",monte_mean_bi_float,140,150);
  RooRealVar sigma_L_bi("sigma_L_bi","sigma_L_bi",monte_sigma_L_bi_float,0.01,5);
  RooRealVar sigma_R_bi("sigma_R_bi","sigma_R_bi",monte_sigma_R_bi_float,0.01,5);
  RooBifurGauss gauss_del_bi("gauss_del_bi","gauss_del_bi",*deltam,mean_bi,sigma_L_bi, sigma_R_bi);

  // Total delta signal (two gaussians + bifurcated gaussian)
  RooRealVar cdel_1("cdel_1", "", 0.5, 0, 1) ;
  RooRealVar cdel_2("cdel_2", "", 0.5, 0, 1) ;
  RooAddPdf del_sig("del_sig","", RooArgList(gauss_del,gauss_del2,gauss_del_bi), RooArgList( cdel_1, cdel_2), true);


  ////BACKGROUND: SHAPE PARAMETER\\\\

  // Deltam combinatorial background (modified exponential)
  RooRealVar M_th_del("M_th_del","M_th_del",139,170) ;
  M_th_del.setVal(139.57);   // pi plus/minus mass as threshold value
  M_th_del.setConstant(); 
  RooRealVar a_del("a_del","a_del",20,0,60) ;
  RooRealVar b_del("b_del","b_del", -20,-100,100) ;
  //RooArgList(*deltam,M_th_del,a_del,b_del);
  RooGenericPdf bkg_del("bkg_del","Combinatorial Background plot","(1-exp(-(deltam-M_th_del)/a_del))*(deltam/M_th_del)^b_del",RooArgList(*deltam,M_th_del,a_del,b_del)) ;


  ////COMBINING: SHAPE PARAMETER\\\\

  // Total signal (deltam: 2 gauss + bifurcated , D0: 2 gauss)
  RooProdPdf sig_pdf("sig_pdf","sig_pdf",RooArgSet(del_sig));

  // Total Combinatorial background (deltam: modified exponential , D0: straight line)
  RooProdPdf com_bkgpdf("com_bkgpdf","com_bkgpdf",RooArgSet(bkg_del));

  // Total Random pi background (deltam: modified exponential, D0: same gaussian as signal
  RooProdPdf ranpi_bkgpdf("ranpi_bkgpdf","ranpi_bkgpdf",RooArgSet(bkg_del));

  //Here the bkg_del pdf is used for both the random pi and combinatorial background
  cout << "Finished creating shape parameters" <<endl;

  ////FINAL  PDF: SHAPE PARAMETER\\\\
  
  auto num_data_mass = data_mass->numEntries();
  RooRealVar nsignal("nsignal", "", num_data_mass*0.3, 0, num_data_mass*1.5);
  RooRealVar n_com_bkg("n_com_bkg","", num_data_mass*0.35,0,num_data_mass*1.5);
  RooRealVar n_ranpi_bkg("n_ranpi_bkg","",num_data_mass*0.35, 0,num_data_mass*1.5);

  // Sum of signal, random pi and combinatorial background components
  RooAddPdf totalpdf("totalpdf", "", RooArgList(sig_pdf, com_bkgpdf, ranpi_bkgpdf), RooArgList(nsignal, n_com_bkg,n_ranpi_bkg));

 cout << "Finished creating total shape parameters" << endl;

  ////FITTING AND PLOTTING: SHAPE PARAMETER\\\\

  //Fitting and printing plots
  totalpdf.fitTo(*data_mass);
  TFile* totalpdf_fit = new TFile("totalpdf_fit.root","RECREATE");
  totalpdf.Write();
}
