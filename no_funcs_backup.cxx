/* Using the reweighted D0Bar dataset, this code reruns the fits
 * and produces time dependent asymmetry plots. These results can
 * then be compared to the unweighted plots to verify improvements
 * in reducing the pseudo-Agamma value
 *
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

void magup16_reweighted_fit(){

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
  //RooArgSet* monte_vars = (RooArgSet*)&(*monte_TObject);
//  RooRealVar* monte_cdel_1 = (RooRealVar*)&(*monte_vars)["monte_cdel_1"];
 
// Getting real data Delta Mass and D0 Mass
 /*
  TFile *_file0 = TFile::Open("/home/ppe/l/ldickson/MSci/2016_datasets/D0_M_1825_1905_magup_Kpipi_2016.root ");
  RooDataSet* data_mass = (RooDataSet*)_file0->Get("data");
  const RooArgSet* temp_vars = (RooArgSet*)data_mass->get(0);
  RooRealVar* deltam = (RooRealVar*)&(*temp_vars)["deltam"];
*/
// Loading MC data for seeding fits
 /*
  TFile *_file_monte_var = TFile::Open("/home/ppe/l/ldickson/MSci/improving_fits_code/datasets/monte_vars_RAS.root");
  TObject * monte_TObject = _file_monte_var->Get("monte_var_values");
  RooArgSet* monte_vars = (RooArgSet*)&(*monte_TObject);
  RooRealVar* monte_cdel_1 = (RooRealVar*)&(*monte_vars)["monte_cdel_1"];
 */
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
 /*
  totalpdf.fitTo(*data_mass);
  auto c2 = new TCanvas();
  auto del_frame = RooPlot_make(deltam,data_mass,totalpdf);
*/
//  RooPlot* del_frame = deltam->frame();
//  data_mass->plotOn(del_frame);
//  totalpdf.plotOn(del_frame);
  /*
  cout << "Finished Fitting with MC seed" <<endl;
  totalpdf.Print("t"); //Print out each component of totalpdf
  totalpdf.plotOn(del_frame, Components("sig_pdf"), LineColor(kBlue), LineStyle(kDashed)) ;
  totalpdf.plotOn(del_frame, Components("gauss_del_bi"), LineColor(kCyan), LineStyle(kDashed)) ;
  totalpdf.plotOn(del_frame, Components("gauss_del"), LineColor(kGreen), LineStyle(kDashed)) ;
  totalpdf.plotOn(del_frame, Components("gauss_del2"), LineColor(kOrange), LineStyle(kDashed)) ;
  totalpdf.plotOn(del_frame, Components("com_bkgpdf"), LineColor(kMagenta), LineStyle(kDashed)) ;
  totalpdf.plotOn(del_frame, Components("ranpi_bkgpdf"), LineColor(kRed), LineStyle(kDashed)) ;
  del_frame->Draw();
*/
/*
  const char* file_name = "deltam_fits_magup_2016.root";
  const char* save_choice = "RECREATE"; 
  void save_file(file_name,del_frame);
*/

  ////GETTING COMBINED DATA  + CREATING VARIABLE BINS\\\\

  // Getting combined data for D0
  
  auto data_D0_comb = get_rootfile<RooDataSet*>("/home/ppe/l/ldickson/MSci/2016_datasets/D0_M_1825_1905_magup_Kpipi_2016.root","data");
  //auto data_D0_comb = *data_D0_comb_load.get();  
  auto ctau_D0 = get_RRV(data_D0_comb,"ctau");
  auto deltam_D0 = get_RRV(data_D0_comb,"deltam");
  float ctau_float;

  // Defining dataset and variables D0_M, deltam and ctau for D0
  /*
  TFile *_file1 = TFile::Open("/home/ppe/l/ldickson/MSci/2016_datasets/D0_M_1825_1905_magup_Kpipi_2016.root ");
  RooDataSet* data_D0_comb = (RooDataSet*) _file1->Get("data");
  const RooArgSet* variables_D0 = data_D0_comb->get(0);
  RooRealVar* ctau_D0 = (RooRealVar*)&(*variables_D0)["ctau"];
  RooRealVar* deltam_D0 = (RooRealVar*)&(*variables_D0)["deltam"];
  cout << "D0 done"  <<endl;
   */
  // Getting combined data for D0_bar

  auto data_D0BAR_comb = get_rootfile<RooDataSet*>("/home/ppe/l/ldickson/MSci/detector_assymmetry_code/code_2_0_1_6/ProjectCode/D0Bar_weightapplied_magup16.root", "data");
 // auto data_D0BAR_comb = *weighted_data_D0BAR_comb_load.get();
  auto ctau_D0BAR = get_RRV(data_D0BAR_comb,"ctau");
  auto deltam_D0BAR = get_RRV(data_D0BAR_comb,"deltam");
  float ctaubar_double;

  // Defining dataset and variables D0_M, deltam and ctau for D0_bar - WEIGHTED 
  /*
  TFile *_file2 = TFile::Open("/home/ppe/l/ldickson/MSci/detector_assymmetry_code/code_2_0_1_6/ProjectCode/D0Bar_weightapplied_magup16.root");
  RooDataSet* data_D0BAR_comb = (RooDataSet*) _file2->Get("data");
  const RooArgSet* variables_D0BAR = data_D0BAR_comb->get(0);
  RooRealVar* ctau_D0BAR = (RooRealVar*)&(*variables_D0BAR)["ctau"];
  RooRealVar* deltam_D0BAR = (RooRealVar*)&(*variables_D0BAR)["deltam"];
  cout << "loaded D0 and D0bar(weighted) data" << endl; 
  */
  // Combining D0 and D0bar and sorting bins
  unsigned int entries_D0 =  data_D0_comb->numEntries() ;
  unsigned  int entries_D0bar =  data_D0BAR_comb->numEntries();
  vector<float> decaytimes_D0 ;
  vector<float> decaytimes_D0bar ;
  vector<float> decaytimes;
  float decaytime_D0 ;
  float decaytime_D0bar ;

  for(unsigned int i = 0 ; i < entries_D0 ; ++i){
    decaytime_D0 = ((RooRealVar*)&(*data_D0_comb->get(i))["ctau"])->getVal() ;
    decaytimes_D0.push_back(decaytime_D0) ;
  }
  
  for(unsigned int q = 0 ; q < entries_D0bar ; ++q){
    decaytime_D0bar = ((RooRealVar*)&(*data_D0BAR_comb->get(q))["ctau"])->getVal() ;
    decaytimes_D0bar.push_back(decaytime_D0bar) ;
  } 

  decaytimes.reserve( decaytimes_D0.size() + decaytimes_D0bar.size() ); // preallocate memory
  decaytimes.insert( decaytimes.end(), decaytimes_D0.begin(), decaytimes_D0.end() );
  decaytimes.insert( decaytimes.end(), decaytimes_D0bar.begin(), decaytimes_D0bar.end() );
//  sort(decaytimes.begin(), decaytimes.end()) ;
cout << "combined datasets with size: " << decaytimes.size();

  ////CREATING VARIABLE BIN EDGES\\\\

/*  
  unsigned  int total_entries = decaytimes.size(); 
  const unsigned  int number_bins = 10;
  unsigned  int number_per_bin = total_entries/number_bins;
  const vector<float> bin_edges_float = get_bin_edges(number_bins,decaytimes);
  cout << "bin edges float:" << bin_edges_float << endl;
  vector<string> bin_edges;
  string temp;
  for (unsigned i = 0; i <bin_edges.size(); i++){
    temp = to_string(bin_edges_float[i]);
    bin_edges.push_back(temp);
  }
*/

 // cout << "bin edges string:" << bin_edges << endl;

unsigned  int total_entries = decaytimes.size(); 
  const unsigned  int number_bins = 10;
  unsigned  int number_per_bin = total_entries/number_bins;

  string decaytimes_str_D0;
  string decaytimes_str_D0bar;
  string decaytimes_str;
  vector<string> bin_edges;
  vector<float> bin_edges_float;
  bin_edges.push_back(to_string(0.));
  //checks the number of data points within a bound as calculated above and produces ctau lower edge
  for (unsigned int j = 1 ; j < total_entries+1 ; j++){
    if (j % number_per_bin == 0.) { 
      decaytimes_str =to_string(decaytimes[j]);
      bin_edges_float.push_back(decaytimes[j]);
      bin_edges.push_back(decaytimes_str);
    }
    else{ continue; }
  }


  ////COMBINING/SORTING VECTORS AND CREATING VARIABLE BINS\\\\

  //Finding the signal strength between bins for BDTvalues_D0
//  auto histdata_D0 = bin_data(bin_edges, BDTvalues_D0);
//  cout << "histdata_D0=" << histdata_D0 << endl;
  // Finding the signal strength between bins for BDTvalues_D0bar
//  auto histdata_D0bar = bin_data(bin_edges, BDTvalues_D0bar);
//  cout << "histdata_D0bar=" << histdata_D0bar << endl;

  //Creating strings for D0 reduced dataset cuts
  vector<string> bin_string_vect;
  string bin_str1;
  string bin_str2;
  string bin_str3;
  string bin_str4;
  string bin_str5;
  string total_bin_str;
  //D0

  for(unsigned int m = 0; m<bin_edges.size() -1; m++){
    bin_str1 =  "ctau>= ";
    bin_str2 =  bin_edges[m] ; 
    bin_str3 = " && " ;
    bin_str4 = "ctau <" ;
    bin_str5 =bin_edges[m+1] ;
    total_bin_str = bin_str1 + bin_str2 + bin_str3 + bin_str4 + bin_str5;
    bin_string_vect.push_back(total_bin_str) ;
  }
cout << "here" << endl;

cout <<"bin string vect"<< bin_string_vect << endl;
  //D0bar
/*
  vector<string>bin_string_vect_bar;
  for(unsigned int m = 0; m<bin_edges.size() -1; m++){
    bin_str1 =  "ctau>= ";
    bin_str2 =  bin_edges[m] ; 
    bin_str3 = " && " ;
    bin_str4 = "ctau <" ;
    bin_str5 =bin_edges[m+1] ;
    total_bin_str = bin_str1 + bin_str2 + bin_str3 + bin_str4 + bin_str5;
    bin_string_vect_bar.push_back(total_bin_str) ;
  }
*/
  //D0 reduce data sets

  vector<RooDataSet*> dataset_d0_vect;
  const char* bin_char;  
 
  for(unsigned int kl = 0; kl<number_bins ; kl++){
    bin_char = bin_string_vect[kl].c_str();
    RooDataSet* d1 = (RooDataSet*) data_D0_comb->reduce(RooArgSet(*deltam_D0, *ctau_D0),bin_char);
    dataset_d0_vect.push_back(d1);
  }

  //D0bar reduced datasets

  vector<RooDataSet*> dataset_d0bar_vect;
  for(unsigned int ky = 0; ky<number_bins  ; ky++){
    bin_char = bin_string_vect[ky].c_str();
    RooDataSet* d1 = (RooDataSet*) data_D0BAR_comb->reduce(RooArgSet(*deltam_D0BAR, *ctau_D0BAR),bin_char);
    dataset_d0bar_vect.push_back(d1);
  }

  ////SETTING SHAPE PARAMETERS CONSTANT\\\\
  // loading variable's data created from fit
  RooArgSet* totalpdf_vars = totalpdf.getParameters(*data_mass);

  // assinging variable values to constant. Signal/background variables left to vary

  // -- Signal Shape Delta mass--
  set_const(totalpdf_vars, "mean_del");
  set_const(totalpdf_vars,"mean_del");
  set_const(totalpdf_vars,"mean_del2");
  set_const(totalpdf_vars,"mean_bi");
  set_const(totalpdf_vars,"sigma_del");
  set_const(totalpdf_vars,"sigma_del2");
  set_const(totalpdf_vars,"sigma_L_bi");
  set_const(totalpdf_vars,"sigma_R_bi");
  set_const(totalpdf_vars,"cdel_1");
  set_const(totalpdf_vars,"cdel_2");
  


  // --Background Shape -- 
  // ((RooRealVar*)&(*totalpdf_vars)["M_th_del"])->setConstant() ; already set constant earlier
  set_const(totalpdf_vars,"a_del");
  set_const(totalpdf_vars,"b_del");

  // -- Signal Value --
 // total_nsignal

  // -- Background Value --
 //total_com_bkg
 //total_ranpi_bkg


////FITTING NEW TIME-DEPENDANT SIGNAL VALUES AND PUTTING INTO VECTOR\\\\
  // D0 decay time itteration over subdatasets from vector
  int k;
  int data_d0_size=dataset_d0_vect.size(); 
  RooRealVar *nsignal_vals_D0 ;
  vector<double> nsignal_D0_vect_double;
  vector<double> nsignal_D0_error_vect_double;
  RooAbsArg* coefficient;
  double nsignal_D0_double;
  double nsignal_D0_error_double;
  RooArgSet* totalpdf_signalvars;
  string title_str;
  const char* title_str_char;
  string loop_str; 
  string file_str;
  const char* file_str_char;
  
// Creating file for saving time dependent fits for D0  
TFile* td_D0_up_16 = new TFile("test_reweighted_time_dependent_D0_deltam_magup_2016.root","RECREATE");

for (k=0; k<data_d0_size; k++) {
    totalpdf.fitTo(*dataset_d0_vect[k]);    
    loop_str = to_string(k);
    title_str = "time dependent D0 deltam" + loop_str;
    title_str_char = title_str.c_str();
    file_str = "D0_time_dependence_deltam" + loop_str;
    file_str_char = file_str.c_str();
    
    // auto canvas_D0 = new TCanvas();
    RooPlot* del_D0_td_frame = deltam_D0->frame(Title(file_str_char));
    dataset_d0_vect[k]->plotOn(del_D0_td_frame);
    totalpdf.plotOn(del_D0_td_frame);
    del_D0_td_frame->Write(file_str_char);

    totalpdf_signalvars = totalpdf.getParameters(*dataset_d0bar_vect[k]);
    nsignal_D0_double = nsignal.getVal();
    nsignal_D0_error_double = nsignal.getError();
    nsignal_D0_vect_double.push_back(nsignal_D0_double);
    nsignal_D0_error_vect_double.push_back(nsignal_D0_error_double);

    //   canvas_D0->Close();
 }
td_D0_up_16->Close();
delete td_D0_up_16;


 //D0_bar decay time itteration over subdatasets from vector
  int p; 
  int data_D0bar_size = dataset_d0bar_vect.size();
  RooRealVar *nsignal_vals_D0bar ;
  vector<double> nsignal_D0bar_vect_double;
  vector<double> nsignal_D0bar_error_vect_double;
  double nsignal_D0bar_double;
  double nsignal_D0bar_error_double;
  RooArgSet* totalpdf_signalvars_BAR;
  string title_str_BAR;
  const char* title_str_char_BAR;
  string loop_str_BAR;
  string file_str_BAR;
  const char* file_str_char_BAR;

// Creating file for saving time dependent fits for D0BAR
TFile* td_D0BAR_up_16 = new TFile("test_time_dependent_D0BAR_deltam_magup_2016.root","RECREATE");

  for (p=0; p<data_D0bar_size; p++) {
    totalpdf.fitTo(*dataset_d0bar_vect[p], SumW2Error(kTRUE));
    loop_str_BAR = to_string(p);
    title_str_BAR = "time dependent D0BAR deltam" + loop_str; 
    title_str_char_BAR = title_str_BAR.c_str();
    file_str_BAR = "D0BAR_time_dependence_deltam" + loop_str;
    file_str_char_BAR = file_str_BAR.c_str();
    
    // auto canvas_D0BAR = new TCanvas();
    RooPlot* del_D0BAR_td_frame = deltam_D0BAR->frame(Title(file_str_char_BAR));
    dataset_d0bar_vect[p]->plotOn(del_D0BAR_td_frame);
    totalpdf.plotOn(del_D0BAR_td_frame);
    del_D0BAR_td_frame->Write(file_str_char_BAR);

    totalpdf_signalvars_BAR = totalpdf.getParameters(*dataset_d0bar_vect[p]);
    nsignal_D0bar_double = nsignal.getVal();
    nsignal_D0bar_error_double = nsignal.getError();
    nsignal_D0bar_vect_double.push_back(nsignal_D0bar_double);
    nsignal_D0bar_error_vect_double.push_back(nsignal_D0bar_error_double);
    
    //  canvas_D0->Close();
  }
td_D0BAR_up_16->Close();
delete td_D0BAR_up_16;


  //                 PLOTTING DATA ON HIST              //
 
//  D0  and D0bar  //
/*  
   	vector <Float_t> edge_vector;
	Double_t edge;
	Float_t edge_float;
	const char* bin_edge_char;
	for (unsigned int z=0; z<bin_edges.size();z++){
		bin_edge_char = bin_edges[z].c_str();
		edge = atof(bin_edge_char);
		edge_float = (const Float_t) edge;
		edge_vector.push_back(edge);
	}
*/

//  TH1D* D0_sig_hist = new TH1D("D0_sig_hist","D0 signal hist", number_bins, &edge_vector[0]); 
//  TH1D* D0bar_sig_hist = new TH1D("D0bar_sig_hist","D0 BAR signal hist", number_bins, &edge_vector[0]);
  
  TH1D* D0_sig_hist = new TH1D("D0_sig_hist","D0 signal hist", number_bins, &bin_edges_float[0]);
  TH1D* D0bar_sig_hist = new TH1D("D0bar_sig_hist","D0 BAR signal hist", number_bins, &bin_edges_float[0]);

  D0_sig_hist->Sumw2();
  D0bar_sig_hist->Sumw2();

  for (unsigned int q=0; q<number_bins; q++) {
    // D0 //
    nsignal_D0_double = nsignal_D0_vect_double[q];
    D0_sig_hist->SetBinContent(q+1,nsignal_D0_double);

    // D0 BAR //
    nsignal_D0bar_double = nsignal_D0bar_vect_double[q];
    D0bar_sig_hist->SetBinContent(q+1,nsignal_D0bar_double);		     
  }  

  auto c6 =  new TCanvas();
  D0_sig_hist->Draw();

  auto c7 = new TCanvas();
  D0bar_sig_hist->Draw();

  TH1D *h3 = (TH1D*)D0_sig_hist->Clone("h3");
  TH1D *h4 = (TH1D*)D0_sig_hist->Clone("h4");
  h3->Add(D0bar_sig_hist,-1);
  h4->Add(D0bar_sig_hist);
  h3->Divide(h4);


  //     CALCULATING ERRORS     //
  

//err = sqrt( errD0^2*(dD0(f))^2 + errD0bar^2*(dD0bar(f))^2 )

  int u;  
  double dD0_denom;
  double dD0_numer;
  double dD0_differential;
  double dD0_term1;
  double dD0bar_denom;
  double dD0bar_numer;
  double dD0bar_differential;
  double dD0bar_term2;
  double error_square_sum;
  vector<double> error_vect;
  for (u=0; u<data_d0_size; u++) {
    // D0 partial differntial squared  
    dD0_numer = 4*pow(nsignal_D0bar_vect_double[u],2);
    dD0_denom = pow(nsignal_D0_vect_double[u]+ nsignal_D0bar_vect_double[u],4);
    dD0_differential = dD0_numer/dD0_denom;
    dD0_term1 = dD0_differential*pow(nsignal_D0_error_vect_double[u],2);

    //D0bar partial differential squared
    dD0bar_numer = 4*pow(nsignal_D0_vect_double[u],2);
    dD0bar_denom = pow(nsignal_D0_vect_double[u]+ nsignal_D0bar_vect_double[u],4);
    dD0bar_differential = dD0bar_numer/dD0bar_denom;
    dD0bar_term2 =  dD0bar_differential*pow(nsignal_D0bar_error_vect_double[u],2);

    // Total error
    error_square_sum = sqrt(dD0_term1+dD0bar_term2);
    error_vect.push_back(error_square_sum);
    h3->SetBinError(u+1 , error_square_sum);
  }

  auto c8 = new TCanvas();
  h3->Draw();

 TFile* hist_file = new TFile("test_testing_wght_app_hist_magup_2016.root","RECREATE");
  h3->Write();
  hist_file->Close();
  delete hist_file;


}
