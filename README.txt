This folder, ProjectCode, contains all datasets and code for 
reweighting the D0BAR data to remove the pseudo-Agamma value.

1) Use the weight_calculation_magup_2016.cxx script for calculating the weights.This script uses the BDT only data sets which were created using the addmva.py
script from manymoo git hub. This applied the TMVA training variables to the 
whole RooDataSet and recorded only the BDT values for each point.
2) This script outputs a root file, asymm_hist_magup_2016.root, which contains the asymmetry between the number of each data points with a given
BDT value as a histogram. The ratio in these bins of BDT values is then the weight which is applied to the D0bar RooDataSet using...
3) applying_weight.cxx script. This uses a method suggested by Dr Alexander to use a histogram method for binning as a histogram. This selects the bin within the histogram that a given BDT value places it in and extracts the ratio from within this bin as the weight for the data point. 
4) This then creates a new RooDataSet that is weighted - D0Bar_weightapplied_magup16.root. Using the magup16_reweighted_fit.cxx script the orginal fitting can be reapplied using a much faster and cleaner templated and functioned code. This first fits the deltam data extracted from the unweighted D0 RooDataSet. The signal value is then extracted from this time integrated value. This is fitted to the time dependent sections for the D0 and D0BAR but now the D0BAR data is reweighted to remove the asymmetry.

note: the 2015 data is located at /nfs/lhcb/d2hh01/hhpi0/data/2015/magup/unweighted_magup15_D0.root or thereabouts 
