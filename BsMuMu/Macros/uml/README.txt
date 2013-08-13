The code in this folder can be used to run the UML fit on the final mass.
Both 1D-BDT and the binned BDT fits can be run.

The fitter requires root version 5.34.09.

Some pointers to the code:
- fit_pdf contains the pdf for the fit
- *_gau are constraints pdf (syst)
- main workspace is read from root files
- define_total_extended combines the pdfs into the fit pdf

How to run:
1D fit 
./main_fitData.o -ws output/ws_simul4_bdt_BF2_PEE_keys_new.root -simul 4 -BF 2 -pee -y all -syst -rare_constr  -final -i input/realdata.txt -cuts_file a -print  -sig 1 -berns

Fit bdt cats
./main_fitData.o -ws output/ws_simul4_bdt_simulAll_BF2_PEE_keys_new.root -simul 4 -BF 2 -pee -y all -syst -simul_all 12   -final -i input/realdata.txt -rare_constr  -sig 1 -berns

Toy for 1D
./main_toyMC.o -ws output/ws_simul4_bdt_BF2_PEE_keys_new.root -simul 4 -BF 2 -pee -y all -syst -nexp 10 -rare_constr -sig 0 -final -berns

Toy for bdt cats
./main_toyMC.o -ws output/ws_simul4_bdt_simulAll_BF2_PEE_keys_new.root -e input/estimates.txt -simul 4 -BF 2 -pee -y all -syst -simul_all 12 -nexp 10 -rare_constr -sig 0 -final -berns
