//most of them are useless now, just keep the same structure as what it has before
const Int_t kNbCharge = 4; //[OS,MM,PP,LS]
Char_t *chargeName[kNbCharge] = {"OS", "MM", "PP", "LS"};
const Int_t kNbCath = 6; //[gl-gl,gl-tr,tr-tr,gl-calo,tr-calo,calo-calo]
Char_t *oniaCatName[kNbCath] = {"glgl", "gltr", "trtr", "glcalo", "trcalo", "calocalo"};
const Int_t kNbCuts = 7;
//Char_t *muCutName[kNbCuts] = {"noCuts","dist", "chi2", "comp", "vertexProb"};
Char_t *muCutName[kNbCuts] = {"Cut0","Cut1", "Cut2", "Cut3", "Cut4", "Cut5", "Cut6"};
const Int_t kNbMuCath = 3;
Char_t *muCatName[kNbMuCath] = {"gl", "tr", "calo"};
const Int_t kNbDimuonSet = 3; //[0]..3.0 < M < 3.2; [1]... M > 2 GeV, [2]... M < 4 GeV
Double_t massMIN[kNbDimuonSet] = {7.0, 2.,  7.};
Double_t massMAX[kNbDimuonSet] = {12.0, 15., 12.};

//cuts to be applied:
Int_t MIN_nhits_trk = 11;//12;
Double_t MAX_normchi2_cal = 3.0;
Double_t MAX_normchi2_trk = 5.0;
Double_t MAX_normchi2_glb = 20.0;
Int_t MIN_nhits_pixel = 1;//2;
Double_t MAX_d0_trk = 5.0;
Double_t MAX_dz_trk = 20.0;
Double_t MIN_vtxprob_jpsi = 0.001;
Double_t MIN_caloComp = 0.9;
//MC matching parameter:
Double_t MAX_deltaR = 0.03;

