#include "CallLibraries.h"  // call libraries from ROOT and C++

//Event information
TH1D *mHistNEvents = new TH1D("mHistNEvents","",2,0.,2.); mHistNEvents->Sumw2();
TH1D *mHistNTrks = new TH1D("mHistNTrks","",100,0.,100.); mHistNTrks->Sumw2();
TH1D *mHistMultiplicity = new TH1D("mHistMultiplicity","",100,0.,100.); mHistMultiplicity->Sumw2();
TH1D *mHistEMultiplicity = new TH1D("mHistEMultiplicity","",100,0.,100.); mHistEMultiplicity->Sumw2();
TH1D *mHistNTrksMC = new TH1D("mHistNTrksMC","",100,0.,100.); mHistNTrksMC->Sumw2();
TH1D *mHistMultiplicityMC = new TH1D("mHistMultiplicityMC","",100,0.,100.); mHistMultiplicityMC->Sumw2();
TH1D *mHistEMultiplicityMC = new TH1D("mHistEMultiplicityMC","",100,0.,100.); mHistEMultiplicityMC->Sumw2();
TH1D *mHistNTrksTRUTHMC = new TH1D("mHistNTrksTRUTHMC","",100,0.,100.); mHistNTrksTRUTHMC->Sumw2();
TH1D *mHistMultiplicityTRUTHMC = new TH1D("mHistMultiplicityTRUTHMC","",100,0.,100.); mHistMultiplicityTRUTHMC->Sumw2();
TH1D *mHistEMultiplicityTRUTHMC = new TH1D("mHistEMultiplicityTRUTHMC","",100,0.,100.); mHistEMultiplicityTRUTHMC->Sumw2();
TH1D *mHistNEvMix = new TH1D("mHistNEvMix","",30,0.,30.); mHistNEvMix->Sumw2();

//B0 tracker hits
TH2D* mHistB0OccLayer0 = new TH2D("mHistB0OccLayer0", "mHistB0OccLayer0", 100, -400, 0, 100, -170, 170); mHistB0OccLayer0->Sumw2();
TH2D* mHistB0OccLayer1 = new TH2D("mHistB0OccLayer1", "mHistB0OccLayer1", 100, -400, 0, 100, -170, 170); mHistB0OccLayer1->Sumw2();
TH2D* mHistB0OccLayer2 = new TH2D("mHistB0OccLayer2", "mHistB0OccLayer2", 100, -400, 0, 100, -170, 170); mHistB0OccLayer2->Sumw2();
TH2D* mHistB0OccLayer3 = new TH2D("mHistB0OccLayer3", "mHistB0OccLayer3", 100, -400, 0, 100, -170, 170); mHistB0OccLayer3->Sumw2();
TH1D* mHistB0HitEDep = new TH1D("mHistB0HitEDep", ";Deposited Energy [keV]", 100, 0.0, 500.0); mHistB0HitEDep->Sumw2();
	
//B0 EMCAL clusters
TH2D* mHistB0Ecal = new TH2D("mHistB0Ecal", ";hit x [mm];hit y [mm]", 100, -400, 0, 100, -170, 170); mHistB0Ecal->Sumw2();
TH1D* mHistB0EcalEDep = new TH1D("mHistB0EcalEDep", ";Cluster Energy [GeV]", 100, 0.0, 100.0); mHistB0EcalEDep->Sumw2();
TH1D* mHistB0EPTrk = new TH1D("mHistB0EPTrk", "E/p", 1000, 0.0, 10.0); mHistB0EPTrk->Sumw2();
TH2D* mHistB0EP2DTrk = new TH2D("mHistB0EP2DTrk", "mHistB0EP2DTrk", 1000, 0.0, 100.0, 1000.0, 0.0, 100.0); mHistB0EP2DTrk->Sumw2();
TH1D* mHistB0DRTrk = new TH1D("mHistB0DRTrk", "mHistB0DRTrk", 500, 0.0, 5.0); mHistB0DRTrk->Sumw2();
TH1D* mHistB0EPTwr = new TH1D("mHistB0EPTwr", "E/p", 1000, 0.0, 10.0); mHistB0EPTwr->Sumw2();
TH2D* mHistB0EP2DTwr = new TH2D("mHistB0EP2DTwr", "mHistB0EP2DTwr", 1000, 0.0, 100.0, 1000.0, 0.0, 100.0); mHistB0EP2DTwr->Sumw2();
TH1D* mHistB0DRTwr = new TH1D("mHistB0DRTwr", "mHistB0DRTwr", 500, 0.0, 5.0); mHistB0DRTwr->Sumw2();

// Define a 5D histogram: η, φ, pT, -t, multiplicity
int binsTrk[5] = {40, 64, 200, 200, 20};
double etaRangeTrk[2] = {4.0, 6.0};
double phiRangeTrk[2] = {-M_PI, M_PI};
double ptRangeTrk[2] = {0, 5};
double mtRangeTrk[2] = {0, 20};
double multRangeTrk[2] = {0, 20};

//Reconstructed tracks (for usage with B0 too!!)
// All Tracks
TH1D* mHistB0TrkPx = new TH1D("mHistB0TrkPx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistB0TrkPx->Sumw2();
TH1D* mHistB0TrkPy = new TH1D("mHistB0TrkPy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistB0TrkPy->Sumw2();
TH1D* mHistB0TrkPt = new TH1D("mHistB0TrkPt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistB0TrkPt->Sumw2();
TH1D* mHistB0TrkPz = new TH1D("mHistB0TrkPz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistB0TrkPz->Sumw2();
TH1D* mHistB0TrkTheta = new TH1D("mHistB0TrkTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistB0TrkTheta->Sumw2();
THnSparseD *mHistB04Trk = new THnSparseD("mHistB04Trk","mHistB04Trk;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistB04Trk->Sumw2();

// Hadrons
TH1D* mHistB0HadPx = new TH1D("mHistB0HadPx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistB0HadPx->Sumw2();
TH1D* mHistB0HadPy = new TH1D("mHistB0HadPy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistB0HadPy->Sumw2();
TH1D* mHistB0HadPt = new TH1D("mHistB0HadPt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistB0HadPt->Sumw2();
TH1D* mHistB0HadPz = new TH1D("mHistB0HadPz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistB0HadPz->Sumw2();
TH1D* mHistB0HadTheta = new TH1D("mHistB0HadTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistB0HadTheta->Sumw2();
THnSparseD *mHistB04Had = new THnSparseD("mHistB04Had","mHistB04Had;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistB04Had->Sumw2();

// Electrons
TH1D* mHistB0ElePx = new TH1D("mHistB0ElePx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistB0ElePx->Sumw2();
TH1D* mHistB0ElePy = new TH1D("mHistB0ElePy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistB0ElePy->Sumw2();
TH1D* mHistB0ElePt = new TH1D("mHistB0ElePt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistB0ElePt->Sumw2();
TH1D* mHistB0ElePz = new TH1D("mHistB0ElePz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistB0ElePz->Sumw2();
TH1D* mHistB0EleTheta = new TH1D("mHistB0EleTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistB0EleTheta->Sumw2();
THnSparseD *mHistB04Ele = new THnSparseD("mHistB04Ele","mHistB04Ele;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistB04Ele->Sumw2();

//MC tracks
// All Tracks
TH1D* mHistMCTrkPx = new TH1D("mHistMCTrkPx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistMCTrkPx->Sumw2();
TH1D* mHistMCTrkPy = new TH1D("mHistMCTrkPy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistMCTrkPy->Sumw2();
TH1D* mHistMCTrkPt = new TH1D("mHistMCTrkPt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistMCTrkPt->Sumw2();
TH1D* mHistMCTrkPz = new TH1D("mHistMCTrkPz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistMCTrkPz->Sumw2();
TH1D* mHistMCTrkTheta = new TH1D("mHistMCTrkTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistMCTrkTheta->Sumw2();
THnSparseD *mHistMC4Trk = new THnSparseD("mHistMC4Trk","mHistMC4Trk;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistMC4Trk->Sumw2();

// Hadrons
TH1D* mHistMCHadPx = new TH1D("mHistMCHadPx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistMCHadPx->Sumw2();
TH1D* mHistMCHadPy = new TH1D("mHistMCHadPy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistMCHadPy->Sumw2();
TH1D* mHistMCHadPt = new TH1D("mHistMCHadPt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistMCHadPt->Sumw2();
TH1D* mHistMCHadPz = new TH1D("mHistMCHadPz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistMCHadPz->Sumw2();
TH1D* mHistMCHadTheta = new TH1D("mHistMCHadTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistMCHadTheta->Sumw2();
THnSparseD *mHistMC4Had = new THnSparseD("mHistMC4Had","mHistMC4Had;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistMC4Had->Sumw2();

// Electrons
TH1D* mHistMCElePx = new TH1D("mHistMCElePx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistMCElePx->Sumw2();
TH1D* mHistMCElePy = new TH1D("mHistMCElePy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistMCElePy->Sumw2();
TH1D* mHistMCElePt = new TH1D("mHistMCElePt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistMCElePt->Sumw2();
TH1D* mHistMCElePz = new TH1D("mHistMCElePz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistMCElePz->Sumw2();
TH1D* mHistMCEleTheta = new TH1D("mHistMCEleTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistMCEleTheta->Sumw2();
THnSparseD *mHistMC4Ele = new THnSparseD("mHistMC4Ele","mHistMC4Ele;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistMC4Ele->Sumw2();


//MCTRUTH tracks
// All Tracks
TH1D* mHistMCTRUTHTrkPx = new TH1D("mHistMCTRUTHTrkPx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistMCTRUTHTrkPx->Sumw2();
TH1D* mHistMCTRUTHTrkPy = new TH1D("mHistMCTRUTHTrkPy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistMCTRUTHTrkPy->Sumw2();
TH1D* mHistMCTRUTHTrkPt = new TH1D("mHistMCTRUTHTrkPt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistMCTRUTHTrkPt->Sumw2();
TH1D* mHistMCTRUTHTrkPz = new TH1D("mHistMCTRUTHTrkPz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistMCTRUTHTrkPz->Sumw2();
TH1D* mHistMCTRUTHTrkTheta = new TH1D("mHistMCTRUTHTrkTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistMCTRUTHTrkTheta->Sumw2();
THnSparseD *mHistMCTRUTH4Trk = new THnSparseD("mHistMCTRUTH4Trk","mHistMCTRUTH4Trk;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistMCTRUTH4Trk->Sumw2();

// Hadrons
TH1D* mHistMCTRUTHHadPx = new TH1D("mHistMCTRUTHHadPx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistMCTRUTHHadPx->Sumw2();
TH1D* mHistMCTRUTHHadPy = new TH1D("mHistMCTRUTHHadPy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistMCTRUTHHadPy->Sumw2();
TH1D* mHistMCTRUTHHadPt = new TH1D("mHistMCTRUTHHadPt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistMCTRUTHHadPt->Sumw2();
TH1D* mHistMCTRUTHHadPz = new TH1D("mHistMCTRUTHHadPz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistMCTRUTHHadPz->Sumw2();
TH1D* mHistMCTRUTHHadTheta = new TH1D("mHistMCTRUTHHadTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistMCTRUTHHadTheta->Sumw2();
THnSparseD *mHistMCTRUTH4Had = new THnSparseD("mHistMCTRUTH4Had","mHistMCTRUTH4Had;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistMCTRUTH4Had->Sumw2();

// Electrons
TH1D* mHistMCTRUTHElePx = new TH1D("mHistMCTRUTHElePx", ";p_{x} [GeV/c]", 200, -10.0, 10.0); mHistMCTRUTHElePx->Sumw2();
TH1D* mHistMCTRUTHElePy = new TH1D("mHistMCTRUTHElePy", ";p_{y} [GeV/c]", 200, -10.0, 10.0); mHistMCTRUTHElePy->Sumw2();
TH1D* mHistMCTRUTHElePt = new TH1D("mHistMCTRUTHElePt", ";p_{T} [GeV/c]", 200, 0.0, 5.0); mHistMCTRUTHElePt->Sumw2();
TH1D* mHistMCTRUTHElePz = new TH1D("mHistMCTRUTHElePz", ";p_{z} [GeV/c]", 200, 0.0, 320.0); mHistMCTRUTHElePz->Sumw2();
TH1D* mHistMCTRUTHEleTheta = new TH1D("mHistMCTRUTHEleTheta", ";#theta [mrad]", 200, 0.0, 25.0); mHistMCTRUTHEleTheta->Sumw2();
THnSparseD *mHistMCTRUTH4Ele = new THnSparseD("mHistMCTRUTH4Ele","mHistMCTRUTH4Ele;#eta;#phi;p_{T};Momentum Transfer, -t [GeV^{2}];Multiplicity",5, binsTrk, new double[5]{etaRangeTrk[0], phiRangeTrk[0], ptRangeTrk[0], mtRangeTrk[0], multRangeTrk[0]}, new double[5]{etaRangeTrk[1], phiRangeTrk[1], ptRangeTrk[1], mtRangeTrk[1], multRangeTrk[1]}); mHistMCTRUTH4Ele->Sumw2();

// Efficiency and Fake
TH2D* mHistB0TrkMatch = new TH2D("mHistB0TrkMatch", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0TrkMatch->Sumw2();
TH2D* mHistB0HadMatch = new TH2D("mHistB0HadMatch", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0HadMatch->Sumw2();
TH2D* mHistB0EleMatch = new TH2D("mHistB0EleMatch", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0EleMatch->Sumw2();
TH2D* mHistB0TrkUnMatch = new TH2D("mHistB0TrkUnMatch", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0TrkUnMatch->Sumw2();
TH2D* mHistB0HadUnMatch = new TH2D("mHistB0HadUnMatch", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0HadUnMatch->Sumw2();
TH2D* mHistB0EleUnMatch = new TH2D("mHistB0EleUnMatch", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0EleUnMatch->Sumw2();
TH2D* mHistB0TrkMult = new TH2D("mHistB0TrkMult", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0TrkMult->Sumw2();
TH2D* mHistB0HadMult = new TH2D("mHistB0HadMult", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0HadMult->Sumw2();
TH2D* mHistB0EleMult = new TH2D("mHistB0EleMult", "; p_{T} [GeV/c]", 200, 0.0, 5.0, 20, 0.0, 20.0); mHistB0EleMult->Sumw2();

TH1D* mHistB0DeltaPt = new TH1D("mHistB0DeltaPt", ";#Delta p_{T} [GeV/c]", 200, -2.0, 2.0); mHistB0DeltaPt->Sumw2();
TH2D* mHistB0DeltaPtRes = new TH2D("mHistB0DeltaPtRes", ";P_{T, MC} [GeV/c]; #Delta p_{T}/p_{T, MC} [percent/100]", 100, 0.0, 5.0, 100, -1.0, 1.0); mHistB0DeltaPtRes->Sumw2(); 	
TH2D* mHistB0DeltaPRes = new TH2D("mHistB0DeltaPRes", ";Three-Momentum,  p_{MC} [GeV/c]; #Delta p/p_{MC} [percent/100]", 100, 0.0, 200.0, 100, -1.0, 1.0); mHistB0DeltaPRes->Sumw2(); 
TH2D* mHistB0DeltaPxRes = new TH2D("mHistB0DeltaPxRes", ";P_{x, MC} [GeV/c]; #Delta p_{x}/p_{x, MC} [percent/100]", 100, -10.0, 10.0, 100, -1.0, 1.0); mHistB0DeltaPxRes->Sumw2(); 
TH2D* mHistB0DeltaPyRes = new TH2D("mHistB0DeltaPyRes", ";P_{y, MC} [GeV/c]; #Delta p_{y}/p_{y, MC} [percent/100]", 100, -10.0, 10.0, 100, -1.0, 1.0); mHistB0DeltaPyRes->Sumw2(); 
TH2D* mHistB0DeltaPzRes = new TH2D("mHistB0DeltaPzRes", ";P_{z, MC} [GeV/c]; #Delta p_{z}/p_{z, MC} [percent/100]", 100, -10.0, 10.0, 100, -1.0, 1.0); mHistB0DeltaPzRes->Sumw2(); 

// Particle decays
TH1D *mPhotonConversion = new TH1D("mPhotonConversion","mPhotonConversion",500,0.,1.); mPhotonConversion->Sumw2();
TH1D* mHistK0sSS  = new TH1D("mHistK0sSS","K^{0}_{S} invariant mass;M(GeV);counts",400,0.4,0.6); mHistK0sSS->Sumw2();
TH1D* mHistK0sOS  = new TH1D("mHistK0sOS","K^{0}_{S} invariant mass;M(GeV);counts",400,0.4,0.6); mHistK0sOS->Sumw2(); 
TH1D* mHistRhoSS  = new TH1D("mHistRhoSS","#rho^{0} invariant mass;M(GeV);counts",2000,0.5,1.5); mHistRhoSS->Sumw2(); 
TH1D* mHistRhoOS  = new TH1D("mHistRhoOS","#rho^{0} invariant mass;M(GeV);counts",2000,0.5,1.5); mHistRhoOS->Sumw2(); 
TH1D* mHistPhiSS  = new TH1D("mHistPhiSS","#phi invariant mass;M(GeV);counts",800,0.8,1.2); mHistPhiSS->Sumw2(); 
TH1D* mHistPhiOS  = new TH1D("mHistPhiOS","#phi invariant mass;M(GeV);counts",800,0.8,1.2); mHistPhiOS->Sumw2(); 
TH1D* mHistLambdaSS = new TH1D("mHistLambdaSS","#Lambda invariant mass;M(GeV);counts",400,1.0,1.2); mHistLambdaSS->Sumw2(); 
TH1D* mHistLambdaOS = new TH1D("mHistLambdaOS","#Lambda invariant mass;M(GeV);counts",400,1.0,1.2); mHistLambdaOS->Sumw2(); 
TH1D* mHistALambdaSS = new TH1D("mHistALambdaSS","#bar{#Lambda} invariant mass;M(GeV);counts",400,1.0,1.2); mHistALambdaSS->Sumw2(); 
TH1D* mHistALambdaOS = new TH1D("mHistALambdaOS","#bar{#Lambda} invariant mass;M(GeV);counts",400,1.0,1.2); mHistALambdaOS->Sumw2(); 
TH1D* mHistJpsiSS = new TH1D("mHistJpsiSS","J/#Psi invariant mass;M(GeV);counts",1400,2.7,3.4); mHistJpsiSS->Sumw2(); 
TH1D* mHistJpsiOS = new TH1D("mHistJpsiOS","J/#Psi invariant mass;M(GeV);counts",1400,2.7,3.4); mHistJpsiOS->Sumw2(); 
TH1D* mHistD0SS   = new TH1D("mHistD0SS","D^{0} invariant mass;M(GeV);counts",600,1.7,2.0); mHistD0SS->Sumw2(); 
TH1D* mHistD0OS   = new TH1D("mHistD0OS","D^{0} invariant mass;M(GeV);counts",600,1.7,2.0); mHistD0OS->Sumw2(); 
TH1D* mHistAD0SS   = new TH1D("mHistAD0SS","D^{0} invariant mass;M(GeV);counts",600,1.7,2.0); mHistAD0SS->Sumw2(); 
TH1D* mHistAD0OS   = new TH1D("mHistAD0OS","D^{0} invariant mass;M(GeV);counts",600,1.7,2.0); mHistAD0OS->Sumw2(); 

TH2D *APPlot_K0s=new TH2D("APPlot_K0s","APPlot_K0s",400,-1.,1.,400,0.,1.); APPlot_K0s->Sumw2();
TH2D *APPlot_Lam=new TH2D("APPlot_Lam","APPlot_Lam",400,-1.,1.,400,0.,1.); APPlot_Lam->Sumw2();
TH1D *LamAsK0sMassHypothesis=new TH1D("LamAsK0sMassHypothesis","LamAsK0sMassHypothesis",400,0.4,0.6); LamAsK0sMassHypothesis->Sumw2();
TH1D *K0sAsLamMassHypothesis=new TH1D("K0sAsLamMassHypothesis","K0sAsLamMassHypothesis",400,1.0,1.2); K0sAsLamMassHypothesis->Sumw2();

TH1D* mHistK0sSSMC  = new TH1D("mHistK0sSSMC","K^{0}_{S} invariant mass;M(GeV);counts",400,0.4,0.6); mHistK0sSSMC->Sumw2(); 
TH1D* mHistK0sOSMC  = new TH1D("mHistK0sOSMC","K^{0}_{S} invariant mass;M(GeV);counts",400,0.4,0.6); mHistK0sOSMC->Sumw2(); 
TH1D* mHistRhoSSMC  = new TH1D("mHistRhoSSMC","#rho^{0} invariant mass;M(GeV);counts",1000,0.5,1.0); mHistRhoSSMC->Sumw2(); 
TH1D* mHistRhoOSMC  = new TH1D("mHistRhoOSMC","#rho^{0} invariant mass;M(GeV);counts",1000,0.5,1.0); mHistRhoOSMC->Sumw2(); 
TH1D* mHistPhiSSMC  = new TH1D("mHistPhiSSMC","#phi invariant mass;M(GeV);counts",800,0.8,1.2); mHistPhiSSMC->Sumw2(); 
TH1D* mHistPhiOSMC  = new TH1D("mHistPhiOSMC","#phi invariant mass;M(GeV);counts",800,0.8,1.2); mHistPhiOSMC->Sumw2(); 
TH1D* mHistLambdaSSMC  = new TH1D("mHistLambdaSSMC","#Lambda invariant mass;M(GeV);counts",400,1.0,1.2); mHistLambdaSSMC->Sumw2(); 
TH1D* mHistLambdaOSMC  = new TH1D("mHistLambdaOSMC","#Lambda invariant mass;M(GeV);counts",400,1.0,1.2); mHistLambdaOSMC->Sumw2(); 
TH1D* mHistALambdaSSMC  = new TH1D("mHistALambdaSSMC","#bar{#Lambda} invariant mass;M(GeV);counts",400,1.0,1.2); mHistALambdaSSMC->Sumw2(); 
TH1D* mHistALambdaOSMC  = new TH1D("mHistALambdaOSMC","#bar{#Lambda} invariant mass;M(GeV);counts",400,1.0,1.2); mHistALambdaOSMC->Sumw2(); 
TH1D* mHistJpsiSSMC  = new TH1D("mHistJpsiSSMC","J/#Psi invariant mass;M(GeV);counts",1400,2.7,3.4); mHistJpsiSSMC->Sumw2(); 
TH1D* mHistJpsiOSMC  = new TH1D("mHistJpsiOSMC","J/#Psi invariant mass;M(GeV);counts",1400,2.7,3.4); mHistJpsiOSMC->Sumw2(); 
TH1D* mHistD0SSMC  = new TH1D("mHistD0SSMC","D^{0} invariant mass;M(GeV);counts",600,1.7,2.0); mHistD0SSMC->Sumw2(); 
TH1D* mHistD0OSMC  = new TH1D("mHistD0OSMC","D^{0} invariant mass;M(GeV);counts",600,1.7,2.0); mHistD0OSMC->Sumw2(); 
TH1D* mHistAD0SSMC   = new TH1D("mHistAD0SSMC","D^{0} invariant mass;M(GeV);counts",600,1.7,2.0); mHistAD0SSMC->Sumw2(); 
TH1D* mHistAD0OSMC   = new TH1D("mHistAD0OSMC","D^{0} invariant mass;M(GeV);counts",600,1.7,2.0); mHistAD0OSMC->Sumw2(); 

// For 2PC
// Define a 5D histogram: Δη, Δφ, pT1, pT2, multiplicity
int bins[5] = {40, 64, 100, 100, 20};
double etaRange[2] = {-2.0, 2.0};
double phiRange[2] = {-M_PI/2.0, 3.0*M_PI/2.0};
double ptRange[2] = {0, 5};
double multRange[2] = {0, 20};
THnSparseD *mHist2PC = new THnSparseD("mHist2PC","Two-Particle Correlations;#Delta#eta;#Delta#phi;p_{T1};p_{T2};Multiplicity",5, bins, new double[5]{etaRange[0], phiRange[0], ptRange[0], ptRange[0], multRange[0]}, new double[5]{etaRange[1], phiRange[1], ptRange[1], ptRange[1], multRange[1]}); mHist2PC->Sumw2();
THnSparseD *mHist2PCMix = new THnSparseD("mHist2PCMix","Two-Particle Correlations Mix;#Delta#eta;#Delta#phi;p_{T1};p_{T2};Multiplicity",5, bins, new double[5]{etaRange[0], phiRange[0], ptRange[0], ptRange[0], multRange[0]}, new double[5]{etaRange[1], phiRange[1], ptRange[1], ptRange[1], multRange[1]}); mHist2PCMix->Sumw2();

TH1D* hTrigger_0_5_1_0 = new TH1D("hTrigger_0_5_1_0", "nTrigger for 0.5 < pT < 1.0 GeV/c", 2, 0, 2); hTrigger_0_5_1_0->Sumw2();
TH1D* hTrigger_1_0_1_5 = new TH1D("hTrigger_1_0_1_5", "nTrigger for 1.0 < pT < 1.5 GeV/c", 2, 0, 2); hTrigger_1_0_1_5->Sumw2();
TH1D* hTrigger_1_5_2_0 = new TH1D("hTrigger_1_5_2_0", "nTrigger for 1.5 < pT < 2.0 GeV/c", 2, 0, 2); hTrigger_1_5_2_0->Sumw2();
TH1D* hTrigger_2_0_2_5 = new TH1D("hTrigger_2_0_2_5", "nTrigger for 2.0 < pT < 2.5 GeV/c", 2, 0, 2); hTrigger_2_0_2_5->Sumw2();
TH1D* hTrigger_2_5_3_0 = new TH1D("hTrigger_2_5_3_0", "nTrigger for 2.5 < pT < 3.0 GeV/c", 2, 0, 2); hTrigger_2_5_3_0->Sumw2();
TH1D* hTrigger_3_0_4_0 = new TH1D("hTrigger_3_0_4_0", "nTrigger for 3.0 < pT < 4.0 GeV/c", 2, 0, 2); hTrigger_3_0_4_0->Sumw2();
TH1D* hTriggerAbove4 = new TH1D("hTriggerAbove4", "nTrigger for pT > 4.0 GeV/c", 2, 0, 2); hTriggerAbove4->Sumw2();

void WriteHistosEvtDet(){

	// Some event quantities
	mHistNEvents->Write();
	mHistNTrks->Write();
	mHistMultiplicity->Write();	
	mHistEMultiplicity->Write();	
	mHistNTrksMC->Write();
	mHistMultiplicityMC->Write();	
	mHistEMultiplicityMC->Write();	
	mHistNTrksTRUTHMC->Write();
	mHistMultiplicityTRUTHMC->Write();	
	mHistEMultiplicityTRUTHMC->Write();	
	mHistNEvMix->Write();
	
	// Detector quantities
	mHistB0OccLayer0->Write();
	mHistB0OccLayer1->Write();
	mHistB0OccLayer2->Write();
	mHistB0OccLayer3->Write();
	mHistB0HitEDep->Write();
	mHistB0Ecal->Write();
	mHistB0EcalEDep->Write();	
	mHistB0EPTrk->Write();
	mHistB0EP2DTrk->Write();
	mHistB0DRTrk->Write();
	mHistB0EPTwr->Write();
	mHistB0EP2DTwr->Write();
	mHistB0DRTwr->Write();
	
}
	
void WriteHistosTracks(){

	// Track quantities	
	mHistB0TrkPx->Write();
	mHistB0TrkPy->Write();
	mHistB0TrkPt->Write();
	mHistB0TrkPz->Write();
	mHistB0TrkTheta->Write();
	mHistB04Trk->Write();
	// Hadrons
	mHistB0HadPx->Write();
	mHistB0HadPy->Write();
	mHistB0HadPt->Write();
	mHistB0HadPz->Write();
	mHistB0HadTheta->Write();
	mHistB04Had->Write();
	// Electrons
	mHistB0ElePx->Write();
	mHistB0ElePy->Write();
	mHistB0ElePt->Write();
	mHistB0ElePz->Write();
	mHistB0EleTheta->Write();
	mHistB04Ele->Write();	

	// For MC
	// Track quantities	
	mHistMCTrkPx->Write();
	mHistMCTrkPy->Write();
	mHistMCTrkPt->Write();
	mHistMCTrkPz->Write();
	mHistMCTrkTheta->Write();
	mHistMC4Trk->Write();
	// Hadrons
	mHistMCHadPx->Write();
	mHistMCHadPy->Write();
	mHistMCHadPt->Write();
	mHistMCHadPz->Write();
	mHistMCHadTheta->Write();
	mHistMC4Had->Write();
	// Electrons
	mHistMCElePx->Write();
	mHistMCElePy->Write();
	mHistMCElePt->Write();
	mHistMCElePz->Write();
	mHistMCEleTheta->Write();
	mHistMC4Ele->Write();	
	
	// For MCTRUTH
	// Track quantities	
	mHistMCTRUTHTrkPx->Write();
	mHistMCTRUTHTrkPy->Write();
	mHistMCTRUTHTrkPt->Write();
	mHistMCTRUTHTrkPz->Write();
	mHistMCTRUTHTrkTheta->Write();
	mHistMCTRUTH4Trk->Write();
	// Hadrons
	mHistMCTRUTHHadPx->Write();
	mHistMCTRUTHHadPy->Write();
	mHistMCTRUTHHadPt->Write();
	mHistMCTRUTHHadPz->Write();
	mHistMCTRUTHHadTheta->Write();
	mHistMCTRUTH4Had->Write();
	// Electrons
	mHistMCTRUTHElePx->Write();
	mHistMCTRUTHElePy->Write();
	mHistMCTRUTHElePt->Write();
	mHistMCTRUTHElePz->Write();
	mHistMCTRUTHEleTheta->Write();
	mHistMCTRUTH4Ele->Write();	

	// For efficiency
	mHistB0TrkMatch->Write();
	mHistB0HadMatch->Write();
	mHistB0EleMatch->Write();
	mHistB0TrkUnMatch->Write();
	mHistB0HadUnMatch->Write();
	mHistB0EleUnMatch->Write();
	mHistB0TrkMult->Write();
	mHistB0HadMult->Write();
	mHistB0EleMult->Write();
	
	// For resolution
	mHistB0DeltaPt->Write();
	mHistB0DeltaPtRes->Write();
	mHistB0DeltaPxRes->Write();
	mHistB0DeltaPyRes->Write();
	mHistB0DeltaPzRes->Write();
	mHistB0DeltaPRes->Write();

}

void WriteHistosDecay(){

	mPhotonConversion->Write();
	mHistJpsiSS->Write();
	mHistJpsiOS->Write();
	mHistRhoSS->Write();
	mHistRhoOS->Write();
	mHistPhiOS->Write();
	mHistPhiSS->Write();
	mHistD0SS->Write();
	mHistD0OS->Write();
	mHistAD0SS->Write();
	mHistAD0OS->Write();
	mHistK0sSS->Write();
	mHistK0sOS->Write();
	APPlot_K0s->Write();
	K0sAsLamMassHypothesis->Write();
	mHistLambdaSS->Write();
	mHistLambdaOS->Write();
	mHistALambdaSS->Write();
	mHistALambdaOS->Write();
	APPlot_Lam->Write();
	LamAsK0sMassHypothesis->Write();

	mHistJpsiSSMC->Write();
	mHistJpsiOSMC->Write();
	mHistRhoSSMC->Write();
	mHistRhoOSMC->Write();
	mHistPhiSSMC->Write();
	mHistPhiOSMC->Write();
	mHistD0SSMC->Write();
	mHistD0OSMC->Write();
	mHistAD0SSMC->Write();
	mHistAD0OSMC->Write();	
	mHistK0sSSMC->Write();
	mHistK0sOSMC->Write();
	mHistLambdaSSMC->Write();
	mHistLambdaOSMC->Write();
	mHistALambdaSSMC->Write();
	mHistALambdaOSMC->Write();

}

void WriteHistos2PC(){

	// 2PC quantities			
	mHist2PC->Write();
	mHist2PCMix->Write();

	hTrigger_0_5_1_0->Write();
	hTrigger_1_0_1_5->Write();
	hTrigger_1_5_2_0->Write();
	hTrigger_2_0_2_5->Write();
	hTrigger_2_5_3_0->Write();
	hTrigger_3_0_4_0->Write();
	hTriggerAbove4->Write();

}
