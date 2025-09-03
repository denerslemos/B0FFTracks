#include "CallLibraries.h"  // call libraries from ROOT and C++

// Reconstructed charged particles (Tracks)
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPx;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPy;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoPz;
std::unique_ptr<TTreeReaderArray<int>> TrkRecoPDG;
std::unique_ptr<TTreeReaderArray<float>> TrkRecoCharge;

// MCParticlesHeadOnFrameNoBeamFX
std::unique_ptr<TTreeReaderArray<double>> TrkMCTRUTHGenPx;
std::unique_ptr<TTreeReaderArray<double>> TrkMCTRUTHGenPy;
std::unique_ptr<TTreeReaderArray<double>> TrkMCTRUTHGenPz;
std::unique_ptr<TTreeReaderArray<double>> TrkMCTRUTHGenE;
std::unique_ptr<TTreeReaderArray<double>> TrkMCTRUTHGenM;
std::unique_ptr<TTreeReaderArray<int>> TrkMCTRUTHGenStatus;
std::unique_ptr<TTreeReaderArray<int>> TrkMCTRUTHGenPDG;
std::unique_ptr<TTreeReaderArray<float>> TrkMCTRUTHGenCharge;

// MCParticles particles
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPx;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPy;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenPz;
std::unique_ptr<TTreeReaderArray<double>> TrkMCGenM;
std::unique_ptr<TTreeReaderArray<int>> TrkMCGenPDG;
std::unique_ptr<TTreeReaderArray<int>> TrkMCGenStatus;
std::unique_ptr<TTreeReaderArray<float>> TrkMCGenCharge;

// Other MC info
std::unique_ptr<TTreeReaderArray<unsigned int>> AssTSCPRec;
std::unique_ptr<TTreeReaderArray<unsigned int>> AssTSCPSim;
std::unique_ptr<TTreeReaderArray<float>> AssTSCPWgt;

// B0 Tracker quantities
std::unique_ptr<TTreeReaderArray<float>> B0HitsX;
std::unique_ptr<TTreeReaderArray<float>> B0HitsY;
std::unique_ptr<TTreeReaderArray<float>> B0HitsZ;
std::unique_ptr<TTreeReaderArray<float>> B0HitsEDep;
std::unique_ptr<TTreeReaderArray<float>> B0MCMomX;
std::unique_ptr<TTreeReaderArray<float>> B0MCMomY;
std::unique_ptr<TTreeReaderArray<float>> B0MCMomZ;

// B0 ECAL quantities
std::unique_ptr<TTreeReaderArray<float>> B0ECalX;
std::unique_ptr<TTreeReaderArray<float>> B0ECalY;
std::unique_ptr<TTreeReaderArray<float>> B0ECalZ;
std::unique_ptr<TTreeReaderArray<float>> B0ECalEDep;

// Reads B0 tree variables
/*
Arguments:
chain: TChain of input files
tree_reader: the object for TTreeReader
*/
void B0FFTrackTreeReader(TChain* chain, std::unique_ptr<TTreeReader>& tree_reader) {

    tree_reader = std::make_unique<TTreeReader>(chain);

	// Reconstructed charged particles (Tracks)
    TrkRecoPx = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TrkRecoPy = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TrkRecoPz = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z");
    TrkRecoCharge = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedTruthSeededChargedParticles.charge");
    TrkRecoPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "ReconstructedTruthSeededChargedParticles.PDG");
    
	// Track matching
	AssTSCPRec = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID");
    AssTSCPSim = std::make_unique<TTreeReaderArray<unsigned int>>(*tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID");
    AssTSCPWgt = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.weight");

	// MC Gen -> add GEANT Stuff
	TrkMCGenPx = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.x");
	TrkMCGenPy = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.y");
	TrkMCGenPz = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.momentum.z");
	TrkMCGenM = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticles.mass");
    TrkMCGenCharge = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "MCParticles.charge");
	TrkMCGenPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticles.PDG");
	TrkMCGenStatus = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticles.generatorStatus");

	// Full gen level
	TrkMCTRUTHGenPx = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.x");
	TrkMCTRUTHGenPy = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.y");
	TrkMCTRUTHGenPz = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticlesHeadOnFrameNoBeamFX.momentum.z");
	TrkMCTRUTHGenM = std::make_unique<TTreeReaderArray<double>>(*tree_reader, "MCParticlesHeadOnFrameNoBeamFX.mass");
    TrkMCTRUTHGenCharge = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "MCParticlesHeadOnFrameNoBeamFX.charge");
	TrkMCTRUTHGenPDG = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticlesHeadOnFrameNoBeamFX.PDG");
	TrkMCTRUTHGenStatus = std::make_unique<TTreeReaderArray<int>>(*tree_reader, "MCParticlesHeadOnFrameNoBeamFX.generatorStatus");
					
	//b0 tracker hits
	B0HitsX = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0TrackerRecHits.position.x");
	B0HitsY = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0TrackerRecHits.position.y");
	B0HitsZ = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0TrackerRecHits.position.z");
	B0HitsEDep = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0TrackerRecHits.edep"); //deposited energy per hit
	B0MCMomX = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0TrackerHits.momentum.x");
	B0MCMomY = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0TrackerHits.momentum.y");
	B0MCMomZ = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0TrackerHits.momentum.z");

	//b0 EMCAL
	B0ECalX = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0ECalClusters.position.x");
   	B0ECalY = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0ECalClusters.position.y");
   	B0ECalZ = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0ECalClusters.position.z");
	B0ECalEDep = std::make_unique<TTreeReaderArray<float>>(*tree_reader, "B0ECalClusters.energy"); //deposited energy in cluster
			
}
