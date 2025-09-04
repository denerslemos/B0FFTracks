/*

Simple analysis code for tracks at B0 far-forward detector	 

Author: Dener De Souza Lemos

Last update: Sept 4th 2025

*/
#include "B0FFTrackTreeReader.h"  // call tree reader
#include "B0FFTrackHistogramDefinition.h" // call for histograms
#include "B0FFTrackFunctions.h" // call functions

void B0FFTrackAnalyzer(TString InputFileList, TString OutputFile, bool DoMixing){

	// Read the list of input file(s)
	fstream FileList;
	FileList.open(Form("%s",InputFileList.Data()), ios::in);
	if(!FileList.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << InputFileList.Data() << endl;}
	// Make a chain and a vector of file names
	std::vector<TString> FileListVector;
	string FileChain;
	while(getline(FileList, FileChain)){FileListVector.push_back(FileChain.c_str());}
	FileList.close();	
	TChain *mychain = new TChain("events");
	for (std::vector<TString>::iterator listIterator = FileListVector.begin(); listIterator != FileListVector.end(); listIterator++){
		TFile *testfile = TFile::Open(*listIterator,"READ");
		if(testfile && !testfile->IsZombie() && !testfile->TestBit(TFile::kRecovered)){ // safety against corrupted files
			cout << "Adding file " << *listIterator << " to the chains" << endl; // adding files to the chains for each step
			mychain->Add(*listIterator);
		}else{cout << "File: " << *listIterator << " failed!" << endl;}
	}

	std::unique_ptr<TTreeReader> tree_reader;
	B0FFTrackTreeReader(mychain, tree_reader); 
	
	// Output
	TFile *OutFile = TFile::Open(Form("%s",OutputFile.Data()),"RECREATE");	
	
	// Loop Through Events
	int NEVENTS = 0;

	vector < vector< TVector3 > > B0HadrMix; 
	vector < int > B0MultMix; 

	while(tree_reader->Next()) {	
	
		mHistNEvents->Fill(1);
	    if(NEVENTS%10000 == 0) cout << "Events Processed: " << NEVENTS << endl;

		float HitX = -9999.;
		float HitY = -9999.;
		float HitZ = -9999.;
		float HitE = -9999.;

		//b0 clusters
		for(int ib0cls = 0; ib0cls < B0ECalX->GetSize(); ib0cls++){	
			HitX = (*B0ECalX)[ib0cls];
			HitY = (*B0ECalX)[ib0cls];
			HitZ = (*B0ECalX)[ib0cls];
			HitE = (*B0ECalEDep)[ib0cls]*1.246; //poor man's calibration constant, for now							
			mHistB0Ecal->Fill(HitX, HitY);
			mHistB0EcalEDep->Fill(HitE);				
		}

		//b0 tracker hits
		for(int ib0hit = 0; ib0hit < B0HitsX->GetSize(); ib0hit++){
			HitX = (*B0HitsX)[ib0hit];
			HitY = (*B0HitsY)[ib0hit];
			HitZ = (*B0HitsZ)[ib0hit];
			HitE = (*B0HitsEDep)[ib0hit]*1e6; //convert GeV --> keV
			mHistB0HitEDep->Fill(HitE);		
			//ITS3 Layout
			//if(hit_z > 5800 && hit_z < 6000){ h_B0_occupancy_map_layer_0->Fill(hit_x, hit_y); }
			//if(hit_z > 6000 && hit_z < 6200){ h_B0_occupancy_map_layer_1->Fill(hit_x, hit_y); }
			//if(hit_z > 6200 && hit_z < 6400){ h_B0_occupancy_map_layer_2->Fill(hit_x, hit_y); }
			//if(hit_z > 6400 && hit_z < 6600){ h_B0_occupancy_map_layer_3->Fill(hit_x, hit_y); }		
			//ACLGAD layout
			if(HitZ > 5700 && HitZ < 5990){ mHistB0OccLayer0->Fill(HitX, HitY); }
			if(HitZ > 6100 && HitZ < 6200){ mHistB0OccLayer1->Fill(HitX, HitY); }
			if(HitZ > 6400 && HitZ < 6500){ mHistB0OccLayer2->Fill(HitX, HitY); }
			if(HitZ > 6700 && HitZ < 6750){ mHistB0OccLayer3->Fill(HitX, HitY); }			
		}
		
		// Tracks start here
		// Find matching reco and gen for fake rate
        vector<bool> isMatched(TrkRecoPx->GetSize(), false);  
        vector<int> numRecperGen(TrkMCGenPx->GetSize(), 0);

        for(size_t iassoc = 0; iassoc < AssTSCPWgt->GetSize(); iassoc++){
            if( (*AssTSCPWgt)[iassoc] > 0.8 ) {
                if ( (*AssTSCPSim)[iassoc] < numRecperGen.size() ) { numRecperGen[(*AssTSCPSim)[iassoc]]++; }
                if ( (*AssTSCPRec)[iassoc] < isMatched.size() ) { isMatched[(*AssTSCPRec)[iassoc]] = true; } 
            }
        }

		//Loop of gen level tracks: MCParticles (following Alex code)
        std::vector<bool> isGenFound(TrkMCGenStatus->GetSize(), false);
		int gntrk = 0;
		int gmult = 0;
		int gemult = 0;
    	vector< TVector3 > MCTrk; 
    	vector< float > MCTrkCharge; 
    	vector< int > MCPDGId; 
    	vector< TVector3 > MCHadr; 
    	vector< float > MCHadrCharge; 
    	vector< TVector3 > MCElec; 
    	vector< float > MCElecCharge; 
    	
		for(size_t imc = 0; imc < TrkMCGenStatus->GetSize(); imc++){

			if ( (*TrkMCGenStatus)[imc] != 1 ) continue; // only stable ones
			if ( fabs((*TrkMCGenCharge)[imc]) != 1 ) continue; // only charged ones
			TVector3 TrkMCGen((*TrkMCGenPx)[imc], (*TrkMCGenPy)[imc], (*TrkMCGenPz)[imc]); 
			TrkMCGen.RotateY(0.025);
    		double GenTrkEta = TrkMCGen.Eta();
    		double GenTrkPhi = TrkMCGen.Phi();
	    	double GenTrkPt  = TrkMCGen.Perp();
	    	double GenTrkTheta = 1000 * TrkMCGen.Theta();
	        double GenTrkP = TrkMCGen.Mag();
			// B0 acceptance
    	    if ( GenTrkP <= 0 ) continue;
			if ( GenTrkPt <= 0 ) continue;    
			if ( GenTrkTheta < 5.5 || GenTrkTheta > 20.0 ) continue;
			
			gntrk = gntrk + 1;
			MCTrk.push_back(TrkMCGen);
			MCTrkCharge.push_back((*TrkMCGenCharge)[imc]);
			MCPDGId.push_back((*TrkMCGenPDG)[imc]);
			
			if( fabs((*TrkMCGenPDG)[imc]) == 11 ) { // electrons
				gemult = gemult + 1;
				MCElec.push_back(TrkMCGen);
				MCElecCharge.push_back((*TrkMCGenCharge)[imc]);
			}

			if( fabs((*TrkMCGenPDG)[imc]) == 211 || fabs((*TrkMCGenPDG)[imc]) == 321 || fabs((*TrkMCGenPDG)[imc]) == 2212 ) { // charged hadrons: pi, kaon, proton
				gmult = gmult + 1;		
				MCHadr.push_back(TrkMCGen);
				MCHadrCharge.push_back((*TrkMCGenCharge)[imc]);
			}
			
			if (numRecperGen[imc] > 0) { 
				// Find a high-quality match for resolution, but only once per gen particle
				if(!isGenFound[imc]){
					for(size_t iassoc = 0; iassoc < AssTSCPWgt->GetSize(); iassoc++){
						if( (*AssTSCPSim)[iassoc] == imc && (*AssTSCPWgt)[iassoc] > 0.8){
							auto recIndex = (*AssTSCPRec)[iassoc];
							if (recIndex < TrkRecoPx->GetSize()) {
								TVector3 RecVec((*TrkRecoPx)[recIndex], (*TrkRecoPy)[recIndex], (*TrkRecoPz)[recIndex]);
								RecVec.RotateY(0.025);	
								double RecTheta = 1000 * RecVec.Theta();			
    	   							if ( RecVec.Mag()  <= 0 ) continue;
								if ( RecVec.Perp() <= 0 ) continue;    
								if ( fabs((*TrkRecoCharge)[recIndex]) != 1 ) continue;
								if (  RecTheta < 5.5 || RecTheta > 20.0 ) continue;
								double delPt = RecVec.Perp() - TrkMCGen.Perp();
								double delP = RecVec.Mag() - TrkMCGen.Mag();
                                				double delPx = RecVec.Px() - TrkMCGen.Px();
                                				double delPy = RecVec.Py() - TrkMCGen.Py();
                                				double delPz = RecVec.Pz() - TrkMCGen.Pz();

                                				double delPt_percent = delPt/TrkMCGen.Perp();
                                				double delP_percent = delP/TrkMCGen.Mag();
                                				double delPx_percent = delPx/TrkMCGen.Px();
                                				double delPy_percent = delPy/TrkMCGen.Py();
                                				double delPz_percent = delPz/TrkMCGen.Pz();

 								mHistB0DeltaPt->Fill(delPt_percent);
								mHistB0DeltaPtRes->Fill(TrkMCGen.Perp(), delPt_percent);
								mHistB0DeltaPxRes->Fill(TrkMCGen.Px(), delPx_percent);
								mHistB0DeltaPyRes->Fill(TrkMCGen.Py(), delPy_percent);
								mHistB0DeltaPzRes->Fill(TrkMCGen.Pz(), delPz_percent);
								mHistB0DeltaPRes->Fill(TrkMCGen.Mag(), delP_percent);
								isGenFound[imc] = true;
							}
							break;
						}
					}
				}
			} 
        } 
  
		//Loop of truth level tracks: MCParticlesHeadOnFrameNoBeamFX (following Alex code)
		int gtntrk = 0;
		int gtmult = 0;
		int gtemult = 0;
    	vector< TVector3 > MCtTrk; 
    	vector< TVector3 > MCtHadr; 
    	vector< float > MCtHadrCharge; 
    	vector< TVector3 > MCtElec; 
    	vector< float > MCtElecCharge; 
    	
		for(size_t imct = 0; imct < TrkMCTRUTHGenStatus->GetSize(); imct++){

			if ( (*TrkMCTRUTHGenStatus)[imct] != 1 ) continue; // only stable ones
			if ( fabs((*TrkMCTRUTHGenCharge)[imct]) != 1 ) continue; // only charged ones
			TVector3 TrkMCTRUTHGen((*TrkMCTRUTHGenPx)[imct], (*TrkMCTRUTHGenPy)[imct], (*TrkMCTRUTHGenPz)[imct]); 
    		double GenTrkEta = TrkMCTRUTHGen.Eta();
    		double GenTrkPhi = TrkMCTRUTHGen.Phi();
	    	double GenTrkPt  = TrkMCTRUTHGen.Perp();
	    	double GenTrkTheta = 1000 * TrkMCTRUTHGen.Theta();
	        double GenTrkP = TrkMCTRUTHGen.Mag();
			// B0 acceptance
    	    if ( GenTrkP <= 0 ) continue;
			if ( GenTrkPt <= 0 ) continue;    
			if ( GenTrkTheta < 5.5 || GenTrkTheta > 20.0 ) continue;
			
			gtntrk = gtntrk + 1;
			MCtTrk.push_back(TrkMCTRUTHGen);
			
			if( fabs((*TrkMCGenPDG)[imct]) == 11 ) { // electrons
				gtemult = gtemult + 1;
				MCtElec.push_back(TrkMCTRUTHGen);
				MCtElecCharge.push_back((*TrkMCTRUTHGenCharge)[imct]);
			}

			if( fabs((*TrkMCGenPDG)[imct]) == 211 || fabs((*TrkMCGenPDG)[imct]) == 321 || fabs((*TrkMCGenPDG)[imct]) == 2212 ) { // charged hadrons: pi, kaon, proton
				gtmult = gtmult + 1;		
				MCtHadr.push_back(TrkMCTRUTHGen);
				MCtHadrCharge.push_back((*TrkMCTRUTHGenCharge)[imct]);
			}
			
        } 		
        
        
		//Loop of reconstructed tracks
		int ntrk = 0;
		int mult = 0;
		int emult = 0;
		std::vector<int> towerAssigned(B0ECalX->GetSize(), -1);
		std::vector<float> towerBestDR(B0ECalX->GetSize(), 999.9); // large initial ΔR
	    // find best match tower
	    float bestDR = 999.9;
    	int bestCalo = -1;
    	vector< TVector3 > B0Trk; 
    	vector< bool > B0TrkMatch; 
    	vector< bool > B0TrkMult; 
    	vector< int > B0TrkFlavour; 
    	vector< float > B0TrkCharge; 
    	vector< TVector3 > B0Hadr; 
    	vector< float > B0HadrCharge; 
    	vector< bool > B0HadrMatch; 
    	vector< bool > B0HadrMult; 
    	vector< TVector3 > B0Elec; 
    	vector< float > B0ElecCharge; 
    	vector< bool > B0ElecMatch; 
    	vector< bool > B0ElecMult; 

		for ( size_t itrk = 0; itrk < TrkRecoPx->GetSize(); itrk++ ){
			TVector3 TrkReco((*TrkRecoPx)[itrk], (*TrkRecoPy)[itrk], (*TrkRecoPz)[itrk]);
			TrkReco.RotateY(0.025);
    		double RecoTrkEta = TrkReco.Eta();
    		double RecoTrkPhi = TrkReco.Phi();
	    	double RecoTrkPt  = TrkReco.Perp();
	    	double RecoTrkTheta = 1000 * TrkReco.Theta();
	        double RecoTrkP = TrkReco.Mag();
			// B0 acceptance
    	    if ( RecoTrkP <= 0 ) continue;
			if ( RecoTrkPt <= 0 ) continue;    
			if ( fabs((*TrkRecoCharge)[itrk]) != 1 ) continue;
			if ( RecoTrkTheta < 5.5 || RecoTrkTheta > 20.0 ) continue;
			for ( size_t jb0cls = 0; jb0cls < B0ECalX->GetSize(); jb0cls++ ) {				
				float ECalEta = computeEta((*B0ECalX)[jb0cls], (*B0ECalY)[jb0cls], (*B0ECalZ)[jb0cls]);
				float ECalPhi = computePhi((*B0ECalX)[jb0cls], (*B0ECalY)[jb0cls]);
				float dR = deltaR(RecoTrkEta, RecoTrkPhi, ECalEta, ECalPhi);
				if (dR < bestDR) { 
					towerAssigned[jb0cls] = itrk; 
            		towerBestDR[jb0cls]   = dR;
		            bestDR = dR;
        		    bestCalo = jb0cls;					
				}
			}
			
			// Compute E/p for track method
			bool isHadron = false;
			bool isElectron = false;
			if ( bestCalo >= 0 ) {	
				double eOverP = (*B0ECalEDep)[bestCalo] / RecoTrkP;
				mHistB0EPTrk->Fill( eOverP ); 
				mHistB0DRTrk->Fill(bestDR); 
				mHistB0EP2DTrk->Fill(RecoTrkP, (*B0ECalEDep)[bestCalo]);
				if( eOverP < 0.8 ) isHadron = true;
				if( eOverP > 0.8 && eOverP < 1.2 ) isElectron = true; 
			} else { isHadron = true; }

			ntrk = ntrk + 1;				
			B0Trk.push_back(TrkReco);
			B0TrkCharge.push_back((*TrkRecoCharge)[itrk]);
			if( isHadron ) { B0TrkFlavour.push_back(1); } else if( isElectron ) { B0TrkFlavour.push_back(2); } else { B0TrkFlavour.push_back(0); }
			B0TrkMatch.push_back(isMatched[itrk]);
            int genId = -1;
            bool isMultipleRec = false;
			for(size_t iassoc = 0; iassoc < AssTSCPWgt->GetSize(); iassoc++){
				if((*AssTSCPRec)[iassoc] == itrk && (*AssTSCPWgt)[iassoc] > 0.8) {
					genId = (*AssTSCPSim)[iassoc];
					if(genId != -1 && numRecperGen[genId] > 1){ isMultipleRec = true; } // If this gen particle was reconstructed more than once, this reco track is a duplicate.
					break;
				}
			}
			B0TrkMult.push_back(isMultipleRec);

			if( isHadron ){
				mult = mult + 1;
				B0Hadr.push_back(TrkReco);
				B0HadrCharge.push_back((*TrkRecoCharge)[itrk]);
				B0HadrMatch.push_back(isMatched[itrk]);
				B0HadrMult.push_back(isMultipleRec);
			}
			
			if( isElectron ){
				B0Elec.push_back(TrkReco);
				B0ElecCharge.push_back((*TrkRecoCharge)[itrk]);			
				B0ElecMatch.push_back(isMatched[itrk]);
				B0ElecMult.push_back(isMultipleRec);
			}

		}

		// Now loop again to compute E/p for tower method
		for (size_t jb0cls = 0; jb0cls < B0ECalX->GetSize(); jb0cls++) {
			int trkIdx = towerAssigned[jb0cls];
			if (trkIdx == -1) continue; // no matched track
			TVector3 TrkReco((*TrkRecoPx)[trkIdx], (*TrkRecoPy)[trkIdx], (*TrkRecoPz)[trkIdx]);
			float matchedE = (*B0ECalEDep)[jb0cls];
			mHistB0EPTwr->Fill(matchedE / TrkReco.Mag());
			mHistB0EP2DTwr->Fill(TrkReco.Mag(), matchedE);
			mHistB0DRTwr->Fill(towerBestDR[jb0cls]);
		}
		
		// lets loop over all particles to add some histograms
		FillBasicHistos(B0Trk, B0TrkMatch, B0TrkMult, mHistB0TrkPx, mHistB0TrkPy, mHistB0TrkPz, mHistB0TrkPt, mHistB0TrkTheta, mHistB04Trk, mHistB0TrkUnMatch, mHistB0TrkMatch, mHistB0TrkMult);
		FillBasicHistos(B0Hadr, B0HadrMatch, B0HadrMult, mHistB0HadPx, mHistB0HadPy, mHistB0HadPz, mHistB0HadPt, mHistB0HadTheta, mHistB04Had, mHistB0HadUnMatch, mHistB0HadMatch, mHistB0HadMult);
		FillBasicHistos(B0Elec, B0ElecMatch, B0ElecMult, mHistB0ElePx, mHistB0ElePy, mHistB0ElePz, mHistB0ElePt, mHistB0EleTheta, mHistB04Ele, mHistB0EleUnMatch, mHistB0EleMatch, mHistB0EleMult);
		// For MC
		FillBasicMCHistos(MCTrk, mHistMCTrkPx, mHistMCTrkPy, mHistMCTrkPz, mHistMCTrkPt, mHistMCTrkTheta, mHistMC4Trk);
		FillBasicMCHistos(MCHadr, mHistMCHadPx, mHistMCHadPy, mHistMCHadPz, mHistMCHadPt, mHistMCHadTheta, mHistMC4Had);
		FillBasicMCHistos(MCElec, mHistMCElePx, mHistMCElePy, mHistMCElePz, mHistMCElePt, mHistMCEleTheta, mHistMC4Ele);
		// For MCTRUTH
		FillBasicMCHistos(MCtTrk, mHistMCTRUTHTrkPx, mHistMCTRUTHTrkPy, mHistMCTRUTHTrkPz, mHistMCTRUTHTrkPt, mHistMCTRUTHTrkTheta, mHistMCTRUTH4Trk);
		FillBasicMCHistos(MCtHadr, mHistMCTRUTHHadPx, mHistMCTRUTHHadPy, mHistMCTRUTHHadPz, mHistMCTRUTHHadPt, mHistMCTRUTHHadTheta, mHistMCTRUTH4Had);
		FillBasicMCHistos(MCtElec, mHistMCTRUTHElePx, mHistMCTRUTHElePy, mHistMCTRUTHElePz, mHistMCTRUTHElePt, mHistMCTRUTHEleTheta, mHistMCTRUTH4Ele);

		// Reconstructed
		InvMassPhys( B0Trk, B0TrkCharge, B0TrkFlavour, mPhotonConversion, mHistJpsiOS,  mHistJpsiSS, mHistRhoOS, mHistRhoSS, mHistPhiOS, mHistPhiSS, mHistD0OS, mHistD0SS, mHistAD0OS, mHistAD0OS, mHistK0sOS, mHistK0sSS, mHistLambdaOS, mHistLambdaSS, mHistALambdaOS, mHistALambdaSS, APPlot_K0s, K0sAsLamMassHypothesis, APPlot_Lam, LamAsK0sMassHypothesis);

		// MC
		InvMassMC( MCTrk , MCTrkCharge, MCPDGId, mHistJpsiOSMC, mHistJpsiSSMC, mHistRhoOSMC, mHistRhoSSMC, mHistPhiOSMC, mHistPhiSSMC, mHistD0OSMC, mHistD0SSMC, mHistAD0OSMC, mHistAD0SSMC, mHistK0sOSMC, mHistK0sSSMC, mHistLambdaOSMC, mHistLambdaSSMC, mHistALambdaOSMC, mHistALambdaSSMC);
						
		// 2PC -- to be moved as a function
		std::vector<int> indices(B0Hadr.size());
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(), [&](int a, int b) { return B0Hadr[a].Perp() > B0Hadr[b].Perp(); }); // sort by pT	
		std::vector<TVector3> B0Hadr_sorted;
		std::vector<int> B0HadrCharge_sorted;
		for (size_t idx : indices) { B0Hadr_sorted.push_back(B0Hadr[idx]); B0HadrCharge_sorted.push_back(B0HadrCharge[idx]); }
		if( B0Hadr_sorted.size() > 1 ){
			TwoPC(B0Hadr_sorted, mHist2PC, mult);
			if ( DoMixing ) B0HadrMix.push_back(B0Hadr_sorted); B0MultMix.push_back(mult); 
		}
				
		for (const auto &trk : B0Hadr) { 
			if ( trk.Perp() >= 0.5 && trk.Perp() < 1.0 ) hTrigger_0_5_1_0->Fill(1);
			if ( trk.Perp() >= 1.0 && trk.Perp() < 1.5 ) hTrigger_1_0_1_5->Fill(1);
			if ( trk.Perp() >= 1.5 && trk.Perp() < 2.0 ) hTrigger_1_5_2_0->Fill(1);
			if ( trk.Perp() >= 2.0 && trk.Perp() < 2.5 ) hTrigger_2_0_2_5->Fill(1);
			if ( trk.Perp() >= 2.5 && trk.Perp() < 3.0 ) hTrigger_2_5_3_0->Fill(1);
			if ( trk.Perp() >= 3.0 && trk.Perp() < 3.5 )  hTrigger_3_0_4_0->Fill(1);
			if ( trk.Perp() >= 4.0 ) hTriggerAbove4->Fill(1);
		} // Loop over tracks in the event

		mHistNTrks->Fill(ntrk);
		mHistMultiplicity->Fill(mult);
		mHistEMultiplicity->Fill(emult);
		mHistNTrksMC->Fill(gntrk);
		mHistMultiplicityMC->Fill(gmult);
		mHistEMultiplicityMC->Fill(gemult);
		mHistNTrksTRUTHMC->Fill(gtntrk);
		mHistMultiplicityTRUTHMC->Fill(gtmult);		
		mHistEMultiplicityTRUTHMC->Fill(gtemult);				

		NEVENTS++;
		
	}
	cout << "Total number of events: " << NEVENTS << endl;

	if( DoMixing ) Mixing(B0MultMix, B0HadrMix, 10, mHist2PCMix, mHistNEvMix, 100);

	OutFile->mkdir("EvtDetHistos");
	OutFile->cd("EvtDetHistos");		
	WriteHistosEvtDet();
	OutFile->mkdir("TrkHistos");
	OutFile->cd("TrkHistos");		
	WriteHistosTracks();
	OutFile->mkdir("DecayHistos");
	OutFile->cd("DecayHistos");		
	WriteHistosDecay();
	OutFile->mkdir("TwoPCHistos");
	OutFile->cd("TwoPCHistos");		
	WriteHistos2PC();
	OutFile->Close();

}

