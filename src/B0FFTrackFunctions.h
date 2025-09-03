static constexpr double mPi   	   = 0.139570;
static constexpr double mP    	   = 0.938272;
static constexpr double mK    	   = 0.493677;
static constexpr double mE    	   = 0.000511;
static constexpr double mK0smin    = 0.4;
static constexpr double mK0smax    = 0.6;
static constexpr double mRhomin    = 0.5;
static constexpr double mRhomax    = 1.5;
static constexpr double mPhimin    = 0.8;
static constexpr double mPhimax    = 1.2;
static constexpr double mLammin    = 1.0;
static constexpr double mLammax    = 1.3;
static constexpr double mD0mmin    = 1.7;
static constexpr double mD0mmax    = 2.0;
static constexpr double mJPsmin    = 2.7;
static constexpr double mJPsmax    = 3.4;

// Function to compute pseudorapidity from position or momentum
float computeEta(float x, float y, float z) {
    float r = std::sqrt(x*x + y*y + z*z);
    float theta = std::acos(z / r);
    return -std::log(std::tan(theta / 2.0));
}

// Function to compute phi
float computePhi(float x, float y) {
    return std::atan2(y, x);
}

// DeltaR
float deltaR(float eta1, float phi1, float eta2, float phi2) {
    float dphi = std::fabs(phi1 - phi2);
    if (dphi > M_PI) dphi = 2*M_PI - dphi;
    float deta = eta1 - eta2;
    return std::sqrt(deta*deta + dphi*dphi);
}

// DeltaPhi
double DeltaPhi(double phi1, double phi2) {
    double dphi = phi1 - phi2;
	while(dphi < -M_PI/2.0) dphi += 2.0 * M_PI;
    while(dphi > 3.0*M_PI/2.0) dphi -= 2.0 * M_PI;
    return dphi;
}

void FillBasicHistos( vector<TVector3> B0Info, vector<bool> B0Match, vector<bool> B0Mult, TH1* px, TH1* py, TH1* pz, TH1* pt, TH1* theta, THnSparse* vec4, TH2* unmatch, TH2* match, TH2* mult){
	for (size_t ib0 = 0; ib0 < B0Info.size(); ib0++) {
		px->Fill(B0Info[ib0].Px());
		py->Fill(B0Info[ib0].Py());
		pz->Fill(B0Info[ib0].Pz());
		pt->Fill(B0Info[ib0].Perp());
		theta->Fill(1000 * B0Info[ib0].Theta());
		double Vec4Trk[5] = {B0Info[ib0].Eta(), B0Info[ib0].Phi(), B0Info[ib0].Perp(), B0Info[ib0].Perp()*B0Info[ib0].Perp(), (double)B0Info.size()};
		vec4->Fill(Vec4Trk);
		if( !B0Match[ib0] ){ unmatch->Fill(B0Info[ib0].Perp(), (double)B0Info.size()); } else { match->Fill(B0Info[ib0].Perp(), (double)B0Info.size()); }
		if( B0Mult[ib0] ){ mult->Fill(B0Info[ib0].Perp(), (double)B0Info.size()); }
	}
}

void FillBasicMCHistos( vector<TVector3> MCInfo, TH1* px, TH1* py, TH1* pz, TH1* pt, TH1* theta, THnSparse* vec4){
	for (size_t ib0 = 0; ib0 < MCInfo.size(); ib0++) {
		px->Fill(MCInfo[ib0].Px());
		py->Fill(MCInfo[ib0].Py());
		pz->Fill(MCInfo[ib0].Pz());
		pt->Fill(MCInfo[ib0].Perp());
		theta->Fill(1000 * MCInfo[ib0].Theta());
		double Vec4Trk[5] = {MCInfo[ib0].Eta(), MCInfo[ib0].Phi(), MCInfo[ib0].Perp(), MCInfo[ib0].Perp()*MCInfo[ib0].Perp(), (double)MCInfo.size()};
		vec4->Fill(Vec4Trk);
	}
}

void InvMassPhys( vector<TVector3> B0Info, vector<float> B0Charge, vector<int> B0Flavour, TH1* PhotonConversion, TH1* InvMassOSJPsi,  TH1* InvMassSSJPsi, TH1* InvMassOSRho,  TH1* InvMassSSRho, TH1* InvMassOSPhi,  TH1* InvMassSSPhi, TH1* InvMassOSD0,  TH1* InvMassSSD0, TH1* InvMassOSAD0,  TH1* InvMassSSAD0, TH1* InvMassOSK0s,  TH1* InvMassSSK0s, TH1* InvMassOSLam,  TH1* InvMassSSLam, TH1* InvMassOSALam,  TH1* InvMassSSALam, TH2* APPlotK0s, TH1* K0sAsLamMassHyp, TH2* APPlotLam, TH1* LamAsK0sMassHyp){
    size_t nTracks = B0Info.size();
    if (nTracks < 2) return; // no pairs possible, double check
	for( size_t i = 0; i < nTracks; i++ ){ // 2-loops in 1
        const TVector3& track_i = B0Info[i];
        const int& flavor_i = B0Flavour[i];
        const float& charge_i = B0Charge[i];
		for( size_t j = i+1; j < nTracks; j++ ){ // 2-loops in 1
        	const TVector3& track_j = B0Info[j];
	        const int& flavor_j = B0Flavour[j];
    	    const float& charge_j = B0Charge[j];
			TLorentzVector track_i_4v_el;
			TLorentzVector track_j_4v_el;
			track_i_4v_el.SetVectM(track_i, mE); 
			track_j_4v_el.SetVectM(track_j, mE); 
			TLorentzVector decayelel = track_i_4v_el + track_j_4v_el;
			// check photon conversion to ee
			if ( (charge_i + charge_j == 0) && decayelel.M() < 1.0 ) { PhotonConversion->Fill(decayelel.M()); }
			// J/Psi
			if ( decayelel.M() > mJPsmin && decayelel.M() < mJPsmax ){ if ( charge_i + charge_j == 0 ) { InvMassOSJPsi->Fill(decayelel.M()); } else { InvMassSSJPsi->Fill(decayelel.M()); } }

			TLorentzVector track_i_4v_pi;
			TLorentzVector track_j_4v_pi;
			track_i_4v_pi.SetVectM(track_i, mPi); 
			track_j_4v_pi.SetVectM(track_j, mPi); 
			TLorentzVector track_i_4v_ka;
			TLorentzVector track_j_4v_ka;
			track_i_4v_ka.SetVectM(track_i, mK); 
			track_j_4v_ka.SetVectM(track_j, mK); 
			TLorentzVector track_i_4v_pr;
			TLorentzVector track_j_4v_pr;
			track_i_4v_pr.SetVectM(track_i, mP); 
			track_j_4v_pr.SetVectM(track_j, mP); 

			TLorentzVector decaypipi = track_i_4v_pi + track_j_4v_pi;
			TLorentzVector decaypika = track_i_4v_pi + track_j_4v_ka;
			TLorentzVector decaypipr = track_i_4v_pi + track_j_4v_pr;
			TLorentzVector decaykapi = track_i_4v_ka + track_j_4v_pi;
			TLorentzVector decayprpi = track_i_4v_pr + track_j_4v_pi;
			TLorentzVector decaykaka = track_i_4v_ka + track_j_4v_ka;
			
			// Rho meson
			if ( decaypipi.M() > mRhomin && decaypipi.M() < mRhomax ){ if ( charge_i + charge_j == 0 ) { InvMassOSRho->Fill(decaypipi.M()); } else { InvMassSSRho->Fill(decaypipi.M()); } }
			// D0 meson
			if ( decaypika.M() > mD0mmin && decaypika.M() < mD0mmax ){ if ( charge_i + charge_j == 0 ) { InvMassOSD0->Fill(decaypika.M()); } else { InvMassSSD0->Fill(decaypika.M()); } }
			// D0bar meson
			if ( decaykapi.M() > mD0mmin && decaykapi.M() < mD0mmax ){ if ( charge_i + charge_j == 0 ) { InvMassOSAD0->Fill(decaykapi.M()); } else { InvMassSSAD0->Fill(decaykapi.M()); } }
			// Phi meson
			if ( decaykaka.M() > mPhimin && decaykaka.M() < mPhimax ){ if ( charge_i + charge_j == 0 ) { InvMassOSPhi->Fill(decaykaka.M()); } else { InvMassSSPhi->Fill(decaykaka.M()); } }
			// K0s
			if ( decaypipi.M() > mK0smin && decaypipi.M() < mK0smax ){ 
				if ( charge_i + charge_j == 0 ) { 

					float Pxp = 0.0, Pyp = 0.0, Pzp = 0.0;
					float Pxn = 0.0, Pyn = 0.0, Pzn = 0.0;
					if ( charge_i > 0 ){ 
						Pxp = track_i.Px(); Pyp = track_i.Py(); Pzp = track_i.Pz();
						Pxn = track_j.Px(); Pyn = track_j.Py(); Pzn = track_j.Pz();
					} else {
						Pxp = track_j.Px(); Pyp = track_j.Py(); Pzp = track_j.Pz();
						Pxn = track_i.Px(); Pyn = track_i.Py(); Pzn = track_i.Pz();
					}
					// Armenteros-Podolanski
 					float Pp=sqrt(Pxp*Pxp+Pyp*Pyp+Pzp*Pzp);
 					float Pn=sqrt(Pxn*Pxn+Pyn*Pyn+Pzn*Pzn);
					float PN=sqrt((Pxp+Pxn)*(Pxp+Pxn)+(Pyp+Pyn)*(Pyp+Pyn)+(Pzp+Pzn)*(Pzp+Pzn));
					float cosx1=(Pxp*(Pxp+Pxn)+Pyp*(Pyp+Pyn)+Pzp*(Pzp+Pzn))/(Pp*PN);
					float cosx2=(Pxn*(Pxp+Pxn)+Pyn*(Pyp+Pyn)+Pzn*(Pzp+Pzn))/(Pn*PN);
					float QT = Pp*sqrt(1-cosx1*cosx1);
					float Alpha;
					Alpha = (Pp*cosx1-Pn*cosx2)/(Pp*cosx1+Pn*cosx2);
					APPlotK0s->Fill(Alpha,QT);
					// K0s <-> Lam: Mass hypotesis
					TLorentzVector track_i_4v_hyp_pi;
					TLorentzVector track_i_4v_hyp_pr;
					TLorentzVector track_j_4v_hyp_pr;
					TLorentzVector track_j_4v_hyp_pi;
					track_i_4v_hyp_pi.SetVectM(track_i, mPi); 
					track_i_4v_hyp_pr.SetVectM(track_i, mP); 
					track_j_4v_hyp_pi.SetVectM(track_j, mPi); 
					track_j_4v_hyp_pr.SetVectM(track_j, mP); 

					TLorentzVector masshyppipr = track_i_4v_hyp_pi + track_j_4v_hyp_pr;	
					TLorentzVector masshypprpi = track_i_4v_hyp_pr + track_j_4v_hyp_pi;	
					
					K0sAsLamMassHyp->Fill(masshyppipr.M());
					K0sAsLamMassHyp->Fill(masshypprpi.M());					

					InvMassOSK0s->Fill(decaypipi.M());

				} else { InvMassSSK0s->Fill(decaypipi.M()); } 
			}

			// Lambda and Anti-Lambda			
			TLorentzVector decayforL;
			if( track_i.Mag() > track_j.Mag() ){ decayforL = decayprpi; } else { decayforL = decaypipr; }
			if ( decayforL.M() > mLammin && decayforL.M() < mLammax ){ 
				if ( charge_i + charge_j == 0 ) { 
					float Pxp = 0.0, Pyp = 0.0, Pzp = 0.0;
					float Pxn = 0.0, Pyn = 0.0, Pzn = 0.0;
					if ( charge_i > 0 ){ 
						Pxp = track_i.Px(); Pyp = track_i.Py(); Pzp = track_i.Pz();
						Pxn = track_j.Px(); Pyn = track_j.Py(); Pzn = track_j.Pz();
					} else {
						Pxp = track_j.Px(); Pyp = track_j.Py(); Pzp = track_j.Pz();
						Pxn = track_i.Px(); Pyn = track_i.Py(); Pzn = track_i.Pz();
					}
					// Armenteros-Podolanski
 					float Pp=sqrt(Pxp*Pxp+Pyp*Pyp+Pzp*Pzp);
 					float Pn=sqrt(Pxn*Pxn+Pyn*Pyn+Pzn*Pzn);
					float PN=sqrt((Pxp+Pxn)*(Pxp+Pxn)+(Pyp+Pyn)*(Pyp+Pyn)+(Pzp+Pzn)*(Pzp+Pzn));
					float cosx1=(Pxp*(Pxp+Pxn)+Pyp*(Pyp+Pyn)+Pzp*(Pzp+Pzn))/(Pp*PN);
					float cosx2=(Pxn*(Pxp+Pxn)+Pyn*(Pyp+Pyn)+Pzn*(Pzp+Pzn))/(Pn*PN);
					float QT = Pp*sqrt(1-cosx1*cosx1);
					float Alpha;
					Alpha = (Pp*cosx1-Pn*cosx2)/(Pp*cosx1+Pn*cosx2);
					APPlotLam->Fill(Alpha,QT);
					// K0s <-> Lam: Mass hypotesis
					TLorentzVector track_i_4v_hyp_pi;
					TLorentzVector track_j_4v_hyp_pi;
					track_i_4v_hyp_pi.SetVectM(track_i, mPi); 
					track_j_4v_hyp_pi.SetVectM(track_j, mPi); 

					TLorentzVector masshyppipi = track_i_4v_hyp_pi + track_j_4v_hyp_pi;	
					LamAsK0sMassHyp->Fill(masshyppipi.M());

					bool isproton_i = ( ( track_i.Mag() > track_j.Mag() ) && charge_i > 0 ); 
					bool isproton_j = ( ( track_j.Mag() > track_i.Mag() ) && charge_j > 0 ); 
					if ( isproton_i || isproton_j ){ InvMassOSLam->Fill(decayforL.M());	} else { InvMassOSALam->Fill(decayforL.M()); }				

				} else {
					bool isproton_i = ( ( track_i.Mag() > track_j.Mag() ) && charge_i > 0 ); 
					bool isproton_j = ( ( track_j.Mag() > track_i.Mag() ) && charge_j > 0 ); 
					if ( isproton_i || isproton_j ){ InvMassSSLam->Fill(decayforL.M());	} else { InvMassSSALam->Fill(decayforL.M()); }		
				}
				
			}

			// Lambdabar

		}
	}
}

void InvMassMC( vector<TVector3> B0Info, vector<float> B0Charge, vector<int> B0PDGID, TH1* InvMassOSJPsi,  TH1* InvMassSSJPsi, TH1* InvMassOSRho,  TH1* InvMassSSRho, TH1* InvMassOSPhi,  TH1* InvMassSSPhi, TH1* InvMassOSD0,  TH1* InvMassSSD0, TH1* InvMassOSAD0,  TH1* InvMassSSAD0, TH1* InvMassOSK0s,  TH1* InvMassSSK0s, TH1* InvMassOSLam,  TH1* InvMassSSLam, TH1* InvMassOSALam,  TH1* InvMassSSALam){
    size_t nTracks = B0Info.size();
    if (nTracks < 2) return; // no pairs possible, double check
	for( size_t i = 0; i < nTracks; i++ ){ 
        const TVector3& track_i = B0Info[i];
        const float& charge_i = B0Charge[i];
        const int& pdgid_i = B0PDGID[i];
        if( charge_i == 0 ) continue;
		float mass1 = 0.0;
		if ( fabs( pdgid_i ) == 11 )    mass1 = mE;
		if ( fabs( pdgid_i ) == 211 )   mass1 = mPi;
		if ( fabs( pdgid_i ) == 321 )   mass1 = mK;
		if ( fabs( pdgid_i ) == 2212 )  mass1 = mP;

		for( size_t j = i+1; j < nTracks; j++ ){ 
        	const TVector3& track_j = B0Info[j];
    	    const float& charge_j = B0Charge[j];
	        const int& pdgid_j = B0PDGID[j];
	        if( charge_j == 0 ) continue;
			float mass2 = 0.0;
			if ( fabs( pdgid_j ) == 11 )    mass2 = mE;
			if ( fabs( pdgid_j ) == 211 )   mass2 = mPi;
			if ( fabs( pdgid_j ) == 321 )   mass2 = mK;
			if ( fabs( pdgid_j ) == 2212 )  mass2 = mP;
			
			TLorentzVector track_i_4v;
			TLorentzVector track_j_4v;
			track_i_4v.SetVectM(track_i, mass1); 
			track_j_4v.SetVectM(track_j, mass2); 
			TLorentzVector decay = track_i_4v + track_j_4v;			
			// J/Psi			
			if ( fabs(pdgid_i) == 11 && fabs(pdgid_j) == 11 && ( decay.M() > mJPsmin && decay.M() < mJPsmax ) ) { if ( charge_i + charge_j == 0 ) { InvMassOSJPsi->Fill(decay.M()); } else { InvMassSSJPsi->Fill(decay.M()); } } 
			// Rho
			if ( fabs(pdgid_i) == 211 && fabs(pdgid_j) == 211 && ( decay.M() > mRhomin && decay.M() < mRhomax ) ) { if ( charge_i + charge_j == 0 ) { InvMassOSRho->Fill(decay.M()); } else { InvMassSSRho->Fill(decay.M()); } } 
			// Rho
			if ( fabs(pdgid_i) == 321 && fabs(pdgid_j) == 321 && ( decay.M() > mPhimin && decay.M() < mPhimax ) ) { if ( charge_i + charge_j == 0 ) { InvMassOSPhi->Fill(decay.M()); } else { InvMassSSPhi->Fill(decay.M()); } } 
			// D0s
			if ( fabs(pdgid_i) == 211 && fabs(pdgid_j) == 321 && ( decay.M() > mD0mmin && decay.M() < mD0mmax ) ) { 
				if ( charge_i + charge_j == 0 ) { if( charge_i > 0) { InvMassOSD0->Fill(decay.M()); } else { InvMassOSAD0->Fill(decay.M()); } 
				} else { if( charge_i > 0) { InvMassSSD0->Fill(decay.M()); } else { InvMassSSAD0->Fill(decay.M()); } }
			} 
			if ( fabs(pdgid_i) == 321 && fabs(pdgid_j) == 211 && ( decay.M() > mD0mmin && decay.M() < mD0mmax ) ) { 
				if ( charge_i + charge_j == 0 ) { if( charge_j > 0) { InvMassOSD0->Fill(decay.M()); } else { InvMassOSAD0->Fill(decay.M()); } 
				} else { if( charge_j > 0) { InvMassSSD0->Fill(decay.M()); } else { InvMassSSAD0->Fill(decay.M()); } }
			} 
			// Lambdas
			if ( fabs(pdgid_i) == 211 && fabs(pdgid_j) == 2122 && ( decay.M() > mLammin && decay.M() < mLammax ) ) { 
				if ( charge_i + charge_j == 0 ) { if( charge_i < 0) { InvMassOSLam->Fill(decay.M()); } else { InvMassOSALam->Fill(decay.M()); } 
				} else { if( charge_i < 0) { InvMassSSLam->Fill(decay.M()); } else { InvMassSSALam->Fill(decay.M()); } }
			} 
			if ( fabs(pdgid_i) == 2122 && fabs(pdgid_j) == 211 && ( decay.M() > mLammin && decay.M() < mLammax ) ) { 
				if ( charge_i + charge_j == 0 ) { if( charge_j < 0) { InvMassOSLam->Fill(decay.M()); } else { InvMassOSALam->Fill(decay.M()); } 
				} else { if( charge_j < 0) { InvMassSSLam->Fill(decay.M()); } else { InvMassSSALam->Fill(decay.M()); } }
			} 
			// K0s
			if ( fabs(pdgid_i) == 211 && fabs(pdgid_j) == 211 && ( decay.M() > mK0smin && decay.M() < mK0smax ) ) { if ( charge_i + charge_j == 0 ) { InvMassOSK0s->Fill(decay.M()); } else { InvMassSSK0s->Fill(decay.M()); } } 

        }
    }
}

void TwoPC(const std::vector<TVector3>& B0Info, THnSparse* vec4, int multiplicity) {

    size_t nTracks = B0Info.size();
    if (nTracks < 2) return; // no pairs possible, double check

    for (size_t i = 0; i < nTracks; i++) {
        const TVector3& track_i = B0Info[i];
        for (size_t j = 0; j < nTracks; j++) {
            if (i == j) continue; // skip self-pairs
            const TVector3& track_j = B0Info[j];
            double deta = track_i.Eta() - track_j.Eta();
            double dphi = DeltaPhi(track_i.Phi(), track_j.Phi());
            double pT1  = track_i.Perp();
            double pT2  = track_j.Perp();
            double Vec2PC[5] = {deta, dphi, pT1, pT2, static_cast<double>(multiplicity)};
            vec4->Fill(Vec2PC);
        }
    }
}

void Mixing(std::vector<int> multiplicity, std::vector<std::vector<TVector3>> B0Trks, int NEventtoMix, THnSparse* vec4, TH1* NEvMix, int multWindow) {

    int aux_n_evts = (int)multiplicity.size();
    std::cout << "Running Mixing... " << std::endl;
    std::cout << "Total # of events: " << aux_n_evts << std::endl;

    std::mt19937 rng(12345); // fixed seed for reproducibility; replace with std::random_device{}() for true randomness

    for (int nevt_trg = 0; nevt_trg < aux_n_evts; nevt_trg++) {

        if (nevt_trg != 0 && (nevt_trg % (aux_n_evts / 10)) == 0) std::cout << nevt_trg << " out of " << aux_n_evts << std::endl;

        auto &TrkTrgEvtVec = B0Trks[nevt_trg];
        int nMix_nevt_trg = TrkTrgEvtVec.size();
        if (nMix_nevt_trg == 0) continue;

        int mult_trg = multiplicity[nevt_trg];

        // Build list of candidate events with similar multiplicity
        std::vector<int> assocCandidates;
        for (int iev = 0; iev < aux_n_evts; iev++) {
            if (iev == nevt_trg) continue;
				if (std::abs(multiplicity[iev] - mult_trg) <= multWindow) assocCandidates.push_back(iev);
        }
        // If no candidates, skip
        if (assocCandidates.empty()) continue;
        // Shuffle candidates
        std::shuffle(assocCandidates.begin(), assocCandidates.end(), rng);
        // Take first NEventtoMix candidates
        int nAssocToUse = std::min(NEventtoMix, (int)assocCandidates.size());
        NEvMix->Fill(nAssocToUse);
        for (int iass = 0; iass < nAssocToUse; iass++) {
            int nevt_assoc = assocCandidates[iass];
            auto &TrkAssEvtVec = B0Trks[nevt_assoc];
            int nMix_nevt_ass = (int)TrkAssEvtVec.size();
            if (nMix_nevt_ass == 0) continue;
            for (int imix = 0; imix < nMix_nevt_trg; imix++) {
                for (int iimix = 0; iimix < nMix_nevt_ass; iimix++) {
                    double deta = TrkTrgEvtVec[imix].Eta() - TrkAssEvtVec[iimix].Eta();
                    double dphi = DeltaPhi(TrkTrgEvtVec[imix].Phi(), TrkAssEvtVec[iimix].Phi());
                    double pT1 = TrkTrgEvtVec[imix].Perp();
                    double pT2 = TrkAssEvtVec[iimix].Perp();
                    double Vec2PC[5] = {deta, dphi, pT1, pT2, (double)mult_trg};
                    vec4->Fill(Vec2PC);
                }
            }
        }
    }
}

