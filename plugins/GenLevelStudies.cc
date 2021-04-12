// -*- C++ -*-
//
// Package:    GenTutorial/GenLevelStudies
// Class:      GenLevelStudies
//
/**\class GenLevelStudies GenLevelStudies.cc GenTutorial/GenLevelStudies/plugins/GenLevelStudies.cc
Description: [one line class summary]
Implementation:
[Notes on implementation]
*/
//
// Original Author:  Sandro Fonseca De Souza
//         Created:  Wed, 05 Feb 2020 09:25:53 GMT
//
//
// CMSSW folder: /afs/cern.ch/work/m/mabarros/CMSSW_10_6_4/src/GenTutorial/GenLevelStudies/plugins
// scp command: scp DemoAnalyzer.cc mabarros@lxplus7.cern.ch:/afs/cern.ch/work/m/mabarros/public/CMSSW_10_6_4/src/GenTutorial/GenLevelStudies/plugins

// Interesting directory: /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_10_6_4/src/SimGeneral/HepPDTRecord/interface

// system include files



#include <memory>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include "TH1.h"
#include "TH3.h"
//#include "Pythia8/Pythia.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::GenParticleCollection;
using namespace std;
using namespace reco;
using namespace edm;
//using namespace Pythia8;

class GenLevelStudies : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit GenLevelStudies(const edm::ParameterSet&);
		~GenLevelStudies();

	private:
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void beginJob() override ;
        virtual void endJob() override ;
        void initialize();
		// ----------member data ---------------------------
		edm::EDGetTokenT<GenParticleCollection> genParticlesToken_;
        
		// members for ParticleListDrawer
		std::string getParticleName(int id) const;

		edm::InputTag src_;
		edm::EDGetTokenT<reco::CandidateView> srcToken_;
		edm::ESHandle<ParticleDataTable> pdt_;

		int maxEventsToPrint_; // Must be signed, because -1 is used for no limit.
		unsigned int nEventAnalyzed_;
		bool printOnlyHardInteraction_;
		bool printVertex_;
		bool printFlags_;
		bool useMessageLogger_;
	
		
		TTree *mc;
        TH1D *histjp, *histDstar, *histpTMuonN, *histPhiMuonN, *histetaMuonN;
		TH1D *jpsivertexX, *dstarvertexX, *jpsivertexY, *dstarvertexY, *jpsivertexZ, *dstarvertexZ, *hist_vertexZ_diff;
		TH3D *histAsso;
		vector<double> genpt;
        vector<double> geneta;
		int runNumber=0; int eventNumber=0;   

		// Variables for vertices;
		vector<double> jpsiVertexX;
		vector<double> jpsiVertexY;
		vector<double> jpsiVertexZ;

		vector<double> DstarVertexX;
		vector<double> DstarVertexY;
		vector<double> DstarVertexZ;

		vector<double> jpsiAssoToDstarVertexX;
		vector<double> jpsiAssoToDstarVertexY;
		vector<double> jpsiAssoToDstarVertexZ;

		vector<double> DstarAssoToJpsiVertexX;
		vector<double> DstarAssoToJpsiVertexY;
		vector<double> DstarAssoToJpsiVertexZ;

		vector<double> jpsimass;
		vector<double> Dstarmass;

		vector<double> jpsiAssoToDstarmass;
		vector<double> DstarAssoToJpsimass;

		          
};

////////////////////////////////////////////////////////////////////////////
GenLevelStudies::GenLevelStudies(const edm::ParameterSet& iConfig)
	:genParticlesToken_(consumes<GenParticleCollection>(edm::InputTag{"genParticles"})),
	 src_(iConfig.getParameter<InputTag>("src")),
	 srcToken_(consumes<reco::CandidateView>(src_)),
	 maxEventsToPrint_(iConfig.getUntrackedParameter<int>("maxEventsToPrint",1)),
	 nEventAnalyzed_(0),
	 printOnlyHardInteraction_(iConfig.getUntrackedParameter<bool>("PrintOnlHardInteraction", false)),
	 printVertex_(iConfig.getUntrackedParameter<bool>("printVertex", false)),
	 printFlags_(iConfig.getUntrackedParameter<bool>("printFlags", false)),
	 useMessageLogger_(iConfig.getUntrackedParameter<bool>("useMessageLogger", false)) 

{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs;
	mc = fs->make<TTree>("mc","mc");
    TFileDirectory Mreconst = fs->mkdir("MrecoHist");
	TFileDirectory MuonDir = fs->mkdir("MuonsHist");
	TFileDirectory jpsiWithDstar = fs->mkdir("AssosParticles");

	histjp = Mreconst.make<TH1D>("histjp", "histjp", 200, 3.0964, 3.0977);
	histDstar = Mreconst.make<TH1D>("histDstar","histDstar", 200, 1.8646, 1.8656);
    
	jpsivertexX = Mreconst.make<TH1D>("jpsix","jpsix", 100, -0.1, 0.1);
	jpsivertexY = Mreconst.make<TH1D>("jpsiy","jpsiy", 100, -0.1, 0.1);
	jpsivertexZ = Mreconst.make<TH1D>("jpsiz","jpsiz", 100, -5, 5);

	dstarvertexX = Mreconst.make<TH1D>("dstarx","dstarx", 100, -0.1, 0.1);
	dstarvertexY = Mreconst.make<TH1D>("dstary","dstary", 100, -0.1, 0.1);
	dstarvertexZ = Mreconst.make<TH1D>("dstarz","dstarz", 100, -5, 5);
	hist_vertexZ_diff = Mreconst.make<TH1D>("deltaVertexZ","deltaVertexZ", 100, -5, 5);

	histAsso = jpsiWithDstar.make<TH3D>("histAssosParticle", "histAssosParticle", 200, 3.0964, 3.0977, 200, 1.8646, 1.8656143, 143, 0, 143);

	
	histpTMuonN = MuonDir.make<TH1D>("pTmuonN", "pTmuonN", 100, 0., 50.);  
	histPhiMuonN = MuonDir.make<TH1D>("phimuonN","phimuonN", 100, -2*M_PI, 2*M_PI);
	histetaMuonN = MuonDir.make<TH1D>("pseudorapMuonN","pseudorapMuonN", 100, -2.5, 2.5);
	
}

std::string GenLevelStudies::getParticleName(int id) const{
	const ParticleData * pd = pdt_->particle(id);
	if (!pd) { // !0 =1 (true)
		std::ostringstream ss;
		ss << "p" << id;
		return ss.str();
	}
	else{
		return pd->name();
	}
} 

GenLevelStudies::~GenLevelStudies()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}
//
// member functions
//

// ------------ method called  ------------
void GenLevelStudies::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;
    //to clear vectors
    initialize();
    //run, event, lumi section
	runNumber= iEvent.id().run();
	eventNumber= iEvent.id().event();

	// Initialization for ParticleDrawerList;
	Handle <reco::CandidateView> particles;
	iEvent.getByToken (srcToken_, particles);
	iSetup.getData(pdt_);

	Handle<GenParticleCollection> genParticles;
	iEvent.getByToken(genParticlesToken_, genParticles);
	cout << "Running " << endl;
 // ------------------------------ Particle name ----------------------------------------------- 
    /* if(maxEventsToPrint_ < 0 || nEventAnalyzed_ < static_cast <unsigned int>(maxEventsToPrint_)){
		ostringstream out;
		char buf[256];

		out << endl
		    << "[ParticleListDrawer] analysing particle collection " << src_.label() << endl;
	
	    snprintf(buf, 256, " idx  |   ID  -   Name |Stat|  Mo1  Mo2  Da1  Da2 |nMo nDa|    pt       eta     phi   |     px         py         pz        m     |");
	    out << buf;
	    if (printVertex_){
			snprintf(buf, 256, "        vx       vy        vz     |");
             out << buf;
        }
        out << endl;
	
        int idx  = -1;
        int iMo1 = -1;
        int iMo2 = -1;
        int iDa1 = -1;
        int iDa2 = -1; 

		vector<const reco::Candidate *> cands;
		vector<const Candidate *>::const_iterator found = cands.begin();
        
		for(CandidateView::const_iterator p = particles->begin(); p != particles->end(); ++p){
			cands.push_back(&*p);
			if (printOnlyHardInteraction_ && p->status() != 3) continue;

			// Particle Name
			int id = p->pdgId();
			string ParticleName = getParticleName(id);

            // Particle Index
			int idx  = -1;
			idx = p - particles->begin();

			// Particles Mothers and Daughters
			
            int iMo1 = -1;
            int iMo2 = -1;
            int iDa1 = -1;
            int iDa2 = -1;

			int nMo = p->numberOfMothers();
			int nDa = p->numberOfDaughters();

			found = find(cands.begin(), cands.end(), p->mother(0));
			if (found != cands.end()) iMo1 = found - cands.begin();

			found = find(cands.begin(), cands.end(), p->mother(nMo-1));
			if (found != cands.end()) iMo2 = found - cands.begin();

			found = find(cands.begin(), cands.end(), p->daughter(0));
			if (found != cands.end()) iDa1 = found - cands.begin();

			found = find(cands.begin(), cands.end(), p->daughter(nDa-1));
			if (found != cands.end()) iDa2 = found - cands.begin();

            char buf[256];
			snprintf(buf, 256,
	         " %4d | %5d - %10s | %2d | %4d %4d %4d %4d | %2d %2d | %7.3f %10.3f %6.3f | %10.3f %10.3f %10.3f %8.3f |",
                       idx,
                       p->pdgId(),
                       ParticleName.c_str(),
                       p->status(),
                       iMo1,iMo2,iDa1,iDa2,nMo,nDa,
                       p->pt(),
                       p->eta(),
                       p->phi(),
                       p->px(),
                       p->py(),
                       p->pz(),
                       p->mass()
                      );
            out << buf;


			if (printVertex_) {
				snprintf(buf, 256, " %10.3f %10.3f %10.3f |",
				         p->vertex().x(),
						 p->vertex().y(),
						 p->vertex().z());
			    out << buf;
			}
			nEventAnalyzed_++;

			if (useMessageLogger_)
			    LogVerbatim("ParticleListDrawer") << out.str();
			else
			   cout << out.str();


		}


	} */
    int cjp = 0;
	int cdstar = 0;
	for(const auto& genParticles : iEvent.get(genParticlesToken_) ){
        
		//cout << " PDG_ID: " << genParticles.pdgId() << " status: "<< genParticles.status() << " pT: "  << genParticles.pt() << "|eta|: " << abs(genParticles.eta()) << " phi: "<< genParticles.phi() << endl;
        		 
        //double pt = genParticles.pt();
        //double eta = abs(genParticles.eta()); 
		//Pythia8::Vec4 muonP, muonN, sum; OBS: I founded a problem by substituting
		// the genParticle.p() into Vec4. Remember: p = 14.0, vec4 = 14.0, 14.0 (maybe p is just three momentum)

		// If the particle is J/Psi and if J/psi -> u+u-
		
		if (genParticles.pdgId() == 443) {
			jpsivertexX->Fill(genParticles.vertex().x());
			jpsivertexY->Fill(genParticles.vertex().y());
			jpsivertexZ->Fill(genParticles.vertex().z());
			cjp += 1;
		}
			
			
					
			
		// If the particle is D*+
		
		if (genParticles.pdgId() == 413 && genParticles.numberOfDaughters()==2){
			if ((genParticles.daughter(0)->pdgId() == 421) || (genParticles.daughter(1)->pdgId() == 421)){
				//std::cout << genParticles.daughter(0)->daughter(0)->pdgId() << std::endl;
				if ((genParticles.daughter(0)->daughter(0)->pdgId() == -321) || (genParticles.daughter(0)->daughter(1)->pdgId() == -321) ){
					dstarvertexX->Fill(genParticles.vertex().x());
					dstarvertexY->Fill(genParticles.vertex().y());
					dstarvertexZ->Fill(genParticles.vertex().z());
					cdstar += 1;

				}
				
			}
	
		// If the particle is D*-
		}
		else if (genParticles.pdgId() == -413 && genParticles.numberOfDaughters()==2){
			if ((genParticles.daughter(0)->pdgId() == -421) || (genParticles.daughter(1)->pdgId() == -421)){
				if ((genParticles.daughter(0)->daughter(0)->pdgId() == 321) || (genParticles.daughter(0)->daughter(1)->pdgId() == 321) ){
					dstarvertexX->Fill(genParticles.vertex().x());
					dstarvertexY->Fill(genParticles.vertex().y());
					dstarvertexZ->Fill(genParticles.vertex().z());
					cdstar += 1;

				}
				
			}
	
	
	
		}	
		
	}
	std::cout << "num jpsi: " << cjp << std::endl;
	std::cout << "num dstar: " << cdstar << std::endl;
	
	hist_vertexZ_diff = new TH1D(*jpsivertexZ); // TH1F or TH1D, same as h_hit_hitsperchannel1
	hist_vertexZ_diff->SetNameTitle("hist_vertexZ_diff", "hist_vertexZ_diff;#Delta Z [cm]; Events");
	if (!(hist_vertexZ_diff->GetSumw2N() > 0)) hist_vertexZ_diff->Sumw2(kTRUE); // ensure proper error propagation
	hist_vertexZ_diff->Add(dstarvertexZ, -1.0);


	mc->Fill();
}

void GenLevelStudies::initialize( )
{
        runNumber=0; eventNumber=0;
        genpt.clear();
        geneta.clear();

}
///////
//++++++++++++++++++
void GenLevelStudies::endJob(){

    cout <<"######################################################################"<<endl;
    cout << "Number of Events: " << eventNumber << " Run Number: " << runNumber << endl;
    
}


/////////////////////
     void
GenLevelStudies::beginJob()
{
    mc->Branch("runNumber",&runNumber,"runNumber/I");
	mc->Branch("eventNumber",&eventNumber,"eventNumber/I");
    mc->Branch("genpt",&genpt);
    mc->Branch("geneta",&geneta);
}



//////////////////////////////////////////
//define this as a plug-in
DEFINE_FWK_MODULE(GenLevelStudies);