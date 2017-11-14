// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {

  


  /// @brief MC validation analysis for single photon events
  class MCGAM : public Analysis {
  public:

    /// Default constructor
    MCGAM()
      : Analysis("MCGAM")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      // General FS
      FinalState fs(-5.0, 5.0);
      declare(fs, "FS");

      // Get leading photon
      LeadingParticlesFinalState photonfs(FinalState(-5.0, 5.0, 10.0*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");


      //test: identified final state
      IdentifiedFinalState photons_test(Cuts::abseta < 5.0 && Cuts::pT > 10*GeV);
      photons_test.acceptId(PID::PHOTON);
      declare(photons_test, "Photon_Test");

      // FS for isolation excludes the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "JetFS");

      _h_photon_pT = bookHisto1D("photon_pT", logspace(50, 10.0, 0.5*(sqrtS()>0.?sqrtS():14000.)));
      _h_photon_pT_lin = bookHisto1D("photon_pT_lin", 50, 0.0, 70.0);
      _h_photon_y = bookHisto1D("photon_y", 50, -5.0, 5.0);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      // Get the photon
      const Particles photons = apply<FinalState>(e, "LeadingPhoton").particles();
      if (photons.size() < 1) vetoEvent;

      const Particles photontest = apply<FinalState>(e, "Photon_Test").particlesByPt();
      /*if (photontest.size()!=photons.size()){
          MSG_INFO("  ################ Unterschied: ###########################  ");
          MSG_INFO(photons);
          MSG_INFO(photontest);
      }*/

      // find four vector of phootn from ME
      FourMomentum photon_me;
      const GenEvent* evt = e.genEvent();
      for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p) {
        if ((*p)->status()!=3) continue;
        if ( (*p)->pdg_id()!=22) continue;
        photon_me = (*p)->momentum();

      }


      // perform dr matching
      double dr(9999999);
      Particle selected_photon;
      foreach (const Particle & test_photon, photontest) {
         if  (deltaR(photon_me,test_photon)< dr){
           selected_photon = test_photon;
           dr = deltaR(photon_me,test_photon);
         }
        }



    /*  if (photons.size() != 1) {
        MSG_INFO("strange: not exactly one leading photon");
        MSG_INFO(photons);
      //  vetoEvent;
      } */
      const FourMomentum photon = photons.front().momentum();

       MSG_INFO("number of photons" << photontest.size() );
       MSG_INFO("hardest photon: " << photon);
       MSG_INFO("matched photon: " << selected_photon);
       MSG_INFO("dr = " << dr << "\n");

      // Get all charged particles
      const FinalState& fs = apply<FinalState>(e, "JetFS");
      if (fs.empty()) {
        vetoEvent;
      }

      // Passed cuts, so get the weight
      const double weight = e.weight();

      _h_photon_pT->fill(photon.pT(),weight);
      _h_photon_pT_lin->fill(photon.pT(),weight);
      _h_photon_y->fill(photon.rapidity(),weight);
    }


    // Finalize
    void finalize() {
      scale(_h_photon_pT, crossSectionPerEvent());
      scale(_h_photon_pT_lin, crossSectionPerEvent());
      scale(_h_photon_y, crossSectionPerEvent());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_photon_pT;
    Histo1DPtr _h_photon_pT_lin;
    Histo1DPtr _h_photon_y;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MCGAM);

}
