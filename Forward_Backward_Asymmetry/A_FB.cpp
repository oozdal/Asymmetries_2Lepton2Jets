#include "SampleAnalyzer/User/Analyzer/A_FB.h"
#include <cmath>
#include <iostream>
#include "SampleAnalyzer/Commons/Service/PDGService.h"

using namespace MA5;
using namespace std;

MAbool A_FB::Initialize(const MA5::Configuration& cfg,
                      const std::map<std::string,std::string>& parameters)
{
  // Initializing PhysicsService for MC
  PHYSICS->mcConfig().Reset();

  // definition of the multiparticle "hadronic"
  PHYSICS->mcConfig().AddHadronicId(-5);
  PHYSICS->mcConfig().AddHadronicId(-4);
  PHYSICS->mcConfig().AddHadronicId(-3);
  PHYSICS->mcConfig().AddHadronicId(-2);
  PHYSICS->mcConfig().AddHadronicId(-1);
  PHYSICS->mcConfig().AddHadronicId(1);
  PHYSICS->mcConfig().AddHadronicId(2);
  PHYSICS->mcConfig().AddHadronicId(3);
  PHYSICS->mcConfig().AddHadronicId(4);
  PHYSICS->mcConfig().AddHadronicId(5);
  PHYSICS->mcConfig().AddHadronicId(21);

  // definition of the multiparticle "invisible"
  PHYSICS->mcConfig().AddInvisibleId(12);

  // ===== Signal region ===== //
  Manager()->AddRegionSelection("Rap_No");
  Manager()->AddRegionSelection("Rap_00");
  Manager()->AddRegionSelection("Rap_01");
  Manager()->AddRegionSelection("Rap_02");  
  Manager()->AddRegionSelection("Rap_04");
  Manager()->AddRegionSelection("Rap_06");

  // ===== Selections ===== //
  Manager()->AddCut("$N_l=2$");
  Manager()->AddCut("$N_j=2$");
//  Manager()->AddCut("OSSF");
  Manager()->AddCut("$|y_{ll}|>0.0$", "Rap_00");
  Manager()->AddCut("$|y_{ll}|>0.1$", "Rap_01");  
  Manager()->AddCut("$|y_{ll}|>0.2$", "Rap_02");  
  Manager()->AddCut("$|y_{ll}|>0.4$", "Rap_04");
  Manager()->AddCut("$|y_{ll}|>0.6$", "Rap_06");

  // ===== Histograms ===== //

//  Manager()->AddHisto("AFB_pos_Rap_No",   30,50.,200., "Rap_No");
//  Manager()->AddHisto("AFB_neg_Rap_No",   30,50.,200., "Rap_No");

//  Manager()->AddHisto("AFB_pos_Rap_00",   15,50.,200., "Rap_00");
//  Manager()->AddHisto("AFB_neg_Rap_00",   15,50.,200., "Rap_00");

//  Manager()->AddHisto("AFB_pos_Rap_01",  15,50.,200., "Rap_01");
//  Manager()->AddHisto("AFB_neg_Rap_01",  15,50.,200., "Rap_01");

//  Manager()->AddHisto("AFB_pos_Rap_02",  15,50.,200., "Rap_02");
//  Manager()->AddHisto("AFB_neg_Rap_02",  15,50.,200., "Rap_02");

//  Manager()->AddHisto("AFB_pos_Rap_04",  15,50.,200., "Rap_04");
//  Manager()->AddHisto("AFB_neg_Rap_04",  15,50.,200., "Rap_04");

//  Manager()->AddHisto("AFB_pos_Rap_06",  15,50.,200., "Rap_06");
//  Manager()->AddHisto("AFB_neg_Rap_06",  15,50.,200., "Rap_06");


  Manager()->AddHisto("AFB_pos_Rap_No",  13,2400.,3600., "Rap_No");
  Manager()->AddHisto("AFB_neg_Rap_No",  13,2400.,3600., "Rap_No");  

  Manager()->AddHisto("AFB_pos_Rap_00",  13,2400.,3600., "Rap_00");
  Manager()->AddHisto("AFB_neg_Rap_00",  13,2400.,3600., "Rap_00");

  Manager()->AddHisto("AFB_pos_Rap_01",  13,2400.,3600., "Rap_01");
  Manager()->AddHisto("AFB_neg_Rap_01",  13,2400.,3600., "Rap_01");  

  Manager()->AddHisto("AFB_pos_Rap_02",  13,2400.,3600., "Rap_02");
  Manager()->AddHisto("AFB_neg_Rap_02",  13,2400.,3600., "Rap_02");

  Manager()->AddHisto("AFB_pos_Rap_04",  6,2400.,3600., "Rap_04");
  Manager()->AddHisto("AFB_neg_Rap_04",  6,2400.,3600., "Rap_04");

  Manager()->AddHisto("AFB_pos_Rap_06",  6,2400.,3600., "Rap_06");
  Manager()->AddHisto("AFB_neg_Rap_06",  6,2400.,3600., "Rap_06");


  Manager()->AddHisto("Eta_lp_Rap_No",  16,-4.,4., "Rap_No");
  Manager()->AddHisto("Eta_lm_Rap_No",  16,-4.,4., "Rap_No");
  
  Manager()->AddHisto("Eta_lp_Rap_00",  16,-4.,4., "Rap_00");
  Manager()->AddHisto("Eta_lm_Rap_00",  16,-4.,4., "Rap_00");

  Manager()->AddHisto("Eta_lp_Rap_01",  16,-4.,4., "Rap_01");
  Manager()->AddHisto("Eta_lm_Rap_01",  16,-4.,4., "Rap_01");

  Manager()->AddHisto("Eta_lp_Rap_02",  16,-4.,4., "Rap_02");
  Manager()->AddHisto("Eta_lm_Rap_02",  16,-4.,4., "Rap_02");  

  Manager()->AddHisto("Eta_lp_Rap_04",  16,-4.,4., "Rap_04");
  Manager()->AddHisto("Eta_lm_Rap_04",  16,-4.,4., "Rap_04");  

  Manager()->AddHisto("Eta_lp_Rap_06",  16,-4.,4., "Rap_06");
  Manager()->AddHisto("Eta_lm_Rap_06",  16,-4.,4., "Rap_06");  

//  Manager()->AddHisto("Mll_Rap_No",   15,50.,200., "Rap_No");
//  Manager()->AddHisto("Mll_Rap_00",   15,50.,200., "Rap_00");
//  Manager()->AddHisto("Mll_Rap_01",   15,50.,200., "Rap_01");
//  Manager()->AddHisto("Mll_Rap_02",   15,50.,200., "Rap_02");
//  Manager()->AddHisto("Mll_Rap_04",   15,50.,200., "Rap_04");
//  Manager()->AddHisto("Mll_Rap_06",   15,50.,200., "Rap_06");

//  Manager()->AddHisto("Mll_Rap_No",   13,2400.,3600., "Rap_No");  
//  Manager()->AddHisto("Mll_Rap_00",   13,2400.,3600., "Rap_00");
//  Manager()->AddHisto("Mll_Rap_01",   13,2400.,3600., "Rap_01");  
//  Manager()->AddHisto("Mll_Rap_02",   13,2400.,3600., "Rap_02");  
//  Manager()->AddHisto("Mll_Rap_04",    6,2400.,3600., "Rap_04");
//  Manager()->AddHisto("Mll_Rap_06",    6,2400.,3600., "Rap_06");


  // No problem during initialization

  myfile.open("A_FB.csv");//, std::ios_base::app);
  myfile << "Mll" << ","
	 << "Mll_boosted" << ","
	 << "Mjj" << ","
	 << "Mlljj" << ","
	 << "MET" << ","
         << "MT"  << ","
	 << "MHT"  << ","
         << "Angle_CM_true" << ","
         << "Angle_CM_rec" << ","
         << "Eta_Boosted_lp" << ","
         << "Eta_Boosted_lm" << "," 
	 << "Eta_lp" << ","
	 << "Eta_lm" << ","
         << "Quark Charge" << ","
	 << "Lepton Charge" << ","
         << "yll_lab" << ","
         << "yll_CM"  << ","
	 << "hasOS" << ","
         << "ml1jj" << ","
         << "ml2jj" << ","
         << "mlpjj" << ","
         << "mlmjj" << ","
         << "angle_CM_true_wrtl1" << ","
	 << "angle_CM_rec_wrtl1" << ","
	 << "EvWeight" << "\n";
 
  return true;
}

MAbool A_FB::Execute(SampleFormat& sample, const EventFormat& event)
{
  MAfloat32 EvWeight = 1.0;
  if (weighted_events_ && event.mc()!=0) EvWeight *= event.mc()->weight();

  //EvWeight *= sample.mc()->xsection(); //* 1000. * lumi;
  if (sample.mc()!=0) sample.mc()->addWeightedEvents(EvWeight);
  Manager()->InitializeForNewEvent(EvWeight);


  // Clearing particle containers
  leptons.clear();
  jets.clear();
  initialquarks.clear();
//  lepton1.clear();
//  lepton2.clear();

  // Filling particle containers
  for (MAuint32 i=0;i<event.mc()->particles().size();i++)
  {	  
    if (is_lepton(&(event.mc()->particles()[i])))
      leptons.push_back(&(event.mc()->particles()[i]));
    if (is_jet(&(event.mc()->particles()[i])))
      jets.push_back(&(event.mc()->particles()[i]));
    if (is_initialquark(&(event.mc()->particles()[i])))
      initialquarks.push_back(&(event.mc()->particles()[i])); 
  }

  // Sorting particles
  // Sorting particle collection according to PTordering 
  SORTER->sort(leptons,PTordering);
  SORTER->sort(jets,PTordering);
  // for getting 2th particle
//  lepton2=SORTER->rankFilter(leptons,2,PTordering);
  // for getting 1th particle
//  lepton1=SORTER->rankFilter(leptons,1,PTordering);


  // DiLepton reconstruction
  MAuint32 hasOS   = 0;
  MAuint32 hasOSSF = 0;
  for (MAuint32 lid = 0; lid < leptons.size(); lid++)
  {
    for (MAuint32 sid = 0; (sid < leptons.size()) && (sid != lid); sid++)
    {
        int sumflvr = leptons[lid]->pdgid() + leptons[sid]->pdgid();
        // OS: 0, 2,  OSSF : 0
        // SSOF : 24, SSSF : 22, 26
        // SSOF : 24, OSOF : 2
        if (std::abs(sumflvr) == 0) {hasOSSF++; hasOS++;}
    }
  }

  if(!Manager()->ApplyCut(leptons.size()==2, "$N_l=2$")) return true;
  if(!Manager()->ApplyCut(jets.size()==2, "$N_j=2$")) return true;
//  if(!Manager()->ApplyCut(hasOSSF>0, "OSSF")) return true;
  
  // Useful quantities
  double mll    = (leptons[0]->momentum() + leptons[1]->momentum()).M();
  double mjj    = (jets[0]->momentum() + jets[1]->momentum()).M();
  double mlljj  = (leptons[0]->momentum() + leptons[1]->momentum() + jets[0]->momentum() + jets[1]->momentum()).M(); 

  double ml1jj  = (leptons[0]->momentum() + jets[0]->momentum() + jets[1]->momentum()).M();
  double ml2jj  = (leptons[1]->momentum() + jets[0]->momentum() + jets[1]->momentum()).M();

  // MET & MT
  double MET_double = PHYSICS->Transverse->EventMET(event.mc());
  MALorentzVector MET = event.mc()->MET().momentum();
  double MtMiss = leptons[0]->mt_met(MET);

  // HT
  // MAdouble64 THT = event.mc()->THT();
  // double THT = PHYSICS->Transverse->EventTHT(event.mc());
  double MHT = PHYSICS->Transverse->EventMHT(event.mc());

  double eta_l1 = leptons[0]->eta();
  double eta_l2 = leptons[1]->eta();
//  double angle_lab  = 0;
  MAdouble64 angle_CM_true = 0;
  MAdouble64 angle_CM_rec  = 0;  
  MAdouble64 cos_angle_CM_star = 0;

  MAdouble64 angle_CM_true_wrtl1 = 0;
  MAdouble64 angle_CM_rec_wrtl1  = 0;

  double etalp  = 0;
  double etalm  = 0;

  MALorentzVector lplepton;
  MALorentzVector lmlepton;
  MALorentzVector l1lepton;
  MALorentzVector l2lepton;

  for (MAuint32 i=0; i<leptons.size(); i++)
  {
  	if (leptons[i]->pdgid()<0)
	{
		lplepton.SetPxPyPzE(leptons[i]->px(),leptons[i]->py(),leptons[i]->pz(),leptons[i]->e());    
	}	
        if (leptons[i]->pdgid()>0)
        {
                lmlepton.SetPxPyPzE(leptons[i]->px(),leptons[i]->py(),leptons[i]->pz(),leptons[i]->e());
        } 	
  }

  l1lepton.SetPxPyPzE(leptons[0]->px(),leptons[0]->py(),leptons[0]->pz(),leptons[0]->e());
  l2lepton.SetPxPyPzE(leptons[1]->px(),leptons[1]->py(),leptons[1]->pz(),leptons[1]->e());

 MALorentzVector lplm_true = lplepton + lmlepton;
 MABoost recoil_true;
 recoil_true.setBoostVector(-lplm_true.Px()/lplm_true.E(), -lplm_true.Py()/lplm_true.E(), -lplm_true.Pz()/lplm_true.E());
 MALorentzVector boosted_lp = lplepton;
 MALorentzVector boosted_lm = lmlepton; 
// recoil_true.boost(boosted_lp);
 recoil_true.boost(boosted_lm);


 MALorentzVector l1l2_true = l1lepton + l2lepton;
 MABoost recoil_l1l2;
 recoil_l1l2.setBoostVector(-l1l2_true.Px()/l1l2_true.E(), -l1l2_true.Py()/l1l2_true.E(), -l1l2_true.Pz()/l1l2_true.E());
 MALorentzVector boosted_l1 = l1lepton;
 MALorentzVector boosted_l2 = l2lepton;
 recoil_l1l2.boost(boosted_l1);
// recoil_true.boost(boosted_l2);

 double mll_boosted = (boosted_lp + boosted_lm).M();
 double eta_boosted_lp = 0;
 double eta_boosted_lm = 0;
 double quarkcharge  = 0;
 double leptoncharge = 0;

 double mlpjj = 0;
 double mlmjj = 0;

 MALorentzVector pos_quarks;
 MALorentzVector neg_quarks;

 MAdouble64 yll_lab    = (leptons[0]->momentum() + leptons[1]->momentum()).Rapidity();
 MAdouble64 yll_CM     = (boosted_lp + boosted_lm).Rapidity();

 signed int signofyll_lab = yll_lab/std::fabs(yll_lab);
 signed int signofyll_CM  = yll_CM/std::fabs(yll_CM);

 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.0), "$|y_{ll}|>0.0$")) return true; 
 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.1), "$|y_{ll}|>0.1$")) return true;
 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.2), "$|y_{ll}|>0.2$")) return true; 
 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.4), "$|y_{ll}|>0.4$")) return true;
 if(!Manager()->ApplyCut(std::fabs(yll_lab>0.6), "$|y_{ll}|>0.6$")) return true;

 for (MAuint32 i=0; i<initialquarks.size(); i++)
 {
	quarkcharge = (PDG->GetCharge(initialquarks[i]->pdgid()))/3.; 
	if ( initialquarks[i]->pdgid() > 0 )
	{	
		pos_quarks.SetPxPyPzE(initialquarks[i]->px(),initialquarks[i]->py(),initialquarks[i]->pz(),initialquarks[i]->e());
	}
        if ( initialquarks[i]->pdgid() < 0 )
	{
                neg_quarks.SetPxPyPzE(initialquarks[i]->px(),initialquarks[i]->py(),initialquarks[i]->pz(),initialquarks[i]->e());
       	}	
 }

 angle_CM_true_wrtl1  = boosted_l1.Angle(pos_quarks);
 angle_CM_rec_wrtl1   = boosted_l1.Angle(boosted_l1 + boosted_l2);

 for (MAuint32 i=0; i<leptons.size(); i++)
 {      
	leptoncharge =  (PDG->GetCharge(leptons[i]->pdgid()))/3.; 
        if ( leptoncharge > 0. )
        {
//      angle_lab = leptons[i]->angle(pos_quarks);
//	angle_CM_true  = boosted_lp.Angle(pos_quarks);
//	angle_CM_rec   = boosted_lp.Angle(boosted_lp + boosted_lm);
//	cos_angle_CM_star = cos(angle_CM_rec)*( yll_lab/std::fabs(yll_lab) );
        etalp = leptons[i]->eta();
//      eta_boosted_lp = boosted_lp.Eta();	

        mlpjj  = (leptons[i]->momentum() + jets[0]->momentum() + jets[1]->momentum()).M();

        Manager()->FillHisto("Eta_lp_Rap_No", etalp);
        Manager()->FillHisto("Eta_lp_Rap_00", etalp);
        Manager()->FillHisto("Eta_lp_Rap_01", etalp);	
        Manager()->FillHisto("Eta_lp_Rap_02", etalp);
        Manager()->FillHisto("Eta_lp_Rap_04", etalp);
        Manager()->FillHisto("Eta_lp_Rap_06", etalp);

//	Manager()->FillHisto("Eta_lp_Rap_No", eta_boosted_lp);
//      Manager()->FillHisto("Eta_lp_Rap_00", eta_boosted_lp);
//      Manager()->FillHisto("Eta_lp_Rap_02", eta_boosted_lp);	
//      Manager()->FillHisto("Eta_lp_Rap_04", eta_boosted_lp);
//      Manager()->FillHisto("Eta_lp_Rap_06", eta_boosted_lp);
	}

	if ( leptoncharge < 0. )
        {
//      angle_lab = leptons[i]->angle(pos_quarks);
        angle_CM_true  = boosted_lm.Angle(pos_quarks);
        angle_CM_rec   = boosted_lm.Angle(boosted_lp + boosted_lm);
//	angle_lab_rec   = boosted_lm.Angle(lplepton + lmlepton);
	cos_angle_CM_star = cos(angle_CM_rec)*( yll_CM/std::fabs(yll_CM) );
        etalm = leptons[i]->eta();
        eta_boosted_lm = boosted_lm.Eta();
 
        mlmjj  = (leptons[i]->momentum() + jets[0]->momentum() + jets[1]->momentum()).M();

//	double test_eta = -log(tan(angle_CM_true/2.)); 
//      cout << angle_CM_true << " " << test_eta  << " " << etalm << " " << eta_boosted_lm <<"\n";

//      double test_eta = -log(tan(angle_CM_rec/2.));
//      cout << cos(angle_CM_rec) << " " << test_eta  << " " << etalm << " " << eta_boosted_lm <<"\n";

//      double test_eta = -log(tan(cos_angle_CM_star/2.));
//      cout << cos(angle_CM_rec) << " " << test_eta  << " " << etalm << " " << eta_boosted_lm <<"\n";

        Manager()->FillHisto("Eta_lm_Rap_No", etalm);
        Manager()->FillHisto("Eta_lm_Rap_00", etalm);
        Manager()->FillHisto("Eta_lm_Rap_01", etalm);	
        Manager()->FillHisto("Eta_lm_Rap_02", etalm);
        Manager()->FillHisto("Eta_lm_Rap_04", etalm);
        Manager()->FillHisto("Eta_lm_Rap_06", etalm);

//      Manager()->FillHisto("Eta_lm_Rap_No", eta_boosted_lm);
//      Manager()->FillHisto("Eta_lm_Rap_00", eta_boosted_lm);
//      Manager()->FillHisto("Eta_lm_Rap_02", eta_boosted_lm);	
//      Manager()->FillHisto("Eta_lm_Rap_04", eta_boosted_lm);
//      Manager()->FillHisto("Eta_lm_Rap_06", eta_boosted_lm);
        }

        if ( leptoncharge < 0 and cos(angle_CM_true)>0 )
		
        {
        Manager()->FillHisto("AFB_pos_Rap_No", mll);
	}
        
	if ( leptoncharge < 0 and cos(angle_CM_true)<0 )
        {
        Manager()->FillHisto("AFB_neg_Rap_No", mll);
	}
  }	

//if ( (signofyll_lab>0 and etalm<0) or (signofyll_lab<0 and etalm>0) ) 
  if ( (signofyll_lab>0 and eta_boosted_lm>0) or (signofyll_lab<0 and eta_boosted_lm<0) )
  {
  Manager()->FillHisto("AFB_pos_Rap_00", mll);
  Manager()->FillHisto("AFB_pos_Rap_01", mll);
  Manager()->FillHisto("AFB_pos_Rap_02", mll);
  Manager()->FillHisto("AFB_pos_Rap_04", mll);
  Manager()->FillHisto("AFB_pos_Rap_06", mll);
  }

//if ( (signofyll_lab>0 and etalm>0) or (signofyll_lab<0 and etalm<0) )
  if ( (signofyll_lab>0 and eta_boosted_lm<0) or (signofyll_lab<0 and eta_boosted_lm>0) )	  
  {
  Manager()->FillHisto("AFB_neg_Rap_00", mll);
  Manager()->FillHisto("AFB_neg_Rap_01", mll);
  Manager()->FillHisto("AFB_neg_Rap_02", mll);
  Manager()->FillHisto("AFB_neg_Rap_04", mll);
  Manager()->FillHisto("AFB_neg_Rap_06", mll);
  }


//  Manager()->FillHisto("Mll_Rap_No",  mll);
//  Manager()->FillHisto("Mll_Rap_00",  mll);
//  Manager()->FillHisto("Mll_Rap_01",  mll);  
//  Manager()->FillHisto("Mll_Rap_02",  mll);  
//  Manager()->FillHisto("Mll_Rap_04",  mll);
//  Manager()->FillHisto("Mll_Rap_06",  mll);
    
  myfile  << mll << ","
	  << mll_boosted << ","
	  << mjj << ","
	  << mlljj << ","
	  << MET_double << ","
          << MtMiss << ","
	  << MHT << ","
          << angle_CM_true << ","
          << angle_CM_rec << "," 
          << eta_boosted_lp << ","
	  << eta_boosted_lm << "," 
	  << etalp << ","
	  << etalm << ","
	  << quarkcharge << ","
	  << leptoncharge << ","
          << yll_lab << ","
	  << yll_CM << ","
	  << hasOS << ","
	  << ml1jj << ","
	  << ml2jj << ","
	  << mlpjj << ","
	  << mlmjj << ","
          << angle_CM_true_wrtl1 << ","
          << angle_CM_rec_wrtl1 << ","
          << EvWeight << "\n";

  return true;
}

void A_FB::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
	myfile.close();
}
