// cuts: (CSVL, CSVM, CSVT) = (0.605, 0.890, 0.970)
namespace btagCuts {
  float CSVL=0.605;
  float CSVM=0.890;
  float CSVT=0.970;
};

void parameterizeEfficiency()
{

  float csvCut = btagCuts::CSVM;

  std::string filename("baby_TT_uncorrected.root");
  TFile *file = TFile::Open(filename.c_str());
  TTree *tree = static_cast<TTree*>(file->Get("tree"));

  // These are the cuts that would be needed in a full parameterization
  //  std::vector<float> etaCuts = { 0.0, 1.2, 2.4, 5.0 };
  // std::vector<float> ptCuts = { 30, 50, 70, 100, 140, 
  // 				200, 300, 670, 1e4};
  std::vector<float> etaCuts = { 0.0, 1.2, 2.4};
  std::vector<float> ptCuts = { 30, 50, 70, 100, 140, 
				200, 1e4};
  std::vector<int> flavorCuts = {0, 4, 5};
  std::vector<float> flavorCutsFloat = {0, 4, 5};

  TH1F *hCSV[etaCuts.size()-1][ptCuts.size()-1][flavorCuts.size()];
  TH3F *btagEfficiency = new TH3F("btagEfficiency", "btagEfficiency", 
			       etaCuts.size()-1, &etaCuts[0],
			       ptCuts.size()-1, &ptCuts[0],
			       flavorCuts.size(), &flavorCutsFloat[0]);

  for(int ieta(0); ieta<etaCuts.size()-1; ieta++ ) {
    for(int ipt(0); ipt<ptCuts.size()-1; ipt++ ) {
      for(int iflavor(0); iflavor<flavorCuts.size(); iflavor++) {
	TString histname(Form("csv_%d_%d_%d",ieta, ipt, iflavor));
	hCSV[ieta][ipt][iflavor] = new TH1F(histname, histname, 2, -0.5, 1.5);
	TString cut(Form("jets_hadronFlavour==%d && %f<=jets_eta && jets_eta <%f && %f <= jets_pt && jets_pt < %f", 
			 flavorCuts.at(iflavor), etaCuts.at(ieta), etaCuts.at(ieta+1), ptCuts.at(ipt), ptCuts.at(ipt+1)));
	TString cutWithCSV=cut;
	cutWithCSV+=" && jets_csv> ";
	cutWithCSV+=csvCut;
	std::cout << cut << std::endl;      
	std::cout << cutWithCSV << std::endl;      
	tree->Project(histname.Data(), cutWithCSV.Data(), cut.Data());
	std::cout << "efficiency: " << hCSV[ieta][ipt][iflavor]->GetMean() << "\n" << std::endl;
	btagEfficiency->SetBinContent(ieta+1, ipt+1, iflavor+1, hCSV[ieta][ipt][iflavor]->GetMean());
      }
    }
  }

  TFile *output = new TFile("bmaker/data/btagEfficiency.root", "recreate");  
  btagEfficiency->Write();
  TNamed documentation("Documentation: this file contains a parameterization in (eta, pt, flavor)", "");
  documentation.Write();
  output->Close();

}
