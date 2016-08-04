#! /usr/bin/env python

import argparse
import ROOT
import array

def ParameterizeEfficiency(out_file_path, ttbar_only, no_cuts):
    c = ROOT.TChain("tree", "tree")
    subdir = ("unskimmed" if no_cuts else "merged_stdnj5")
    if ttbar_only:
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_TTJets_Tune*.root")
    else:
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_TTJets*Lept*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_TTJets_HT*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_WJetsToLNu*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_ST_*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_TTWJets*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_TTZ*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*DYJetsToLL*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*QCD_HT*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_ZJet*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_ttHJetTobb*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_TTGJets*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_TTTT*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_WH_HToBB*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_ZH_HToBB*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_WWTo*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_WZ*.root")
        c.Add("/net/cms2/cms2r0/babymaker/babies/2016_06_14/mc/"+subdir+"/*_ZZ_*.root")

    btags = [("loose", 0.460),
             ("medium", 0.800),
             ("tight", 0.935)]

    eta_cuts = [0., 1.2, 2.4]
    pt_cuts = [30., 50., 70., 100., 140., 200., 300., 670., 1.e4]
    flavor_cuts = [-0.5, 3.5, 4.5, 5.5]

    numerators = [ROOT.TH3D("btagEfficiency_"+btag[0], "btagEfficiency_"+btag[0],
                            len(eta_cuts)-1, array.array('d', eta_cuts),
                            len(pt_cuts)-1, array.array('d', pt_cuts),
                            len(flavor_cuts)-1, array.array('d', flavor_cuts))
                  for btag in btags]
    denominator = numerators[0].Clone("btagEfficiency_denominator")

    entry = 0
    num_entries = c.GetEntries()
    for event in c:
        if entry % 1000 == 0:
            print str(100.*entry/num_entries)+"%"
        entry = entry + 1
        if not (c.stitch and getattr(c,"pass")): continue
        pass_cut = (True if no_cuts else (c.nleps>=1 and c.ht>500. and c.met>200. and c.njets>=5))
        if not pass_cut: continue
        for ijet in xrange(len(c.jets_csv)):
            flavor = abs(c.jets_hflavor[ijet])
            pt = c.jets_pt[ijet]
            eta = abs(c.jets_eta[ijet])
            weight = c.weight
            csv = c.jets_csv[ijet]

            denominator.Fill(eta, pt, flavor, weight)
            for inum in xrange(len(btags)):
                if csv > btags[inum][1]:
                    numerators[inum].Fill(eta, pt, flavor, weight)

    for n in numerators:
        n.Divide(denominator)
                
    out_file = ROOT.TFile(out_file_path, "recreate")
    for n in numerators:
        n.Write()
    doc = ROOT.TNamed("Documentation: this file contains a parameterization in (eta, pt, flavor)", "")
    doc.Write()
    out_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes b-tagging efficiency as a function of pT, eta, and flavor",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t", "--ttbar_only", action="store_true",
                        help="Compute efficiency using only inclusive ttbar")
    parser.add_argument("-o", "--out_file", metavar="OUT_FILE", default="bmaker/data/btagEfficiency.root",
                        help="Save efficiences to %(metavar)s")
    parser.add_argument("-i", "--inclusive", action="store_true",
                        help="Use all available events, applying only basic filters")
    args = parser.parse_args()

    ParameterizeEfficiency(args.out_file, args.ttbar_only, args.inclusive)
