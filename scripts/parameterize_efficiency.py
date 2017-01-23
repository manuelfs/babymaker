#! /usr/bin/env python

import argparse
import ROOT
import array

def ParameterizeEfficiency(out_file_path, docuts):
    ROOT.TH1.SetDefaultSumw2()
    c = ROOT.TChain("tree", "tree")
    c.Add("/net/cms2/cms2r0/babymaker/babies/2017_01_21/mc/unprocessed/*_TTJets_*.root")

    btags = [("loose", 0.5426),
             ("medium", 0.8484),
             ("tight", 0.9535)]

    btags_deep = [("loose", 0.2219),
                  ("medium", 0.6324),
                  ("tight", 0.8958)]

    eta_cuts = [0., 1.2, 2.4]
    pt_cuts = [30., 50., 70., 100., 140., 200., 300., 670., 1.e4]
    flavor_cuts = [-0.5, 3.5, 4.5, 5.5]

    numerators = [ROOT.TH3D("btagEfficiency_"+btag[0], "btagEfficiency_"+btag[0],
                            len(eta_cuts)-1, array.array('d', eta_cuts),
                            len(pt_cuts)-1, array.array('d', pt_cuts),
                            len(flavor_cuts)-1, array.array('d', flavor_cuts))
                  for btag in btags]
    numerators_deep = [ROOT.TH3D("btagEfficiency_deep_"+btag[0], "btagEfficiency_deep_"+btag[0],
                            len(eta_cuts)-1, array.array('d', eta_cuts),
                            len(pt_cuts)-1, array.array('d', pt_cuts),
                            len(flavor_cuts)-1, array.array('d', flavor_cuts))
                  for btag in btags_deep]

    denominator = numerators[0].Clone("btagEfficiency_denominator")

    entry = 0
    num_entries = c.GetEntries()
    for event in c:
        if entry % 10000 == 0:
            print "Completed", '{:.2f}'.format(100.*entry/num_entries)+"%"
        entry = entry + 1
        # if not (c.stitch and getattr(c,"pass")): continue # if using HT bins
        if not (getattr(c,"pass")): continue
        if docuts:
            if (c.nleps<1 or c.st<=500. or c.met<=200. or c.njets<5): 
                continue
        for ijet in xrange(len(c.jets_csv)):
            if (c.jets_islep[ijet]): continue
            flavor = abs(c.jets_hflavor[ijet])
            pt = c.jets_pt[ijet]
            eta = abs(c.jets_eta[ijet])
            csv = c.jets_csv[ijet]
            csvd = c.jets_csvd[ijet]

            denominator.Fill(eta, pt, flavor)
            for inum in xrange(len(btags)):
                if csv > btags[inum][1]:
                    numerators[inum].Fill(eta, pt, flavor)
            for inum in xrange(len(btags_deep)):
                if csvd > btags_deep[inum][1]:
                    numerators_deep[inum].Fill(eta, pt, flavor)

    for n in numerators:
        n.Divide(denominator)
    for n in numerators_deep:
        n.Divide(denominator)
                
    out_file = ROOT.TFile(out_file_path, "recreate")
    for n in numerators:
        n.Write()
    doc = ROOT.TNamed("Documentation: this file contains a parameterization in (eta, pt, flavor) for CSVv2 b-tagger", "")
    doc.Write()
    out_file.Close()

    out_file = ROOT.TFile(out_file_path.replace(".root","_deep.root"), "recreate")
    for n in numerators_deep:
        n.Write()
    doc = ROOT.TNamed("Documentation: this file contains a parameterization in (eta, pt, flavor) for the deep CSV b-tagger", "")
    doc.Write()
    out_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes b-tagging efficiency as a function of pT, eta, and flavor",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--out_file", metavar="OUT_FILE", default="btagEfficiency.root",
                        help="Save efficiences to %(metavar)s")
    parser.add_argument("-c", "--docuts", action="store_true", 
                        help="Use all available events, applying only basic filters")
    args = parser.parse_args()

    ParameterizeEfficiency(args.out_file, args.docuts)
