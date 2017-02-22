#! /usr/bin/env python

from __future__ import print_function

import argparse
import os

import ROOT

import utilities

def expandCut(cut):
    pass_ra2b="globalTightHalo2016Filter==1&&HBHENoiseFilter==1&&HBHEIsoNoiseFilter==1&&eeBadScFilter==1"
    if cut=="standard":
        return "nleps>=1&&st>500&&met>100"
    elif cut=="stdnj5":
        return "nleps>=1&&(ht+Sum$(leps_pt))>500&&met>150&&njets>=5"
    elif cut=="met150":
        return "nleps>=1&&st>500&&met>150&&met<=200&&njets>=5"
    elif cut=="higtight":
        return "met>200&&nvleps==0&&njets>=4&&njets<=5&&nbt>=2&&!low_dphi&&hig_drmax<2.2"
    elif cut=="abcd":
        return "nleps==1&&st>500&&max(met,met_tru)>200&&njets>=6&&nbm>=1&&mj14>250&&nveto==0"
    elif cut=="baseline":
        return "nleps==1&&st>500&&met>200&&njets>=6&&nbm>=1"
    elif cut=="sys_abcd":
        return "nleps==1&&max(ht,Max$(sys_ht))>500&&max(met,Max$(sys_met))>200&&max(njets,Max$(sys_njets))>=6&&max(nbm,Max$(sys_nbm))>=1&&max(mj,Max$(sys_mj))>250"
    elif cut=="zisr":
        return "nvleps==2&&nleps>=1&&Max$(leps_pt)>30&&((elelv_m>80&&elelv_m<100)||(mumuv_m>80&&mumuv_m<100))"
    elif cut=="dy_ht300":
        return "nvleps==2&&nleps>=1&&Max$(leps_pt)>30&&((elelv_m>80&&elelv_m<100)||(mumuv_m>80&&mumuv_m<100))&&ht>300"
    elif cut=="ttisr":
        return "nvleps==2&&nleps>=1&&max(Max$(mus_pt*(mus_tight&&mus_reliso<.1)),Max$(els_pt*(els_tight&&els_reliso<.1)))>30&&nbm==2"
    elif cut=="wisr":
        return "met>100&&max(Max$(mus_pt*(mus_tight&&mus_reliso<.1)),Max$(els_pt*(els_tight&&els_reliso<.1)))>30&&nbl==0"
    elif cut=="wisrht200":
        return "ht>200&&met>100&&max(Max$(mus_pt*(mus_tight&&mus_reliso<.1)),Max$(els_pt*(els_tight&&els_reliso<.1)))>30&&nbl==0"
    elif cut=="ttdilep_ht300":
        return "nels==1&&nmus==1&&Max$(leps_pt)>30&&ht>300&&met>100&&nbm>=1"
    elif cut=="qcd":
        return "ht>1000&&met<50&&(nvmus+nvels)==0"
    elif cut=="qcd_njet10":
        return "ht>1000&&met<50&&(nvmus+nvels)==0&&njets>=10"
    elif cut=="mm_std":
        return "Sum$(mm_nleps>=1&&mm_ht>500.&&mm_met>200.)>0"
    elif cut=="mm_std_nj5mj250":
        return "Sum$(mm_nleps>=1&&mm_ht>500&&mm_met>200&&mm_njets>=5&&mm_mj14_lep>250)>0||Sum$(mm_nleps>=1&&mm_ht>500&&mm_met>200&&mm_njets>=5&&mm_mj14_nolep>250)>0"
    elif cut=="ra2_qcd":
        return pass_ra2b+"&&(@Electrons.size()+@Muons.size())==0&&NJets>=3"
    elif cut=="ra2_ht300":
        return pass_ra2b+"&&HT>300"
    elif cut=="ra2_eht300":
        return pass_ra2b+"&&Max$(Electrons.Pt()*(abs(Electrons.Eta())<2))>35&&HT>300"
    elif cut=="ra2_zmht200":
        return pass_ra2b+"&&@ZCandidates.size()>=1&&MHT>200"
    else:
        return cut

def skimFiles(in_files, out_file, cut, keep_existing):
    in_files = [ utilities.fullPath(in_file) for in_file in in_files ]
    out_file = utilities.fullPath(out_file)

    utilities.ensureDir(os.path.dirname(out_file))

    cut = expandCut(cut)

    print("INPUT FILES:",in_files,"\n")
    print("OUTPUT FILE:",out_file,"\n")
    print("        CUT:",cut,"\n")

    if keep_existing and os.path.exists(out_file):
        print("Keeping pre-existing "+out_file+"\n")
        return

    in_tree = ROOT.TChain("tree", "tree")
    for in_file in in_files:
        in_tree.Add(in_file)

    with utilities.ROOTFile(out_file, "recreate") as out:
        out_tree = in_tree.CopyTree(cut)
        out_tree.Write()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Skims non-SMS ntuples.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("cut", help="Skim cut to apply.")
    parser.add_argument("out_file", help="File in which to save skim.")
    parser.add_argument("in_files", nargs="+", help="Files to skim.")
    parser.add_argument("-k","--keep_existing", action="store_true",
                        help="Do not overwrite output file if it already exists.")
    args = parser.parse_args()

    skimFiles(args.in_files, args.out_file, args.cut, args.keep_existing)
