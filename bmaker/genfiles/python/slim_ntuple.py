#! /usr/bin/env python

import argparse
import fnmatch
import ROOT

def GetBranchNames(tree):
    return [ branch.GetName() for branch in tree.GetListOfBranches() ]

def GetRules(slim_file_name):
    rules = [ line.rstrip('\n').split() for line in open(slim_file_name) ]
    good_rules = [ rule for rule in rules if len(rule)==2 and (rule[0]=="keep" or rule[0]=="drop") ]
    bad_rules = [ rule for rule in rules if rule not in good_rules ]
    for rule in bad_rules:
        print "Invalid rule:",rule
    return good_rules

def PassRules(branch, rules):
    matched_rules =  [ rule for rule in rules if fnmatch.fnmatch(branch, rule[1]) ]
    return len(matched_rules)==0 or matched_rules[-1][0] == "keep"

def SlimNtuple(slim_file_name, output_file_name, input_file_names, test_mode):
    print "     INPUT FILES:",input_file_names
    print "     OUTPUT FILE:",output_file_name
    print "      RULES FILE:",slim_file_name

    output_file = ROOT.TFile(output_file_name, "recreate")

    in_tree = ROOT.TChain("tree", "tree")
    in_treeglobal = ROOT.TChain("treeglobal", "treeglobal")
    for input_file_name in input_file_names:
        in_tree.Add(input_file_name)
        in_treeglobal.Add(input_file_name)

    branch_names = GetBranchNames(in_tree)
    rules = GetRules(slim_file_name)
    kept_branches = [ branch for branch in branch_names if PassRules(branch, rules) ]
    kept_branches.sort()
    dropped_branches = [ branch for branch in branch_names if branch not in kept_branches ]
    dropped_branches.sort()

    print "DROPPED BRANCHES:",dropped_branches
    print "   KEPT BRANCHES:",kept_branches
    if test_mode: return

    for branch in branch_names:
        if branch in kept_branches: in_tree.SetBranchStatus(branch, True)
        else:                       in_tree.SetBranchStatus(branch, False)

    out_tree = in_tree.CloneTree()
    out_treeglobal = in_treeglobal.CloneTree()

    out_tree.Write()
    out_treeglobal.Write()
    output_file.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prunes branches from an ntuple",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-t", "--test", action="store_true",
                        help="Run in test mode, quickly diplaying the list of kept and dropped branchs without actually copying the trees.")
    parser.add_argument("slim_file",
                        help="File containing rules for pruning branches (one rule per line). Rules are are the form \"keep XXX\" or \"drop YYY\". Unix shell-style wildcards (e.g., '*') allow pattern matching. Branches are kept by default if no matching rule is found for the branch. If multiple rules match, the last takes precedence.")
    parser.add_argument("output_file",
                        help="File in which to save the slimmed and merged ntuple.")
    parser.add_argument("input_files", nargs="+",
                        help="Files containing ntuples to be slimmed and merged.")
    args = parser.parse_args()

    SlimNtuple(args.slim_file, args.output_file, args.input_files, args.test)
