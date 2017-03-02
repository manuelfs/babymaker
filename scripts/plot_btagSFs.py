#! /usr/bin/env python
import argparse
from array import array
from sys import argv
#Batch mode hack
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Plot the btag SFs')
    parser.add_argument('infile', help='Path to CSV SF file to read in')
    parser.add_argument('type', help='Measurement type: incl, mujets, or comb')
    parser.add_argument('jetflavor', help='0,1 for b,c (incl/mujets). 2 for light (incl).')
    parser.add_argument('-wp','--workingpoint', help='Working point to plot. Default = 1 (medium)')
    args = parser.parse_args()

    # Read in file
    with open(args.infile,'r') as file:
        info = [line.split(', ') for line in file]

    #Relevant information
    wp = args.workingpoint if args.workingpoint else '1'
    mtype = args.type
    jetflav = args.jetflavor
    systype = ['central','up','down']
    ptBins = ['20'] if mtype=='incl' else ['20','30','50','70','100','140','200','300','600','1000']
    
    #Make a map between the systype and pt to the corresponding SF function
    mapping = dict()
    for stype in systype:
        mapping[stype] = dict()
        for pt in ptBins:
            # Central only has one bin
            if stype=='central' and pt!='20': continue
            for line in info:
                if line[0]==wp and line[1]==mtype and line[2]==stype and line[3]==jetflav and line[6]==pt:
                    mapping[stype][pt] = line[10].strip().strip('"')

    #Setup ROOT plots
    xpts = [20., 25., 30., 40., 50., 60., 70., 85., 100., 120., 140., 170., 200., 250., 300., 450., 600., 800., 1000., 1500.]
    xerr = [0., 5., 0., 10., 0., 10., 0., 15., 0., 20., 0., 30., 0., 50., 0., 150., 0., 200., 0., 500.]
    yerr = [0.]*len(xpts)

    # Get central, up and down values
    ypts = []
    for stype in systype:
        if mtype=='incl':
            #Easier since only one pt bin for all systypes
            ytmp = [eval(mapping[stype]['20'].replace('x',str(x))) for x in xpts]
            ypts.append(ytmp)

        else:
            if stype=='central':
                ytmp = [eval(mapping['central']['20'].replace('x',str(x))) for x in xpts]
            else:
                #Use the function of the lower bin to handle the points in the middle of bins
                ytmp = [eval(mapping[stype][str(int(xpts[i]))].replace('x',str(xpts[i]))) if i%2==0 else eval(mapping[stype][str(int(xpts[i-1]))].replace('x',str(xpts[i]))) 
                        for i in range(len(xpts)-2)]
                
                #Handle special cases for pt=1000,1500 since these go beyond the pt binning of the SFs
                ytmp.append(eval(mapping[stype]['600'].replace('x',str(xpts[len(xpts)-2]))))
                ytmp.append(eval(mapping[stype]['600'].replace('x',str(xpts[len(xpts)-1]))))
        
            ypts.append(ytmp)

    yerrlo = [abs(i-j) for i,j in zip(ypts[0],ypts[2])]
    yerrhi = [abs(i-j) for i,j in zip(ypts[0],ypts[1])]

    #Make the plots
    c = ROOT.TCanvas()

    graph = ROOT.TGraphAsymmErrors(len(xpts),array('d',xpts),array('d',ypts[0]),array('d',xerr),array('d',xerr),array('d',yerrlo),array('d',yerrhi))
             
    #Prety plot
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1)
    graph.SetMarkerColor(2)
    graph.SetLineColor(2)
    graph.SetLineWidth(3)
    graph.SetFillColor(4)
    graph.SetFillStyle(3013)

    title = 'b-flavor SFs' if jetflav=='0' else 'c-flavor' if jetflav=='1' else 'Light-flavor SFs'
    graph.SetTitle(title+' for '+mtype+' measurement')
    graph.SetName('graph_mtype'+mtype+'_jetflav'+jetflav)
    graph.SetMinimum(0.8)
    graph.SetMaximum(1.2)
    graph.GetXaxis().SetRangeUser(30,1500)
    graph.GetXaxis().SetMoreLogLabels(1)
    graph.Draw('A02CP')
    
    c.SetGridy(1)
    title = 'sfs_'+mtype+'_jetflav'+jetflav+'.png'
    c.SetLogx(1)
    c.SaveAs(title)

# If you want to save the graph to a root file
#    f = ROOT.TFile('btagSFs.root','update')
#    graph.Write()
#    f.Write()
