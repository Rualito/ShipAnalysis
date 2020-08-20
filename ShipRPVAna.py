B# example for accessing smeared hits and fitted tracks
from __future__ import division
import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from decorators import *
import shipRoot_conf
shipRoot_conf.configure()

from ROOT import TMath
from ROOT import TLorentzVector
from ROOT import TVector3
from argparse import ArgumentParser
import struct

debug = False
chi2CutOff  = 4.
PDG = ROOT.TDatabasePDG.Instance()
inputFile  = None
geoFile    = None
dy         = None
nEvents    = 9999999
fiducialCut = True
measCutFK = 25
measCutPR = 22
docaCut = 2.
flag0 = 1

h = {}


# Mass_bins = 600
# ut.bookHist(h,'Invariant_Mass','Mass (GeV)',nbinsx=Mass_bins, xmin=0.0, xmax=3.0)
# h['Invariant_Mass'].GetXaxis().SetTitle('Mass (GeV)')
# h['Invariant_Mass'].GetYaxis().SetTitle('Counts')
# LMom_bins = 300
# ut.bookHist(h,'Linear_Momentum','Linear Momentum of the Daughter (GeV)',nbinsx=LMom_bins, xmin=0.0, xmax=300.0)
# h['Linear_Momentum'].GetXaxis().SetTitle('Momentum (GeV)')
# h['Linear_Momentum'].GetYaxis().SetTitle('Counts')
# TotLMom_bins = 300
# ut.bookHist(h,'Total_Linear_Momentum','Linear Momentum of the Mother (GeV)',nbinsx=TotLMom_bins, xmin=0.0, xmax=300.0)
# h['Total_Linear_Momentum'].GetXaxis().SetTitle('Momentum (GeV)')
# h['Total_Linear_Momentum'].GetYaxis().SetTitle('Counts')
# TMom_bins = 300
# ut.bookHist(h,'Transverse_Momentum','Transverse Momentum of the Daughter (GeV)',nbinsx=LMom_bins, xmin=0.0, xmax=5.0)
# h['Transverse_Momentum'].GetXaxis().SetTitle('Momentum (GeV)')
# h['Transverse_Momentum'].GetYaxis().SetTitle('Counts')
# TotTMom_bins = 300
# ut.bookHist(h,'Total_Transverse_Momentum','Transverse Momentum of the Mother (GeV/c)',nbinsx=TotTMom_bins, xmin=0.0, xmax=5.0)
# h['Total_Transverse_Momentum'].GetXaxis().SetTitle('Momentum (GeV)')
# h['Total_Transverse_Momentum'].GetYaxis().SetTitle('Counts')
# TLMom_bins = 300
# ut.bookHist(h,'TransLin_Momentum','Fraction of Transverse Momentum of the Daughters',nbinsx=TLMom_bins, xmin=0.0, xmax=1.0)
# h['TransLin_Momentum'].GetXaxis().SetTitle('Fraction of Transverse Momentum')
# h['TransLin_Momentum'].GetYaxis().SetTitle('Counts')
# TotTLMom_bins = 100
# ut.bookHist(h,'Total_TransLin_Momentum','Fraction of Transverse Momentum of the Mother',nbinsx=TotTLMom_bins, xmin=0.0, xmax=0.1)
# h['Total_TransLin_Momentum'].GetXaxis().SetTitle('Fraction of Transverse Momentum')
# h['Total_TransLin_Momentum'].GetYaxis().SetTitle('Counts')
# OpAng_bins = 300
# ut.bookHist(h,'Opening_Angle','Opening Angle (rad)',nbinsx=OpAng_bins, xmin=0, xmax=TMath.Pi()/4)
# h['Opening_Angle'].GetXaxis().SetTitle('Opening Angle (rad)')
# h['Opening_Angle'].GetYaxis().SetTitle('Counts')
# DecAng_bins = int(TMath.Pi()/0.0001)
# ut.bookHist(h,'Decay_Angle','Decay Angle (rad)', nbinsx=DecAng_bins, xmin=TMath.Pi()/2, xmax=TMath.Pi()*3/2)
# h['Decay_Angle'].GetXaxis().SetTitle('Decay Angle (rad)')
# h['Decay_Angle'].GetYaxis().SetTitle('Counts')
# ImpPar_bins = 3000
# ut.bookHist(h,'Impact_Parameter','Impact Parameter (cm)', nbinsx=ImpPar_bins, xmin=0, xmax=300)
# h['Impact_Parameter'].GetXaxis().SetTitle('Impact Parameter (cm)')
# h['Impact_Parameter'].GetYaxis().SetTitle('Counts')
# ZRes_bins = 1000
# ut.bookHist(h,'Z_at_Target_Resolution','Resolution of Z Reconstruction at the Target (cm)', nbinsx=ZRes_bins, xmin=0, xmax=10)
# h['Z_at_Target_Resolution'].GetXaxis().SetTitle('Resolution of Z (cm)')
# h['Z_at_Target_Resolution'].GetYaxis().SetTitle('Counts')
# ZRec_bins = 5036
# ut.bookHist(h,'Z_Rec','Z Reconstruction in the Vessel (cm)', nbinsx=ZRec_bins, xmin=-2481, xmax=2555)
# h['Z_Rec'].GetXaxis().SetTitle('Z position (cm)')
# h['Z_Rec'].GetYaxis().SetTitle('Counts')
# XRec_bins = 2000
# ut.bookHist(h,'X_Rec','X Reconstruction in the Vessel (cm)', nbinsx=XRec_bins, xmin=-1000, xmax=1000)
# h['X_Rec'].GetXaxis().SetTitle('X Reconstruction (cm)')
# h['X_Rec'].GetYaxis().SetTitle('Counts')
# YRec_bins = 1595
# ut.bookHist(h,'Y_Rec','Y Reconstruction in the Vessel (cm)', nbinsx=YRec_bins, xmin=-1000, xmax=595)
# h['Y_Rec'].GetXaxis().SetTitle('Y Reconstruction (cm)')
# h['Y_Rec'].GetYaxis().SetTitle('Counts')

# ut.bookHist(h,'Impact_vs_TransMom','Impact Parameter vs Transverse Momentum', nbinsx=150, xmin=0, xmax=300, nbinsy=50, ymin=0.0, ymax=5.0)
# h['Impact_vs_TransMom'].GetXaxis().SetTitle('Impact Parameter (cm)')
# h['Impact_vs_TransMom'].GetYaxis().SetTitle('Transverse Momentum (GeV/c)')
# ut.bookHist(h,'Impact_vs_TotLinMom','Impact Parameter vs Total Linear Momentum', nbinsx=150, xmin=0, xmax=300, nbinsy=300, ymin=0.0, ymax=300.0)
# h['Impact_vs_TotLinMom'].GetXaxis().SetTitle('Impact Parameter (cm)')
# h['Impact_vs_TotLinMom'].GetYaxis().SetTitle('Linear Momentum (GeV/c)')
# ut.bookHist(h,'Impact_vs_TotTransMom','Impact Parameter vs Total Transverse Momentum', nbinsx=150, xmin=0, xmax=300, nbinsy=50, ymin=0.0, ymax=5.0)
# h['Impact_vs_TotTransMom'].GetXaxis().SetTitle('Impact Parameter (cm)')
# h['Impact_vs_TotTransMom'].GetYaxis().SetTitle('Transverse Momentum (GeV/c)')
# ut.bookHist(h,'Impact_vs_TotTransLinMom','Impact Parameter vs Fraction of Transverse Momentum', nbinsx=150, xmin=0, xmax=300, nbinsy=100, ymin=0.0, ymax=0.1)
# h['Impact_vs_TotTransLinMom'].GetXaxis().SetTitle('Impact Parameter (cm)')
# h['Impact_vs_TotTransLinMom'].GetYaxis().SetTitle('Fraction of Transverse Momentum')
# ut.bookHist(h,'Impact_vs_Opening','Impact Parameter vs Opening_Angle', nbinsx=150, xmin=0, xmax=300, nbinsy=300, ymin=0.0, ymax=TMath.Pi()/4)
# h['Impact_vs_Opening'].GetXaxis().SetTitle('Impact Parameter (cm)')
# h['Impact_vs_Opening'].GetYaxis().SetTitle('Opening Angle (rad)')
# ut.bookHist(h,'LinMom_vs_Opening','Linear Momentum vs Opening_Angle', nbinsx=300, xmin=0, xmax=300, nbinsy=300, ymin=0.0, ymax=TMath.Pi()/4)
# h['LinMom_vs_Opening'].GetXaxis().SetTitle('Linear Momentum (GeV/c)')
# h['LinMom_vs_Opening'].GetYaxis().SetTitle('Opening Angle (rad)')
# ut.bookHist(h,'X_vs_Y','Reconstruction of X vs Y', nbinsx=XRec_bins, xmin=-1000, xmax=1000, nbinsy=YRec_bins, ymin=-1000, ymax=595)
# h['X_vs_Y'].GetXaxis().SetTitle('X Reconstruction (cm)')
# h['X_vs_Y'].GetYaxis().SetTitle('Y Reconstruction (cm)')
# ut.bookHist(h,'Z_vs_R','Reconstruction of Z vs the Transversal Position', nbinsx=ZRec_bins, xmin=-2481, xmax=2555, nbinsy=int(TMath.Sqrt(2)*1000), ymin=0, ymax=int(TMath.Sqrt(2)*1000))
# h['Z_vs_R'].GetXaxis().SetTitle('Z Reconstruction (cm)')
# h['Z_vs_R'].GetYaxis().SetTitle('Sqrt(X^2 + Y^2) Reconstruction (cm)')
# ut.bookHist(h,'Z_vs_Impact','Reconstruction of Z vs the Impact Parameter', nbinsx=ZRec_bins, xmin=-2481, xmax=2555, nbinsy=600, ymin=0, ymax=300)
# h['Z_vs_Impact'].GetXaxis().SetTitle('Z Reconstruction (cm)')
# h['Z_vs_Impact'].GetYaxis().SetTitle('Impact Parameter (cm)')
# ut.bookHist(h,'Impact_vs_R','Impact Parameter vs the Transversal Position', nbinsx=600, xmin=0, xmax=300, nbinsy=int(TMath.Sqrt(2)*1000), ymin=0, ymax=int(TMath.Sqrt(2)*1000))
# h['Impact_vs_R'].GetXaxis().SetTitle('Impact Parameter (cm)')
# h['Impact_vs_R'].GetYaxis().SetTitle('Sqrt(X^2 + Y^2) Reconstruction (cm)')



k0L_count = 0
n_e_count = 0
K_count = 0
e_count = 0
n_mu_count = 0
mu_count = 0

k0L_pdg = 130 
e_pdg = 11
n_e_pdg = 12
 
n_mu_pdg = 14
mu_pdg = 13

countDict = dict()
# chargeDict = dict()

def isInFiducial(X,Y,Z):
   
   if Z > ShipGeo.TrackStation1.z : return False
   if Z < ShipGeo.tauMudet.zMudetC + ShipGeo.tauMudet.Ztot/2. : return False #ShipGeo.vetoStation.z+100.*u.cm : return False
   if dist2InnerWall(X,Y,Z)<5*u.cm: return False
   
   return True 

from array import array
def dist2InnerWall(X,Y,Z):
  dist = 0
 # return distance to inner wall perpendicular to z-axis, if outside decayVolume return 0.
  node = sGeo.FindNode(X,Y,Z)
  # If X, Y or Z are out of bounds we want to avoid checking the name of something that doesn't exist, so we just return it as being out of the range of interest
  try:
    if ShipGeo.tankDesign < 5:
       if not 'cave' in node.GetName(): return dist  # TP 
    else:
       if not 'block' in node.GetName(): return dist
  except:
    return dist
  start = array('d',[X,Y,Z])
  nsteps = 8
  dalpha = 2*ROOT.TMath.Pi()/nsteps
  rsq = X**2+Y**2
  minDistance = 100 *u.m
  for n in range(nsteps):
    alpha = n * dalpha
    sdir  = array('d',[ROOT.TMath.Sin(alpha),ROOT.TMath.Cos(alpha),0.])
    node = sGeo.InitTrack(start, sdir)
    nxt = sGeo.FindNextBoundary()
    if ShipGeo.tankDesign < 5 and nxt.GetName().find('I')<0: return 0    
    distance = sGeo.GetStep()
    if distance < minDistance  : minDistance = distance
  return minDistance


# Need to check why they do this implementation
def ImpactParameter(point,tPos,tMom):
  t = 0
  if hasattr(tMom,'P'): P = tMom.P()
  else:			P = tMom.Mag()

  for i in range(3):	t += tMom(i)/P*(point(i)-tPos(i))
  dist = 0
  for i in range(3):	dist += (point(i)-tPos(i)-t*tMom(i)/P)**2
  dist = TMath.Sqrt(dist)
  return dist


# Z seems reasonable but X and Y are pure garbage (This is the same implemention as in Impact_Parameter)
def ImpactPosition(point,tPos,tMom):
  t = 0
  if hasattr(tMom,'P'):	P = tMom.P()
  else:			P = tMom.Mag()
  for i in range(3):	t += tMom(i)/P*(point(i)-tPos(i))
  ImpPos = TVector3(tPos.X()+t*tMom.X()/P, tPos.Y()+t*tMom.Y()/P, tPos.Z()+t*tMom.Z()/P)
  return ImpPos

def myVertex(t1,t2,PosDir):
  # Calculating closest distance between 2 tracks
  # d = |pq . u x v|/|u x v|
  a = TVector3(PosDir[t1][0](0), PosDir[t1][0](1), PosDir[t1][0](2))
  u = TVector3(PosDir[t1][1](0), PosDir[t1][1](1), PosDir[t1][1](2))
  c = TVector3(PosDir[t2][0](0), PosDir[t2][0](1), PosDir[t2][0](2))
  v = TVector3(PosDir[t2][1](0), PosDir[t2][1](1), PosDir[t2][1](2))
  pq = a-c
  uCrossv = u.Cross(v)
  dist = pq.Dot(uCrossv)/(uCrossv.Mag()+1E-8)
  E = u.Dot(a) - u.Dot(c)
  F = v.Dot(a) - v.Dot(c)
  A,B = u.Mag2(), -u.Dot(v)
  C,D = u.Dot(v), -v.Mag2()
  t = -(C*E-A*F)/(B*C-A*D)
  X = c.x()+v.x()*t
  Y = c.y()+v.y()*t
  Z = c.z()+v.z()*t
  return X,Y,Z,TMath.Abs(dist)

def match2HSP(p):
  matched = False
  hnlKey = []
  for t in [p.GetDaughter(0),p.GetDaughter(1)]:
    mcp = sTree.fitTrack2MC[t]
    while mcp > -0.5:
      mo = sTree.MCTrack[mcp]
      if TMath.Abs(mo.GetPdgCode()) == 9900015 or TMath.Abs(mo.GetPdgCode()) == 4900023:
        hnlKey.append(mcp)
        break
      mcp = mo.GetMotherId()
  if len(hnlKey) == 2:
    if hnlKey[0]==hnlKey[1]: matched = True
    return matched


def myEventLoop(n):
   
  rc = sTree.GetEntry(n)
 
  #    if not isInFiducial(HNLMom.X(),HNLMom.Y(),HNLMom.Z()): continue 

  global Events_w_more_than_2_hits
  global Total_rec
  global Rec_within_decay_vessel 
  global Rec_with_2_part
  
  global countDict

  if not hasattr(sTree, 'Particles'): return 0
  
  genDaughters = []
  for n, particle in enumerate(sTree.MCTrack):
     if n!=0 and sTree.MCTrack[particle.GetMotherId()].GetPdgCode() == 9900015 and n<10:
        genDaughters.append(particle.GetPdgCode())

  genDaughtersID = ""
  for id in sorted(genDaughters):
     genDaughtersID+=str(PDG.GetParticle(id).GetName())+" "

  if genDaughtersID not in countDict:
     countDict[genDaughtersID]=dict()
     
  l = 0
  for phit in sTree.TimeDetPoint:
    l = l + 1
    if l > 1: break

  if l > 1:
    
    # Here there is either 1 or 0 HNL particles reconstructed, but might be different in other cases

    Events_w_more_than_2_hits += 1
    for HNL in sTree.Particles:

      Total_rec += 1
 
      #TParticle.ProductionVertex defines a TLorentzVector as the point where TParticle decays
      HNLPos = TLorentzVector()
      HNL.ProductionVertex(HNLPos)
      HNLMom = TLorentzVector()
      HNL.Momentum(HNLMom)

      # If the HNL (or reconstructed particle) has decayed before the Muon Detector, after the TrackStation1 or within 5 cm of the Walls, we skip the event
      if not isInFiducial(HNLPos.X(),HNLPos.Y(),HNLPos.Z()): continue

      Rec_within_decay_vessel += 1

      # Trying to find a Pi0. If there is one we say the reconstructed particle has pi0's momentum added to it
      Ispi0 = False
      if hasattr(sTree, 'EcalReconstructed'):
        possiblepi0 = []
        for pi0 in pi0Reco.findPi0(sTree,HNLPos):
          if TMath.Abs(pi0.M()-0.135)>0.02: continue #0.02
          Ispi0 = True
          possiblepi0.append(pi0)
          
      if Ispi0:
        maxmom = 0
        num = 0
        for index, pi0 in enumerate(possiblepi0):
          if pi0.P() > maxmom:
            maxmom = pi0.P()
            num = index
        HNLwithPi0 = HNLMom + possiblepi0[num]


      t1,t2 = HNL.GetDaughter(0),HNL.GetDaughter(1)
      #Drawing the daughter related histograms      
      if t2 - t1 == 1 or True:
        Rec_with_2_part += 1
        # Drawing the invariant mass histogram
        if not Ispi0:
          inv_mass = HNL.GetMass()
          tot_Pt = HNLMom.Pt()
          tot_P = HNLMom.P()
        else:
          inv_mass = HNLwithPi0.M()
          tot_Pt = HNLwithPi0.Pt()
          tot_P = HNLwithPi0.P()
        
        #Daughters = []
        daughtersID = ""
        daughtersIDarr = []
        for i in range(0, t2-t1):
           currentDaughter = sTree.FitTracks[t1+i].getFittedState() 
           #Daughters.append(currendDaughter)
           daughtersIDarr.append(currentDaughter.getPDG())

        if Ispi0:
           daughtersIDarr.append(111)
        
        # print(sorted(daughtersIDarr))
        for idS in sorted(daughtersIDarr):
           daughtersID+=str(PDG.GetParticle(idS).GetName()) + " "
      
        Daughter1 = sTree.FitTracks[t1].getFittedState()
        Daughter2 = sTree.FitTracks[t2].getFittedState()
        #print(Daughter1.getPDG())
	#print(Daughter2.getPDG())
        
        #if Daughter1.getPDG() not in countDict:
           #countDict[Daughter1.getPDG()] = 0
        #if Daughter2.getPDG() not in countDict:
           #countDict[Daughter2.getPDG()] = 0
      
        #countDict[Daughter1.getPDG()]+=1
        #countDict[Daughter2.getPDG()]+=1
        

        if daughtersID not in countDict[genDaughtersID]:
           countDict[genDaughtersID][daughtersID] = 0
           
        countDict[genDaughtersID][daughtersID]+=1




###############################################################
###############################################################
########## READING THE FILES WITH THE SIMULATIONS #############
###############################################################
###############################################################

 
import os

Events_w_more_than_2_hits = 0
Total_rec = 0
Rec_within_decay_vessel = 0
Rec_with_2_part = 0

rootdir = '/eos/home-g/gmachado/gmachado/VACUUM/HNL/HNLtoPie/Mass=1.0/Coupling=5e-10/e/'
parser = ArgumentParser()
parser.add_argument("-r", "--rootdir", dest="rootdir", help="Root directory from which the macro starts searching \'_rec.root\' files", required=False, default=False)

options = parser.parse_args()

if options.rootdir:
  rootdir = options.rootdir


flag_Pi = False

for subdir, dirs, files in os.walk(rootdir):
  for file in files:
    if file.startswith("ship.conical") and file.endswith("_rec.root"):
      inputFile =  os.path.join(subdir, file)

      f = ROOT.TFile(inputFile)
      if not hasattr(f,'cbmsim'): continue
      sTree = f.cbmsim
      print(" INPUT FILE: ", inputFile, sTree.GetEntries())
      nEvents = 9999999

      if flag0:
        
        if not geoFile:
          geoFile = inputFile.replace('ship.','geofile_full.').replace('_rec.','.')
        fgeo = ROOT.TFile(geoFile)

        # new geofile, load Shipgeo dictionary written by run_simScript.py
        upkl    = Unpickler(fgeo)
        ShipGeo = upkl.load('ShipGeo')
        ecalGeoFile = ShipGeo.ecal.File
        dy = ShipGeo.Yheight/u.m

        # -----Create geometry----------------------------------------------
        import shipDet_conf
        run = ROOT.FairRunSim()
        run.SetName("TGeant4")  # Transport engine
        run.SetOutputFile(ROOT.TMemFile('output', 'recreate'))  # Output file
        run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for the mag field
        vrtdb = run.GetRuntimeDb()
        # -----Create geometry----------------------------------------------
        print(ShipGeo)
        modules = shipDet_conf.configure(run,ShipGeo)

	

        import geomGeant4
        if hasattr(ShipGeo.Bfield,"fieldMap"):
          fieldMaker = geomGeant4.addVMCFields(ShipGeo, '', True, withVirtualMC = False)
        else:
          print("no fieldmap given, geofile too old, not anymore support")
          exit(-1)
        sGeo   = fgeo.FAIRGeom
        geoMat =  ROOT.genfit.TGeoMaterialInterface()
        ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
        bfield = ROOT.genfit.FairShipFields()
        bfield.setField(fieldMaker.getGlobalField())
        fM = ROOT.genfit.FieldManager.getInstance()
        fM.init(bfield)

        flag0 = 0

        if not flag_Pi:
          import pi0Reco
          flag_Pi = True

      nEvents = min(sTree.GetEntries(),nEvents)

      #for n in range(10):
      for n in range(nEvents):
         myEventLoop(n)


	    

if inputFile == None:
   print("there is no \"ship.conical\" file in the desired directory")


# output histograms
hfile = 'tree_original.root'
ROOT.gROOT.cd()
ut.writeHists(h,hfile)

# reconstructed_events.write("Events with more than 2 hits on a detector = " + str(Events_w_more_than_2_hits) + "\n")
# reconstructed_events.write("Total number of reconstructions = " + str(Total_rec) + "\n")
# reconstructed_events.write("Total number of reconstructions within the decay vessel = " + str(Rec_within_decay_vessel) + "\n")
# reconstructed_events.write("Events reconstructed with 2 daughter particles that decayed within the decay vessel = " + str(Rec_with_2_part) + "\n")


# close some files
# reconstructed_events.close()
# data_file.close()

for genPdgcode in countDict:
   for dauPdgcode in countDict[genPdgcode]: 
      print(genPdgcode,":", dauPdgcode, ":",  countDict[genPdgcode][dauPdgcode])
