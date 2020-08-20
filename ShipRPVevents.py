# example for accessing smeared hits and fitted tracks
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
 
  
  global Events_w_more_than_2_hits
  global Total_rec
  global Rec_within_decay_vessel 
  global Rec_with_2_part
  
  global countDict
  global event_file 
  global current_file
  global reconStats

  if not hasattr(sTree, 'Particles'): return 0

  # event_file.write(current_file + ", ") #
  event_file.write(str(n)+", ")
  
  genLeptons = 0
  genDaughters = []
  for n, particle in enumerate(sTree.MCTrack):
     if n!=0 and sTree.MCTrack[particle.GetMotherId()].GetPdgCode() == 9900015 and n<10:
        genDaughters.append(particle.GetPdgCode())
        if 11 <= abs(genDaughters[-1]) <= 16:
           genLeptons += (genDaughters[-1]>0)*2 -1

  genDaughtersID = ""
  for id in sorted(genDaughters):
     genDaughtersID+=str(PDG.GetParticle(id).GetName())+" "

  if genDaughtersID not in countDict:
     countDict[genDaughtersID]=dict()
  
  event_file.write(genDaughtersID+", ")
  

  # if there's at least one point in ...
  l = 0
  for phit in sTree.TimeDetPoint:
    l = l + 1
    if l > 1: break

  if l > 1:
   
    # Here there is either 1 or 0 HNL particles reconstructed, but might be different in other cases

    Events_w_more_than_2_hits += 1
    
    for HNL in sTree.Particles: # all of the reconstructed particles?

      Total_rec += 1
      
      #TParticle.ProductionVertex defines a TLorentzVector as the point where TParticle decays
      HNLPos = TLorentzVector()
      HNL.ProductionVertex(HNLPos)
      HNLMom = TLorentzVector()
      HNL.Momentum(HNLMom)

      # If the HNL (or reconstructed particle) has decayed before the Muon Detector, after the TrackStation1 or within 5 cm of the Walls, we skip the event
      # if it decays inside the decay vessel, with some tolerance
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
        chargeTotal = 0
        leptonCount = 0
        for i in range(0, t2-t1+1):
           currentDaughter = sTree.FitTracks[t1+i].getFittedState() 
           #Daughters.append(currendDaughter)
           daughtersIDarr.append(currentDaughter.getPDG())
           
           if 11 <= abs(daughtersIDarr[-1]) <= 16:
              leptonCount += 2*(daughtersIDarr[-1]>0)-1 

           chargeTotal+=PDG.GetParticle(daughtersIDarr[-1]).Charge()

        if Ispi0:
           daughtersIDarr.append(111)
        
        # print(sorted(daughtersIDarr))
        for idS in sorted(daughtersIDarr):
           daughtersID+=str(PDG.GetParticle(idS).GetName()) + " "

        if daughtersID not in countDict[genDaughtersID]:
           countDict[genDaughtersID][daughtersID] = 0
        countDict[genDaughtersID][daughtersID]+=1

        if chargeTotal not in reconStats[current_file]["ChargeSum"]:
           reconStats[current_file]["ChargeSum"][chargeTotal]=0
        reconStats[current_file]["ChargeSum"][chargeTotal]+=1
        
        if leptonCount not in reconStats[current_file]["LeptonNumber"]:
           reconStats[current_file]["LeptonNumber"][leptonCount]=0
        reconStats[current_file]["LeptonNumber"][leptonCount]+=1
        
        if chargeTotal == 0 and abs(leptonCount-genLeptons)<=1:
           reconStats[current_file]["ValidReconstruct"]+=1

        event_file.write(daughtersID)
  event_file.write("\n")
      
        

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

event_file_name = "event_decay_list.dat"
event_file = open(event_file_name, "w")

event_file.write("File, ")
event_file.write("Event N, ")
event_file.write("Mother Decay, ")
event_file.write("Recon Decay ")
event_file.write("\n")

current_file = ""

reconStats = dict()

for subdir, dirs, files in os.walk(rootdir):
  for file in files:
    if file.startswith("ship.conical") and file.endswith("_rec.root"):
      inputFile =  os.path.join(subdir, file)

      f = ROOT.TFile(inputFile)
      if not hasattr(f,'cbmsim'): continue
      sTree = f.cbmsim

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

      print(" INPUT FILE: ", inputFile, sTree.GetEntries())
      nEvents = 9999999
      nEvents = min(sTree.GetEntries(),nEvents)
      current_file = str(inputFile)
      reconStats[current_file] = dict() 
      reconStats[current_file]["nEvents"] = nEvents
      reconStats[current_file]["ChargeSum"] = dict()
      reconStats[current_file]["LeptonNumber"] = dict()
      reconStats[current_file]["ValidReconstruct"] = 0
      #for n in range(10):
      
      for n in range(nEvents):
         myEventLoop(n)

event_file.close()
	    

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
      
# for fileName in reconStats:
#    print("\nIn ", fileName)
#    for chargeVal in reconStats[fileName]["ChargeSum"]:
      
#       if chargeVal % 3 !=0:
#          chStr = str(int(chargeVal))+"/3"
#       else:
#          chStr = str(int(chargeVal/3))
#       print("\t charge {0:>2}: {1:>4}".format(chStr, reconStats[fileName]["ChargeSum"][chargeVal]) )


valid_recons = open("valid_reconstructions_stats.dat", "w")

for filename in reconStats:
   percentage = 100*(float(reconStats[filename]["ValidReconstruct"])/(reconStats[filename]["nEvents"]*1.0))
   print("\n  In ", filename)
   print("\t Valid reconstructions: {0:>4}%".format(percentage))
   valid_recons.write("{0} {1}\n".format(filename, round(percentage, 4)))
valid_recons.close()

# # do table
# for filename in reconStats: 
#    print("\n  In ", filename)
#    if 0 in reconStats[filename]["ChargeSum"]: 
#       print("\t Percentage of reconstructed events with charge 0: " + str(round(100*(reconStats[filename]["ChargeSum"][0]/reconStats[filename]["nEvents"]), 3)))
# print("\n")



