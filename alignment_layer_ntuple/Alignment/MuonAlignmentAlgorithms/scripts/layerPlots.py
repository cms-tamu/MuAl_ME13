from signConventions import *
import sys, os, time, math
from array import array

# List of conventions used here (many will parallel pre-existing conventions)
#   - endcap is 1 for + endcap for CSC, 2 for - endcap for CSC, 0 for all DT
#   - ringwheel is ring number (if CSC) or wheel number (DT)
#   - isCSC is boolean (true if describing CSC, of course)
#   - chambers are uniquely specified via a tuple of the form:
#       - if CSC, ("CSC", endcap#, station#, ring#, chamber#)
#       - if DT, ("DT", 0, station#, wheel#, chamber#)

if len(sys.argv) < 2:
    print "Filename must be first parameter"
    sys.exit()

filename = sys.argv[1]

print ">>> Making layer plots from", filename

sys.argv.append('-b') # hack to start it in batch mode
import ROOT as r

nSigma = 1.5
fh = r.TFile(filename,"read")
tt = fh.Get("csc_layer_ttree")
# tt.Draw("(2*(charge==1)-1)/pt","endcap==1 && station==1 && ring ==3 && chamber==17 && nlayers > 5 && select")


def chamberRadius(type, station, ringwheel):
    if(type == "CSC"):
        return signConventions[("CSC", 1, station, ringwheel, 1)][3]
    else:
        return signConventions[("DT", station, ringwheel, 1)][3]

def plotRanges(station, ring):
    if(station is 1):
        if(ring is 1): return 30, 80
        if(ring is 2): return 60, 100
        if(ring is 3): return 60, 100
    if(station is 2):
        if(ring is 1): return 80, 120
        if(ring is 2): return 80, 180
    if(station is 3):
        if(ring is 1): return 80, 150
        if(ring is 2): return 70, 110
    if(station is 4):
        if(ring is 1): return 80, 150
        if(ring is 2): return 80, 150
    
    return

def chamberPrettyString(type, endcap, station, ringwheel, chamber):
    if(type == "CSC"):
        prefix = "+" if endcap is 1 else "-"
        return "ME%s%d/%d/%d" % (prefix, station, ringwheel, chamber)
    else:
        prefix = "+" if station > 0 else "-"
        return "MB%s%d/%d/%d" % (prefix, abs(station), ringwheel, chamber)

def chamberToDirectoryPath(type, endcap, station, ringwheel, chamber):
    dir = "layer_plots/"
    if(type == "CSC"):
        prefix = "+" if endcap is 1 else "-"
        return "%sME%s/%d/%d/%d/" % (dir,prefix, station, ringwheel, chamber)
    else:
        prefix = "+" if station > 0 else "-"
        return "%sMB/MB%s%d/%d/%d/" % (dir,prefix, abs(station), ringwheel, chamber)

def binData(indexPt, Nbins, minPt):
    d = {}

    #minPt -= 0.01

    binWidth = 1.0/abs(minPt)/Nbins
    #binWidth = 0.0019226
    #print binWidth, minPt
    
    for i in range(Nbins+1):
        d[i] = []
        
    for elem in indexPt:
        idx = elem[0]
        pt = elem[1]
        binIdx = int(math.floor(1./abs(pt)/binWidth));
        # print binIdx, pt
        d[binIdx].append(idx)
    
    return d
    
def getGoodMuonIndices(posIndexPt, negIndexPt, Nbins, minPt):
    dP = binData(posIndexPt, Nbins, minPt)
    dN = binData(negIndexPt, Nbins, minPt)
    indicesToConsider = []
    for bin in range(Nbins):
        minMuonsInBin = min(len(dP[bin]), len(dN[bin]))
        indicesToConsider.extend(dP[bin][:minMuonsInBin])
        indicesToConsider.extend(dN[bin][:minMuonsInBin])
        
        # print bin, len(dP[bin]), len(dN[bin]), minMuonsInBin #, d[bin]
    return indicesToConsider
    
            
def equalizeCharges(chamberKey):
    posIndexPt, negIndexPt, minPt = [], [], 999.0
    nPosMu, nNegMu = 0,0
    for e in dChambers[chamberKey]:
        index, charge, pt = e
        if(charge is 1):
            posIndexPt.append([index,pt])
            nPosMu += 1
        else:
            negIndexPt.append([index,pt])
            nNegMu += 1
        if(pt < minPt): minPt = pt

    goodIndices = getGoodMuonIndices(posIndexPt, negIndexPt, 17, minPt)
    prettyStr = chamberPrettyString(*chamberKey)
    # print ">>> [%s] Found %i positive muons and %i negative muons" % (prettyStr, nPosMu, nNegMu)
    nMuonsToConsider = len(goodIndices) 
    print ">>> [%s] %i pos muons and %i neg muons => %d muons after equalization" % (prettyStr, nPosMu, nNegMu, nMuonsToConsider)
    # print ">>> [%s] Considering %d total muons after {q/pT,q/pz} equalization" % (prettyStr, nMuonsToConsider)
    
    # overwrite indices with the good indices we just found
    dChambers[chamberKey] = goodIndices
    return nMuonsToConsider
    
def fitCut(hist, sigmas, opts):
    lower, upper = hist.GetMean()-sigmas*hist.GetRMS(), hist.GetMean()+sigmas*hist.GetRMS()
    hist.Fit("gaus",opts, "", lower, upper)

def getFitParamsGauss(hist):
    fit = hist.GetFunction("gaus")
    
    p0 = fit.GetParameter(0) # const
    p0e = fit.GetParError(0) # const
    
    p1 = fit.GetParameter(1) # mean
    p1e = fit.GetParError(1) # mean error
    
    return p0, p0e, p1, p1e

def makePlotsCSC(chamberKey):
    h2D_cnt_actual = []
    h2D_cnt_tracks = []
    h2D_res_x_actual = []
    h2D_res_x_tracks = []
    h1D_res_x_rproj = []
    h1D_res_x_rphiproj = []
    h1D_res_x = []

    h1D_trans_layers = r.TH1F("h1D_trans_layers", "x offset vs layer;layer number;x offset (microns)",  6,0.5,6.5  )
    
    nXbins, nYbins = 50, 50
    XMax, YMax = 60., 100.
    station = chamberKey[2]
    ring = chamberKey[3]
    prettyStr = chamberPrettyString(*chamberKey)
    try:
        XMax, YMax = plotRanges(station, ring)
        print ">>> [%s] Detected ME%d/%d, so using XMax = %d and YMax = %d" % (prettyStr, station, ring, XMax, YMax)
    except: pass

    for i in range(6):
        laypfx = "L" + str(i+1) + " "
        # print laypfx
        h2D_cnt_actual.append(  r.TH2F(laypfx + "h2D_cnt_actual", laypfx + "actual hit locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
        h2D_cnt_tracks.append(  r.TH2F(laypfx + "h2D_cnt_tracks", laypfx + "track hit locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

        h2D_res_x_actual.append(  r.TProfile2D(laypfx + "h2D_res_x_actual", laypfx + "avg (per bin) x residuals at actual hit positions;local x (cm);local y (cm);cm",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
        h2D_res_x_tracks.append(  r.TProfile2D(laypfx + "h2D_res_x_tracks", laypfx + "avg (per bin) x residuals at track positions;local x (cm);local y (cm);cm",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

        # h2D_pull_tracks.append(  r.TH2F(laypfx + "h2D_pull_tracks", laypfx + "sum of x residuals (pull)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
        # h2D_statsig_tracks.append(  r.TH2F(laypfx + "h2D_statsig_tracks", laypfx + "avg res_x / sigma (significance)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

        h1D_res_x.append(  r.TH1F(laypfx + "h1D_res_x", laypfx + "x residuals;cm;counts",  100,-8.0,8.0)  )

        h1D_res_x_rproj.append(  r.TProfile(laypfx + "h1D_res_x_rproj", laypfx + "x residual vs r;cm;cm",  nYbins,-YMax,YMax)  )
        # h1D_res_x_xproj.append(  r.TProfile(laypfx + "h1D_res_x_xproj", laypfx + "x residual vs x",  nXbins,-100.0,100.0)  )
        h1D_res_x_rphiproj.append(  r.TProfile(laypfx + "h1D_res_x_rphiproj", laypfx + "x residual vs rphi;cm;cm",  nXbins,-YMax,YMax)  )

        # there are 2.530335 cm between each of the 48 wgs, so to get a binsize of 2.530335, we need (80-(-80))/2.530335=63.2=64 bins
        # h1D_res_rphi_yproj.append(  r.TProfile(laypfx + "h1D_res_rphi_yproj", laypfx + "rphi residual vs local y",  64,-80.0,80.0)  )
    
    print ">>> [%s] Booked histos" % (prettyStr)
    
    # convert the list of good indices into a dictionary for speedup
    # since we'll be checking it every time while looping through the ttree
    goodMuons = { }
    for idx in dChambers[chamberKey]:
        goodMuons[idx] = 1
        
    
    pinRadius = chamberRadius("CSC", chamberKey[2], chamberKey[3])
    # print pinRadius

    
    for idx, muon in enumerate(tt):
        endcap, station, ring, chamber = ord(muon.endcap), ord(muon.station), ord(muon.ring), ord(muon.chamber)

            
        if(idx not in goodMuons): continue
        
        for i in range(6):
            actual_x, actual_y = muon.hit_x[i], muon.hit_y[i]
            res_x, res_y = muon.res_x[i], muon.res_y[i]
            track_x = actual_x + res_x
            track_y = actual_y + res_y
            angle = math.atan(1.0 * actual_x / (pinRadius + actual_y))
            angle_track = math.atan(1.0 * track_x / (pinRadius + track_y))
            angle_ytrack = math.atan(1.0 * actual_x / (pinRadius + track_y))
            rphi_track = (pinRadius)*math.atan(track_x / (pinRadius + actual_y))
            rphi_ytrack = (pinRadius)*math.atan(track_x / (pinRadius + track_y))
            
            

            h1D_res_x[i].Fill(res_x)
            
            h2D_cnt_actual[i].Fill(actual_x, actual_y)
            h2D_cnt_tracks[i].Fill(track_x, track_y)

            h2D_res_x_actual[i].Fill(actual_x, actual_y, res_x)
            h2D_res_x_tracks[i].Fill(track_x, track_y, res_x)

            
            h1D_res_x_rproj[i].Fill(actual_y, res_x)
            h1D_res_x_rphiproj[i].Fill(rphi_track, res_x)
            
    print ">>> [%s] Filled histos" % (prettyStr)
    
    c1 = r.TCanvas("Canvas1", "Alignment Visualizations")
    c1.SetRightMargin(0.32);

    c1.SetGridx()
    c1.SetGridy()
    
    prefix = chamberToDirectoryPath(*chamberKey)
    os.system("mkdir -p %s" % (prefix))
    os.system("cp -pv" + " indexbase.php " + prefix+"_index.php") #
   
    layerTranslation = []
    for i in range(6):

        r.gStyle.SetPalette(1)
        r.gStyle.SetOptFit(0)

        suffix = "_LAY" + str(i+1) + ".png"     

        h2D_cnt_actual[i].Draw("colz")
        c1.SaveAs(prefix + "h2D_cnt_actual" + suffix)

        h2D_cnt_tracks[i].Draw("colz")
        c1.SaveAs(prefix + "h2D_cnt_tracks" + suffix)


        levels = [0.00, 0.35, 0.58, 0.78, 1.00]
        red    = [0.00, 0.00, 0.87, 1.00, 0.51]
        green  = [0.00, 0.00, 0.00, 0.00, 0.00]
        blue   = [0.51, 1.00, 0.12, 0.00, 0.00]
        levels = array('d', levels)
        ncontours = 15

        r.TColor.CreateGradientColorTable(len(levels), levels, array('d', red), array('d', green), array('d', blue), ncontours)
        r.gStyle.SetNumberContours(ncontours)
        h2D_res_x_actual[i].GetZaxis().SetRangeUser(-3.5, 3.5)
        h2D_res_x_actual[i].Draw("colz")
        c1.SaveAs(prefix + "h2D_res_x_actual" + suffix)

        h2D_res_x_tracks[i].GetZaxis().SetRangeUser(-3.5, 3.5)
        h2D_res_x_tracks[i].Draw("colz")
        c1.SaveAs(prefix + "h2D_res_x_tracks" + suffix)

        r.gStyle.SetOptFit(1) # display fitting parameters

        h1D_res_x_rproj[i].GetYaxis().SetRangeUser(-1.4, 1.4)
        h1D_res_x_rproj[i].Draw("E0")
        h1D_res_x_rproj[i].Fit("pol1","QC")
        c1.SaveAs(prefix + "h1D_res_x_rproj" + suffix)

       
        h1D_res_x_rphiproj[i].GetYaxis().SetRangeUser(-2.0, 2.0)
        h1D_res_x_rphiproj[i].Draw("E0")
        h1D_res_x_rphiproj[i].Fit("pol1","QC")
        c1.SaveAs(prefix + "h1D_res_x_rphiproj" + suffix)
        
        h1D_res_x[i].Draw()
        fitCut(h1D_res_x[i], nSigma, "QC")
        c1.SaveAs(prefix + "h1D_res_x" + suffix)

        fitparamsGauss = getFitParamsGauss(h1D_res_x[i])
        layerTranslation.append([fitparamsGauss[2], fitparamsGauss[3]]) #mean,meanerror


    for i, val in enumerate(layerTranslation):
        # convert cm to microns
        h1D_trans_layers.SetBinContent(i+1, 1.0e4*val[0])
        h1D_trans_layers.SetBinError(i+1, 1.0e4*val[1])

    h1D_trans_layers.Fit("pol1","QC")
    h1D_trans_layers.Draw("E0")
    c1.SaveAs(prefix + "h1D_trans_layers" + ".png")


    r.gStyle.SetPalette(1)




# mu = muon()

# nlayers = 0
# tt.SetBranchAddress("*",r.AddressOf(mu,"nlayers"))
# muon = tt.GetEntry(300)
# print nlayers

dChambers = { }

count = 0

# ls = []
# residuals = []
# avg, sig = -0.00645, 3.3555
print ">>> Looping through all elements in tree and collecting indices for charge equalization"
for i,muon in enumerate(tt):
    # if(count > 5000): break
    # ls.append(muon)
    endcap, station, ring, chamber = ord(muon.endcap), ord(muon.station), ord(muon.ring), ord(muon.chamber)
    isCSC = True
    # change ring to ringwheel in code
    chamberKey = ("CSC" if isCSC else "DT", 0 if not isCSC else endcap, station, ring, chamber)
    
    if(not muon.select): continue
    if(muon.nlayers < 6): continue
    
    
    
    if(muon.charge == 1): charge = 1
    else: charge = -1
    
    count += 1
    
    if(chamberKey in dChambers):
        # pz if CSC, but pt if DT
        dChambers[chamberKey].append([i,charge,muon.pz])
            
    else:
        dChambers[chamberKey] = [[i,charge,muon.pz]]
        
print ">>> Done looping through all elements in tree"

for chamber in dChambers.keys():
    nMuons = equalizeCharges(chamber)
    if(nMuons > 1000):
        t1 = time.time()
        makePlotsCSC(chamber)
        print ">>> [%s] Took %.2f seconds" % (chamberPrettyString(*chamberKey), time.time() - t1)
    else:
        print ">>> [%s] Low statistics, so skipping plots" % (chamberPrettyString(*chamberKey))

