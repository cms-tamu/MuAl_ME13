import ROOT as r
import struct, math, os, sys
from array import array

selectedCSC = (1, 2, 2, 1)  # endcap (p-1 m-0), station, wheel, chamber
specialSuffix = ""

# suppress "Info in..." messages, as there will be a lot
r.gErrorIgnoreLevel = r.kWarning


nSigma = 1.5 # number of sigma to do gaussian fit with
ME13pinR = 595.1500244141 # radius of ME13 ring alignment pin
isMC = False


print ">>> Args:", sys.argv
if(len(sys.argv) > 1):
    filename = sys.argv[1]
    cscStr = sys.argv[1].split("ME")[1]
    cscInfo = cscStr.split("_")
    if(cscInfo[0] == "p"): cscInfo[0] = 1
    else: cscInfo[0] = 0
    selectedCSC = ( cscInfo[0], int(cscInfo[1]), int(cscInfo[2]), int(cscInfo[3]) )
    if(cscInfo[4] == "MC"):
        specialSuffix = cscInfo[5]
        isMC = True
    else:
        specialSuffix = cscInfo[4]

    print ">>> Found argument with filename:", filename
    print ">>> Setting selectedCSC to", selectedCSC
else:
    print ">>> Using default file:", filename
    print ">>> Setting selectedCSC to default:", selectedCSC

if(specialSuffix is not ""): print ">>> Special suffix is",specialSuffix
    
prefix = "ME%s_%d_%d_%d" % ("p" if selectedCSC[0] == 1 else "m", selectedCSC[1], selectedCSC[2], selectedCSC[3])
if(isMC): prefix += "_MC"
prefix += "_%s_plots/" % (specialSuffix)
print ">>> Output into folder",prefix


def rProjection(localr, angle):
    return (localr + ME13pinR) / math.cos(angle) - ME13pinR


def fitCut(hist, sigmas, opts):
    lower, upper = hist.GetMean()-sigmas*hist.GetRMS(), hist.GetMean()+sigmas*hist.GetRMS()
    hist.Fit("gaus",opts, "", lower, upper)

def getFitParams(hist):
    fit = hist.GetFunction("pol1")
    
    p0 = fit.GetParameter(0) # offset from x axis
    p0e = fit.GetParError(0) # offset from x axis
    
    p1 = fit.GetParameter(1) # slope
    p1e = fit.GetParError(1) # slope
    
    return p0, p0e, p1, p1e

def evalInCenter(hist):
    fit = hist.GetFunction("pol1")
    p0 = fit.GetParameter(0) # offset from x axis
    p0e = fit.GetParError(0) # offset from x axis
    p1 = fit.GetParameter(1) # slope
    p1e = fit.GetParError(1) # slope

    # evaluate at layer 3.5
    val = p0+3.5*p1
    err = math.sqrt(p0e**2 + 3.5**2 * p1e**2)
    
    return val, err


def getFitParamsGauss(hist):
    fit = hist.GetFunction("gaus")
    
    p0 = fit.GetParameter(0) # const
    p0e = fit.GetParError(0) # const
    
    p1 = fit.GetParameter(1) # mean
    p1e = fit.GetParError(1) # mean error
    
    #print p0, p0e, p1, p1e
    return p0, p0e, p1, p1e

def fixFlooredBins(hist,minZ=-3.5):
    for x in range(hist.GetNbinsX()+1):
        for y in range(hist.GetNbinsY()+1):
            mu = hist.GetBinContent(x,y)
            if(mu <= minZ):
                hist.SetBinContent(x,y,minZ)

def zoomXrange(h):
    nonEmptyBins = [(e,h.GetXaxis().GetBinCenter(e)) for e in range(h.GetNbinsX()) if h.GetBinContent(e)>0] 
    minX, maxX = min(nonEmptyBins)[1], max(nonEmptyBins)[1]
    tolerance = (maxX-minX)*0.15
    h.GetXaxis().SetRangeUser(minX-tolerance,maxX+tolerance)

def binData(indexPt, Nbins, minPt):
    d = {}

    binWidth = 1.0/abs(minPt)/Nbins
    
    for i in range(Nbins+1):
        d[i] = []
        
    for elem in indexPt:
        idx = elem[0]
        pt = elem[1]
        binIdx = int(math.floor(1./abs(pt)/binWidth));
        #print binIdx, pt
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
        
        #print bin, len(dP[bin]), len(dN[bin]), minMuonsInBin #, d[bin]
    return indicesToConsider

def applyRotation(x,y,phiz):
    xp = math.cos(phiz)*x - math.sin(phiz)*y
    yp = math.sin(phiz)*x + math.cos(phiz)*y
    return xp,yp

c1 = r.TCanvas("Canvas1", "Alignment Visualizations")

typePrefix = "ME%s%d/%d/%d" % ("+" if selectedCSC[0] == 1 else "-", selectedCSC[1], selectedCSC[2], selectedCSC[3])
if(isMC):
    typePrefix += " #font[2]{#color[4]{MC}} "
else:
    typePrefix += " #font[2]{#color[2]{DATA}} "

nXbins, nYbins = 50, 50
XMax, YMax = 60., 100.
xyRatio = YMax/XMax 

resRange = 8.0
NUMLAYERS = 6

h2D_cnt_actual = []
h2D_cnt_tracks = []
h2D_res_x_actual = []
h2D_res_x_tracks = []
h2D_pull_tracks = []
h1D_res_x_rproj = []
#h1D_res_x_xproj = []
h1D_res_x_rphiproj = []
h1D_res_rphi_yproj = []
h1D_res_x = []
h2D_statsig_tracks = []

#h1D_res_x_avg = r.TH1F("h1D_res_x_avg", "avg x residuals", 100,-8.0,8.0) 
h1D_actual_localy = r.TH1F("h1D_actual_localy", typePrefix+"distribution of actual y hits (on L3);cm;counts",  48*2,-90,90)
h1D_tracks_localy = r.TH1F("h1D_tracks_localy", typePrefix+"distribution of track y positions (on L3);cm;counts",  nYbins,-90,90)
h1D_actual_angle = r.TH1F("h1D_actual_angle", typePrefix+"angular distribution of hits (on L3);phi;counts",  nYbins,-0.08,0.08)
h1D_tracks_angle = r.TH1F("h1D_tracks_angle", typePrefix+"angular distribution of tracks (on L3);phi;counts",  nYbins,-0.08,0.08)
h1D_rot_dxdr_layers = r.TH1F("h1D_rot_dxdr_layers", typePrefix+"phiz (dx/dr) vs layer;layer;dx/dr (urad)",  6,0.5,6.5  )
h1D_trans_layers = r.TH1F("h1D_trans_layers", typePrefix+"x offset vs layer;layer; x offset (microns)",  6,0.5,6.5  )
h2D_nlayers_hit = r.TProfile2D("h2D_nlayers_hit", typePrefix+"num layers hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
h2D_nDT_hit = r.TProfile2D("h2D_nDT_hit", typePrefix+"num DTs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
h2D_nCSC_hit = r.TProfile2D("h2D_nCSC_hit", typePrefix+"num CSCs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
h2D_nDTCSC_hit = r.TProfile2D("h2D_nDTCSC_hit", typePrefix+"num DTs+CSCs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
h2D_nTracker_hit = r.TProfile2D("h2D_nTracker_hit", typePrefix+"num tracker hits;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
h1D_pt = r.TH1F("h1D_pt", typePrefix+"p_{T};p_{T} (GeV/c);counts",  100, 20, 210)
h1D_pz = r.TH1F("h1D_pz", typePrefix+"p_{z};p_{z} (GeV/c);counts",  100, 20, 210)
h1D_p = r.TH1F("h1D_p", typePrefix+"p;p (GeV/c);counts",  100, 20, 210) 
h1D_eta = r.TH1F("h1D_eta", typePrefix+"#eta;#eta;counts", 800, -2.2, 2.2) 
h1D_phi = r.TH1F("h1D_phi", typePrefix+"#phi;#phi;counts",  800, -math.pi, math.pi)


for i in range(NUMLAYERS):
    laypfx = " #font[2]{L" + str(i+1) + "} "
    print ">>> Booking histos for layer",i+1
    h2D_cnt_actual.append(  r.TH2F(laypfx + "h2D_cnt_actual", typePrefix+laypfx + "actual hit locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
    h2D_cnt_tracks.append(  r.TH2F(laypfx + "h2D_cnt_tracks", typePrefix+laypfx + "track locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )

    h2D_res_x_actual.append(  r.TProfile2D(laypfx + "h2D_res_x_actual", typePrefix+laypfx + "avg x res at actual hit positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
    h2D_res_x_tracks.append(  r.TProfile2D(laypfx + "h2D_res_x_tracks", typePrefix+laypfx + "avg x res at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )

    h2D_pull_tracks.append(  r.TProfile2D(laypfx + "h2D_pull_tracks", typePrefix+laypfx + "pull at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
    h2D_statsig_tracks.append(  r.TH2F(laypfx + "h2D_statsig_tracks", typePrefix+laypfx + "significance at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )

    h1D_res_x.append(  r.TH1F(laypfx + "h1D_res_x", typePrefix+laypfx + "x residuals;x residual (cm);counts",  100,-resRange,resRange)  )

    h1D_res_x_rproj.append(  r.TProfile(laypfx + "h1D_res_x_rproj", typePrefix+laypfx + "x residual vs r;r (cm);x residual (cm)",  26,-80.0,80.0)  )
    #h1D_res_x_xproj.append(  r.TProfile(laypfx + "h1D_res_x_xproj", laypfx + "x residual vs x",  nXbins,-100.0,100.0)  )
    # rphiproj has 50 bins, but rproj has 26 bins. This is so that the number of bins for PHYSICAL regions is the same
    h1D_res_x_rphiproj.append(  r.TProfile(laypfx + "h1D_res_x_rphiproj", typePrefix+laypfx + "x residual vs scaled-phi;scaled-phi (cm);x residual (cm)",  50,-80.0,80.0)  )

    # there are 2.530335 cm between each of the 48 wgs, so to get a binsize of 2.530335, we need (80-(-80))/2.530335=63.2=64 bins
    h1D_res_rphi_yproj.append(  r.TProfile(laypfx + "h1D_res_rphi_yproj", typePrefix+laypfx + "rphi residual vs local y;local y (cm);rphi residual (cm)",  64,-80.0,80.0)  )

    #h1D_res_x[i].StatOverflows(True)
    #h1D_res_x_rproj[i].StatOverflows(True)
    #h1D_res_x_rphiproj[i].StatOverflows(True)

fh = r.TFile(filename)
tt = r.gDirectory.Get("csc_layer_ttree")

nPosMu = 0
nNegMu = 0
nMuons = 0
posIndexPt = []
negIndexPt = []
minPt = 9999.0


for i,muon in enumerate(tt):
    nMuons += 1
    pt = muon.pz
    if(not muon.select or muon.nlayers < 6): continue
    endcap, station, ring, chamber = ord(muon.endcap), ord(muon.station), ord(muon.ring), ord(muon.chamber)
    if (not (selectedCSC == (endcap, station, ring, chamber))): continue

            
    if muon.charge == 1:
        nPosMu += 1
        posIndexPt.append([i,pt])
    else:
        nNegMu += 1
        negIndexPt.append([i,pt])
    if(pt < minPt): minPt = pt

# print nMuons
goodIndices = getGoodMuonIndices(posIndexPt, negIndexPt, 17, minPt)
print ">>> Found %i positive muons and %i negative muons" % (nPosMu, nNegMu)
nMuonsToConsider = len(goodIndices) 
print ">>> Considering ",nMuonsToConsider,"total muons after {q/pT,q/pz} equalization"

goodMuons = {}
for i in range(nMuons): goodMuons[i] = 0
for e in goodIndices: goodMuons[e] = 1

count = 0
for idx, muon in enumerate(tt):
    count += 1
    if(count % 2000 == 0): print ">>>",count 

    #if(count > 5000): break
    if(not muon.select or muon.nlayers < 6): continue
    
    endcap, station, ring, chamber = ord(muon.endcap), ord(muon.station), ord(muon.ring), ord(muon.chamber)

    if(selectedCSC == (endcap, station, ring, chamber)):

        if(goodMuons[idx] != 1): continue

        h1D_pt.Fill(muon.pt)
        h1D_pz.Fill(muon.pz)
        h1D_p.Fill(math.sqrt(muon.pt*muon.pt+muon.pz*muon.pz))
        h1D_eta.Fill(muon.eta)
        phi = muon.phi
        #if(muon.phi < 0): phi += 2*3.14159265
        h1D_phi.Fill(phi)

        for i in range(NUMLAYERS):
            try:
                actual_x, actual_y = muon.hit_x[i], muon.hit_y[i]
                res_x, res_y = muon.res_x[i], muon.res_y[i]
            except:
                layStr = "muon.lay%i_" % (i+1)
                actual_x, actual_y = eval(layStr+"x"), eval(layStr+"y")
                res_x, res_y = eval(layStr+"res_x"), eval(layStr+"res_y")

            # roughly discard tails
            if(abs(res_x) > resRange): continue 

            # propagated track = actual + residual
            track_x = actual_x + res_x
            track_y = actual_y + res_y

            # do translations or rotations

            #actual_x, actual_y = 0.99*actual_x, 0.99*actual_y # make chamber 1% smaller wrt to x,y=0,0
            #actual_x, actual_y = actual_x, actual_y+0.30 # shift y up by 3 mm
            #actual_x, actual_y = applyRotation(actual_x, actual_y, 0.001) # rotate CCW by 1 mrad
            #actual_x, actual_y = actual_x + 0.03, actual_y # shift x up by 300 microns
            # recalculate residuals
            #res_x, res_y = track_x - actual_x, track_y - actual_y

            ######

            angle = math.atan(1.0 * actual_x / (ME13pinR + actual_y))
            angle_track = math.atan(1.0 * track_x / (ME13pinR + track_y))
            rphi_track = (ME13pinR)*math.atan(track_x / (ME13pinR + track_y))
            res_rphi = math.cos(angle)*res_x + math.sin(angle)*res_y
        
            if(i==2): # if layer 3 ("center" of chamber)
                h1D_actual_localy.Fill(actual_y)
                h1D_tracks_localy.Fill(track_y)
            
                h1D_actual_angle.Fill(angle)
                h1D_tracks_angle.Fill(angle_track)

                h2D_nlayers_hit.Fill(actual_x, actual_y, muon.nlayers)
                h2D_nDT_hit.Fill(actual_x, actual_y, muon.nDT)
                h2D_nCSC_hit.Fill(actual_x, actual_y, muon.nCSC)
                h2D_nDTCSC_hit.Fill(actual_x, actual_y, muon.nCSC+muon.nDT)
                h2D_nTracker_hit.Fill(actual_x, actual_y, muon.nTracker)


            h1D_res_x[i].Fill(res_x)
            
            h2D_cnt_actual[i].Fill(actual_x, actual_y)
            h2D_cnt_tracks[i].Fill(track_x, track_y)

            h2D_res_x_actual[i].Fill(actual_x, actual_y, res_x)
            h2D_res_x_tracks[i].Fill(track_x, track_y, res_x)

            h2D_pull_tracks[i].Fill(track_x, track_y, res_x)

            h1D_res_x_rproj[i].Fill(rProjection(actual_y,angle), res_x)
            h1D_res_x_rphiproj[i].Fill(rphi_track, res_x)
            h1D_res_rphi_yproj[i].Fill(actual_y, res_rphi)
       

c1.SetRightMargin(0.32);
c1.SetGridx()
c1.SetGridy()
c1.SetCanvasSize(696,472)


os.system("mkdir -p " + prefix)
os.system("cp -p " + " indexbase.php " + prefix+"_index.php") #

r.gStyle.SetOptStat("rme")

layerRotationR = []
layerTranslation = []

for i in range(NUMLAYERS):

    r.gStyle.SetPalette(1)
    r.gStyle.SetOptFit(0)

    suffix = "_LAY" + str(i+1) + ".png"
    
    if(i == 0):
        h1D_pt.Draw()
        c1.SaveAs(prefix + "h1D_pt" + suffix)
        
        h1D_pz.Draw()
        c1.SaveAs(prefix + "h1D_pz" + suffix)

        h1D_p.Draw()
        c1.SaveAs(prefix + "h1D_p" + suffix)

        h1D_phi.Draw()
        zoomXrange(h1D_phi)
        c1.SaveAs(prefix + "h1D_phi" + suffix)

        zoomXrange(h1D_eta)
        h1D_eta.Draw()
        c1.SaveAs(prefix + "h1D_eta" + suffix)

        h1D_actual_angle.Draw()
        c1.SaveAs(prefix + "h1D_angle_actual" + suffix)

        h1D_tracks_angle.Draw()
        c1.SaveAs(prefix + "h1D_angle_tracks" + suffix)


        h1D_actual_localy.Draw()
        c1.SaveAs(prefix + "h1D_localy_actual" + suffix)

        h1D_tracks_localy.Draw()
        c1.SaveAs(prefix + "h1D_localy_tracks" + suffix)

        
        c1.SetCanvasSize(555,672)

        h2D_nlayers_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nlayers_hit" + suffix)

        h2D_nDT_hit.GetZaxis().SetRangeUser(0, 5)
        h2D_nDT_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nDT_hit" + suffix)

        h2D_nCSC_hit.GetZaxis().SetRangeUser(0, 5)
        h2D_nCSC_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nCSC_hit" + suffix)

        h2D_nDTCSC_hit.GetZaxis().SetRangeUser(0, 5)
        h2D_nDTCSC_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nDTCSC_hit" + suffix)

        h2D_nTracker_hit.GetZaxis().SetRangeUser(0, 21)
        h2D_nTracker_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nTracker_hit" + suffix)
        
        c1.SetCanvasSize(696,472)


    c1.SetCanvasSize(555,672)

    for x in range(h2D_res_x_tracks[i].GetNbinsX()): 
        for y in range(h2D_res_x_tracks[i].GetNbinsY()): 
            mu = h2D_res_x_tracks[i].GetBinContent(x,y)
            sig = h2D_res_x_tracks[i].GetBinError(x,y)
            ibin = h2D_res_x_tracks[i].GetBin(x,y)
            entries = h2D_res_x_tracks[i].GetBinEntries(ibin)

            if(sig > 0.001 and entries > 2):
                h2D_statsig_tracks[i].SetBinContent(x,y,abs(mu)/sig)
            else:
                h2D_statsig_tracks[i].SetBinContent(x,y,0.0)

    h2D_statsig_tracks[i].GetZaxis().SetRangeUser(0, 7)
    h2D_statsig_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_statsig_tracks" + suffix)


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

    # multiply avg resx by hit occupancy to get "force" of pulling
    for x in range(h2D_pull_tracks[i].GetNbinsX()): 
        for y in range(h2D_pull_tracks[i].GetNbinsY()): 
            ibin = h2D_pull_tracks[i].GetBin(x,y)
            mu = h2D_pull_tracks[i].GetBinContent(ibin)
            entries = h2D_pull_tracks[i].GetBinEntries(ibin)
            # setbincontent sets the sum of the bin, so bin mean ends up 
            # being sum/entries. we want bin to be mu*entries, so we need
            # to set it to mu*entries*entries
            if(entries > 0): h2D_pull_tracks[i].SetBinContent(ibin,mu*entries*entries)
    h2D_pull_tracks[i].GetZaxis().SetRangeUser(-80.0,80.0)
    h2D_pull_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_pull_tracks" + suffix)

    h2D_res_x_actual[i].GetZaxis().SetRangeUser(-3.5, 3.5)
    fixFlooredBins(h2D_res_x_actual[i], minZ=-3.5)
    h2D_res_x_actual[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_actual" + suffix)

    h2D_res_x_tracks[i].GetZaxis().SetRangeUser(-3.5, 3.5)
    fixFlooredBins(h2D_res_x_tracks[i], minZ=-3.5)
    h2D_res_x_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_tracks" + suffix)

    c1.SetCanvasSize(696,472)


    r.gStyle.SetOptFit(1) # display fitting parameters

    h1D_res_x_rproj[i].GetYaxis().SetRangeUser(-2.0, 2.0)
    # h1D_res_x_rproj[i].BuildOptions(0,0,"")
    h1D_res_x_rproj[i].Draw("E0")
    h1D_res_x_rproj[i].Fit("pol1","QC")
    c1.SaveAs(prefix + "h1D_res_x_rproj" + suffix)


    fitparams = getFitParams(h1D_res_x_rproj[i])
    layerRotationR.append([fitparams[2], fitparams[3]]) #p1, p1error


    h1D_res_rphi_yproj[i].GetYaxis().SetRangeUser(-2.0, 2.0)
    #h1D_res_rphi_yproj[i].GetXaxis().SetRangeUser(-0.1, 0.1)
    h1D_res_rphi_yproj[i].Draw("E0")
    h1D_res_rphi_yproj[i].Fit("pol1","QC")
    c1.SaveAs(prefix + "h1D_res_rphi_yproj" + suffix)
   
    h1D_res_x_rphiproj[i].GetYaxis().SetRangeUser(-2.0, 2.0)
    # draw horizontal line at mean_y value for the tprofile
    f1 = r.TF1("f1",str(h1D_res_x_rphiproj[i].GetMean(2)),h1D_res_x_rphiproj[i].GetXaxis().GetXmin(),h1D_res_x_rphiproj[i].GetXaxis().GetXmax())
    leg = r.TLegend(0.70,0.15,0.95,0.30)
    leg.AddEntry(f1,"Mean y = %.5f" % h1D_res_x_rphiproj[i].GetMean(2),"l")
    h1D_res_x_rphiproj[i].Draw("E0")
    leg.Draw("same")
    f1.Draw("same")
    #h1D_res_x_rphiproj[i].Fit("pol0","QCWW")
    c1.SaveAs(prefix + "h1D_res_x_rphiproj" + suffix)
 
    h1D_res_x[i].Draw()
    fitCut(h1D_res_x[i], nSigma, "QC")
    c1.SaveAs(prefix + "h1D_res_x" + suffix)

    mu = h1D_res_x[i].GetMean()
    muerr = h1D_res_x[i].GetMeanError()
    print ">>> Unfitted mean x residual for layer %i: %f +/- %f" % (i, mu, muerr)
    print ">>> Saved plots for layer %d" % (i+1)

    fitparamsGauss = getFitParamsGauss(h1D_res_x[i])
    layerTranslation.append([fitparamsGauss[2], fitparamsGauss[3]]) #mean,meanerror


for i, val in enumerate(layerRotationR):
    h1D_rot_dxdr_layers.SetBinContent(i+1, 1.0e6*val[0])
    h1D_rot_dxdr_layers.SetBinError(i+1, 1.0e6*val[1])

for i, val in enumerate(layerTranslation):
    h1D_trans_layers.SetBinContent(i+1, 1.0e4*val[0])
    h1D_trans_layers.SetBinError(i+1, 1.0e4*val[1])


h1D_rot_dxdr_layers.GetYaxis().SetRangeUser(-1500, 1500)
h1D_rot_dxdr_layers.Fit("pol1","QC")
h1D_rot_dxdr_layers.Draw("E0")
centerValue = evalInCenter(h1D_rot_dxdr_layers)
b = centerValue
print ">>> Rotation on layer 3.5: %.1f +/- %.1f urad" % centerValue
leg = r.TLegend(0.70,0.15,0.98,0.30)
leg.AddEntry(h1D_rot_dxdr_layers.GetFunction("pol1"),"#scale[2.0]{L3.5 fit = %.0f #pm%.0f #murad}" % centerValue,"l")
leg.Draw("same")
c1.SaveAs(prefix + "h1D_rot_dxdr_layers" + ".png")

h1D_trans_layers.GetYaxis().SetRangeUser(-700,700)
h1D_trans_layers.Fit("pol1","QC")
h1D_trans_layers.Draw("E0")
centerValue = evalInCenter(h1D_trans_layers)
a = centerValue
print ">>> X offset on layer 3.5: %.1f +/- %.1f microns" % centerValue
leg = r.TLegend(0.70,0.15,0.98,0.30)
leg.AddEntry(h1D_trans_layers.GetFunction("pol1"),"#scale[2.0]{L3.5 fit = %.0f#pm%.0f #mum}" % centerValue,"l")
leg.Draw("same")
c1.SaveAs(prefix + "h1D_res_x_layers" + ".png")

# make it easier to extract various results/parameters from script running by grepping for DBG
# chamber#, specialsuffix, mc/data, nMuons, xoffset_layer3.5, xoffset_layer3.5error, phiz_layer3.5, phiz_layer3.5error
print "DBG %d %s %s %d %.6f %.6f %.6f %.6f" % (selectedCSC[-1], specialSuffix, "MC" if isMC else "DATA", nMuonsToConsider, a[0]*10**(-4), a[1]*10**(-4), b[0]*10**(-6), b[1]*10**(-6))

r.gStyle.SetPalette(1)
