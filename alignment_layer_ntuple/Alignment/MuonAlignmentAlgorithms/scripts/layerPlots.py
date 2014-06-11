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

if len(sys.argv) < 3:
    print "usage: python layerPlots.py [filename] [minNCrossedChambers]"
    sys.exit()

filename = sys.argv[1]
minChambers = int(sys.argv[2])

print ">>> Making layer plots from", filename

sys.argv.append('-b') # hack to start it in batch mode
import ROOT as r

# suppress "Info in..." messages, as there will be a lot
r.gErrorIgnoreLevel = r.kWarning

nSigma = 1.5
fh = r.TFile(filename,"read")
tt = fh.Get("csc_layer_ttree")


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

def rProjection(localr, angle, pinR):
    return (localr + pinR) / math.cos(angle) - pinR

def getFitParams(hist):
    fit = hist.GetFunction("pol1")
    
    p0 = fit.GetParameter(0) # offset from x axis
    p0e = fit.GetParError(0) # offset from x axis
    
    p1 = fit.GetParameter(1) # slope
    p1e = fit.GetParError(1) # slope
    
    return p0, p0e, p1, p1e


def getFitParamsGauss(hist):
    fit = hist.GetFunction("gaus")
    
    p0 = fit.GetParameter(0) # const
    p0e = fit.GetParError(0) # const
    
    p1 = fit.GetParameter(1) # mean
    p1e = fit.GetParError(1) # mean error
    
    return p0, p0e, p1, p1e

def fixFlooredBins(hist,minZ=-3.5):
    for x in range(hist.GetNbinsX()+1):
        for y in range(hist.GetNbinsY()+1):
            mu = hist.GetBinContent(x,y)
            if(mu <= minZ):
                hist.SetBinContent(x,y,minZ)


def makePlotsCSC(chamberKey):
    h2D_cnt_actual = []
    h2D_cnt_tracks = []
    h2D_res_x_actual = []
    h2D_res_x_tracks = []
    h2D_pull_tracks = []
    h1D_res_x_rproj = []
    h1D_res_x_rphiproj = []
    h1D_res_rphi_yproj = []
    h1D_res_x = []
    h2D_statsig_tracks = []


    nXbins, nYbins = 50, 50
    XMax, YMax = 60., 100.
    resRange = 8.0
    station = chamberKey[2]
    ring = chamberKey[3]
    prettyStr = chamberPrettyString(*chamberKey)
    try:
        XMax, YMax = plotRanges(station, ring)
        print ">>> [%s] Detected ME%d/%d, so using XMax = %d and YMax = %d" % (prettyStr, station, ring, XMax, YMax)
    except: pass

    typePrefix = "#font[2]{%s} " % prettyStr

    h1D_pt = r.TH1F("h1D_pt", typePrefix+"p_T;p_T (GeV/c);counts",  100, 20, 210)
    h1D_pz = r.TH1F("h1D_pz", typePrefix+"p_z;p_z (GeV/c);counts",  100, 20, 210)
    h1D_actual_localy = r.TH1F("h1D_actual_localy", typePrefix+"distribution of actual y hits (on L3);cm;counts",  2*nYbins,-YMax,YMax)
    h1D_tracks_localy = r.TH1F("h1D_tracks_localy", typePrefix+"distribution of track y positions (on L3);cm;counts",  nYbins,-YMax,YMax)
    h1D_actual_angle = r.TH1F("h1D_actual_angle", typePrefix+"angular distribution of hits (on L3);phi;counts",  nYbins,-0.15,0.15)
    h1D_tracks_angle = r.TH1F("h1D_tracks_angle", typePrefix+"angular distribution of tracks (on L3);phi;counts",  nYbins,-0.15,0.15)
    h1D_rot_dxdr_layers = r.TH1F("h1D_rot_dxdr_layers", typePrefix+"phiz (dx/dr) vs layer;layer;dx/dr (urad)",  6,0.5,6.5  )
    h1D_trans_layers = r.TH1F("h1D_trans_layers", typePrefix+"x offset vs layer;layer; x offset (microns)",  6,0.5,6.5  )
    h2D_nlayers_hit = r.TProfile2D("h2D_nlayers_hit", typePrefix+"num layers hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
    h2D_nDT_hit = r.TProfile2D("h2D_nDT_hit", typePrefix+"num DTs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
    h2D_nCSC_hit = r.TProfile2D("h2D_nCSC_hit", typePrefix+"num CSCs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
    h2D_nTracker_hit = r.TProfile2D("h2D_nTracker_hit", typePrefix+"num tracker hits;actual hit x;actual hit y",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)

    for i in range(6):
        laypfx = " #font[2]{L" + str(i+1) + "} "
        h2D_cnt_actual.append(  r.TH2F(laypfx + "h2D_cnt_actual", typePrefix+laypfx + "actual hit locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
        h2D_cnt_tracks.append(  r.TH2F(laypfx + "h2D_cnt_tracks", typePrefix+laypfx + "track hit locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

        h2D_res_x_actual.append(  r.TProfile2D(laypfx + "h2D_res_x_actual", typePrefix+laypfx + "avg (per bin) x residuals at actual hit positions;local x (cm);local y (cm);cm",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
        h2D_res_x_tracks.append(  r.TProfile2D(laypfx + "h2D_res_x_tracks", typePrefix+laypfx + "avg (per bin) x residuals at track positions;local x (cm);local y (cm);cm",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

        h2D_pull_tracks.append(  r.TProfile2D(laypfx + "h2D_pull_tracks", typePrefix+laypfx + "sum of x residuals (pull) at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
        h2D_statsig_tracks.append(  r.TH2F(laypfx + "h2D_statsig_tracks", typePrefix+laypfx + "avg res_x / sigma (significance) at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

        h1D_res_x.append(  r.TH1F(laypfx + "h1D_res_x", typePrefix+laypfx + "x residuals;cm;counts",  100,-resRange,resRange)  )

        h1D_res_x_rproj.append(  r.TProfile(laypfx + "h1D_res_x_rproj", typePrefix+laypfx + "x residual vs r;r (cm);x residual (cm)",  nYbins,-YMax,YMax)  )
        h1D_res_x_rphiproj.append(  r.TProfile(laypfx + "h1D_res_x_rphiproj", typePrefix+laypfx + "x residual vs scaled-phi;scaled-phi (cm);x residual (cm)",  nXbins,-YMax,YMax)  )

        # there are 2.530335 cm between each of the 48 wgs, so to get a binsize of 2.530335, we need (80-(-80))/2.530335=63.2=64 bins
        h1D_res_rphi_yproj.append(  r.TProfile(laypfx + "h1D_res_rphi_yproj", typePrefix+laypfx + "rphi residual vs local y;local y (cm);rphi residual (cm)",  nYbins,-YMax,YMax)  )
    
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

        h1D_pt.Fill(muon.pt)
        h1D_pz.Fill(muon.pz)
        
        for i in range(6):

            try:
                actual_x, actual_y = muon.hit_x[i], muon.hit_y[i]
                res_x, res_y = muon.res_x[i], muon.res_y[i]
            except:
                layStr = "muon.lay%i_" % (i+1)
                actual_x, actual_y = eval(layStr+"x"), eval(layStr+"y")
                res_x, res_y = eval(layStr+"res_x"), eval(layStr+"res_y")

            track_x = actual_x + res_x
            track_y = actual_y + res_y
            angle = math.atan(1.0 * actual_x / (pinRadius + actual_y))
            angle_track = math.atan(1.0 * track_x / (pinRadius + track_y))
            rphi_track = (pinRadius)*math.atan(track_x / (pinRadius + track_y))
            res_rphi = math.cos(angle)*res_x + math.sin(angle)*res_y
            
            if(i==2): # if layer 3 ("center" of chamber)
                h1D_actual_localy.Fill(actual_y)
                h1D_tracks_localy.Fill(track_y)
            
                h1D_actual_angle.Fill(angle)
                h1D_tracks_angle.Fill(angle_track)

            h1D_res_x[i].Fill(res_x)
            
            h2D_cnt_actual[i].Fill(actual_x, actual_y)
            h2D_cnt_tracks[i].Fill(track_x, track_y)

            h2D_res_x_actual[i].Fill(actual_x, actual_y, res_x)
            h2D_res_x_tracks[i].Fill(track_x, track_y, res_x)

            h2D_pull_tracks[i].Fill(track_x, track_y, res_x)

            h2D_nlayers_hit.Fill(actual_x, actual_y, muon.nlayers)
            h2D_nDT_hit.Fill(actual_x, actual_y, muon.nDT)
            h2D_nCSC_hit.Fill(actual_x, actual_y, muon.nCSC)
            h2D_nTracker_hit.Fill(actual_x, actual_y, muon.nTracker)

            if(abs(res_x) <= resRange):
                h1D_res_x_rproj[i].Fill(rProjection(actual_y,angle,pinRadius), res_x)
                h1D_res_x_rphiproj[i].Fill(rphi_track, res_x)
                h1D_res_rphi_yproj[i].Fill(actual_y, res_rphi)
            
            
    print ">>> [%s] Filled histos" % (prettyStr)
    
    c1 = r.TCanvas("Canvas1", "Alignment Visualizations")
    c1.SetRightMargin(0.32);

    c1.SetGridx()
    c1.SetGridy()
    
    prefix = chamberToDirectoryPath(*chamberKey)
    os.system("mkdir -p %s" % (prefix))
    os.system("cp -pv" + " indexbase.php " + prefix+"_index.php") #
   
    r.gStyle.SetOptStat("rme")

    layerTranslation = []
    layerRotationR = []
    for i in range(6):

        r.gStyle.SetPalette(1)
        r.gStyle.SetOptFit(0)

        suffix = "_LAY" + str(i+1) + ".png"     

        if(i == 0):
            h1D_pt.Draw()
            c1.SaveAs(prefix + "h1D_pt" + suffix)
            
            h1D_pz.Draw()
            c1.SaveAs(prefix + "h1D_pz" + suffix)
            
            h1D_actual_angle.Draw()
            c1.SaveAs(prefix + "h1D_angle_actual" + suffix)

            h1D_tracks_angle.Draw()
            c1.SaveAs(prefix + "h1D_angle_tracks" + suffix)


            h1D_actual_localy.Draw()
            c1.SaveAs(prefix + "h1D_localy_actual" + suffix)

            h1D_tracks_localy.Draw()
            c1.SaveAs(prefix + "h1D_localy_tracks" + suffix)

            
            h2D_nlayers_hit.Draw("colz")
            c1.SaveAs(prefix + "h2D_nlayers_hit" + suffix)

            h2D_nDT_hit.GetZaxis().SetRangeUser(0, 4)
            h2D_nDT_hit.Draw("colz")
            c1.SaveAs(prefix + "h2D_nDT_hit" + suffix)

            h2D_nCSC_hit.GetZaxis().SetRangeUser(0, 4)
            h2D_nCSC_hit.Draw("colz")
            c1.SaveAs(prefix + "h2D_nCSC_hit" + suffix)

            h2D_nTracker_hit.GetZaxis().SetRangeUser(0, 21)
            h2D_nTracker_hit.Draw("colz")
            c1.SaveAs(prefix + "h2D_nTracker_hit" + suffix)



        for x in range(nXbins):
            for y in range(nYbins):
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
        for x in range(nXbins):
            for y in range(nYbins):
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

        r.gStyle.SetOptFit(1) # display fitting parameters

        h1D_res_x_rproj[i].GetYaxis().SetRangeUser(-2.0, 2.0)
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
        print ">>> [%s] Unfitted mean x residual for layer %i: %f +/- %f" % (prettyStr, i, mu, muerr)
        print ">>> [%s] Saved plots for layer %d" % (prettyStr, i+1)

        fitparamsGauss = getFitParamsGauss(h1D_res_x[i])
        layerTranslation.append([fitparamsGauss[2], fitparamsGauss[3]]) #mean,meanerror



    for i, val in enumerate(layerRotationR):
        h1D_rot_dxdr_layers.SetBinContent(i+1, 1.0e6*val[0])
        h1D_rot_dxdr_layers.SetBinError(i+1, 1.0e6*val[1])
    for i, val in enumerate(layerTranslation):
        # convert cm to microns
        h1D_trans_layers.SetBinContent(i+1, 1.0e4*val[0])
        h1D_trans_layers.SetBinError(i+1, 1.0e4*val[1])

    #h1D_rot_dxdr_layers.GetYaxis().SetRangeUser(-1500, 1500)
    h1D_rot_dxdr_layers.Fit("pol1","QC")
    h1D_rot_dxdr_layers.Draw("E0")
    centerValue = evalInCenter(h1D_rot_dxdr_layers)
    print ">>> [%s] Rotation on layer 3.5: %.1f +/- %.1f urad" % (prettyStr, centerValue[0], centerValue[1])
    leg = r.TLegend(0.70,0.15,0.98,0.30)
    leg.AddEntry(h1D_rot_dxdr_layers.GetFunction("pol1"),"#scale[2.0]{L3.5 fit = %.0f #pm%.0f #murad}" % centerValue,"l")
    leg.Draw("same")
    c1.SaveAs(prefix + "h1D_rot_dxdr_layers" + ".png")

    #h1D_trans_layers.GetYaxis().SetRangeUser(-700,700)
    h1D_trans_layers.Fit("pol1","QC")
    h1D_trans_layers.Draw("E0")
    centerValue = evalInCenter(h1D_trans_layers)
    print ">>> [%s] X offset on layer 3.5: %.1f +/- %.1f microns" % (prettyStr, centerValue[0], centerValue[1])
    leg = r.TLegend(0.70,0.15,0.98,0.30)
    leg.AddEntry(h1D_trans_layers.GetFunction("pol1"),"#scale[2.0]{L3.5 fit = %.0f#pm%.0f #mum}" % centerValue,"l")
    leg.Draw("same")
    c1.SaveAs(prefix + "h1D_res_x_layers" + ".png")

    r.gStyle.SetPalette(1)



dChambers = { }

count = 0

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
    if(muon.nDT + muon.nCSC < minChambers): continue
    
    
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
    if(nMuons > 200):
        t1 = time.time()
        makePlotsCSC(chamber)
        print ">>> [%s] Took %.2f seconds" % (chamberPrettyString(*chamberKey), time.time() - t1)
    else:
        print ">>> [%s] Low statistics, so skipping plots" % (chamberPrettyString(*chamberKey))

