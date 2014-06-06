import ROOT as r
import struct, math, os, sys
from array import array

selectedCSC = (1, 2, 2, 1)  # endcap (p-1 m-0), station, wheel, chamber
specialSuffix = ""

# suppress "Info in..." messages, as there will be a lot
r.gErrorIgnoreLevel = r.kWarning


doChargeCut = (specialSuffix is "_POSMU" or specialSuffix is "_NEGMU")
chargeToConsider = 1 if specialSuffix is "_POSMU" else -1

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


#def setAxisTitles(hist, xtitle, ytitle):
    #hist.GetXaxis().SetTitle(xtitle);
    #hist.GetYaxis().SetTitle(ytitle);


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

c1 = r.TCanvas("Canvas1", "Alignment Visualizations")

typePrefix = "ME%s%d/%d/%d" % ("+" if selectedCSC[0] == 1 else "-", selectedCSC[1], selectedCSC[2], selectedCSC[3])
if(isMC):
    typePrefix += " #font[2]{#color[4]{MC}} "
else:
    typePrefix += " #font[2]{#color[2]{DATA}} "

nXbins, nYbins = 50, 50
XMax, YMax = 60., 100.
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
h2D_nlayers_hit = r.TProfile2D("h2D_nlayers_hit", typePrefix+"num layers hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
h2D_nDT_hit = r.TProfile2D("h2D_nDT_hit", typePrefix+"num DTs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
h2D_nCSC_hit = r.TProfile2D("h2D_nCSC_hit", typePrefix+"num CSCs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
h2D_nTracker_hit = r.TProfile2D("h2D_nTracker_hit", typePrefix+"num tracker hits;actual hit x;actual hit y",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)

for i in range(NUMLAYERS):
    laypfx = " #font[2]{L" + str(i+1) + "} "
    print ">>> Booking histos for layer",i+1
    h2D_cnt_actual.append(  r.TH2F(laypfx + "h2D_cnt_actual", typePrefix+laypfx + "actual hit locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
    h2D_cnt_tracks.append(  r.TH2F(laypfx + "h2D_cnt_tracks", typePrefix+laypfx + "track locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

    h2D_res_x_actual.append(  r.TProfile2D(laypfx + "h2D_res_x_actual", typePrefix+laypfx + "avg (per bin) x residuals at actual hit positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
    h2D_res_x_tracks.append(  r.TProfile2D(laypfx + "h2D_res_x_tracks", typePrefix+laypfx + "avg (per bin) x residuals at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

    h2D_pull_tracks.append(  r.TH2F(laypfx + "h2D_pull_tracks", typePrefix+laypfx + "sum of x residuals (pull) at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
    h2D_statsig_tracks.append(  r.TH2F(laypfx + "h2D_statsig_tracks", typePrefix+laypfx + "avg res_x / sigma (significance) at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

    h1D_res_x.append(  r.TH1F(laypfx + "h1D_res_x", typePrefix+laypfx + "x residuals;x residual (cm);counts",  100,-8.0,8.0)  )

    h1D_res_x_rproj.append(  r.TProfile(laypfx + "h1D_res_x_rproj", typePrefix+laypfx + "x residual vs r;x residual (cm);r (cm)",  50,-75.0,75.0)  )
    #h1D_res_x_xproj.append(  r.TProfile(laypfx + "h1D_res_x_xproj", laypfx + "x residual vs x",  nXbins,-100.0,100.0)  )
    h1D_res_x_rphiproj.append(  r.TProfile(laypfx + "h1D_res_x_rphiproj", typePrefix+laypfx + "x residual vs scaled-phi;x residual (cm);scaled-phi (cm)",  50,-75.0,75.0)  )

    # there are 2.530335 cm between each of the 48 wgs, so to get a binsize of 2.530335, we need (80-(-80))/2.530335=63.2=64 bins
    h1D_res_rphi_yproj.append(  r.TProfile(laypfx + "h1D_res_rphi_yproj", typePrefix+laypfx + "scaled-phi residual vs local y;rphi residual (cm);local y (cm)",  64,-80.0,80.0)  )

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

    # if(count > 5000): break
    if(not muon.select or muon.nlayers < 6): continue
    
    endcap, station, ring, chamber = ord(muon.endcap), ord(muon.station), ord(muon.ring), ord(muon.chamber)

    if(doChargeCut and not(charge == chargeToConsider)): continue

    if(selectedCSC == (endcap, station, ring, chamber)):

        
        if(goodMuons[idx] != 1): continue


        for i in range(NUMLAYERS):
        
          

            try:
                actual_x, actual_y = muon.hit_x[i], muon.hit_y[i]
                res_x, res_y = muon.res_x[i], muon.res_y[i]
            except:
                layStr = "muon.lay%i_" % (i+1)
                actual_x, actual_y = eval(layStr+"x"), eval(layStr+"y")
                res_x, res_y = eval(layStr+"res_x"), eval(layStr+"res_y")

            # propagated track = actual + residual
            track_x = actual_x + res_x
            track_y = actual_y + res_y

            if( actual_x < -998.0 or actual_y < -998.0): continue
            if( res_x < -998.0 or res_y < -998.0): continue
            
            angle = math.atan(1.0 * actual_x / (ME13pinR + actual_y))
            angle_track = math.atan(1.0 * track_x / (ME13pinR + track_y))
            angle_ytrack = math.atan(1.0 * actual_x / (ME13pinR + track_y))

            rphi_track = (ME13pinR)*math.atan(track_x / (ME13pinR + track_y))

            res_rphi = math.cos(angle)*res_x + math.sin(angle)*res_y

        
            if(i==2): # if layer 3 ("center" of chamber)
                h1D_actual_localy.Fill(actual_y)
                h1D_tracks_localy.Fill(track_y)
            
            #sumx += res_x
            #sumr += rProjection(track_y, angle_track)
            #sumxcnt += 1

            h1D_actual_angle.Fill(angle)
            h1D_tracks_angle.Fill(angle_track)

            h1D_res_x[i].Fill(res_x)
            
            h2D_cnt_actual[i].Fill(actual_x, actual_y)
            h2D_cnt_tracks[i].Fill(track_x, track_y)

            h2D_res_x_actual[i].Fill(actual_x, actual_y, res_x)
            h2D_res_x_tracks[i].Fill(track_x, track_y, res_x)

            h2D_pull_tracks[i].Fill(track_x, track_y, res_x)

            h2D_nlayers_hit.Fill(actual_x, actual_y, muon.nlayers)
            try:
                h2D_nDT_hit.Fill(actual_x, actual_y, muon.nDT)
                h2D_nCSC_hit.Fill(actual_x, actual_y, muon.nCSC)
                h2D_nTracker_hit.Fill(actual_x, actual_y, muon.nTracker)
            except: pass

            h1D_res_x_rproj[i].Fill(rProjection(actual_y,angle), res_x)

            h1D_res_x_rphiproj[i].Fill(rphi_track, res_x)
            #h1D_res_x_rphiproj[i].Fill(rphi_ytrack, res_x)

            h1D_res_rphi_yproj[i].Fill(actual_y, res_rphi)
       
        # if(sumxcnt > 0):
            
            # TODO make this into a linear fit and evaluate at layer 3.5
            #avgres_x = sumx / sumxcnt 
            #h1D_res_x_avg.Fill(avgres_x)

            #avgr = sumr / sumxcnt
            #h1D_res_x_rproj_avg.Fill(avgr, avgres_x)
            

c1.SetRightMargin(0.32);
c1.SetGridx()
c1.SetGridy()


os.system("mkdir -p " + prefix)
os.system("cp -p " + " indexbase.php " + prefix+"_index.php") #

r.gStyle.SetOptStat("sirmen") # add skewness 

layerRotationR = []
# layerRotationRphi = []
layerTranslation = []

for i in range(NUMLAYERS):

    r.gStyle.SetPalette(1)
    r.gStyle.SetOptFit(0)

    suffix = "_LAY" + str(i+1) + ".png"
    
    if(i == 0):
        r.gStyle.SetOptFit(1)
        #setAxisTitles(h1D_res_x_avg, "cm", "counts")
        #h1D_res_x_avg.Draw()
        #fitCut(h1D_res_x_avg, nSigma, "QC")
        #c1.SaveAs(prefix + "h1D_res_x_avg" + suffix)
        
        #setAxisTitles(h1D_res_x_rproj_avg, "cm", "cm")
        #h1D_res_x_rproj_avg.GetYaxis().SetRangeUser(-1.4, 1.4)
        #h1D_res_x_rproj_avg.Draw("E0")
        #h1D_res_x_rproj_avg.Fit("pol1","QC")
        #c1.SaveAs(prefix + "h1D_res_x_rproj_avg" + suffix)
       
        r.gStyle.SetOptFit(0)
        
        h1D_actual_angle.Draw()
        c1.SaveAs(prefix + "h1D_angle_actual" + suffix)

        h1D_tracks_angle.Draw()
        c1.SaveAs(prefix + "h1D_angle_tracks" + suffix)


        h1D_actual_localy.Draw()
        c1.SaveAs(prefix + "h1D_localy_actual" + suffix)

        h1D_tracks_localy.Draw()
        c1.SaveAs(prefix + "h1D_loclay_tracks" + suffix)

        
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

    h2D_res_x_actual[i].GetZaxis().SetRangeUser(-3.5, 3.5)
    fixFlooredBins(h2D_res_x_actual[i], minZ=-3.5)
    h2D_res_x_actual[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_actual" + suffix)

    h2D_res_x_tracks[i].GetZaxis().SetRangeUser(-3.5, 3.5)
    fixFlooredBins(h2D_res_x_tracks[i], minZ=-3.5)
    h2D_res_x_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_tracks" + suffix)

    # multiply resx heatmap by hit occupancy to get "force" of pulling
    # h2D_pull_tracks[i].Divide(h2D_cnt_tracks[i])
    h2D_pull_tracks[i].GetZaxis().SetRangeUser(-80.0,80.0)
    h2D_pull_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_pull_tracks" + suffix)


    r.gStyle.SetOptFit(1) # display fitting parameters

    #setAxisTitles(h1D_res_x_rproj[i], "cm", "cm")
    h1D_res_x_rproj[i].GetYaxis().SetRangeUser(-1.4, 1.4)
    # h1D_res_x_rproj[i].BuildOptions(0,0,"")
    h1D_res_x_rproj[i].Draw("E0")
    h1D_res_x_rproj[i].Fit("pol1","QC")
    c1.SaveAs(prefix + "h1D_res_x_rproj" + suffix)

    fitparams = getFitParams(h1D_res_x_rproj[i])
    layerRotationR.append([fitparams[2], fitparams[3]]) #p1, p1error


    #setAxisTitles(h1D_res_rphi_yproj[i], "cm", "cm")
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
    
    #setAxisTitles(h1D_res_x[i], "cm", "counts")
    h1D_res_x[i].Draw()
    fitCut(h1D_res_x[i], nSigma, "QC")
    c1.SaveAs(prefix + "h1D_res_x" + suffix)

    # mu = h1D_res_x[i].GetMean()
    # muerr = h1D_res_x[i].GetMeanError()
    # print ">>> Unfitted mean x residual for layer %i: %f +/- %f" % (i, mu, muerr)
    print ">>> Saved plots for layer %d" % (i+1)

    fitparamsGauss = getFitParamsGauss(h1D_res_x[i])
    layerTranslation.append([fitparamsGauss[2], fitparamsGauss[3]]) #mean,meanerror


for i, val in enumerate(layerRotationR):
    h1D_rot_dxdr_layers.SetBinContent(i+1, 1.0e6*val[0])
    h1D_rot_dxdr_layers.SetBinError(i+1, 1.0e6*val[1])

for i, val in enumerate(layerTranslation):
    h1D_trans_layers.SetBinContent(i+1, 1.0e4*val[0])
    h1D_trans_layers.SetBinError(i+1, 1.0e4*val[1])


#setAxisTitles(h1D_rot_dxdr_layers, "layer", "dx/dr (urad)")
h1D_rot_dxdr_layers.GetYaxis().SetRangeUser(-1500, 1500)
h1D_rot_dxdr_layers.Fit("pol1","QC")
h1D_rot_dxdr_layers.Draw("E0")
centerValue = evalInCenter(h1D_rot_dxdr_layers)
print ">>> Rotation on layer 3.5: %.1f +/- %.1f urad" % centerValue
leg = r.TLegend(0.70,0.15,0.98,0.30)
leg.AddEntry(h1D_rot_dxdr_layers.GetFunction("pol1"),"#scale[2.0]{L3.5 fit = %.0f #pm%.0f #murad}" % centerValue,"l")
leg.Draw("same")
c1.SaveAs(prefix + "h1D_rot_dxdr_layers" + ".png")

#setAxisTitles(h1D_trans_layers, "layer", "microns")
h1D_trans_layers.GetYaxis().SetRangeUser(-700,700)
h1D_trans_layers.Fit("pol1","QC")
h1D_trans_layers.Draw("E0")
centerValue = evalInCenter(h1D_trans_layers)
print ">>> X offset on layer 3.5: %.1f +/- %.1f urad" % centerValue
leg = r.TLegend(0.70,0.15,0.98,0.30)
leg.AddEntry(h1D_trans_layers.GetFunction("pol1"),"#scale[2.0]{L3.5 fit = %.0f#pm%.0f #mum}" % centerValue,"l")
leg.Draw("same")
c1.SaveAs(prefix + "h1D_res_x_layers" + ".png")

r.gStyle.SetPalette(1)
