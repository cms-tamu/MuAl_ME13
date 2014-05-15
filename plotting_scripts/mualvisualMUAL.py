import ROOT as r
import struct, math, os, sys
from array import array

selectedCSC = (1, 2, 2, 1)  # endcap (p-1 m-0), station, wheel, chamber
#selectedCSC = (1, 1, 3, 1)
#print selectedCSC
specialSuffix = ""
#specialSuffix = "POSNEG"
#specialSuffix = "_TEST"
#specialSuffix = "_EDGECUTSTEST"
#specialSuffix = "_EDGECUTSBOTH"
#specialSuffix = "_EDGECUTS"
#specialSuffix = "_BOXCUTS"
#specialSuffix = "_BOXFIDCUTS"
#specialSuffix = "_SIDECUTS"
#specialSuffix = "_POSMU"
#specialSuffix = "_NEGMU"
#specialSuffix = "_TESTRES"
#specialSuffix = "_BINS"


doEdgeCuts = (specialSuffix is "_EDGECUTS" or specialSuffix is "_EDGECUTSTEST")
doBoxCuts = specialSuffix is "_BOXCUTS"
doBoxFidCuts = specialSuffix is "_BOXFIDCUTS"
doChargeCut = (specialSuffix is "_POSMU" or specialSuffix is "_NEGMU")
#doRotate, doTranslate = True, True
chargeToConsider = 1 if specialSuffix is "_POSMU" else -1
nSigma = 1.5 # number of sigma to do gaussian fit with
ME13pinR = 595.1500244141 # radius of ME13 ring alignment pin


print ">>> Args:", sys.argv
if(len(sys.argv) > 1):
    filename = sys.argv[1]
    if("layers" in filename):
        cscStr = sys.argv[1].split("_layers")[0]
    else:
        cscStr = sys.argv[1].split("_ALI")[0]
    #print cscStr
    cscStr = cscStr.split("ME")[1]
    cscInfo = cscStr.split("_")
    if(cscInfo[0] == "p"): cscInfo[0] = 1
    else: cscInfo[0] = 0
    selectedCSC = ( cscInfo[0], int(cscInfo[1]), int(cscInfo[2]), int(cscInfo[3]) )

    print ">>> Found argument with filename:", filename
    print ">>> Setting selectedCSC to", selectedCSC
else:
    print ">>> Using default file:", filename
    print ">>> Setting selectedCSC to default:", selectedCSC

if(specialSuffix is not ""): print ">>> Special suffix is",specialSuffix
    

def setAxisTitles(hist, xtitle, ytitle):
    hist.GetXaxis().SetTitle(xtitle);
    hist.GetYaxis().SetTitle(ytitle);


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

def getFitParamsGauss(hist):
    fit = hist.GetFunction("gaus")
    
    p0 = fit.GetParameter(0) # const
    p0e = fit.GetParError(0) # const
    
    p1 = fit.GetParameter(1) # mean
    p1e = fit.GetParError(1) # mean error
    
    #print p0, p0e, p1, p1e
    return p0, p0e, p1, p1e

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

nXbins, nYbins = 50, 50
XMax, YMax = 60., 100.
NUMLAYERS = 6

h2D_cnt_actual = []
h2D_cnt_tracks = []
h2D_res_x_actual = []
h2D_res_x_tracks = []
h2D_pull_tracks = []
h1D_res_x_rproj = []
h1D_res_x_xproj = []
h1D_res_x_rphiproj = []
h1D_res_rphi_yproj = []
h1D_res_x = []
h2D_statsig_tracks = []

h1D_res_x_avg = r.TH1F("h1D_res_x_avg", "avg x residuals", 100,-8.0,8.0) 
h1D_actual_localy = r.TH1F("h1D_actual_localy", "distribution of actual y hits",  nYbins,-75,75)
h1D_tracks_localy = r.TH1F("h1D_tracks_localy", "distribution of track y hits",  nYbins,-75,75)
h1D_actual_angle = r.TH1F("h1D_actual_angle", "angular distribution of hits",  nYbins,-0.08,0.08)
h1D_tracks_angle = r.TH1F("h1D_tracks_angle", "angular distribution of tracks",  nYbins,-0.08,0.08)
h1D_rot_dxdr_layers = r.TH1F("h1D_rot_dxdr_layers", "phiz (dx/dr) vs layer",  6,0.5,6.5  )
# h1D_rot_dxdrphi_layers = r.TH1F("h1D_rot_dxdrphi_layers", "phiz (dx/drphi) vs layer",  6,0.5,6.5  )
h1D_trans_layers = r.TH1F("h1D_trans_layers", "x offset vs layer",  6,0.5,6.5  )
h1D_res_x_rproj_avg = r.TProfile("h1D_res_x_rproj_avg", "avg x residual vs r",  nYbins,-75.0,75.0)
h2D_nlayers_hit = r.TProfile2D("h2D_nlayers_hit", "num layers hit",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
h2D_nDT_hit = r.TProfile2D("h2D_nDT_hit", "num DTs hit",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
h2D_nCSC_hit = r.TProfile2D("h2D_nCSC_hit", "num CSCs hit",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)
h2D_nTracker_hit = r.TProfile2D("h2D_nTracker_hit", "num tracker hits",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)

for i in range(NUMLAYERS):
    laypfx = "L" + str(i+1) + " "
    print ">>> Booking histos for layer",i+1
    # print laypfx
    h2D_cnt_actual.append(  r.TH2F(laypfx + "h2D_cnt_actual", laypfx + "actual hit locations",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
    h2D_cnt_tracks.append(  r.TH2F(laypfx + "h2D_cnt_tracks", laypfx + "track hit locations",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

    h2D_res_x_actual.append(  r.TProfile2D(laypfx + "h2D_res_x_actual", laypfx + "avg (per bin) x residuals at actual hit positions",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
    h2D_res_x_tracks.append(  r.TProfile2D(laypfx + "h2D_res_x_tracks", laypfx + "avg (per bin) x residuals at track positions",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

    h2D_pull_tracks.append(  r.TH2F(laypfx + "h2D_pull_tracks", laypfx + "sum of x residuals (pull)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )
    h2D_statsig_tracks.append(  r.TH2F(laypfx + "h2D_statsig_tracks", laypfx + "avg res_x / sigma (significance)",  nXbins,-XMax,XMax,  nYbins, -YMax,YMax)  )

    h1D_res_x.append(  r.TH1F(laypfx + "h1D_res_x", laypfx + "x residuals",  100,-8.0,8.0)  )

    h1D_res_x_rproj.append(  r.TProfile(laypfx + "h1D_res_x_rproj", laypfx + "x residual vs r",  nYbins,-75.0,75.0)  )
    h1D_res_x_xproj.append(  r.TProfile(laypfx + "h1D_res_x_xproj", laypfx + "x residual vs x",  nXbins,-100.0,100.0)  )
    h1D_res_x_rphiproj.append(  r.TProfile(laypfx + "h1D_res_x_rphiproj", laypfx + "x residual vs rphi",  nXbins,-70.0,70.0)  )

    # there are 2.530335 cm between each of the 48 wgs, so to get a binsize of 2.530335, we need (80-(-80))/2.530335=63.2=64 bins
    h1D_res_rphi_yproj.append(  r.TProfile(laypfx + "h1D_res_rphi_yproj", laypfx + "rphi residual vs local y",  64,-80.0,80.0)  )

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
    pt = 45

    ## q/pt equalization actually uses pz
    try: 
        if(abs(muon.pz) > 30): pt = abs(muon.pz)
        #print muon.pz
    except:
        try:
            if(muon.pt > 30): pt = muon.pt
        except: pass
        
    endcap, station, ring, chamber = struct.unpack("b", muon.endcap)[0], struct.unpack("b", muon.station)[0], struct.unpack("b", muon.ring)[0], struct.unpack("b", muon.chamber)[0]
    if (not (selectedCSC == (endcap, station, ring, chamber))): continue

            
    if muon.charge == 1:
        nPosMu += 1
        posIndexPt.append([i,pt])
    else:
        nNegMu += 1
        negIndexPt.append([i,pt])
    if(pt < minPt):
        #print pt
        minPt = pt

# print nMuons
goodIndices = getGoodMuonIndices(posIndexPt, negIndexPt, 17, minPt)
print ">>> Found %i positive muons and %i negative muons" % (nPosMu, nNegMu)
nMuonsToConsider = len(goodIndices) 
print ">>> Considering ",nMuonsToConsider,"total muons after {q/pT,q/pz} equalization"

goodMuons = {}
for i in range(nMuons): goodMuons[i] = 0
for e in goodIndices: goodMuons[e] = 1

nPosMu = 0
nNegMu = 0
count = 0
for idx, muon in enumerate(tt):
    count += 1
    if(count % 10000 == 0): print count 

    # if(count > 5000): break
    
    endcap, station, ring, chamber = ord(muon.endcap), ord(muon.station), ord(muon.ring), ord(muon.chamber)
    if muon.charge == 1:
        charge = 1
        nPosMu += 1
    else:
        charge = -1
        nNegMu += 1

    if(doChargeCut and not(charge == chargeToConsider)): continue

    if(selectedCSC == (endcap, station, ring, chamber)):

        #h1D_eta_before.Fill(muon.eta)
        
           # FIXME ANA
        if(goodMuons[idx] != 1): continue
        #h1D_eta_after.Fill(muon.eta)

        sumx = 0.0
        sumr = 0.0
        sumxcnt = 0

        # fid cuts FIXME (fix later)
        skipTrack = False
        for i in range(NUMLAYERS):
            layStr = "muon.lay%i_" % (i+1)
            actual_x, actual_y = eval(layStr+"x"), eval(layStr+"y")
            res_x, res_y = eval(layStr+"res_x"), eval(layStr+"res_y")
            track_x = actual_x + res_x
            track_y = actual_y + res_y
            angle = math.atan(1.0 * actual_x / (ME13pinR + actual_y))
            angle_track = math.atan(1.0 * track_x / (ME13pinR + track_y))
            angle_ytrack = math.atan(1.0 * actual_x / (ME13pinR + track_y))
            rphi_track = (ME13pinR)*math.atan(track_x / (ME13pinR + actual_y))
            rphi_ytrack = (ME13pinR)*math.atan(track_x / (ME13pinR + track_y))

            if(doBoxCuts):
                # MC removes box
                if(actual_x > 18.0 and (20.0 < actual_y < 65.0)): 
                    skipTrack = True
                    break

            if(doBoxFidCuts):
                # data
                #if(rphi_track > 14 and (20.0 < track_y < 65.0)): continue

                # MC removes box
                # if(actual_x > 18.0 and (20.0 < actual_y < 65.0)):
                #     skipTrack = True
                #     break
                # MC removes fiducial tracks around box
                if(rphi_track > 12 and (20.0 < track_y < 63.0)):
                # if(rphi_track < -14 and (20.0 < track_y < 63.0)):
                    skipTrack = True
                    break


            if(doEdgeCuts):
                # for ch 17, this makes data look good (x offset)
                #if(abs(rphi_track) > 16): continue
                if(abs(rphi_track) > 30):
                    skipTrack = True
                    break

        # if fiducial loop deemed this a bad track, skip it
        if(skipTrack): continue

        
        for i in range(NUMLAYERS):
        
            layStr = "muon.lay%i_" % (i+1)
           
          

            actual_x, actual_y = eval(layStr+"x"), eval(layStr+"y")
            res_x, res_y = eval(layStr+"res_x"), eval(layStr+"res_y")

            # propagated track = actual + residual
            track_x = actual_x + res_x
            track_y = actual_y + res_y
            #track_x = (eval(layStr+"x") + eval(layStr+"res_x"))
            #track_y = (eval(layStr+"y") + eval(layStr+"res_y"))

            if( actual_x < -998.0 or actual_y < -998.0): continue
            if( res_x < -998.0 or res_y < -998.0): continue
            #if( (eval(layStr+"x") == -999) or (eval(layStr+"y") == -999) ): continue
            #if( (eval(layStr+"res_x") == -999) or (eval(layStr+"res_y") == -999) ): continue
            
            
            angle = math.atan(1.0 * actual_x / (ME13pinR + actual_y))
            angle_track = math.atan(1.0 * track_x / (ME13pinR + track_y))
            angle_ytrack = math.atan(1.0 * actual_x / (ME13pinR + track_y))

            rphi_track = (ME13pinR)*math.atan(track_x / (ME13pinR + actual_y))
            rphi_ytrack = (ME13pinR)*math.atan(track_x / (ME13pinR + track_y))

            # TODO we're using rphi residuals now
            res_rphi = math.cos(angle)*res_x + math.sin(angle)*res_y

        
            if(i==2): # if layer 3 ("center" of chamber)
                h1D_actual_localy.Fill(actual_y)
                h1D_tracks_localy.Fill(track_y)
            
            sumx += res_x
            sumr += rProjection(track_y, angle_track)
            sumxcnt += 1

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

            #h1D_res_x_rproj[i].Fill(rProjection(track_y, angle_track), res_x)
            h1D_res_x_rproj[i].Fill(actual_y, res_x)

            #h1D_res_x_rproj[i].Fill(actual_y, res_x)
            #h1D_res_x_rproj[i].Fill(rProjection(actual_y, angle), res_x)
            
            #h1D_res_x_xproj[i].Fill(actual_x, res_x)
            #h1D_res_x_rphiproj[i].Fill((rProjection(track_y, angle_ytrack)+ME13pinR)*angle_ytrack, res_x)

            h1D_res_x_rphiproj[i].Fill(rphi_track, res_x)
            #h1D_res_x_rphiproj[i].Fill(track_x, res_x)
            #h1D_res_x_rphiproj[i].Fill(angle_track, res_x)

            h1D_res_rphi_yproj[i].Fill(actual_y, res_rphi)
       
        if(sumxcnt > 0):
            avgres_x = sumx / sumxcnt

            #if(abs(avgres_x - meanresx) < nSigma*sigmaresx):
            avgr = sumr / sumxcnt
            h1D_res_x_avg.Fill(avgres_x)
            h1D_res_x_rproj_avg.Fill(avgr, avgres_x)
            

print ">>> analyzed %i positive muons and %i negative muons" % (nPosMu, nNegMu)


c1.SetRightMargin(0.32);

c1.SetGridx()
c1.SetGridy()

prefix = filename.replace(".root",specialSuffix+"_plots/")

os.system("mkdir " + prefix)
os.system("cp -pv" + " indexbase.php " + prefix+"_index.php") #

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
        setAxisTitles(h1D_res_x_avg, "cm", "counts")
        h1D_res_x_avg.Draw()
        fitCut(h1D_res_x_avg, nSigma, "QC")
        c1.SaveAs(prefix + "h1D_res_x_avg" + suffix)
        
        setAxisTitles(h1D_res_x_rproj_avg, "cm", "cm")
        h1D_res_x_rproj_avg.GetYaxis().SetRangeUser(-1.4, 1.4)
        h1D_res_x_rproj_avg.Draw("E0")
        h1D_res_x_rproj_avg.Fit("pol1","QC")
        c1.SaveAs(prefix + "h1D_res_x_rproj_avg" + suffix)
       
        r.gStyle.SetOptFit(0)
        
        setAxisTitles(h1D_actual_angle, "phi", "counts")
        h1D_actual_angle.Draw()
        c1.SaveAs(prefix + "h1D_actual_angle" + suffix)

        setAxisTitles(h1D_tracks_angle, "phi", "counts")
        h1D_tracks_angle.Draw()
        c1.SaveAs(prefix + "h1D_tracks_angle" + suffix)


        setAxisTitles(h1D_actual_localy, "cm", "counts")
        h1D_actual_localy.Draw()
        c1.SaveAs(prefix + "h1D_actual_localy" + suffix)

        setAxisTitles(h1D_tracks_localy, "cm", "counts")
        h1D_tracks_localy.Draw()
        c1.SaveAs(prefix + "h1D_tracks_localy" + suffix)

        
        h2D_nlayers_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nlayers_hit" + suffix)

        h2D_nDT_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nDT_hit" + suffix)

        h2D_nCSC_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nCSC_hit" + suffix)

        h2D_nTracker_hit.Draw("colz")
        c1.SaveAs(prefix + "h2D_nTracker_hit" + suffix)

        # h1D_eta_before.Draw()
        # c1.SaveAs(prefix + "h1D_eta_before" + suffix)

        # h1D_eta_after.Draw()
        # c1.SaveAs(prefix + "h1D_eta_after" + suffix)



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
    h2D_res_x_actual[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_actual" + suffix)

    h2D_res_x_tracks[i].GetZaxis().SetRangeUser(-3.5, 3.5)
    h2D_res_x_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_tracks" + suffix)

    # multiply resx heatmap by hit occupancy to get "force" of pulling
    # h2D_pull_tracks[i].Divide(h2D_cnt_tracks[i])
    h2D_pull_tracks[i].GetZaxis().SetRangeUser(-80.0,80.0)
    h2D_pull_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_pull_tracks" + suffix)


    r.gStyle.SetOptFit(1) # display fitting parameters

    setAxisTitles(h1D_res_x_rproj[i], "cm", "cm")
    h1D_res_x_rproj[i].GetYaxis().SetRangeUser(-1.4, 1.4)
    # h1D_res_x_rproj[i].BuildOptions(0,0,"")
    h1D_res_x_rproj[i].Draw("E0")
    h1D_res_x_rproj[i].Fit("pol1","QC")
    c1.SaveAs(prefix + "h1D_res_x_rproj" + suffix)

    fitparams = getFitParams(h1D_res_x_rproj[i])
    layerRotationR.append([fitparams[2], fitparams[3]]) #p1, p1error


    #setAxisTitles(h1D_res_x_xproj[i], "cm", "cm")
    #h1D_res_x_xproj[i].GetYaxis().SetRangeUser(-0.2, 0.2)
    #h1D_res_x_xproj[i].GetXaxis().SetRangeUser(-30.0, 30.0)
    #h1D_res_x_xproj[i].Draw("E0")
    #h1D_res_x_xproj[i].Fit("pol1","QC")
    #c1.SaveAs(prefix + "h1D_res_x_xproj" + suffix)

    setAxisTitles(h1D_res_rphi_yproj[i], "cm", "cm")
    h1D_res_rphi_yproj[i].GetYaxis().SetRangeUser(-2.0, 2.0)
    #h1D_res_rphi_yproj[i].GetXaxis().SetRangeUser(-0.1, 0.1)
    h1D_res_rphi_yproj[i].Draw("E0")
    h1D_res_rphi_yproj[i].Fit("pol1","QC")
    c1.SaveAs(prefix + "h1D_res_rphi_yproj" + suffix)
   
    setAxisTitles(h1D_res_x_rphiproj[i], "cm", "cm")
    h1D_res_x_rphiproj[i].GetYaxis().SetRangeUser(-2.0, 2.0)
    #h1D_res_x_rphiproj[i].GetXaxis().SetRangeUser(-0.1, 0.1)
    h1D_res_x_rphiproj[i].Draw("E0")
    h1D_res_x_rphiproj[i].Fit("pol1","QC")
    c1.SaveAs(prefix + "h1D_res_x_rphiproj" + suffix)

    # fitparams = getFitParams(h1D_res_x_rphiproj[i])
    # layerRotationRphi.append([fitparams[2], fitparams[3]]) #p1, p1error
    
    setAxisTitles(h1D_res_x[i], "cm", "counts")
    h1D_res_x[i].Draw()
    fitCut(h1D_res_x[i], nSigma, "QC")
    c1.SaveAs(prefix + "h1D_res_x" + suffix)

    mu = h1D_res_x[i].GetMean()
    muerr = h1D_res_x[i].GetMeanError()
    print ">>> Unfitted mean x residual for layer %i: %f +/- %f" % (i, mu, muerr)

    fitparamsGauss = getFitParamsGauss(h1D_res_x[i])
    layerTranslation.append([fitparamsGauss[2], fitparamsGauss[3]]) #mean,meanerror


for i, val in enumerate(layerRotationR):
    h1D_rot_dxdr_layers.SetBinContent(i+1, 1.0e6*val[0])
    h1D_rot_dxdr_layers.SetBinError(i+1, 1.0e6*val[1])

# for i, val in enumerate(layerRotationRphi):
    # h1D_rot_dxdrphi_layers.SetBinContent(i+1, 1.0e6*val[0])
    # h1D_rot_dxdrphi_layers.SetBinError(i+1, 1.0e6*val[1])

for i, val in enumerate(layerTranslation):
    h1D_trans_layers.SetBinContent(i+1, 1.0e4*val[0])
    h1D_trans_layers.SetBinError(i+1, 1.0e4*val[1])


setAxisTitles(h1D_rot_dxdr_layers, "layer", "dx/dr (urad)")
h1D_rot_dxdr_layers.Fit("pol1","QC")
h1D_rot_dxdr_layers.Draw("E0")
c1.SaveAs(prefix + "h1D_rot_dxdr_layers" + ".png")

# setAxisTitles(h1D_rot_dxdrphi_layers, "layer", "dx/drphi (urad)")
# h1D_rot_dxdrphi_layers.Fit("pol1","QC")
# h1D_rot_dxdrphi_layers.Draw("E0")
# c1.SaveAs(prefix + "h1D_rot_dxdrphi_layers" + ".png")

setAxisTitles(h1D_trans_layers, "layer", "microns")
h1D_trans_layers.Fit("pol1","QC")
h1D_trans_layers.Draw("E0")
c1.SaveAs(prefix + "h1D_trans_layers" + ".png")

r.gStyle.SetPalette(1)

# setAxisTitles(h2D_res_x_ntracker, "hits", "cm")
# h2D_res_x_ntracker.Draw("colz")
# c1.SaveAs(prefix + "h2D_res_x_ntracker" + ".png")


