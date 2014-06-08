import ROOT as r
import struct, math, os, sys

c1 = r.TCanvas("Canvas1", "Alignment Visualizations")

h2D_cnt_actual = r.TH2F("h2D_cnt_actual", "MEp13 actual hit locations (L3)",  600,-800,800,  600, -800,800  )
h2D_cnt_tracks = r.TH2F("h2D_cnt_track", "MEp13 track hit locations (L3)",  600,-800,800,  600, -800,800  )
h2D_nDT = r.TProfile2D("h2D_nDT", "MEp13 nDT at actual hit locations",  600,-800,800,  600, -800,800  )
h2D_nCSC = r.TProfile2D("h2D_nCSC", "MEp13 nCSC at actual hit locations",  600,-800,800,  600, -800,800  )
h2D_nDTCSC = r.TProfile2D("h2D_nDTCSC", "MEp13 nDT + nCSC at actual hit locations",  600,-800,800,  600, -800,800  )
h2D_nTracker = r.TProfile2D("h2D_nTracker", "MEp13 nTracker at actual hit locations",  600,-800,800,  600, -800,800  )
h2D_nLayers = r.TProfile2D("h2D_nLayers", "MEp13 nLayers at actual hit locations",  600,-800,800,  600, -800,800  )

h1D_res_x = r.TH1F("h1D_res_x", "MEp13 x residuals (L3)", 100,-10,10)


def etaToR(eta):
    ME13Z = 695.15875
    theta = 2.0 * math.atan(math.exp(-1*eta))
    return ME13Z * math.tan(theta)

def polarToCartesian(r, phi):
    x = r * math.cos(phi)
    y = r * math.sin(phi)
    return x, y

def localToGlobalPolar(x, y, chambernum):
    ME13R = 595.1500244141#695.1587524414
    #r = ME13Z + y
    r = ((ME13R + y)**2 + x**2)**(0.5)
    #print chambernum
    #phi = 2.0 * 3.141592 / 36.0 * chambernum + x / r
    phi = 2.0 * 3.141592 / 36.0 * (chambernum - 1) + math.asin(x / r)
    #print "phi",phi
    return r, phi

pinRadius = 595.1500244   
histRange = pinRadius * 1.5 
nChambers = 36
    
count = 0

print "args:",sys.argv
print "filenames =",sys.argv[1:]
filenames = sys.argv[1:]

for filename in filenames:
    print "#### opened", filename
    fh = r.TFile(filename)
    tt = r.gDirectory.Get("csc_layer_ttree")
    for muon in tt:
        count += 1
        if(count%100000 == 0):
            print count

        #if(count > 200000): break

        if(not muon.select or muon.nlayers < 6): continue
        endcap, station, ring, chamber = ord(muon.endcap), ord(muon.station), ord(muon.ring), ord(muon.chamber)


        if ( endcap == 1 and station == 1 and ring == 3 ):
            # layer 3 values
            i = 3
            actual_x, actual_y = muon.hit_x[i], muon.hit_y[i]
            res_x, res_y = muon.res_x[i], muon.res_y[i]

            #if(actual_y > 80.0): continue
            
            globr, globphi = localToGlobalPolar( actual_x, actual_y, chamber )
            globrTraj, globphiTraj = localToGlobalPolar( actual_x+res_x, actual_y+res_y, chamber )

            xpos, ypos = polarToCartesian( globr, globphi )
            xposTraj, yposTraj = polarToCartesian( globrTraj, globphiTraj )
            h2D_cnt_actual.Fill(-1.0 * xpos, ypos )
            h2D_cnt_tracks.Fill(-1.0 * xposTraj, yposTraj )

            h2D_nCSC.Fill(-1.0 * xpos, ypos, muon.nCSC )
            h2D_nDT.Fill(-1.0 * xpos, ypos, muon.nDT )
            h2D_nTracker.Fill(-1.0 * xpos, ypos, muon.nTracker )
            h2D_nLayers.Fill(-1.0 * xpos, ypos, muon.nlayers )
            h2D_nDTCSC.Fill(-1.0 * xpos, ypos, muon.nCSC+muon.nDT )

            h1D_res_x.Fill( res_x )
            #if(res_x < -900 or res_y < -900): print "oopsie", count
            
        else:
            continue

c1.SetGridx()
c1.SetGridy()
c1.SetCanvasSize(1700,1230)

#prefix = basename.replace(".root","_")
prefix = "ring_plots/"
os.system("mkdir -p %s" % prefix)
os.system("cp indexbase.php %s" % prefix)
suffix = ".png"

c1.SetRightMargin(0.32);

# draw lines separating chambers in phi
lines = []
text = []
for i in range(1,nChambers+1):
    phi = 2.0 * 3.141592 / nChambers * (i - 0.5)
    lines.append( r.TLine( 0,0, 0.8*histRange*math.cos(phi),0.8*histRange*math.sin(phi) ) )
    phi = 2.0 * 3.141592 / nChambers * (i - 1)
    phi = 3.141592 - phi
    label = r.TText( 0.5*histRange*math.cos(phi),0.5*histRange*math.sin(phi), str(i) )
    label.SetTextAlign(22)
    label.SetTextSize(0.02)
    text.append(label)
    
h2D_cnt_actual.Draw("colz")
for line in lines: line.Draw()
for t in text: t.Draw()
c1.SaveAs(prefix + "h2D_cnt_actual" + suffix)

h2D_cnt_tracks.Draw("colz")
for line in lines: line.Draw()
for t in text: t.Draw()
c1.SaveAs(prefix + "h2D_cnt_tracks" + suffix)

h2D_nCSC.Draw("colz")
for line in lines: line.Draw()
for t in text: t.Draw()
c1.SaveAs(prefix + "h2D_nCSC" + suffix)

h2D_nDT.Draw("colz")
for line in lines: line.Draw()
for t in text: t.Draw()
c1.SaveAs(prefix + "h2D_nDT" + suffix)

h2D_nDTCSC.Draw("colz")
for line in lines: line.Draw()
for t in text: t.Draw()
c1.SaveAs(prefix + "h2D_nDTCSC" + suffix)

h2D_nTracker.GetZaxis().SetRangeUser(0, 21)
h2D_nTracker.Draw("colz")
for line in lines: line.Draw()
for t in text: t.Draw()
c1.SaveAs(prefix + "h2D_nTracker" + suffix)

h2D_nLayers.Draw("colz")
for line in lines: line.Draw()
for t in text: t.Draw()
c1.SaveAs(prefix + "h2D_nLayers" + suffix)

h1D_res_x.Draw()
c1.SaveAs(prefix + "h1D_res_x" + suffix)

c1.Update()
#raw_input()


