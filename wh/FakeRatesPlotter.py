from glob import glob
import os
from FinalStateAnalysis.PlotTools.HistToTGRaphErrors import HistStackToTGRaphErrors
from FinalStateAnalysis.PlotTools.RebinView  import RebinView
from FinalStateAnalysis.PlotTools.Plotter    import Plotter
from FinalStateAnalysis.Utilities.floatformatting import smart_float_format
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors
from rootpy.utils import asrootpy
import rootpy.plotting.views as views
import rootpy.plotting as plotting
import logging
import ROOT
import sys
from pdb import set_trace

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )

class FakeRatesPlotter(Plotter):
    def __init__(self, input_dir, output_dir):
        jobid      = os.environ['jobid']
        resultdir  = 'results/%s' % jobid
        lumidir    = 'inputs/%s' % jobid
        self.sqrts = 7 if '7TeV' in jobid else 8
        
        self.output_dir = os.path.join(
            resultdir,
            'plots',
            output_dir
        )
        
        input_files = glob(
            os.path.join(
                resultdir,
                input_dir,
                '*.root'
            )
        )
        lumi_files = glob(
            os.path.join(
                lumidir,
                '*.lumicalc.sum'
            )
        )
        
        
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        super(FakeRatesPlotter, self).__init__(input_files, lumi_files, 
                                            self.output_dir, None)

    def plot_fakes(self, neighbors, idIso, variable, rebin=1, xaxis='', fit=None):
        dirpath = os.path.join(neighbors, 'pass', idIso)
        wz_view = views.SubdirectoryView( 
            RebinView(
                self.get_view('WZJetsTo3LNu*ZToTauTau*'),
                rebin
            ),
            dirpath
            )
        zz_view = views.SubdirectoryView( 
            RebinView(
                self.get_view('ZZJetsTo4L*'),
                rebin
            ),
            dirpath
            )

        data_view = RebinView(
            self.get_view('data'),
            rebin
        )

        fail_view = views.SubdirectoryView( 
            data_view,
            os.path.join(neighbors, 'fail', idIso)
            )
        pass_view = views.SubdirectoryView( 
            data_view,
            os.path.join(neighbors, 'pass', idIso)
            )

        diboson_view = views.TitleView(
            views.SumView(wz_view, zz_view),
            'diboson'
        )
        fakes_view   = views.TitleView(
            views.StyleView(
                fail_view, 
                **remove_name_entry(data_styles['Zjets*'])
            ),
            'Fakes'
        )

        expected_view = views.StackView(
            diboson_view,
            fakes_view
        )
        
        data_hist      = pass_view.Get(variable)
        expected_hist  = expected_view.Get(variable)
        estimate_error = sum( expected_hist.GetHists() )
        estimate_error.SetFillStyle(3013)
        estimate_error.SetMarkerSize(0)
        estimate_error.SetFillColor(ROOT.EColor.kBlack)

        expected_hist.SetMaximum(2 * max(
            data_hist.GetMaximum(), 
            expected_hist.GetMaximum())
        )
        expected_hist.Draw()
        data_hist.Draw('same')
        estimate_error.Draw('PE2 same')

        if xaxis:
            expected_hist.GetXaxis().SetTitle(xaxis)
        
        legend  = self.add_legend(expected_hist, leftside=False, entries=4)
        legend.AddEntry(estimate_error)
        legend.AddEntry(data_hist)
        legend.Draw()
        ratio_plot = self.add_ratio_plot(data_hist, expected_hist, ratio_range=2.)

        if fit:
            #set_trace()
            fitrange = fit.get('range', False)
            if not fitrange:
                nbins = ratio_plot.GetNbinsX()
                fitrange = x_range if x_range else [ ratio_plot.GetBinLowEdge(1), 
                   ratio_plot.GetBinLowEdge(nbins)+ratio_plot.GetBinWidth(nbins)]
            self.lower_pad.cd()
            function = self.fit_shape(ratio_plot, fit['model'], fitrange, fit.get('options','IRMENS'))
            toprint  = '#chi^{2} / DoF = %.2f / %i\n' % (function.GetChisquare() , function.GetNDF())
            for i in range(function.GetNpar()):
                name  = function.GetParName(i) 
                value = function.GetParameter(i)
                error = function.GetParError(i)
                toprint += '%s = %s\n' % (name, smart_float_format((value, error))) #%f #pm %f
        
            stat_box = self.make_text_box(toprint[:-1],fit.get('stat position','bottom-left'))
            stat_box.Draw()
            self.keep.append(stat_box)
            #print toprint
            self.pad.cd()


        
if __name__ == '__main__':
    fitmodel = { 
        'model' : ("slope*x + offset", "offset[0, -1, 1], slope[0,-1,1]"), 
        'range' : [50, 500], 
        'options' : 'WLQIRMENS', 
        'stat position' : 'bottom-left'
    }
    #
    # MM Fakes
    #
    plotter = FakeRatesPlotter('FakeRatesAnalyzerMM','MMFakes')
    sqrts   = plotter.sqrts

    for i in ['k50', 'k20']:
        plotter.plot_fakes(i, 'h2taucuts020', 'tagMuonJetMass', rebin=20, 
                           xaxis='m_{tag #mu Jet} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('tagMuonJetMass-%s-h2taucuts020' % i)

        plotter.plot_fakes(i, 'h2taucuts020', 'probeMuonJetMass', rebin=20, 
                           xaxis='m_{probe #mu Jet} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('probeMuonJetMass-%s-h2taucuts020' % i)

        plotter.plot_fakes(i, 'h2taucuts020', 'tagMuonJetMass', rebin=200, 
                           xaxis='') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('yields2-%s-h2taucuts020' % i)

        plotter.plot_fakes(i, 'h2taucuts020', 'probeMuonJetMass', rebin=200, 
                           xaxis='') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('yields-%s-h2taucuts020' % i)

        plotter.plot_fakes(i, 'h2taucuts020', 'LT', rebin=5, 
                           xaxis='LT (GeV)', fit=fitmodel ) #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('LT-%s-h2taucuts020' % i)

        plotter.plot_fakes(i, 'h2taucuts020', 'muonPt', rebin=10, 
                           xaxis='p_{T #mu} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('pt-%s-h2taucuts020' % i)

        plotter.plot_fakes(i, 'h2taucuts020', 'muonJetPt', rebin=10, 
                           xaxis='Jet p_{T #mu} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('jetpt-%s-h2taucuts020' % i)


    ## plotter.plot_fakes(i, 'h2taucuts', 'tagMuonJetMass', rebin=20, 
    ##                    xaxis='m_{tag #mu Jet} (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('tagMuonJetMass-h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'h2taucuts', 'probeMuonJetMass', rebin=20, 
    ##                    xaxis='m_{probe #mu Jet} (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('probeMuonJetMass-h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'h2taucuts', 'probeMuonJetMass', rebin=200, 
    ##                    xaxis='') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('yields-h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'h2taucuts', 'LT', rebin=5, 
    ##                    xaxis='LT (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('LT-h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'h2taucuts', 'muonPt', rebin=10, 
    ##                    xaxis='LT (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('pt-h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'h2taucuts', 'muonJetPt', rebin=10, 
    ##                    xaxis='LT (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('jetpt-h2taucuts')


    #
    # EM Fakes
    #
    plotter = FakeRatesPlotter('FakeRatesAnalyzerEM_muon','EMFakes')
    sqrts   = plotter.sqrts

    for i in ['k50', 'k20']:
        plotter.plot_fakes(i, 'h2taucuts', 'electronJetMass', rebin=20, 
                           xaxis='m_{tag e Jet} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('tagElectronJetMass-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'h2taucuts', 'muonJetMass', rebin=20, 
                           xaxis='m_{probe #mu Jet} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('probeMuonJetMass-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'h2taucuts', 'muonJetMass', rebin=200, 
                           xaxis='') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('muon-yields-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'h2taucuts', 'LT', rebin=5, 
                           xaxis='LT (GeV)', fit=fitmodel) #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('LT-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'h2taucuts', 'muonPt', rebin=10, 
                           xaxis='p_{T #mu} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('muonpt-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'h2taucuts', 'muonJetPt', rebin=10, 
                           xaxis='Jet p_{T #mu} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('muonjetpt-%s-h2taucuts' % i)

    plotter = FakeRatesPlotter('FakeRatesAnalyzerEM_electron','EMFakes')
    sqrts   = plotter.sqrts

    for i in ['k50', 'k20']:
        plotter.plot_fakes(i, 'eid12Medium_h2taucuts', 'electronJetMass', rebin=20, 
                           xaxis='m_{probe e Jet} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('probeElectronJetMass-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts', 'electronJetMass', rebin=200, 
                           xaxis='') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('electron-yields-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts', 'muonJetMass', rebin=20, 
                           xaxis='m_{tag #mu Jet} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('tagMuonJetMass-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts', 'LT', rebin=5, 
                           xaxis='LT (GeV)', fit=fitmodel) #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('LT-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts', 'electronPt', rebin=10, 
                           xaxis='p_{T e} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('elept-%s-h2taucuts' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts', 'electronJetPt', rebin=10, 
                           xaxis='Jet p_{T e} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('elejetpt-%s-h2taucuts' % i)

    #
    # EE Fakes
    #
    plotter = FakeRatesPlotter('FakeRatesAnalyzerEE','EEFakes')
    sqrts   = plotter.sqrts

    for i in ['k50', 'k20']:
        plotter.plot_fakes(i, 'eid12Medium_h2taucuts020', 'tagElectronJetMass', rebin=20, 
                           xaxis='m_{tag e Jet} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('tagElectronJetMass-%s-eid12Medium_h2taucuts020' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts020', 'tagElectronJetMass', rebin=200, 
                           xaxis='') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('yields-%s-eid12Medium_h2taucuts020' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts020', 'probeElectronJetMass', rebin=20, 
                           xaxis='m_{probe e Jet} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('probeElectronJetMass-%s-eid12Medium_h2taucuts020' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts020', 'LT', rebin=5, 
                           xaxis='LT (GeV)', fit=fitmodel) #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('LT-%s-eid12Medium_h2taucuts020' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts020', 'electronPt', rebin=10, 
                           xaxis='p_{T e} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('pt-%s-h2taucuts020' % i)

        plotter.plot_fakes(i, 'eid12Medium_h2taucuts020', 'electronJetPt', rebin=10, 
                           xaxis='Jet p_{T e} (GeV)') #, leftside=False)
        plotter.add_cms_blurb(sqrts)
        plotter.save('jetpt-%s-h2taucuts020' % i)



    ## plotter.plot_fakes(i, 'eid12Tight_h2taucuts', 'tagElectronJetMass', rebin=20, 
    ##                    xaxis='m_{tag e Jet} (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('tagElectronJetMass-eid12Tight_h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'eid12Tight_h2taucuts', 'tagElectronJetMass', rebin=200, 
    ##                    xaxis='') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('yields-eid12Tight_h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'eid12Tight_h2taucuts', 'probeElectronJetMass', rebin=20, 
    ##                    xaxis='m_{probe e Jet} (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('probeElectronJetMass-eid12Tight_h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'eid12Tight_h2taucuts', 'LT', rebin=5, 
    ##                    xaxis='LT (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('LT-eid12Tight_h2taucuts')
    ## 
    ## plotter.plot_fakes(i, 'eid12Tight_h2taucuts', 'electronPt', rebin=10, 
    ##                    xaxis='LT (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('pt-h2taucuts020')
    ## 
    ## plotter.plot_fakes(i, 'eid12Tight_h2taucuts', 'electronJetPt', rebin=10, 
    ##                    xaxis='LT (GeV)') #, leftside=False)
    ## plotter.add_cms_blurb(sqrts)
    ## plotter.save('jetpt-h2taucuts020')
