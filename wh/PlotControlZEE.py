'''

Make inclusive Z->mumu control plots

'''

import os
import glob
from FinalStateAnalysis.PlotTools.Plotter import Plotter
from FinalStateAnalysis.PlotTools.MedianView import MedianView
from FinalStateAnalysis.PlotTools.DifferentialView import DifferentialView
from FinalStateAnalysis.PlotTools.SubtractionView import SubtractionView
import rootpy.plotting.views as views
import rootpy.plotting as plotting
from FinalStateAnalysis.PlotTools.HistToTGRaphErrors import HistStackToTGRaphErrors
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors
import ROOT
import logging
import sys
logging.basicConfig(stream=sys.stderr, level=logging.INFO)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

class ControlZEEPlotter(Plotter):
    def __init__(self):
        self.jobid = os.environ['jobid']
        self.output_dir = os.path.join('results', self.jobid, 'plots', 'zee')
        samples = [
            'Zjets_M50',
            'WZ*',
            'ZZ*',
            'WW*',
            'TT*',
            'WplusJets*',
            "data_DoubleE*",
            ]
        files = []
        lumifiles = []
        self.sqrts = 7 if '7TeV' in self.jobid else 8
        for x in samples:
            files.extend(glob.glob('results/%s/ControlZEE/%s.root' % (self.jobid, x)))
            lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (self.jobid, x)))
        super(ControlZEEPlotter, self).__init__(files, lumifiles, self.output_dir)
        self.mc_samples = filter(lambda x: 'data' not in x, samples)

    def make_fakes_view(self, data_view, sign, weight_type):
        # View of weighted obj1-fails data
        obj1_view = views.SubdirectoryView(data_view, '%s/f1p2/%s' % (sign, weight_type))
        # View of weighted obj2-fails data
        obj2_view = views.SubdirectoryView(data_view, '%s/p1f2/%s' % (sign, weight_type))
        # View of weighted obj1&2-fails data
        obj12_view = views.SubdirectoryView(data_view, '%s/f1f2/%s' % (sign, weight_type))
        # Give the individual object views nice colors
        subtract_obj12_view = views.ScaleView(obj12_view, -1)
        return obj1_view, obj2_view, subtract_obj12_view
    

    def get_flip_data(self, rebin=1, xaxis='', data_type='data'):
        data_view = self.get_view(data_type)
        data_view = self.rebin_view(data_view, rebin) if rebin != 1 else data_view
        #Get ss/p1p2 views
        ss_p1p2_view = views.SubdirectoryView( data_view, 'ss/p1p2')
        ss_p1p2_view = views.TitleView( views.StyleView(ss_p1p2_view, **data_styles['data*']), 'observed;%s' % xaxis)

        #def make_fakes_view(sign, weight_type):
        #    # View of weighted obj1-fails data
        #    obj1_view = views.SubdirectoryView(data_view, '%s/f1p2/%s' % (sign, weight_type))
        #    # View of weighted obj2-fails data
        #    obj2_view = views.SubdirectoryView(data_view, '%s/p1f2/%s' % (sign, weight_type))
        #    # View of weighted obj1&2-fails data
        #    obj12_view = views.SubdirectoryView(data_view, '%s/f1f2/%s' % (sign, weight_type))
        #    # Give the individual object views nice colors
        #    subtract_obj12_view = views.ScaleView(obj12_view, -1)
        #    return obj1_view, obj2_view, subtract_obj12_view

        #Get fakes according to WJets or QCD
        #ss_f1p2_qcd_view, ss_p1f2_qcd_view, ss_f1f2_qcd_view = self.make_fakes_view(data_view, 'ss','qcd_w')
        ss_f1p2_wje_view, ss_p1f2_wje_view, ss_f1f2_wje_view = self.make_fakes_view(data_view, 'ss','wjet_w')

        ss_fakes_1   = ss_f1p2_wje_view  #MedianView(lowv=ss_f1p2_qcd_view, highv=ss_f1p2_wje_view)
        ss_fakes_2   = ss_p1f2_wje_view  #MedianView(lowv=ss_p1f2_qcd_view, highv=ss_p1f2_wje_view)
        ss_fakes_12  = ss_f1f2_wje_view  #MedianView(lowv=ss_f1f2_qcd_view, highv=ss_f1f2_wje_view)
        ss_fakes_est = views.SumView(ss_fakes_1, ss_fakes_2, ss_fakes_12)
        ss_fakes_est = views.TitleView( views.StyleView(ss_fakes_est, **data_styles['Zjets*']), 'Fakes;%s' % xaxis)

        os_flip_est_up  = views.SubdirectoryView( data_view, 'os/p1p2/charge_weightSysUp')
        os_flip_est     = views.SubdirectoryView( data_view, 'os/p1p2/charge_weight')
        #os_flip_est     = MedianView(highv=os_flip_est_up, centv=os_flip_est)
        os_flip_est_nofake = os_flip_est#views.SumView(os_flip_est, neg_os_fakes)
        os_flip_est_nofake = views.TitleView( views.StyleView(os_flip_est_nofake, **data_styles['WZ*']), 'charge-fakes;%s' % xaxis)
        return ss_p1p2_view, ss_fakes_est, os_flip_est_nofake
            
    def make_charge_flip_control_plot(self, variable, xaxis='', rebin=1, legend_on_the_left=False, 
                                      data_type='data', x_range=None, apply_scale='', show_ratio=False,
                                      differential=False):
        ss_p1p2_view, ss_fakes_est, os_flip_est_nofake = self.get_flip_data(rebin, xaxis, data_type)
        
        if differential:
            ss_p1p2_view       = DifferentialView(ss_p1p2_view      )
            ss_fakes_est       = DifferentialView(ss_fakes_est      )
            os_flip_est_nofake = DifferentialView(os_flip_est_nofake)

        fakes_hist = ss_fakes_est.Get(variable)
        flip_hist  = os_flip_est_nofake.Get(variable)
        if apply_scale:
            flip_hist = MedianView.apply_view(flip_hist, os_flip_est_nofake.Get(variable+apply_scale))
        
        obs_hist       = ss_p1p2_view.Get(variable)
        estimate_hist  = plotting.HistStack()
        estimate_hist.Add(fakes_hist)
        estimate_hist.Add(flip_hist)

        estimate_error = HistStackToTGRaphErrors( estimate_hist )
        estimate_error.SetFillStyle(3013)
        estimate_error.SetFillColor(ROOT.EColor.kBlack)
        estimate_error.SetTitle('Error on estimate')

        #from pdb import set_trace; set_trace()
        sum_stack = sum(estimate_hist.hists)
        print "variable %s: data integral: %.1f (%.1f/%.1f), estimate: %.1f (%.1f/%.1f) (under/overflow)" % (variable, \
            obs_hist.Integral(), obs_hist.GetBinContent(0), obs_hist.GetBinContent(obs_hist.GetNbinsX()+1), \
            sum_stack.Integral(), sum_stack.GetBinContent(0), sum_stack.GetBinContent(sum_stack.GetNbinsX()+1) )
        hmax = max( [ estimate_hist.GetMaximum(), max(list(obs_hist)) ] )
        obs_hist.GetYaxis().SetRangeUser(0,hmax*1.3)
        if x_range:
            obs_hist.GetXaxis().SetRangeUser(x_range[0],x_range[1])
        
        obs_hist.Draw()
        estimate_hist.Draw('same')
        self.canvas.Update()
        estimate_error.Draw('2 same')
        obs_hist.Draw('same')
        self.keep.extend([
            estimate_hist,
            estimate_error,
            obs_hist
            ])

        legend = self.add_legend([obs_hist], leftside=legend_on_the_left, entries=4)
        legend.AddEntry(estimate_hist,'f')
        #legend.AddEntry(estimate_error,'f')        
        legend.Draw()
        self.add_cms_blurb(self.sqrts)
        if show_ratio:
            self.add_ratio_plot(obs_hist, estimate_hist, x_range, ratio_range=0.2)

    
    def plot_zee_control(self, variable, xaxis='', rebin=1, legend_on_the_left=False, 
                         x_range=None, show_ratio=False, logscale=False):
        data_view = self.get_view('data')
        data_view = self.rebin_view(data_view, rebin) if rebin != 1 else data_view
        mc_views  = [self.get_view(i) for i in [ 'ZZ*', 'WZ*', 'WW*', 'TT*', 'Zjets_M50' ] ]
        if rebin != 1:
            mc_views = [self.rebin_view(i, rebin) for i in mc_views]
        
        zee_data  = views.SubdirectoryView( data_view, 'os/p1p2/')
        zee_mcs   = [ views.SubdirectoryView( i, 'os/p1p2/') for i in mc_views]
        
        #os_f1p2_qcd_view, os_p1f2_qcd_view, os_f1f2_qcd_view = self.make_fakes_view(data_view, 'os','qcd_w')
        os_f1p2_wje_view, os_p1f2_wje_view, os_f1f2_wje_view = self.make_fakes_view(data_view, 'os','wjet_w')
        os_fakes_1   = os_f1p2_wje_view  #MedianView(lowv=os_f1p2_qcd_view, highv=os_f1p2_wje_view)
        os_fakes_2   = os_p1f2_wje_view  #MedianView(lowv=os_p1f2_qcd_view, highv=os_p1f2_wje_view)
        os_fakes_12  = os_f1f2_wje_view  #MedianView(lowv=os_f1f2_qcd_view, highv=os_f1f2_wje_view)
        os_fakes_est = views.SumView(os_fakes_1, os_fakes_2, os_fakes_12)
        os_fakes_est = views.TitleView( views.StyleView(os_fakes_est, **data_styles['WplusJets*']), 'Fakes;%s' % xaxis)

        zee_mcs = zee_mcs[:-1]+[os_fakes_est]+zee_mcs[-1:]
        events_estimate = views.StackView( *zee_mcs)
        estimate_hist  = events_estimate.Get(variable)
        obs_hist       = zee_data.Get(variable)
        hmax = max( [ estimate_hist.GetMaximum(), max(list(obs_hist)) ] )
        if logscale:
            obs_hist.GetYaxis().SetRangeUser(10**-2,hmax*10**4)
            self.pad.SetLogy(True)
        else:
            obs_hist.GetYaxis().SetRangeUser(0.,hmax*1.3)

        if x_range:
            obs_hist.GetXaxis().SetRangeUser(x_range[0],x_range[1])
        
        obs_hist.Draw()
        estimate_hist.Draw('same')
        obs_hist.Draw('same')
        self.canvas.Update()
        self.keep.extend([
            estimate_hist,
            obs_hist
            ])

        legend = self.add_legend([obs_hist], leftside=legend_on_the_left, entries=len(zee_mcs)+1)
        legend.AddEntry(estimate_hist,'f')
        #legend.AddEntry(estimate_error,'f')        
        legend.Draw()
        if show_ratio:
            self.add_ratio_plot(obs_hist, estimate_hist, x_range, ratio_range=0.2)
        self.add_cms_blurb(self.sqrts)



plotter = ControlZEEPlotter()

###########################################################################
##  CHARGE FLIP CONTROL PLOTS     #########################################
###########################################################################

#Charge flip control plots
plotter.make_charge_flip_control_plot('TrkMass','Tracker Inv Mass (GeV)', 
                                      [40,50,60,70,75,80]+[80+i for i in range(1,21,2)]+[105,110,120,130,140], 
                                      apply_scale='_high', show_ratio=True, differential=True)
plotter.save('EE_Charge_Flip_xcheck_trk_invMass')
plotter.canvas.SaveAs('test.root')
#plotter.save('EE_Charge_Flip_xcheck_trk_invMass_log')

plotter.make_charge_flip_control_plot('TrkMass','Tracker Inv Mass (GeV)',[40, 150],apply_scale='_high', 
                                      show_ratio=True)
plotter.save('EE_Charge_Flip_xcheck_trk_invMass_counting')

plotter.make_charge_flip_control_plot('TrkMass_NOSCALE','Tracker Inv Mass (GeV)',2)
plotter.save('EE_Charge_Flip_xcheck_trk_invMass_NoScale')

plotter.make_charge_flip_control_plot('ePt','electron p_{T}',2,x_range=[0,200])
plotter.save('EE_Charge_Flip_xcheck_ePt')

plotter.make_charge_flip_control_plot('eAbsEta','electron |#eta|',4, legend_on_the_left=True)
plotter.save('EE_Charge_Flip_xcheck_eAbsEta')

plotter.make_charge_flip_control_plot('SCEnergy','Super cluster energy (GeV)',5)
plotter.save('EE_Charge_Flip_xcheck_eSCEnergy')

plotter.make_charge_flip_control_plot('e1Pt','electron p_{T}',2,x_range=[0,200])
plotter.save('EE_Charge_Flip_xcheck_e1Pt')

plotter.make_charge_flip_control_plot('e1AbsEta','electron |#eta|',10, legend_on_the_left=True)
plotter.save('EE_Charge_Flip_xcheck_e1AbsEta')

plotter.make_charge_flip_control_plot('e2Pt','electron p_{T}',2,x_range=[0,200])
plotter.save('EE_Charge_Flip_xcheck_e2Pt')

plotter.make_charge_flip_control_plot('e2AbsEta','electron |#eta|',10, legend_on_the_left=True)
plotter.save('EE_Charge_Flip_xcheck_e2AbsEta')

plotter.make_charge_flip_control_plot('type1_pfMetEt','electron |#eta|',10, legend_on_the_left=True)
plotter.save('EE_Charge_Flip_xcheck_type1_pfMetEt')


#Check same method against DYJets, to see if we screwed the method
print 'cosure tests'
plotter.make_charge_flip_control_plot('TrkMass','Tracker Inv Mass (GeV)',2, data_type='Zjets_M50')
plotter.save('EE_Charge_Flip_closure_trk_invMass')

plotter.make_charge_flip_control_plot('eAbsEta','electron |#eta|',4, legend_on_the_left=True, data_type='Zjets_M50')
plotter.save('EE_Charge_Flip_closure_absEta')

plotter.make_charge_flip_control_plot('ePt','electron p_{T}',2,x_range=[0,200], data_type='Zjets_M50')
plotter.save('EE_Charge_Flip_closure_ePt')

###########################################################################
##  DATA/MC PLOTS                 #########################################
###########################################################################

#plotter.plot_mc_vs_data('os/p1p2/', 'TrkMass', rebin=2, xaxis='m_{#mu#mu} (GeV)', 
#                        leftside=False, show_ratio=True)
#plotter.add_cms_blurb(plotter.sqrts)
#for obj in plotter.keep:
#    if obj.ClassName().startswith('TH1F'):
#        print obj.GetTitle(), obj.Integral()
#    if obj.ClassName().startswith('THStack'):
#        for histo in obj.GetHists():
#            print histo.GetTitle(), histo.Integral()
#
#plotter.save('mass')
#
#
#plotter.plot_mc_vs_data('os/p1p2/', 'TrkMass', rebin=2, xaxis='m_{#mu#mu} (GeV)', 
#                        leftside=False, show_ratio=True)
#plotter.add_cms_blurb(plotter.sqrts)
#plotter.pad.SetLogy(True)
#plotter.save('mass_log')
#
#plotter.plot_mc_vs_data('os/p1p2/', 'TrkMass', rebin=10, xaxis='m_{#mu#mu} (GeV)', 
#                        leftside=False, show_ratio=True)
#plotter.add_cms_blurb(plotter.sqrts)
#plotter.save('mass_rebin')
#
#plotter.plot_mc_vs_data('os/p1p2/', 'e1Pt', leftside=False, show_ratio=True)
#plotter.save('e1Pt')
#plotter.plot_mc_vs_data('os/p1p2/', 'e1Pt', 5, leftside=False, show_ratio=True)
#plotter.save('e1Pt_rebin')
#plotter.plot_mc_vs_data('os/p1p2/', 'e2Pt', leftside=False, show_ratio=True)
#plotter.save('e2Pt')
#plotter.plot_mc_vs_data('os/p1p2/', 'e2Pt', 5, leftside=False, show_ratio=True)
#plotter.save('e2Pt_rebin')
#
#plotter.plot_mc_vs_data('os/p1p2/', 'e1AbsEta', leftside=False, show_ratio=True)
#plotter.save('e1AbsEta')
#plotter.plot_mc_vs_data('os/p1p2/', 'e2AbsEta', leftside=False, show_ratio=True)
#plotter.save('e2AbsEta')
#
#plotter.plot_mc_vs_data('os/p1p2/', 'e1AbsEta', 5, leftside=False, show_ratio=True)
#plotter.save('e1AbsEta_rebin')
#plotter.plot_mc_vs_data('os/p1p2/', 'e2AbsEta', 5, leftside=False, show_ratio=True)
#plotter.save('e2AbsEta_rebin')
#
#
#plotter.plot_mc_vs_data('os/p1p2/', "trig_weight" , 5, leftside=False, show_ratio=True)
#plotter.save('trig_weight')
#
#plotter.plot_mc_vs_data('os/p1p2/', "PU_weight"   , 5, leftside=False, show_ratio=True)
#plotter.save('PU_weight')
#
#plotter.plot_mc_vs_data('os/p1p2/', "idIso_weight", 5, leftside=False, show_ratio=True)
#plotter.save('idIso_weight')


#plotter.plot_mc_vs_data('os/p1p2/', 'nvtx', show_ratio=True)
#plotter.save('nvtx')

###########################################################################
##  DATA/MC+FAKES PLOTS           #########################################
###########################################################################

plotter.plot_zee_control('TrkMass', rebin=2, xaxis='m_{#mu#mu} (GeV)', 
                        legend_on_the_left=False, show_ratio=True, logscale=True)
plotter.save('mass_wfakes_log')

#print data_styles.keys()
plotter.plot_zee_control('TrkMass', rebin=2, xaxis='m_{#mu#mu} (GeV)', 
                         legend_on_the_left=False, show_ratio=True)
plotter.save('mass_wfakes')

plotter.plot_zee_control('e1Pt', rebin=5, xaxis='p_{T}^{#mu1} (GeV)', 
                         legend_on_the_left=False, show_ratio=True, x_range=[0,200])
plotter.save('e1Pt_wfakes')

plotter.plot_zee_control('e2Pt', rebin=5, xaxis='p_{T}^{#mu2} (GeV)', 
                         legend_on_the_left=False, show_ratio=True, x_range=[0,200])
plotter.save('e2Pt_wfakes')

plotter.plot_zee_control('e1AbsEta', rebin=5, xaxis='|#eta|^{#mu1}', 
                         legend_on_the_left=False, show_ratio=True)
plotter.save('e1AbsEta_wfakes')

plotter.plot_zee_control('e2AbsEta', rebin=5, xaxis='|#eta|^{#mu2}', 
                         legend_on_the_left=False, show_ratio=True)
plotter.save('e2AbsEta_wfakes')


###########################################################################
##  SHAPE FILE FOR FIT            #########################################
###########################################################################

flip_shape_file = ROOT.TFile(
    os.path.join(
        plotter.outputdir, 
        'charge_flip_shapes.root'), 
    'RECREATE'
)
shape_dir = flip_shape_file.mkdir('data')
shape_dir.cd()
ss_p1p2_view, ss_fakes_est, os_flip_est_nofake = plotter.get_flip_data(2, '', 'data')

ss_view = SubtractionView(ss_p1p2_view, ss_fakes_est, restrict_positive=True)

h_ss = ss_view.Get('TrkMass_NOSCALE')      
h_os = os_flip_est_nofake.Get('TrkMass_NOSCALE')

h_ss.SetTitle('ss_trkMass')
h_os.SetTitle('os_trkMass')
h_ss.SetName('ss_trkMass')
h_os.SetName('os_trkMass')

h_ss.Write()      
h_os.Write()

flip_shape_file.Close()
