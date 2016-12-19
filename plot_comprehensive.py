from plot_utils import plot_map
import numpy as np
import matplotlib.pyplot as plt

# - Our input datafile
dataFile = 'xsec_scan.dat'
print('Opening file {}.'.format(dataFile) )
# - Load data into a python
data = np.genfromtxt( dataFile, delimiter=None, skip_header=0, skip_footer=0, names=True )
#################################################################

x = data['mA']
y = data['gt'] 
S = data['PS'][0]
s = data['rs'][0]

#################################################################

z = data['xsec']

log_zmin  = np.log10( z.min() )
cb_levels = np.logspace( log_zmin, 1, 20)

plotSettings = {
'figname'          : './figures/xsec.pdf',
'title'            : r'MWTC - $\sigma( p p \rightarrow Z h)$' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$m_{A}$ [GeV]',
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'cb_label'         : r'$\sigma $ [pb]',
'cb_levels'        : cb_levels,
'cb_logscale'      : True,
'contour_levels'   : [0.5, 1.0, 2.0],
'cb_ticks'         : [ 0.01  , 0.1,   0.5,  1.0, 5.0,  10],
'cb_ticklabels'    : ['0.01' ,'0.1', '0.5', '1','5','10'],
'contour_colormap' : plt.cm.Greens,
'cb_colormap'      : plt.cm.gnuplot
               }

plot_map( plotSettings, x, y, z )

#################################################################

z = data['M1N']

cb_levels = np.linspace( 300, 2700, 25)

plotSettings = {
'figname'          : './figures/M1N.pdf',
'title'            : r'MWTC,  $R^{0}_{1}$ mass' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$m_{R^{0}_1}$ [GeV]',
'cb_levels'        : cb_levels,
'contour_levels'   : [1000.0, 1500.0, 2000.0],
'contour_colormap' : plt.cm.Greens_r,
'cb_colormap'      : plt.cm.Blues
               }

plot_map( plotSettings, x, y, z )

#################################################################

z = data['M2N']

cb_levels = np.linspace( 300, 2700, 25)

plotSettings = {
'figname'          : './figures/M2N.pdf',
'title'            : r'MWTC,  $R^{0}_{2}$ mass' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$m_{R^{0}_1}$ [GeV]',
'cb_levels'        : cb_levels,
'contour_levels'   : [1000.0, 1500.0, 2000.0],
'contour_colormap' : plt.cm.Greens_r,
'cb_colormap'      : plt.cm.Blues
               }

plot_map( plotSettings, x, y, z )


#################################################################

z = data['M2N']-data['M1N']

plotSettings = {
'figname'          : './figures/dMN.pdf',
'title'            : r'MWTC, $R^{0}_{2} - R^{0}_1$ mass difference' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$m_{R^{0}_{2}} - m_{R^{0}_1}$ [GeV]',
'cb_levels'        : 'auto',
'contour_levels'   : [50.0, 100.0, 500.0],
'contour_colormap' : plt.cm.Greens_r,
'cb_colormap'      : plt.cm.Blues
               }

plot_map( plotSettings, x, y, z )

#################################################################

z = data['M2C']-data['M1C']

plotSettings = {
'figname'        : './figures/dMC.pdf',
'title'          : r'MWTC, $R^{\pm}_{2} - R^{\pm}_1$ mass difference' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'         : [ 250.0, 2050.0],
'yrange'         : [ 0.5, 9.5 ],
'xlabel'         : r'$m_{A}$ [GeV]',
'ylabel'         : r'$\widetilde{g}$',
'cb_label'       : r'$m_{R^{\pm}_{2}} - m_{R^{\pm}_1}$ [GeV]',
'cb_levels'      : 'auto',
'contour_levels' : [0.5, 1.0, 2.0],
'cb_colormap'    : plt.cm.Blues
               }

plot_map( plotSettings, x, y, z )

#################################################################

z = data['w1N']

plotSettings = {
'figname'          : './figures/Gamma_R1N.pdf',
'title'            : r'MWTC, $\Gamma_{R^{0}_{1}}$' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$\Gamma_{R^{0}_{1}}$ [GeV]',
'cb_levels'        : 'auto',
'contour_levels'   : [0.5, 1.0, 2.0],
'contour_colormap' : plt.cm.Greens,
'cb_colormap'      : plt.cm.Blues
               }

plot_map( plotSettings, x, y, z )

#################################################################

z = data['w2N']

log_zmin  = np.log10( z.min() )
log_zmax  = np.log10( z.max() )
cb_levels = np.logspace( log_zmin, log_zmax, 20 )

plotSettings = {
'figname'          : './figures/Gamma_R2N.pdf',
'title'            : r'MWTC, $\Gamma_{R^{0}_{2}}$' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$\Gamma_{R^{0}_{2}}$ [GeV]',
'cb_levels'        : cb_levels,
'cb_logscale'      : True,
'contour_levels'   : [10.0, 50.0, 100.0],
'contour_colormap' : plt.cm.Greens,
'cb_colormap'      : plt.cm.Blues,
'cb_ticks'         : [  1.0   ,  10.0,    100.0,    1000.0  ],
'cb_ticklabels'    : [ '1.0'  , '10.0',  '100.0',  '1000.0' ]
               }

plot_map( plotSettings, x, y, z )


#################################################################

z = (data['M2N']-data['M1N']) / ( data['w1N'] + data['w2N'] )

plotSettings = {
'figname'          : './figures/dMN_div_Gamma.pdf',
'title'            : r'MWTC, $R^{0}_{2} - R^{0}_1$ mass degeneracy' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$(m_{R^{0}_{2}} - m_{R^{0}_1})/(\Gamma_{R^{0}_{1}} +  \Gamma_{R^{0}_2} )$',
'cb_levels'        : 'auto',
'contour_levels'   : [0.5, 1.0, 2.0],
'contour_colormap' : plt.cm.Greens,
'cb_colormap'      : plt.cm.Blues
               }

plot_map( plotSettings, x, y, z )

#################################################################

z = (data['M2C']-data['M1C']) / ( data['w1C'] + data['w2C'] )

plotSettings = {
'figname'          : './figures/dMC_div_Gamma.pdf',
'title'            : r'MWTC, $R^{\pm}_{2} - R^{\pm}_1$ mass degeneracy' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$(m_{R^{\pm}_{2}} - m_{R^{\pm}_1})/(\Gamma_{R^{\pm}_{1}} +  \Gamma_{R^{\pm}_2} )$',
'cb_levels'        : 'auto',
'contour_levels'   : [0.5, 1.0, 2.0],
'contour_colormap' : plt.cm.Greens,
'cb_colormap'      : plt.cm.Blues
               }

plot_map( plotSettings, x, y, z )

#################################################################

z = (data['M1C']-data['M1N'])

plotSettings = {
'figname'          : './figures/dM1CN.pdf',
'title'            : r'MWTC, $R^{\pm}_{1} - R^{0}_1$ mass difference' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$m_{R^{\pm}_{1}} - m_{R^{0}_1}$',
'cb_levels'        : 'auto',
'contour_levels'   : [-10.0, -5.0, -2.0],
'contour_colormap' : plt.cm.Greens_r,
'cb_colormap'      : plt.cm.Blues_r
               }

plot_map( plotSettings, x, y, z )

#################################################################

z = (data['M2C']-data['M2N'])

plotSettings = {
'figname'          : './figures/dM2CN.pdf',
'title'            : r'MWTC, $R^{\pm}_{2} - R^{0}_2$ mass difference' + '\n' + r'S = {:.1f}, s = {:.1f}'.format(S, s),
'xrange'           : [ 250.0, 2050.0],
'yrange'           : [ 0.5, 9.5 ],
'xlabel'           : r'$m_{A}$ [GeV]',
'ylabel'           : r'$\widetilde{g}$',
'cb_label'         : r'$m_{R^{\pm}_{2}} - m_{R^{0}_2}$',
'cb_levels'        : 'auto',
'contour_colormap' : plt.cm.Greens_r,
'contour_levels'   : [-2.0,  -1.0],
'cb_colormap'      : plt.cm.Blues_r
               }

plot_map( plotSettings, x, y, z )
