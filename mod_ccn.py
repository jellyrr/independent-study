# mod_ccn
# version : 
# 20200330 main_v0   : setting target and build structure 
# 20200330 main_v0.1 : testing => read calibration data with time stamp
# 20200330 main_v0.2 : testing => read discontinue ccn data

# new	  	: 
# structure : rea_ccn_calib
# target 	: 
# 1. calibration
# 	 (1) S curve fitting
# 	 (2) regression line
# 2. measurement
#	 (1) reproduce Eulian's result

print('\nprogram : mod_ccn')
print('version : \n')
from datetime import datetime as dtm
start = dtm.now().timestamp()
import func as fc


# function class
# reader(start,final,**kwarg)
# default = {'path_ccn'  : 'ccn/', 
#		     'path_dma'  : 'dma/', 
#		     'path_cpc'  : 'cpc/', 
#		     'path_smps' : 'smps/'}
#
## ccn_raw, cpc_raw, smps_raw, dma_raw, modi_ccndata_calib

# calibration(data,date,**kwarg)
# default = {'kappa'   : .61, 
#		  	 'inst_T'  : 299.15, ## lab temperature [K]
#		  	 'rho_w_T' : 997., 
#		  	 'sig_wa'  : .072}
#
## calib_data(SS,get_dc=True,kohler_calib=False)
##
## S_curve(plot_dc=False,**kwarg)
## default = {'color'	   : '#b8dbff', 
##			  'line_color' : '#ff5f60', 
##			  'mfc'	  	   : '#7fdfff', 
##			  'mec'	  	   : '#297db7', 
##			  'order' 	   : [ i-1 for i in range(self.size) ], ## order of axes
##			  'fig_name'   : r'picture/rea_calib_scurve.png'}
##
## calib_line(**kwarg)
## default = {'color'	   : '#008c69', 
##	  	      'line_color' : '#ff5f60', 
##		      'fig_name'   : r'picture/rea_calib_line.png'}

# calibration
## parameter
start_dtm = dtm(2018,12,1,15,0,1)
final_dtm = dtm(2018,12,1,16,59,59)
path_ccn = r'test/msrt/'
path_cpc = r'test/calib/CPC/'
path_dma = r'test/calib/DMA/'

## data
read = fc.reader(start_dtm,final_dtm,path_ccn=path_ccn,path_cpc=path_cpc,path_dma=path_dma)
# data = read.modi_ccndata_calib()
ccn = read.ccn_raw() 
print(ccn)

## plot
# cal = fc.calibration(data,date=start_dtm.strftime('%Y/%m/%d'))
# cal.S_curve(plot_dc=True)
# cal.calib_line()

# measurement




'''
from PIL import Image
Image.open('picture/rea_calib_scurve.png').show()
#'''

 








#=============================================================================
if __name__ == '__main__':
	final = dtm.now().timestamp()
	runtime = (final-start)/60.
	print('\nrunning time = {:3d} min {:6.3f} s'.format(int(runtime),(runtime-int(runtime))*60.))