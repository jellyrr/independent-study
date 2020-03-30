# mod_ccn
# version : 
# 20200330 main_v0 : setting target and build structure 

# start   	: 20200118
# update  	: 20200211
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
import time as tm
start = tm.time()
import datetime as dtm
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
start_dtm = dtm.datetime(2017,12,25,18,20,1)
final_dtm = dtm.datetime(2017,12,25,21,20,0)
path_ccn = r'../171225_cal/ccn/'
path_cpc = r'../171225_cal/cpc/'
path_dma = r'../171225_cal/dma/'

## data
read = fc.reader(start_dtm,final_dtm,path_ccn=path_ccn,path_cpc=path_cpc,path_dma=path_dma)
data = read.modi_ccndata_calib()

## plot
cal = fc.calibration(data,date=start_dtm.strftime('%Y/%m/%d'))
cal.S_curve(plot_dc=True)
cal.calib_line()

# measurement




'''
from PIL import Image
Image.open('picture/rea_calib_scurve.png').show()
#'''










#=============================================================================
if __name__ == '__main__':
	final = tm.time()
	time = (final-start)/60.
	print('\nrunning time = ',int(time),' min','{:6.3f}'.format((time-int(time))*60.),' s')