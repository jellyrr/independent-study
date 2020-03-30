# func of mod_ccn
# version : 
# 20200330 main_v0   : setting target and build structure
# 20200330 main_v0.2 : add file name identify, nan data with discontinue data, data collect from chosen start time

# target : make a CCNc exclusive func with python 3.6
# 0. read file with discontinue data : CCNc, SMPS, CPC, DMS
# 	method : use time stamp -> test by calibraion data
# 1. calibration and measurement data process:
# 	(1) calibration : deal with time due to the data may be unstable in first 4 min by CCNc,
#		then DMA changes diameter every 2 min, however, first and last 30s 
#		may either be unstable
# 	(2) measurement : different deltatime of different SS table
# 2. calibration : S curve, calibration line
# 3. measurement : Eulia's result

import os
import numpy as n
import scipy.stats as st
import matplotlib.pyplot as pl
from datetime import datetime as dtm
from datetime import timedelta as dtmdt

# file reader
class reader:
	def __init__(self,start,final,**kwarg):
		## set data path parameter
		default = {'path_ccn'  : 'ccn/', 
				   'path_dma'  : 'dma/', 
				   'path_cpc'  : 'cpc/', 
				   'path_smps' : 'smps/'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
		default.update(kwarg)
	
		## set class parameter
		self.start = start ## datetime object
		self.final = final ## datetime object
		self.kil_date  = lambda time: dtm.strptime(dtm.strftime(time,'%X'),'%X')
		self.path_ccn  = default['path_ccn']
		self.path_dma  = default['path_dma']
		self.path_cpc  = default['path_cpc']
		self.path_smps = default['path_smps']
		
	## SMPS
	## change time every 5 min
	## keys index : 'Date' : 1 (%D)
	## 				'Start Time' : 2 (%X)
	## 				'diameter' : 4~110
	def smps_raw(self):
		path, fout, switch= self.path_smps, [], True
		for file in os.listdir(path):
			with open(path+file,errors='replace') as f:
				if switch==True:
					[ f.readline() for i in range(15) ]
					header = f.readline()[:-1].split('\t')
					switch = False
				else: [ f.readline() for i in range(16) ]
				[ fout.append(n.array(line[:-1].split('\t'))) for line in f ]
		return dict(zip(header,n.array(fout).T))

	## CCNc	
	## name : CCN 100 data %y%m%d%M0000.csv
	## change time every second
	## keys index : '    Time' : 0 (%X)
	## 				' CCN Number Conc' : 45
	## 
	def ccn_raw(self):
		path, fout, delT, tm_start = self.path_ccn, [], dtmdt(seconds=1.), self.start
		switch, begin = True, True
		for file in os.listdir(path):
			if '.csv' not in file: continue
			with open(path+file,errors='replace') as f:
				## get header and skip instrument information
				if switch:
					[ f.readline() for i in range(4) ]
					header = f.readline()[:-1].split(',')[:-1]
					f.readline()
					switch = False
				else: [ f.readline() for i in range(6) ]
				## collect data from start time, and make nan array for discontinue data
				for line in f:
					if begin: 
						begin = True if line[0:8] != dtm.strftime(tm_start,'%X') else False
						if begin: continue 

					while line[0:8] != dtm.strftime(tm_start,'%X'):
						fout.append(n.array([dtm.strftime(tm_start,'%X')]+[n.nan]*len(header[1::])))
						tm_start += delT
					
					## data colllect
					fout.append(n.array(line[:-1].split(',')))
					tm_start += delT
		return dict(zip(header,n.array(fout).T))
		 
		 
	
	## DMA
	## name : %Y%m%d.txt
	## change time every second
	## keys index : 'Date' : 0 (%Y/%m/%d)
	## 				'Time' : 1 (%X)
	## 				'Diameter' : 3
	def dma_raw(self):
		path, fout = self.path_dma, []
		header = ['Date','Time','Diameter','SPD']
		for file in os.listdir(path):
			if '.txt' not in file: continue
			with open(path+file,errors='replace') as f:
				[ fout.append(n.array(line[:-1].split('\t'))) for line in f ]

		return dict(zip(header,n.array(fout).T))

	## CPC
	## name : %y%m%d_cpc.csv
	## change time every second
	## keys index : 'Time' : 0 (%X)
	## 				'Concentration (#/cm?)' : 1
	def cpc_raw(self):
		path, fout, switch= self.path_cpc, [], True
		for file in os.listdir(path):
			if '.csv' not in file: continue
			with open(path+file,errors='replace') as f:
				if switch==True:
					[ f.readline() for i in range(17) ]
					header = f.readline()[:-1].split(',')[:-1]
					switch = False
				else: [ f.readline() for i in range(18) ]
				[ fout.append(n.array(line[:-1].split(',')[:-1])) for line in f ]

		return dict(zip(header,n.array(fout[:-5]).T))

	## data for CCNc calibration
	## use data of CCNc, DMA, CPC, and modified by time
	def modi_ccndata_calib(self):
		ccn, dma, cpc = self.ccn_raw(), self.dma_raw(), self.cpc_raw()
		time_ccn, SS_ccn = ccn['    Time'], ccn[' Current SS']
		conc_ccn = ccn[' CCN Number Conc'].astype(float)
		conc_cpc = cpc['Concentration (#/cm�)'].astype(float)
		d_dma = dma['Diameter']
		modi_ccndata = {}
	
		start_index = n.where(time_ccn==self.start.strftime('%X'))[0][0]
		final_index = n.where(time_ccn==self.final.strftime('%X'))[0][0]
		
		## for calibration, SS change every 30 mins, and diameter change every 2 mins,
		## however, CCNc is unstable on first 4 mins, then DMA is unstable 30 seconds 
		## after the start and 30 seconds before the end
		## | ---- || -- || -- |.....| -- |	for ch : | - || -- || - |
		## < uns  >< ch >< ch >.....< ch >	  	   	 <uns>< st ><uns>
		## <            30 min           > 	   	   	 <    2 min     >
		## stable(st), unstable(uns), diameter change(ch)
		## 30min = 1800s, 4min 30s = 270s, 4min 90s = 330s, 2min = 120s
		## (30-4)/2 = 13 different diameter
	
		for i in range(0,len(time_ccn),1800):
			i_ccn = start_index+i ## start time may not equal to CCNc's first data
			if i_ccn>final_index: break
		
			conc_ccn_ave, conc_cpc_ave, diameter = [], [], []
			for num in range(13):
				diameter.append(round(float(d_dma[i+270+120*num])))
				conc_ccn_ave.append(n.mean(conc_ccn[i_ccn+270+120*num:i_ccn+330+120*num]))
				conc_cpc_ave.append(n.mean(conc_cpc[i+270+120*num:i+330+120*num]))
			
			modi_ccndata.setdefault(SS_ccn[i_ccn],[n.array(diameter),n.array(conc_ccn_ave),n.array(conc_cpc_ave)])
		return modi_ccndata

class calibration:
	def __init__(self,data,date,**kwarg):
		## set calculating parameter
		default = {'kappa'	 : .61, 
				   'inst_T'  : 299.15, ## lab temperature [K]
				   'rho_w_T' : 997., 
				   'sig_wa'	 : .072}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## set class parameter 
		self.data = data
		self.size = len(data)
		self.fs = 12
		self.date = date
		self.kappa = default['kappa']
		self.coe_A = 4.*default['sig_wa']/(461.*default['inst_T']*default['rho_w_T'])*1e9 ## diameter [nm]
	
	## activation of CCN may be like a S curve, there is a very sharp activation changing which variety
	## by diameter
	def calib_data(self,SS,get_dc=True,kohler_calib=False):
		d_ = self.data[SS][0] ## diameter
		activation = (self.data[SS][1]/self.data[SS][2])*100. ## CCN/CPC
		if get_dc==True: 
			## data for S curve and regression line of activation changing
			line_bottom = n.where((activation-activation.min())<10.)[0][-1]
			line_top = n.where((100.-activation)<10.)[0][0]
			line_coe = st.linregress(d_[line_bottom:line_top+1],activation[line_bottom:line_top+1])
			line_x = n.array([0.,250.])
			line_y = line_coe[0]*line_x+line_coe[1]
			dc = (50.-line_coe[1])/line_coe[0]

			if kohler_calib==False:
				return {'dia' : d_, 
						'acti' : activation,
						'line_data' : [line_x,line_y],
						'dc' : dc}
			else:
				## data for calibration line, which calibrating by kappa Kohler equation
				SS_calib = lambda dc_ : n.exp((4.*self.coe_A**3./(27.*self.kappa*dc_**3.))**.5)*100.-100.
				return n.array([float(SS),SS_calib(dc)])
		else:
			## data for S curve
			return {'dia' : d_, 'acti' : activation}

	## plot S curve to find out criticle diameter
	def S_curve(self,plot_dc=False,**kwarg):
		## set plotting parameter
		default = {'color'	  	: '#b8dbff', 
				   'line_color' : '#ff5f60', 
				   'mfc'	  	: '#7fdfff', 
				   'mec'	  	: '#297db7', 
				   'order' 		: [ i-1 for i in range(self.size) ], ## order of axes
				   'fig_name' 	: r'picture/rea_calib_scurve.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## parameter
		fs = self.fs
		if self.size==1: fig_row = 1
		else: fig_row = 2
		fig_col = round(self.size/2)
		
		## plot
		fig, ax = pl.subplots(fig_row,fig_col,figsize=(8,6),dpi=200,sharex=True,sharey=True)
		ax = ax.reshape(fig_col*fig_row)
		
		for i,SS in zip(default['order'],self.data.keys()):
			calib_dt = self.calib_data(SS,get_dc=plot_dc) ## if plot dc, should get dc first
			ax[i].plot(calib_dt['dia'],calib_dt['acti'],c=default['color'],
					   marker='o',mfc=default['mfc'],mec=default['mec'],label='S curve')
			if plot_dc==True:
				## plot regression line
				dc = calib_dt['dc']
				ax[i].plot(calib_dt['line_data'][0],calib_dt['line_data'][1],
						   c=default['line_color'],label='regression',ls='--')
				ax[i].plot(dc,50.,c=default['line_color'],ms=5.,marker='o')
				ax[i].text(dc+15.,45.,'$d_c$ = {:.2f} nm'.format(dc),fontsize=fs-2)
						
			ax[i].tick_params(direction='in')
			ax[i].spines['top'].set_visible(False)
			ax[i].spines['right'].set_visible(False)       
			ax[i].set(xlim=(-30.,299.),ylim=(0.,119.))
		
			ax[i].set_title('SS = {:4.2f} %'.format(float(SS)),fontsize=fs-2)
			
		ax[1].legend(framealpha=0.,fontsize=fs-3.,loc=4)
		fig.text(.38,.03,'Dry Particle Diameter (nm)',fontsize=fs)
		fig.text(.04,.4,'Activation (%)',rotation=90.,fontsize=fs)
		## for odd axes
		if ((self.size%2==1)&(self.size!=1)): ax[-1].remove() 

		fig.suptitle('Activation of CCN (Date : {:})'.format(self.date),fontsize=fs+3,style='italic')	
		fig.savefig(default['fig_name'])
		pl.close()

	## use kappa kohler eq. and modified data to calculate the SS, then plotting regression line
	## with CCNc's SS
	def calib_line(self,**kwarg):
		## set plotting parameter
		default = {'color'	    : '#008c69', 
				   'line_color'	: '#ff5f60', 
				   'fig_name'   : r'picture/rea_calib_line.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## parameter
		fs, SS_ = self.fs, []

		## plot
		fig, ax = pl.subplots(figsize=(8,6),dpi=200)
		for SS in self.data.keys():
			SS_.append(self.calib_data(SS,kohler_calib=True))
		SS_ = n.array(SS_).T
		SS_.sort()
		
		line_coe = st.linregress(SS_[0],SS_[1])
		line_x = n.array([SS_[0][0],SS_[0][-1]])
		line_y = line_coe[0]*line_x+line_coe[1]
		
		ax.plot(line_x,line_y, c=default['line_color'],ls='--')
		ax.scatter(SS_[0],SS_[1],c=default['color'],marker='o')
		
		ax.tick_params(which='both',direction='in',labelsize=fs,right=True,top=True)
		ax.set(xlim=(.05,.75),ylim=(.05,.95))
		
		ax.set_xlabel('CCNc SS (%)',fontsize=fs)
		ax.set_ylabel('Calibration SS (%)',fontsize=fs)
		ax.set_title('y = {:.1f}x {:+.3f}'.format(line_coe[0],line_coe[1]),fontsize=fs)
			
		fig.suptitle('Calibration Line of Supersaturation (Date : {:})'.format(self.date),
					 fontsize=fs+3,style='italic')	
		fig.savefig(default['fig_name'])
		pl.close()








