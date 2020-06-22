# func of mod_ccn
# version :
# 20200330 main_v0   : setting target and build structure
# 20200330 main_v0.2 : add file name identify, nan data with discontinue data, data collect from chosen start time
# 20200331 main_v0.3 : replace time when same time at diff. row
# 20200331 main_v0.4 : alarm code sensor
# 20200401 main_v0.5 : calib_SS : add table on figure
# 20200402 main_v0.6 : calibration : add figure path
# 20200404 main_v0.7 : read smps data completely, add function of checkout time
# 20200405 main_v0.8 : smps2date fig fix out, correct_time_data : recieve kwarg, use object array, save datetime object instead of datetime string
# 20200405 main_v0.9 : correct_time_data : continue data test rewrite (check out next time with tolerence)
# 20200408 main_v1.0 : cpc_raw : use time correct, correct_time_data : deal with blank line adta and unfortunate header, mdfy_data_mesr, plot_smps2date, plot_cpc2date : build completely
# 20200415 main_v1.1 : data from start to final time, use os.path.join to build path, cpc_raw : start time + 1s to fit its true start time, start building calculate kappa function 
# 20200415 main_v1.2 : correct_time_data : from processing cpc_raw to processing correct_time_data, skip start time test, add nan data to current time, change to the time which has first data
# 20200416 main_v1.3 : mdfy_data_mesr : kappa calculate building, read kappa necessary file done(without test), use accumulate array instead of simps integral to calculate the activate diameter | refresh targets
# 20200417 main_v1.4 : mdfy_data_mesr : kappa calculate complete | measurement kappa calculate function building
# 20200420 main_v1.5 : mdfy_data_mesr : calculate Da before output data | measurement : kappa calculated data split building
# 20200511 main_v1.6 : mdfy_data_mesr : rebuilding calculate kappa function (not yet)
# 20200525 re_v1	 : correct_time_data : rebuilding | kappa calculate finish | output data finish | figure xticklabels at 1200 | approximation of kappa
# 20200525 main_v1.7 : merge with re complete
# 20200601 main_v1.8 : add outDt at mdfy_data_calib
# 20200601 re_v1	 : rewrite the file reading by package panda
# 20200614 re_v1.01	 : cpc data reading
# 20200622 re_v1.02	 : all instrument's data reading | cpc_cor(correct coefficient) manual | plot function completely

# target : make a CCNc exclusive func with python 3.6
# 0. read file with discontinue data : CCNc, SMPS, CPC, DMS
# 	method : use time stamp to correct file
# 1. calibration and measurement data process:
# 	(1) calibration : deal with time due to the data may be unstable in first 4 min by CCNc,
#		then DMA changes diameter every 2 min, however, first and last 30s
#		may either be unstable
# 	(2) measurement : different deltatime of different SS table
# 2. calibration : S curve, calibration table
# 3. measurement : cpc2date, smps2date, kappa2date

from os import listdir
from os.path import join as pth
import numpy as n
from scipy.stats import linregress
from matplotlib.pyplot import subplots, close
from datetime import datetime as dtm
from datetime import timedelta as dtmdt
from h5py import File as hpFile, string_dtype
import pandas as pd
import warnings
warnings.simplefilter('error')
warnings.simplefilter('ignore', category=RuntimeWarning)

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
		self.index		= lambda _freq: pd.date_range(start,final,freq=_freq)
		self.path_ccn   = default['path_ccn']
		self.path_dma   = default['path_dma']
		self.path_cpc   = default['path_cpc']
		self.path_smps  = default['path_smps']

	## SMPS
	## change time every 5 min
	def smps_raw(self):
		fList = []
	
		for file in listdir(self.path_smps):
			if '.txt' not in file: continue
			with open(pth(self.path_smps,file),'r',encoding='utf-8',errors='ignore') as f:
				"""
				rawDt : data, use 'Date' and 'Start Time' as index
				fList : all reading file
				"""
				fList.append(pd.read_csv(f,delimiter='\t',skiprows=17,parse_dates=[['Date','Start Time']])
										 .set_index('Date_Start Time').resample('5T').mean())

		return pd.concat(fList).reindex(self.index('5T'))

	## SMPS_others
	## skiprows and key column different
	## change time every 5 min
	def smpsOthers_raw(self):
		fList = []
	
		for file in listdir(self.path_smps):
			if '.txt' not in file: continue
			with open(pth(self.path_smps,file),'r',encoding='utf-8',errors='ignore') as f:
				"""
				rawDt : data, use 'Date' and 'Start Time' as index
				fList : all reading file
				"""
				fList.append(pd.read_csv(f,delimiter='\t',skiprows=15,parse_dates=[['Date','Start Time']])
										 .set_index('Date_Start Time').resample('5T').mean())

		return pd.concat(fList).reindex(self.index('5T'))

	## CCNc
	## name : CCN 100 data %y%m%d%M0000.csv
	## change time every second
	## keys index : '    Time' : 0 (%X)
	## 				' CCN Number Conc' : 45
	def ccn_raw(self):
		fList = []
	
		for file in listdir(self.path_ccn):
			if '.csv' not in file: continue

			stTime = dtm.strptime(file[15:25],'%m%d%H%M%S').replace(year=self.start.year)
			if float(stTime.strftime('%m%d%H'))<float(self.start.strftime('%m%d%H')): continue
			if float(stTime.strftime('%m%d%H%M%S'))>float(self.final.strftime('%m%d%H%M%S')): break
	
			with open(pth(self.path_ccn,file),'r',encoding='utf-8',errors='ignore') as f:
				"""
				rawDt		: data, which index is '    Time'
				stTime 		: first datetime of each file
				finTime		: end datetime of each file 
				indx		: make new index for rawDt
				fList		: all reading file
				"""
				rawDt = pd.read_csv(f,header=3,index_col=['    Time'])
	
				finTime = dtm.strptime(rawDt.index[-1],'%X').replace(year=stTime.year,day=stTime.day,
									   month=stTime.month)
				fList.append(rawDt.set_index(pd.date_range(stTime,finTime,freq='s')))
	
		return pd.concat(fList).reindex(self.index('s'))

	## DMA
	## change time every second
	def dma_raw(self):
		fList = []
	
		for file in listdir(self.path_dma):
			if '.txt' not in file: continue
			with open(pth(self.path_dma,file),'r',encoding='utf-8',errors='ignore') as f:
				"""
				rawDt		: data, which header is named manual
							('Date','Time','Diameter','SPD')
				fList		: all reading file
				"""
				fList.append(pd.read_csv(f,delimiter='\t',parse_dates=[['Date','Time']],usecols=['Date','Time','Diameter'],
										 names=['Date','Time','Diameter','SPD'],na_values=['OK'])
										 .set_index('Date_Time').astype(float).resample('s').mean())
	
		return pd.concat(fList).reindex(self.index('s'))

	## CPC
	## change time every second
	## each part have 28880 datas for most(8 hrs), and rest for 5 mins (default)
	## mesr_time = 480 mins (default)
	def cpc_raw(self,mesr_time=480):
		fList = []
 
		for file in listdir(self.path_cpc):
			if '.csv' not in file: continue
			with open(pth(self.path_cpc,file),'r',encoding='utf-8',errors='ignore') as f:
				"""
				mesr_time	: time period of each part (minutes)
				rawDt		: data, which header is named manual
							('Time','Concentration','Count','Analog 1','Analog 2')
				startDate 	: first date of each file 
				stTimeIndx  : start time in each part of one file, due to the reason under:
							time is continue in same file, but is not necessarily same in different file
				indx		: fit index for rawDt
				fList		: all reading file
				"""
				rawDt = pd.read_csv(f,index_col=['Time'],names=['Time','Concentration','Count',
									'Analog 1','Analog 2'])
	
				stDate = rawDt['Concentration']['Start Date']
				if type(stDate) != str: stDate = stDate[0]
	
				## drop raws which does not have data by check out the type of index
				rawDt.drop([ _row for _row in rawDt.index if ':' not in str(_row) ],inplace=True)
	
				## cpc output file every 8 hr 5 min(485 mins, 28800 data) in each measurement(default)
				## could manager the output time
				periods = len(rawDt)//(mesr_time*60)+1 if len(rawDt)!=(mesr_time*60) else len(rawDt)//(mesr_time*60)
				stTimeIndx = pd.date_range(f'{stDate} {rawDt.index[0]}',freq=f'{int(mesr_time+5)}T',
										   periods=periods).to_pydatetime()
	
				indx = pd.DatetimeIndex([])
				for _st in stTimeIndx:                           
					indx = indx.append(pd.date_range(_st,_st+dtmdt(minutes=mesr_time),closed='left',freq='s'))
				fList.append(rawDt.set_index(indx[:len(rawDt)]).astype(float))
	
		return pd.concat(fList).reindex(self.index('s'))



	## data for CCNc calibration
	## use data of CCNc, DMA, CPC, and modified by time
	def mdfy_data_calib(self):
		ccn, dma, cpc = self.ccn_raw(), self.dma_raw(), self.cpc_raw()
		time_ccn, SS_ccn = ccn['    Time'], ccn[' Current SS']
		conc_ccn = ccn[' CCN Number Conc'].astype(float)
		conc_cpc = cpc['Concentration (#/cm�)'].astype(float)
		d_dma = dma['Diameter']
		mdfy_data_calib = {}

		final_index = n.where(time_ccn==self.final)[0][0]

		## for calibration, SS change every 30 mins, and diameter change every 2 mins,
		## however, CCNc is unstable on first 4 mins, then DMA is unstable 30 seconds
		## after the start and 30 seconds before the end
		## | ---- || -- || -- |.....| -- |	for ch : | - || -- || - |
		## < uns  >< ch >< ch >.....< ch >	  	   	 <uns>< st ><uns>
		## <           30 mins           > 	   	   	 <    2 mins    >
		## stable(st), unstable(uns), diameter change(ch)
		## 30min = 1800s, 4min 30s = 270s, 4min 90s = 330s, 2min = 120s
		## (30-4)/2 = 13 different diameter

		for i in range(0,len(time_ccn),1800):
			if i>final_index: break

			conc_ccn_ave, conc_cpc_ave, diameter = [], [], []
			for num in range(13):
				diameter.append(round(float(d_dma[i+270+120*num])))
				conc_ccn_ave.append(n.nanmean(conc_ccn[i+270+120*num:i+330+120*num]))
				conc_cpc_ave.append(n.nanmean(conc_cpc[i+270+120*num:i+330+120*num]))
			mdfy_data_calib.setdefault(SS_ccn[i],[n.array(diameter),n.array(conc_ccn_ave),n.array(conc_cpc_ave)])
		return mdfy_data_calib

	def mdfy_data_mesr(self,cpc_cor_slope,cpc_cor_inte,smps_data=True,cpc_data=True,kappa_data=True,
					   smpsOther_data=False,calib_SS_date=None,outDt=False):
		## smps data
		if smps_data:
			smps = self.smps_raw()
			smpsData = {'time' 	   : smps.index,
						'bins' 	   : n.array(smps.keys()[6:113]).astype(float),
						'bin_data' : smps[smps.keys()[6:113]]}
		else: smpsData = None

		## smps other data
		if smpsOther_data:
			smps = self.smps_raw()
			smpsData = {'time' 	   : smps.index,
						'bins' 	   : n.array(smps.keys()[2:109]).astype(float),
						'bin_data' : smps[smps.keys()[2:109]]}
		elif not smps_data: smpsData = None

		## cpc data
		if cpc_data:
			cpc = self.cpc_raw()
			cpcData = {'time' : cpc.index,
					   'conc' : (cpc['Concentration']-cpc_cor_inte)/cpc_cor_slope}
		else: cpcData = None

		## kappa data
		## calculate kappa by cpc, ccn, and smps
		## ccn : data/s, cpc : data/s, smps : data/5 mins
		## ratio of activation = nCCN/nCN = ccnConc/cpcConc (5 mins average)
		## | ----- || ----- || ..... || ----- |
		## < 5 min >< 5 min >  .....  < 5 min >
		## neglect first smps data due to no corresponding with average ratio of activation
		## in addition, unstable data occurs on SS change timing, which we could not use these data  
		## 
		## real nCN  = sum(all bin data)
		## real nCCN = nCCN/nCN * real nCN
		## use accumulate function to find when is sum(bin data) = real nCCN

		if kappa_data:
			## get ccn, smps, cpc data and take necessary information
			## ccn
			ccn = self.ccn_raw()
			ccnConc = ccn[' CCN Number Conc']
			ccnSS	= ccn[' Current SS']

			## get smps data, if smps_data = False, get raw data
			## first smps bin data could not correspond to activate ratio, neglect it
			## resample as 10 min
			try:
				time	  = smpsData['time'][1::2]
				smpsBin   = smpsData['bins']
				smpsBinDt = n.array(smpsData['bin_data'][1:].asfreq('10T'))
			except:
				if smpsOther_data:
					smps = self.smpsOthers_raw()
					smpsBin   = n.array(smps.keys()[2:109]).astype(float)
					smpsBinDt = n.array(smps[smps.keys()[2:109]][1:].asfreq('10T'))
					time	  = smps.index[1::2]
				else:
					smps = self.smps_raw()
					smpsBin = n.array(smps.keys()[6:113]).astype(float)
					smpsBinDt = n.array(smps[smps.keys()[6:113]][1:].asfreq('10T'))
					time	  = smps.index[1::2]
	
			## get cpc data, if cpc_data = False, get raw data
			try:
				cpcConc = cpcData['conc']
			except:
				cpc = self.cpc_raw()
				cpcConc = (cpc['Concentration']-cpc_cor_inte)/cpc_cor_slope

			## due to discountinue bin data, each bin diameter means the average diameter in a range,
			## then we use accumulate function from large bin Dp to smaller to check out the bin Dp 
			## which is most similar to Da(activate Dp)
			## 
			##     | 			 |
			## 1.0 | -----		 |     
			##     |	  \		 | --> nCCN/nCN  		Accumulate plot of bin data to Dp(log scale)
			##     |	   \	 |   					
			## 0.5 |	    \	 |   					each data have been standardization by divide
			##     |		 ----+----					real nCN
			##     |		  	 |    \
			## 0.0 |		  	 |     -----
			##     --------------o----------
			## 	    9	    100	  \	   400	 (nm)(log scale)
			## 				       Da
			##
			## unsable data process
			## 5 mins unstable CCN data and 5 mins mean conc. with last SMPS data
			## | --------- || --------- | SMPS
			## | --------- || --------- | CCN
			## <     us    ><     st    >

			## calculate kappa
			wrongDt = ccnConc>cpcConc
			cpcConc[wrongDt] = n.nan
			ccnConc[wrongDt] = n.nan

			## last ccn and cpc data could not use(closed='left'), and label should be same as smps data([:-1])
			## resample to 10 minute
			actRat = (ccnConc/cpcConc).resample('5T',label='right').mean()[:-1].asfreq('10T')
			ccnSS  = ccnSS.resample('5T',label='right').median()[:-1].asfreq('10T')

			if calib_SS_date is not None:
				with hpFile('output.hdf5','r') as f:
					dset = f['test/Calibration']
					calSS = dset[calib_SS_date]
					for _indSS, _calibSS in zip(calSS['ssSet'],calSS['ssCalib']):
						ccnSS.replace(_indSS,_calibSS,inplace=True)
					ssCalib = calSS['ssCalib'][:]

			## calculate activative data
			## accumulate bins data and normalize to 1, then compare with activative ratio, get the miniumum index
			## normalize bins data with miniumum index divide into activative ratio then multiply 
			## original smps bin diameter with miniumum index to get real activative diameter
			def actDaFunc(_bin_dt,_f_act):
				## nan test
				if (n.isnan(_bin_dt).all()|n.isnan(_f_act)): return n.nan
				"""
				smpsBin   : smps bin value, out of this function
				_accu_bin : reversed accumulate bins data
				_act_indx : index of activate diameter
				"""

				_accu_bin = n.add.accumulate(_bin_dt[::-1])[::-1]/_bin_dt.sum()
				_act_indx = n.abs(_accu_bin-_f_act).argmin()

				try:
					if _accu_bin[_act_indx]>_f_act:
						return n.poly1d(n.polyfit(_accu_bin[_act_indx:_act_indx+2],smpsBin[_act_indx:_act_indx+2],1))(_f_act)
					else:
						return n.poly1d(n.polyfit(_accu_bin[_act_indx-1:_act_indx+1],smpsBin[_act_indx-1:_act_indx+1],1))(_f_act)
				except:
					return _accu_bin[_act_indx]

			actDaDt = n.array([ actDaFunc(bin_dt,f_act) for bin_dt, f_act in zip(smpsBinDt,actRat) ])
			actDaDt[(actDaDt<smpsBin[0])|(actDaDt>smpsBin[-1])] = n.nan
			
			## calculate kappa
			## use simplify if SS < 0.3
			## parameter
			coeA  = 4.*.072/(461.*299.15*997.)*1e9
			limit = ccnSS<.3 ## use primitive function

			## primitive function
			## calculate critical S
			dEq  = lambda _da : 10**n.arange(n.log10(_da),4.+.0001,.0001) ## nm
			criS = lambda _da, _kappa : n.max(((dEq(_da)**3.-_da**3)/(dEq(_da)**3.-_da**3+_kappa*_da**3))*n.exp(coeA/dEq(_da)))

			## kappa approximation
			def kappaApprox(_ss,_da):
				if (n.isnan(_ss)|n.isnan(_da)): return n.nan
				iniKappa = n.arange(.01,1.31,.01) ## a.u.
			
				## approximation function
				## precision from .01 -> .001 -> .0001
				def approxFunc(_precision,kappa_list,output=False):
					for indx, _k in enumerate(kappa_list):
						ssApprox = (criS(_da,_k)-1.)*100.
						if ssApprox < _ss:
							if ((ssApprox<0.)|(indx==0)): return 0. if output else iniKappa
							return n.mean(kappa_list[indx-1:indx+1]) if output else n.arange(kappa_list[indx-1],_k+_precision,_precision)
				try:
					return round(approxFunc(None,approxFunc(.0001,approxFunc(.001,iniKappa)),output=True),5)
				except:
					return n.nan

			## calculate kappa
			## simply function
			kappa = 4.*coeA**3./(27.*actDaDt**3.)/(n.log(ccnSS/100.+1.)**2)
			
			## primitive
			kappa[limit] = n.array([ _k for _k in map(kappaApprox,ccnSS[limit],actDaDt[limit]) ])

			kappaData = {'time' 	  : time,
						 'act_dia'    : actDaDt,
						 'kappa'	  : kappa,
						 'calSS'	  : ssCalib,
						 'mesrSS'	  : ccnSS}

		else: kappaData = None

		mdfy_data_mesr = {'smps'  : smpsData, 
						  'cpc'   : cpcData, 
						  'kappa' : kappaData}

		## save output data
		if (outDt&smps_data&cpc_data&kappa_data):
			with hpFile('output.hdf5','r+') as f:
				print('it is test program, plz change to other key in real')
				mesr = f['test/Measurement']

				dset = mesr.require_group(self.start.strftime('%Y%m%d'))
				dset.clear()

				## save data
				for namDt, procsDt in mdfy_data_mesr.items():
					_nam = dset.create_group(namDt)

					for _key, _val in procsDt.items():
						if 'time' in _key:
							_dtm = _nam.create_dataset(_key,(len(_val),),dtype=string_dtype())
							_dtm[:] = _val
						else:
							_nam[_key] = _val

		return mdfy_data_mesr

# calibration output
class calibration:
	def __init__(self,data,start,output=False,**kwarg):
		## set calculating parameter
		default = {'kappa'	  : .61,
				   'inst_T'   : 299.15, ## lab temperature [K]
				   'rho_w_T'  : 997.,
				   'fig_path' : './',
				   'sig_wa'	  : .072}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## set class parameter
		self.output = output
		self.data = data
		self.size = len(data)
		self.fs	  = 12
		self.start = start
		self.kappa = default['kappa']
		self.coe_A = 4.*default['sig_wa']/(461.*default['inst_T']*default['rho_w_T'])*1e9 ## diameter [nm]
		self.figPath = default['fig_path']

	## activation of CCN may be like a S curve, there is a very sharp activation changing which variety
	## by diameter
	def calib_data(self,SS,get_dc=True,kohler_calib=False):
		d_ = self.data[SS][0] ## diameter
		activation = (self.data[SS][1]/self.data[SS][2])*100. ## CCN/CPC
		if get_dc:
			## data for S curve and regression line of activation changing
			line_bottom = n.where((activation-activation.min())<10.)[0][-1]
			line_top = n.where((100.-activation)<10.)[0][0]
			line_coe = linregress(d_[line_bottom:line_top+1],activation[line_bottom:line_top+1])
			line_x = n.array([0.,250.])
			line_y = line_coe[0]*line_x+line_coe[1]
			dc = (50.-line_coe[1])/line_coe[0]

			if kohler_calib==False:
				return {'dia' : d_,
						'acti' : activation,
						'line_data' : [line_x,line_y],
						'dc' : dc}
			else:
				## data for calibration SS, which calibrating by kappa Kohler equation
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
				   'fig_name' 	: r'calib_Scurve.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## parameter
		fs = self.fs
		if self.size==1: fig_row = 1
		else: fig_row = 2
		fig_col = round(self.size/2)

		## plot
		fig, ax = subplots(fig_row,fig_col,figsize=(8,6),dpi=200,sharex=True,sharey=True)
		ax = ax.reshape(fig_col*fig_row)

		for i, SS in zip(default['order'],self.data.keys()):
			calib_dt = self.calib_data(SS,get_dc=plot_dc) ## if plot dc, should get dc first
			ax[i].plot(calib_dt['dia'],calib_dt['acti'],c=default['color'],
					   marker='o',mfc=default['mfc'],mec=default['mec'],label='S curve')
			if plot_dc:
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

		fig.suptitle('Activation of CCN (Date : {:})'.format(self.start.strftime('%Y/%m/%d')),fontsize=fs+3,style='italic')
		fig.savefig(pth(self.figPath,default['fig_name']))
		close()

	## use kappa kohler eq. and modified data to calculate the SS, then plotting regression line
	## with CCNc's SS
	def calib_SS(self,**kwarg):
		## set plotting parameter
		default = {'mec'	    : '#00b386',
				   'mfc'		: '#ffffff',
				   'ms'			: 6.,
				   'mew'		: 1.8,
				   'fig_name'   : r'calib_CalibTable.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## parameter
		fs, SS_ = self.fs, []
		height = .15

		## plot
		fig, ax = subplots(figsize=(8,6),dpi=200)
		box = ax.get_position()
		ax.set_position([box.x0,box.y0+.03,box.width,box.height*.95])

		for SS in self.data.keys():
			SS_.append(self.calib_data(SS,kohler_calib=True))
		SS_ = n.array(SS_).T
		SS_.sort()

		## save calibration SS to file
		if self.output:
			with hpFile('output.hdf5','r+') as f:
				dset = f['test/Calibration']
				try:
					cur_cal = dset.create_group(self.start.strftime('%Y%m%d'))
				except:
					cur_cal = dset[self.start.strftime('%Y%m%d')]
					cur_cal.clear()
				cur_cal['ssSet'] = SS_[0]
				cur_cal['ssCalib'] = SS_[1]

		## plot 6 calibration data, to align with table, x data use index data
		ax.plot(n.arange(.1,.7,.1),SS_[1],mec=default['mec'],mew=default['mew'],mfc=default['mfc'],ms=default['ms'],marker='o',ls='')
		table = ax.table(cellText=[list(n.round(SS_[1],3))],rowLabels=['Calib SS'],colLabels=list(SS_[0]),
						 cellLoc='center',rowLoc='center',colLoc='center',
						 bbox=[0.,-height,1.,height])
		table.auto_set_font_size(False)
		table.set_fontsize(self.fs)

		ax.tick_params(which='both',direction='in',labelsize=fs,right=True,top=True,labelbottom=False)
		ax.set(xlim=(.05,.65),ylim=(.05,.95))

		# ax.set_xlabel('CCNc SS (%)',fontsize=fs)
		ax.set_ylabel('Calibration SS (%)',fontsize=fs)
		ax.set_title('CCNc SS (%)',fontsize=fs)

		fig.suptitle('Calibration Table of Supersaturation (Date : {:})'.format(self.start.strftime('%Y/%m/%d')),
					 fontsize=fs+3,style='italic')
		fig.savefig(pth(self.figPath,default['fig_name']))
		close()

# measurement output
class measurement:
	def __init__(self,start,final,data=None,dtDate=None,**kwarg):
		## set calculating parameter
		default = {'fig_path' : './'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## func
		index = lambda _freq: pd.date_range(start,final,freq=_freq)
		def ticks(split_hr):
			_tick = index(f'{split_hr}h')
			return _tick, _tick.strftime('%m-%d%n%X')

		## get data
		if dtDate is not None:
			data = {}
			with hpFile('output.hdf5','r') as f:
				for _nam, _dt in f[f'test/Measurement/{dtDate}'].items():
					_data = {}
					for _key, _val in _dt.items():
						if 'time' in _key:
							_data.setdefault(_key,_val[:].astype(n.datetime64).astype(dtm))
						else:
							_data.setdefault(_key,_val[:])
					data.setdefault(_nam,_data)
		elif data is None:
			raise ValueError("Input data to 'data' variable or 'useOutput = True' then set 'dtDate'")
		
		self.smpsData  = data['smps']
		self.cpcData   = data['cpc']
		self.kappaData = data['kappa']

		## set class parameter
		self.fs = 13.
		self.ticks = ticks
		self.index = index
		self.start = start
		self.final = final
		self.figPath = default['fig_path']

	## plot smps with date and set log scale data
	def plot_smps2date(self,plotTogether=None,**kwarg):
		from matplotlib.colors import LogNorm
		## set plot parameter
		default = {'splt_hr'  : 12,
				   'cmap'	  : 'jet',
				   'fig_name' : r'mesr_smps2date.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## set plot variable
		if self.smpsData is not None: smps = self.smpsData 
		else: raise ValueError('SMPS data is None !!!')

		## data set
		time = self.index('5T')
		data = smps['bin_data'].reindex(time).T
		data.replace(0.,1e-5,inplace=True)
		data.replace(n.nan,0.,inplace=True)

		xTick, xTickLab = self.ticks(default['splt_hr'])
		fs = self.fs

		## plot
		if plotTogether is None:
			fig, ax = subplots(figsize=(10.,6.),dpi=150.)
		else:
			fig, ax = plotTogether

		## plot data and necessary setting
		pm = ax.pcolormesh(time,smps['bins'],data,cmap=default['cmap'],norm=LogNorm(vmin=10**.5,vmax=10**5.3))

		box = ax.get_position()
		ax.set_position([box.x0,box.y0+0.02,box.width,box.height])
		cax = fig.add_axes([.92,box.y0+0.02,.015,box.height])

		cb = fig.colorbar(pm,cax=cax)

		ax.tick_params(which='major',length=6.,labelsize=fs-2.)
		ax.tick_params(which='minor',length=3.5)
		cb.ax.tick_params(which='major',length=5.,labelsize=fs-2.)
		cb.ax.tick_params(which='minor',length=2.5)
		ax.set(yscale='log',xticks=xTick)

		## single plot
		if plotTogether is None:
			cb.ax.set_title('number conc.\n(#/$cm^3$/$\Delta log D_p$)',fontsize=fs-2.)
			ax.set_ylabel('Electric modify diameter (nm)',fontsize=fs)
			ax.set_xlabel(f"Time({self.start.year})",fontsize=fs)
			ax.set_xticklabels(xTickLab)

			fig.suptitle(f"SMPS data from ({self.start.strftime('%Y/%m/%d %X')}) to ({self.final.strftime('%Y/%m/%d %X')})",fontsize=fs+2.,style='italic')
			fig.savefig(pth(self.figPath,default['fig_name']))
			close()

		## together plot
		else:
			ax.set_ylabel('Diameter (nm)',fontsize=fs-1.)
			cb.ax.set_title('number conc.\n(#/$cm^3$/$\Delta log D_p$)',fontsize=fs-3.)
			ax.set_xticklabels('')
			ax.set_title(r'SMPS data',fontsize=fs-1.)

			return xTickLab

	## plot cpc with date
	def plot_cpc2date(self,plotTogether=None,**kwarg):
		## set plot parameter
		default = {'splt_hr'  : 12,
				   'fig_name' : r'mesr_cpc2date.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## set plot variable
		if self.cpcData is not None: cpc = self.cpcData 
		else: raise ValueError('CPC data is None !!!')

		## data set
		time = self.index('s')
		data = cpc['conc'].reindex(time)/1e3

		xTick, xTickLab = self.ticks(default['splt_hr'])
		fs = self.fs

		## plot
		if plotTogether is None:
			fig, ax = subplots(figsize=(10.,6.),dpi=150.)
		else:
			fig, ax = plotTogether
		box = ax.get_position()

		## plot data and necessary setting
		ax.plot(time,data,c='#7396ff',lw=1.5)

		ax.tick_params(which='major',direction='in',length=7,labelright=True,right=True,labelsize=fs-2.5)
		ax.set(ylim=(-1.,35.),xticks=xTick)

		## single plot
		if plotTogether is None:
			ax.set_ylabel(r'Number Conc. (#$\times 10^3 /cm^3$)',fontsize=fs)
			ax.set_xlabel(f"Time({self.start.strftime('%Y')})",fontsize=fs)
			ax.set_xticklabels(xTickLab)

			fig.suptitle(f"CPC data from ({self.start.strftime('%Y/%m/%d %X')}) to ({self.final.strftime('%Y/%m/%d %X')})",fontsize=fs+2.,style='italic')
			fig.savefig(pth(self.figPath,default['fig_name']))
			close()

		## together plot
		else:
			ax.set_ylabel(r'Conc. (#$\times 10^3 /cm^3$)',fontsize=fs-1.)
			ax.set_xticklabels('')
			ax.set_title(r'CPC data',fontsize=fs-1.)

			return xTickLab

	## plot cpc with date
	def plot_kappa2date(self,plotTogether=None,**kwarg):
		## set plot parameter
		default = {'splt_hr'  : 12,
				   'fig_name' : r'mesr_kappa2date.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## set plot variable
		if self.kappaData is not None: kappa = self.kappaData 
		else: raise ValueError('Kappa data is None !!!')

		time = self.index('5T')
		data = kappa['kappa'].reindex(time)
		ss   = kappa['mesrSS'].reindex(time)

		xTick, xTickLab = self.ticks(default['splt_hr'])
		fs = self.fs

		color = ['#7617FF','#FF4924','#0A82FF','#efe934','#17FFA4','#FF9C2F']
		
		## plot
		if plotTogether is None:
			fig, ax = subplots(figsize=(10.,6.),dpi=150.)
		else:
			fig, ax = plotTogether
		collectArt = []
		## plot data and necessary setting
		for _ss, _colr in zip(kappa['calSS'],color):
			artist, = ax.plot(time[ss==_ss],data[ss==_ss],mec=_colr,ls='',ms=5.,
							  mfc='#ffffff',mew=2.,marker='o',label=round(_ss,3))
			collectArt.append(artist)

		ax.tick_params(which='major',direction='in',length=7,right=True,labelsize=fs-2.5)
		ax.set(ylim=(-.05,1.),xticks=xTick,xlim=(xTick[0],xTick[-1]))

		ax.set_ylabel(r'$\kappa$ (a.u.)',fontsize=fs)

		## single plot
		if plotTogether is None:
			ax.set_xlabel(f"Time({self.start.strftime('%Y')})",fontsize=fs)
			ax.set_xticklabels(xTickLab)

			ax.legend(framealpha=0,fontsize=fs-2.,title='SS (%)',title_fontsize=fs-3.)
			fig.suptitle(f"$\kappa$ from ({self.start.strftime('%Y/%m/%d %X')}) to ({self.final.strftime('%Y/%m/%d %X')})",fontsize=fs+2.,style='italic')
			fig.savefig(pth(self.figPath,default['fig_name']))
			close()
		
		## together plot
		else:
			box = ax.get_position()
			ax.set_position([box.x0,box.y0,box.width,box.height])
			cax = fig.add_axes([.935,box.y0*1.05,.015,box.height])
			cax.set_axis_off()

			cax.legend(handles=collectArt,framealpha=0,fontsize=fs-2.,title='SS (%)',title_fontsize=fs-3.,
					   loc=10)
			ax.set_xticklabels('')
			ax.set_title(r'$\kappa$',fontsize=fs-1.)

			return xTickLab

	def plot_together(self,**kwarg):
		## set plot parameter
		default = {'splt_hr'  : 12,
				   'order'    : ['smps','cpc','kappa'],
				   'fig_name' : 'mesr_together2date.png'}

		plotFunc = {'smps'  : self.plot_smps2date,
					'cpc'   : self.plot_cpc2date,
					'kappa' : self.plot_kappa2date}

		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		fig, axes = subplots(3,1,figsize=(10.,6.),dpi=150.)

		for _nam, ax in zip(default['order'],axes):
			xTickLab = plotFunc[_nam](plotTogether=(fig,ax),**default)

		ax.set_xticklabels(xTickLab)
		fig.suptitle(f"Date from ({self.start.strftime('%Y/%m/%d %X')}) to ({self.final.strftime('%Y/%m/%d %X')})",fontsize=self.fs+2.,style='italic')

		fig.savefig(pth(self.figPath,default['fig_name']))



