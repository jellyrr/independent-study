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
# 20200623 re_v1.03	 : nan value function | splt_hr cahange to figSet | set low_memory to read_csv of cpc_raw 
# 20200705 re_v1.04	 : read and write hdf5 file(mdfy_data_mesr) | detect path exist or not, then create it | change mesr_data to dataframe and store
# 20200711 main_v1.9 : merge from re | calibration rewrite | mdfy_data_calib | save calibration data in hdf5
# 20200711 main_v2.0 : mdfy_data_calib : ccn=cpc where ccn>cpc, calculating activation reaplace mean by max | cpc_raw : input parameter mesr_time remove, then use self parameter
# 20200716 main_v2.1 : mdfy_data_calib, mdfy_data_mesr : add new storage parameter and calculate activation after mean | build new class : raw_data 

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

from os import listdir, mkdir
from os.path import join as pth, exists
import numpy as n
from datetime import datetime as dtm
from datetime import timedelta as dtmdt
import pandas as pd
from matplotlib.pyplot import subplots, close
import warnings
warnings.simplefilter('ignore', category=RuntimeWarning)

# file reader
class reader:
	def __init__(self,start,final,**kwarg):
		print('\n'+'='*50)
		print(f"Reading file and process data")
		## set data path parameter
		default = {'path_ccn'  : 'ccn/',
				   'path_dma'  : 'dma/',
				   'path_cpc'  : 'cpc/',
				   'path_smps' : 'smps/',
				   'path_output' : './',
				   'cpc_mesr_tm' : 480,
				   'cpc_nan' 	 : False,
				   'smps_nan' 	 : False,
				   'ccn_nan'	 : False,
				   'dma_nan'	 : False}
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
		self.data_prcs	= {'cpc_mesr_tm' : default['cpc_mesr_tm'],
						   'cpc_nan' 	 : default['cpc_nan'],
						   'smps_nan' 	 : default['smps_nan'],
						   'ccn_nan'	 : default['ccn_nan'],
						   'dma_nan'	 : default['dma_nan']}
		self.path_output = default['path_output']
		print(f" from {start.strftime('%Y-%m-%d %X')} to {final.strftime('%Y-%m-%d %X')}")
		print('='*50)
		print(f"{dtm.now().strftime('%m/%d %X')}")

	## SMPS
	## change time every 5 min
	def smps_raw(self):
		print(f"\n	{dtm.now().strftime('%m/%d %X')} : Reading file of SMPS")
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

		fout = pd.concat(fList).reindex(self.index('5T'))

		if self.data_prcs['smps_nan']:
			for nanRang in self.data_prcs['smps_nan']:
				fout.loc[pd.date_range(nanRang[0],nanRang[1],freq='5T')] = n.nan

		return fout

	## SMPS_others
	## skiprows and key column different
	## change time every 5 min
	def smpsOthers_raw(self):
		print(f"\n	{dtm.now().strftime('%m/%d %X')} : Reading file of SMPS")
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

		fout = pd.concat(fList).reindex(self.index('5T'))

		if self.data_prcs['smps_nan']:
			for nanRang in self.data_prcs['smps_nan']:
				fout.loc[pd.date_range(nanRang[0],nanRang[1],freq='5T')] = n.nan

		return fout

	## CCNc
	## name : CCN 100 data %y%m%d%M0000.csv
	## change time every second
	## keys index : '    Time' : 0 (%X)
	## 				' CCN Number Conc' : 45
	def ccn_raw(self):
		print(f"\n	{dtm.now().strftime('%m/%d %X')} : Reading file of CCN")
		fList = []
	
		for file in listdir(self.path_ccn):
			if '.csv' not in file: continue

			stTime = dtm.strptime(file[15:25],'%m%d%H%M%S').replace(year=self.start.year)
			if float(stTime.strftime('%m%d%H'))<float(self.start.strftime('%m%d%H')): continue
			if float(stTime.strftime('%m%d%H%M%S'))>float(self.final.strftime('%m%d%H%M%S')): break
	
			with open(pth(self.path_ccn,file),'r',encoding='utf-8',errors='ignore') as f:
				"""
				rawDt	 : data, which index is '    Time'
				stTime	 : first datetime of each file
				finTime  : end datetime of each file 
				indx	 : make new index for rawDt
				fList	 : all reading file
				"""
				rawDt = pd.read_csv(f,header=3,index_col=['    Time'])
	
				finTime = dtm.strptime(rawDt.index[-1],'%X').replace(year=stTime.year,day=stTime.day,
									   month=stTime.month)
				fList.append(rawDt.set_index(pd.date_range(stTime,finTime,freq='s')))
	
		fout = pd.concat(fList).reindex(self.index('s'))

		if self.data_prcs['ccn_nan']:
			for nanRang in self.data_prcs['ccn_nan']:
				fout.loc[pd.date_range(nanRang[0],nanRang[1],freq='s')] = n.nan
		fout.loc[fout[' Alarm Code']!=0.,' CCN Number Conc'] = n.nan ## set nan data if alarm code != 0.0

		return fout

	## DMA
	## change time every second
	def dma_raw(self):
		print(f"\n	{dtm.now().strftime('%m/%d %X')} : Reading file of DMA")
		fList = []
	
		for file in listdir(self.path_dma):
			if '.txt' not in file: continue
			with open(pth(self.path_dma,file),'r',encoding='utf-8',errors='ignore') as f:
				"""
				rawDt : data, which header is named manual
						('Date','Time','Diameter','SPD')
				fList : all reading file
				"""
				fList.append(pd.read_csv(f,delimiter='\t',parse_dates=[['Date','Time']],usecols=['Date','Time','Diameter'],
										 names=['Date','Time','Diameter','SPD'],na_values=['OK'])
										 .set_index('Date_Time').astype(float).resample('s').mean())
	
		fout = pd.concat(fList).reindex(self.index('s'))

		if self.data_prcs['dma_nan']:
			for nanRang in self.data_prcs['dma_nan']:
				fout.loc[pd.date_range(nanRang[0],nanRang[1],freq='s')] = n.nan

		return fout

	## CPC
	## change time every second
	## each part have 28880 datas for most(8 hrs), and rest for 5 mins (default)
	## mesr_time = 480 mins (default)
	def cpc_raw(self):
		print(f"\n	{dtm.now().strftime('%m/%d %X')} : Reading file of CPC")
		fList = []
		mesr_time = self.data_prcs['cpc_mesr_tm']
 
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
									'Analog 1','Analog 2'],low_memory=False)
	
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
	
		fout = pd.concat(fList).reindex(self.index('s'))

		if self.data_prcs['cpc_nan']:
			for nanRang in self.data_prcs['cpc_nan']:
				fout.loc[pd.date_range(nanRang[0],nanRang[1],freq='s')] = n.nan

		return fout

	## data for CCNc calibration SS
	## use data of CCNc, DMA, CPC, and modified by time
	def mdfy_data_calib(self,outDt=False,**kwarg):
		print('\n'+'-'*50)
		print(f"{dtm.now().strftime('%m/%d %X')} : Modify calibration data")
		print('-'*50)
		from scipy.interpolate import interp1d
		## parameter
		default = {'kappa'		  : .61,
				   'instrument_T' : 299.15, ## lab temperature [K]
				   'rho_w'		  : 997.,	## density of water [kg/m3]
				   'sig_wa'		  : .072}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)
		
		## read file
		ccn, dma, cpc = self.ccn_raw(), self.dma_raw(), self.cpc_raw()
		
		## set data
		SS, ccnConc = ccn[' Current SS'].copy(), ccn[' CCN Number Conc'].copy()
		cpcConc	    = cpc['Concentration'].copy()
		diam		= dma['Diameter'].copy()

		## activation
		# ccnConc, cpcConc = rawCcnConc.copy(), rawCpcConc.copy()
		# wrongDt = ccnConc>cpcConc
		# ccnConc[wrongDt] = cpcConc[wrongDt]
		# actRat = ccnConc/cpcConc
		
		## loop parameter
		stTimeArr = pd.date_range(self.start,self.final,freq='30T',closed='left')
		
		## store parameter
		calSS = pd.DataFrame()
		actTable, diaTable = pd.DataFrame(), pd.DataFrame()
		actDia, SSList = [], []

		## for calibration, SS change every 30 mins, and diameter change every 2 mins,
		## however, CCNc is unstable on first 4 mins, then DMA is unstable 30 seconds
		## after the start and 30 seconds before the end
		## | ---- || -- || -- |.....| -- |	for ch : | - || ---- |
		## < uns  >< ch >< ch >.....< ch >	  	   	 <uns><  st  >
		## <           30 mins           > 	   	   	 <   2 mins  >
		## stable(st), unstable(uns), diameter change(ch)
		## 30min = 1800s, 4min 30s = 270s, 4min 90s = 330s, 2min = 120s
		## (30-4)/2 = 13 different diameter
		## class the data with dma data
		for stTime in stTimeArr:
			## parameter
			_SS = SS[stTime+dtmdt(seconds=270)]
			diaList = []
			actList = []
		
			for changSt in pd.date_range(stTime+dtmdt(seconds=270),stTime+dtmdt(seconds=1799),freq='2T'):
				## mean dianmeter
				diaList.append(diam[changSt:changSt+dtmdt(seconds=90)].mean())

				## mean activation and mean number conc
				meanCcnConc = ccnConc[changSt:changSt+dtmdt(seconds=90)].mean()
				meanCpcConc = cpcConc[changSt:changSt+dtmdt(seconds=90)].mean()
				actList.append(meanCcnConc/meanCpcConc*100.)
		
			## put in DataFrame
			actTable[_SS], diaTable[_SS] = actList, diaList
		
			## calculate activate diameter
			func = interp1d(diaList,actList,kind='cubic')
			interX = n.linspace(diaList[0],diaList[10],2000)
			actDia.append(interX[(func(interX)-50.).__abs__().argmin()])

			## data append
			SSList.append(_SS)

		## calculate SS
		coe_A = 4.*default['sig_wa']/(461.*default['instrument_T']*default['rho_w'])*1e9 ## diameter [nm]
		calSS = n.exp((4.*coe_A**3./(27.*default['kappa']*n.array(actDia)**3.))**.5)*100.-100.

		## data saving
		mdfy_data_calib = {'calib_SS'	: pd.Series(calSS.tolist()).sort_values(),
						   'instr_SS'   : pd.Series(SSList).sort_values(),
						   'act_dia'    : pd.DataFrame(actDia).set_index(pd.Index(SSList)).T,
						   'activation' : actTable,
						   'diameter'   : diaTable,
						   'raw'  		: pd.DataFrame({'ccn' : ccnConc, 'cpc' : cpcConc, 'dma' : diam})}
						   # 'ccn_cpc_coe' : pd.Series(actTable[actTable>80.].mean().mean()/100.)}

		## save output data, use 'w' mode to overwrite
		if outDt:
			_file = pth(self.path_output,f"calibration_{self.start.strftime('%Y%m%d')}.hdf5")
			print(f"\n{dtm.now().strftime('%m/%d %X')} : Saving file -- {_file}")
			with pd.HDFStore(_file,mode='w') as f:
				## save data
				for namDt, data in mdfy_data_calib.items():
					data.to_hdf(f,namDt)

		return mdfy_data_calib

	def mdfy_data_mesr(self,cpc_cor_slope=1.,ccn_cor_slope=1.,smps_data=True,cpc_data=True,
					   kappa_data=True,smpsOther_data=False,calib_SS_date=None,outDt=False):
		print('\n'+'-'*50)
		print(f"{dtm.now().strftime('%m/%d %X')} : Modify measurement data")
		print('-'*50)
		## smps data
		if smps_data:
			print(f"\n{dtm.now().strftime('%m/%d %X')} : Processing SMPS data")
			smps = self.smps_raw()
			smpsData = smps[smps.keys()[6:113]] ## data (DataFrame)
			smpsBin = n.array(smps.keys()[6:113]).astype(float)

		else: smpsData = None

		## other smps data
		if smpsOther_data:
			print(f"\nProcessing SMPS data -- {dtm.now().strftime('%m/%d %X')}")
			smps = self.smps_raw()
			smpsData = smps[smps.keys()[2:109]]
			smpsBin = n.array(smps.keys()[2:109]).astype(float)

		elif not smps_data: smpsData = None

		## cpc data
		if cpc_data:
			print(f"\n{dtm.now().strftime('%m/%d %X')} : Processing CPC data")
			cpc = self.cpc_raw()
			rawCpcData = cpc['Concentration'].copy()
			cpcData = rawCpcData/cpc_cor_slope ## data (Series)
					   
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
			print(f"\n{dtm.now().strftime('%m/%d %X')} : Processing kappa data")
			## get ccn, smps, cpc data and take necessary information
			## ccn
			ccn = self.ccn_raw()
			rawCcnData = ccn[' CCN Number Conc'].copy()
			ccnConc = rawCcnData/ccn_cor_slope
			ccnSS	= ccn[' Current SS'].copy()

			## get smps data, if smps_data = False, get raw data
			## first smps bin data could not correspond to activate ratio, neglect it
			## resample as 10 min
			try:
				smpsBinDt = n.array(smpsData[1:].asfreq('10T'))
			except:
				if smpsOther_data:
					smps = self.smpsOthers_raw()
					smpsBin   = n.array(smps.keys()[2:109]).astype(float)
					smpsBinDt = n.array(smps[smps.keys()[2:109]][1:].asfreq('10T'))
				else:
					smps = self.smps_raw()
					smpsBin = n.array(smps.keys()[6:113]).astype(float)
					smpsBinDt = n.array(smps[smps.keys()[6:113]][1:].asfreq('10T'))
			time = smps.index[1::2]
	
			## get cpc data, if cpc_data = False, get raw data
			try:
				cpcConc = cpcData.copy()
			except:
				cpc = self.cpc_raw()
				rawCpcData = cpc['Concentration'].copy()
				cpcConc = rawCpcData/cpc_cor_slope ## data (Series)

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
			# wrongDt = ccnConc>cpcConc
			# ccnConc[wrongDt] = cpcConc[wrongDt]

			## last ccn and cpc data could not use(closed='left'), and label should be same as smps data([:-1])
			## resample to 10 minute
			meanCcn = ccnConc.resample('5T',label='right').mean()[:-1].asfreq('10T')
			meanCpc = cpcConc.resample('5T',label='right').mean()[:-1].asfreq('10T')
			actRat  = meanCcn/meanCpc
			ccnSS  = ccnSS.resample('5T',label='right').mean()[:-1].round(1).asfreq('10T')

			if calib_SS_date is not None:
				print(f"\n	{dtm.now().strftime('%m/%d %X')} : Use calibration SS")
				calSSfile = pth('calibration',f'{calib_SS_date}',f'calibration_{calib_SS_date}.hdf5')
				calibSS = pd.read_hdf(calSSfile,'calib_SS').values
				instrSS = pd.read_hdf(calSSfile,'instr_SS').values.round(1)

				for _instrSS, _calibSS in zip(instrSS,calibSS):
					ccnSS.replace(_instrSS,_calibSS,inplace=True)

			## calculate activative diameter
			## accumulate bins data and normalize to 1, then compare with activative ratio, get the miniumum index
			## normalize bins data with miniumum index divide into activative ratio then multiply 
			## original smps bin diameter with miniumum index to get real activative diameter
			print(f"\n	{dtm.now().strftime('%m/%d %X')} : Calculating activate diameter")
			from scipy.interpolate import interp1d
			def actDaFunc(_bin_dt,_f_act):
				## nan test
				if (n.isnan(_bin_dt).all()|n.isnan(_f_act)): return n.nan
				_accu_bin = n.add.accumulate(_bin_dt[::-1])[::-1]/_bin_dt.sum()
			
				## inter function
				func = interp1d(smpsBin,_accu_bin,kind='cubic')
				interX = n.logspace(n.log10(smpsBin[0]),n.log10(smpsBin[-1]),5000)[5:-5]
			
				return interX[(func(interX)-_f_act).__abs__().argmin()]

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
			print(f"\n	{dtm.now().strftime('%m/%d %X')} : Calculating kappa")
			kappa = 4.*coeA**3./(27.*actDaDt**3.)/(n.log(ccnSS/100.+1.)**2)
			
			## primitive
			kappa[limit] = n.array([ _k for _k in map(kappaApprox,ccnSS[limit],actDaDt[limit]) ])
			kappaData = pd.DataFrame({'kappa' : kappa, 'mesrSS' : ccnSS}).set_index(time)

		else: kappaData = None

		mdfy_data_mesr = {'smps'  : smpsData,
						  'cpc'   : cpcData, 
						  'kappa' : kappaData,
						  'raw'   : pd.DataFrame({'cpc' : rawCpcData, 'ccn' : rawCcnData, 'act' : actRat*100.}),
						  'SMPS_bins' : pd.Series(smpsBin)}

		## save output data, use 'w' mode to overwrite
		if (outDt&smps_data&cpc_data&kappa_data):
			_file = pth(self.path_output,f"output_{self.start.strftime('%Y%m%d')}.hdf5")
			print(f"\nSaving file : {_file} -- {dtm.now().strftime('%m/%d %X')}")
			with pd.HDFStore(_file,'w') as f:
				## save data
				for namDt, data in mdfy_data_mesr.items():
					data.to_hdf(f,namDt)
				if calib_SS_date: f['Calibration_SS_Date'] = pd.Series(calib_SS_date) ## save calibration date
				else: f['Calibration_SS_Date'] = pd.Series('Without calibration SS')

		return mdfy_data_mesr

# calibration output
class calibration:
	def __init__(self,calib_date,data=None,**kwarg):
		print('\n'+'='*50)
		print(f"Plot calibration data")
		## set calculating parameter
		default = {'fig_path' : pth('./'),
				   'path_input_data' : pth('calibration',f'{calib_date}',f'calibration_{calib_date}.hdf5')}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## detect dir exists or not
		if not exists(default['fig_path']): mkdir(default['fig_path'])

		## get data
		if data is None:
			if exists(default['path_input_data']):
				with pd.HDFStore(default['path_input_data'],'r') as f:
					print(f"\n{dtm.now().strftime('%m/%d %X')} : Loading file -- {default['path_input_data']}")
					self.calibSS = f['calib_SS']
					self.instrSS = f['instr_SS']
					self.actDia  = f['act_dia']
					self.diam	 = f['diameter']
					self.activ	 = f['activation']
					self.raw	 = f['raw']

			else: raise OSError(f"File '{default['path_input_data']}' does not exist !!!")
		else:
			self.calibSS = data['calib_SS']
			self.instrSS = data['instr_SS']
			self.actDia  = data['act_dia']
			self.diam	 = data['diameter']
			self.activ	 = data['activation']
			self.raw	 = data['raw']

		## set class parameter
		self.fs = 13.
		self.calib_date = calib_date
		self.figPath = default['fig_path']
		print(f" calibration date : {calib_date}")
		print('='*50)
		print(f"{dtm.now().strftime('%m/%d %X')}")

	## plot S curve to find out criticle diameter
	def plot_Scurve(self,**kwarg):
		print(f"\n{dtm.now().strftime('%m/%d %X')} : Plotting S curve")
		## set plotting parameter
		default = {'color'	  	: '#b8dbff',
				   'line_color' : '#ff5f60',
				   'mfc'	  	: '#7fdfff',
				   'mec'	  	: '#297db7',
				   'fig_row'	: 2,
				   'fig_col'	: 3,
				   'fig_name' 	: r'calib_Scurve.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## parameter
		fs = self.fs
		size = self.instrSS.size

		## plot
		fig, axes = subplots(default['fig_row'],default['fig_col'],figsize=(8,6),dpi=200,sharex=True,sharey=True)

		for ax, _SS in zip(axes.flatten(),self.instrSS):
			## variables
			critiDia = self.actDia[_SS].loc[0]

			## plot
			ax.plot(self.diam[_SS],self.activ[_SS],c=default['color'],marker='o',mfc=default['mfc'],mec=default['mec'])
			ax.plot(critiDia,50.,c=default['line_color'],ms=5.,marker='o')

			ax.tick_params(which='major',direction='in',length=7,labelsize=fs-2.5)
			[ ax.spines[axis].set_visible(False) for axis in ['right','top'] ]

			ax.set(xlim=(-30.,299.),ylim=(0.,110.))
			ax.set_title(f'SS = {_SS:4.2f} %',fontsize=fs-2)
			ax.text(critiDia+15.,45.,f'$D_c$ = {critiDia:.2f} nm',fontsize=fs-2)

		fig.text(.38,.03,'Dry Particle Diameter (nm)',fontsize=fs)
		fig.text(.04,.4,'Activation (%)',rotation=90.,fontsize=fs)

		## for odd axes
		if ((size%2==1)&(size!=1)): ax[-1].remove()

		fig.suptitle(f"Activation of CCN (Date : {self.calib_date})",fontsize=fs+3,style='italic')
		fig.savefig(pth(self.figPath,default['fig_name']))
		close()

	## use kappa kohler eq. and modified data to calculate the SS, then plotting regression line
	## with CCNc's SS
	def plot_SStable(self,**kwarg):
		print(f"\n{dtm.now().strftime('%m/%d %X')} : Plotting SS calibration table")
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
		fs = self.fs
		height = .15
		size = self.instrSS.size
		xPosition = n.arange(size)

		## plot
		fig, ax = subplots(figsize=(8,6),dpi=200)
		box = ax.get_position()
		ax.set_position([box.x0,box.y0+.03,box.width,box.height*.95])

		## plot 6 calibration data, to align with table, x data use index data
		ax.plot(xPosition,self.calibSS,mec=default['mec'],mew=default['mew'],mfc=default['mfc'],
				ms=default['ms'],marker='o',ls='')
		table = ax.table(cellText=[self.calibSS.round(3).values],rowLabels=['Calib SS'],colLabels=self.instrSS.values,
						 cellLoc='center',rowLoc='center',colLoc='center',
						 bbox=[0.,-height,1.,height])
		table.auto_set_font_size(False)
		table.set_fontsize(fs)

		ax.tick_params(which='both',direction='in',labelsize=fs,right=True,top=True,labelbottom=False)
		ax.set(xlim=(xPosition[0]-.5,xPosition[-1]+.5),ylim=(.05,.95))

		ax.set_ylabel('Calibration SS (%)',fontsize=fs)
		ax.set_title('CCNc SS (%)',fontsize=fs)

		fig.suptitle(f'Calibration Table of Supersaturation (Date : {self.calib_date:})',
					 fontsize=fs+3,style='italic')
		fig.savefig(pth(self.figPath,default['fig_name']))
		close()

# measurement output
class measurement:
	def __init__(self,start,final,data=None,**kwarg):
		print('\n'+'='*50)
		print(f"Plot measurement data")
		## set calculating parameter
		default = {'fig_path' : pth('./'),
				   'splt_hr'  : 12,
				   'path_input_data' : pth('measurement',f"{start.strftime('%Y%m%d')}",
										   f"output_{start.strftime('%Y%m%d')}.hdf5")}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## func
		index = lambda _freq: pd.date_range(start,final,freq=_freq)
		def ticks(split_hr):
			_tick = index(f'{split_hr}h')
			return _tick, _tick.strftime('%m-%d%n%X')

		## detect dir exists or not
		if not exists(default['fig_path']): mkdir(default['fig_path'])

		## get data
		if data is None:
			if exists(default['path_input_data']):
				with pd.HDFStore(default['path_input_data'],'r') as f:
					print(f"\n{dtm.now().strftime('%m/%d %X')} : Loading file -- {default['path_input_data']}")
					print(f"	Calibration SS date : {f['Calibration_SS_Date'][0]}")
					self.smpsBins  = f['SMPS_bins'].values
					self.smpsData  = f['smps']
					self.cpcData   = f['cpc']
					self.kappaData = f['kappa']
					self.rawData   = f['raw']
			else: raise OSError(f"File '{default['path_input_data']}' does not exist !!!")
		else:
			self.smpsBins  = data['SMPS_bins'].values
			self.smpsData  = data['smps']
			self.cpcData   = data['cpc']
			self.kappaData = data['kappa']
			self.rawData   = data['raw']

		## set class parameter
		self.fs = 13.
		self.ticks = ticks
		self.index = index
		self.start = start
		self.final = final
		self.figPath = default['fig_path']
		self.splt_hr = default['splt_hr']
		print(f" from {start.strftime('%Y-%m-%d %X')} to {final.strftime('%Y-%m-%d %X')}")
		print('='*50)
		print(f"{dtm.now().strftime('%m/%d %X')}")

	## plot smps with date and set log scale data
	def plot_smps2date(self,plotTogether=None,**kwarg):
		print(f"\n{dtm.now().strftime('%m/%d %X')} : Plotting SMPS data")
		from matplotlib.colors import LogNorm
		## set plot parameter
		default = {'splt_hr'  : self.splt_hr,
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
		data = smps.reindex(time).T
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
		pm = ax.pcolormesh(time,self.smpsBins,data,cmap=default['cmap'],norm=LogNorm(vmin=10**.5,vmax=10**5.3))

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
		print(f"\n{dtm.now().strftime('%m/%d %X')} : Plotting CPC data")
		## set plot parameter
		default = {'splt_hr'  : self.splt_hr,
				   'fig_name' : r'mesr_cpc2date.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## set plot variable
		if self.cpcData is not None: cpc = self.cpcData 
		else: raise ValueError('CPC data is None !!!')

		## data set
		time = self.index('s')
		data = cpc.reindex(time)/1e3

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
		print(f"\n{dtm.now().strftime('%m/%d %X')} : Plotting Kappa data")
		## set plot parameter
		default = {'splt_hr'  : self.splt_hr,
				   'fig_name' : r'mesr_kappa2date.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## set plot variable
		if self.kappaData is not None: kappa = self.kappaData 
		else: raise ValueError('Kappa data is None !!!')

		## data set
		time  = self.index('5T')
		data = kappa['kappa'].reindex(time).groupby(kappa['mesrSS'])

		xTick, xTickLab = self.ticks(default['splt_hr'])
		fs = self.fs

		color = ['#7617FF','#FF4924','#0A82FF','#efe934','#17FFA4','#FF9C2F','#e699ff','#000000']
		
		## plot
		if plotTogether is None:
			fig, ax = subplots(figsize=(10.,6.),dpi=150.)
		else:
			fig, ax = plotTogether
		collectArt = []
		## plot data and necessary setting
		for (_ss, _kappa), _colr in zip(data,color):
			artist, = ax.plot(_kappa,mec=_colr,ls='',ms=5.,mfc='#ffffff',mew=2.,marker='o',label=round(_ss,3))
			collectArt.append(artist)

		ax.tick_params(which='major',direction='in',length=7,right=True,labelsize=fs-2.5)
		ax.set(xticks=xTick,xlim=(xTick[0],xTick[-1]))
		ax.set_ylim(bottom=-.02)

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
		print(f"\n{dtm.now().strftime('%m/%d %X')} : Plotting CPC, SMPS, Kappa data")
		## set plot parameter
		default = {'splt_hr'  : self.splt_hr,
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

class raw_data:
	def __init__(self,start,final,data=None,**kwarg):
		print('\n'+'='*50)
		print(f"Plot raw data")
		## set calculating parameter
		default = {'fig_path' : pth('./'),
				   'splt_hr'  : 12,
				   'path_input_data' : pth('measurement',f"{start.strftime('%Y%m%d')}",
										   f"output_{start.strftime('%Y%m%d')}.hdf5")}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## func
		index = lambda _freq: pd.date_range(start,final,freq=_freq)
		def ticks(split_hr):
			_tick = index(f'{split_hr}h')
			return _tick, _tick.strftime('%m-%d%n%X')

		## detect dir exists or not
		if not exists(default['fig_path']): mkdir(default['fig_path'])

		## get data
		if data is None:
			if exists(default['path_input_data']):
				with pd.HDFStore(default['path_input_data'],'r') as f:
					print(f"\n{dtm.now().strftime('%m/%d %X')} : Loading file -- {default['path_input_data']}")
					self.rawData = f['raw']
			else: raise OSError(f"File '{default['path_input_data']}' does not exist !!!")
		else:
			self.rawData = data['raw']

		## set class parameter
		self.fs = 13.
		self.ticks = ticks
		self.index = index
		self.start = start
		self.final = final
		self.figPath = default['fig_path']
		self.splt_hr = default['splt_hr']
		print(f" from {start.strftime('%Y-%m-%d %X')} to {final.strftime('%Y-%m-%d %X')}")
		print('='*50)
		print(f"{dtm.now().strftime('%m/%d %X')}")

	def plot_raw(self,plot_act=True,**kwarg):
		print(f"\n{dtm.now().strftime('%m/%d %X')} : Plotting CPC, CCN, time series and satter")
		default = {'splt_hr'  : self.splt_hr,
				   'fig_name_scatter'  : 'raw_scatter.png',
				   'fig_name_tmSeries' : 'raw_tmSeries.png'}

		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## data set
		time = self.index('s')
		rawCcn = self.rawData['ccn'].reindex(time)
		rawCpc = self.rawData['cpc'].reindex(time)

		xTick, xTickLab = self.ticks(default['splt_hr'])
		fs = self.fs
		
		## plot
		## time series
		print(f"\n	{dtm.now().strftime('%m/%d %X')} : plot time series")
		fig, ax = subplots(figsize=(10,6),dpi=150.)

		l1, = ax.plot(rawCcn,c='#ffb973',label='CCN',lw=1.2)
		l2, = ax.plot(rawCpc,c='#5982ff',label='CPC',lw=1.2)
		artist = [l1,l2]

		ax.tick_params(which='major',direction='in',length=7,labelsize=fs-2.5)
		[ ax.spines[axis].set_visible(False) for axis in ['right','top'] ]

		ax.set_xlabel(f"Time({self.start.strftime('%Y')})",fontsize=fs)
		ax.set_ylabel(r'Number Concentration (#$/cm^3$)',fontsize=fs)

		if plot_act:
			try:
				act = self.rawData['act'].dropna()

				axt = ax.twinx()
				axt.set(xlim=ax.get_xlim(),ylim=(-50.,105.))

				l3, = axt.plot(act.index,act,c='#686859',label='activation',lw=1.5,ls='--',marker=',')
				artist.append(l3)

				axt.set(xlim=ax.get_xlim(),ylim=(-50.,105.))

				yTick = axt.get_yticks()
				axt.set_yticks(yTick[(yTick>=0.)&(yTick<=110.)])
				axt.set_xticks(xTick)

				axt.set_ylabel(r'Activation (%)',fontsize=fs)
			except:
				print('	no activation data in raw data')

		ax.set_xticks(xTick)
		ax.set_xticklabels(xTickLab)
		ax.legend(handles=artist,framealpha=0,fontsize=fs-2.)

		fig.suptitle(f"Time series of CPC and CCN from ({self.start.strftime('%Y/%m/%d %X')}) to ({self.final.strftime('%Y/%m/%d %X')})",
					 fontsize=fs+2.,style='italic')
		fig.savefig(pth(self.figPath,default['fig_name_tmSeries']))
		close()

		## scatter
		print(f"\n	{dtm.now().strftime('%m/%d %X')} : plot scatter")
		fig, ax = subplots(figsize=(8,6),dpi=150.)
		
		ax.scatter(rawCcn,rawCpc,ec='#26c9ff',s=35,facecolor='#ffffff')
		
		ax.tick_params(which='major',direction='in',length=7,labelsize=fs-2.5)
		[ ax.spines[axis].set_visible(False) for axis in ['right','top'] ]
		
		
		ax.set_xlabel('CCN (#/$cm^3$)',fontsize=fs)
		ax.set_ylabel('CPC(ver. 2.1.2) (#/$cm^3$) ',fontsize=fs)
		
		# ax.legend(framealpha=0,fontsize=fs-2.)
		#ax.legend(handles=[],framealpha=0,fontsize=fs-2.)

		fig.suptitle(f"Scatter of CCN and CPC \nfrom {self.start.strftime('%Y/%m/%d %X')}) to ({self.final.strftime('%Y/%m/%d %X')})",fontsize=fs+2.,style='italic')
		fig.savefig(pth(self.figPath,default['fig_name_scatter']))
		# show()
		close()












































