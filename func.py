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
import warnings
warnings.simplefilter('error')
warnings.simplefilter('ignore', category=RuntimeWarning)

# file reader
class reader:
	def __init__(self,start,final,td1DtPrces=False,**kwarg):
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
		self.td1DtPrces = td1DtPrces
		self.path_ccn   = default['path_ccn']
		self.path_dma   = default['path_dma']
		self.path_cpc   = default['path_cpc']
		self.path_smps  = default['path_smps']
		# cc=codecs.open('../CCN/CCN 100 data 200415132942.csv', 'r', encoding='utf-8',errors='ignore')

	## data of corrected time
	## start test
	## wrong current time test
	## discountinue data test
	## output data
	def correct_time_data(self,_begin_,_fout_,_line,tm_start,td1_start=False,**_par):
		## raw data
		nan_dt = _par['nanDt']
		_data_ = _par['dtSplt'](_line)

		try:
			## deal with blank line data and unfortunate header
			## skip it
			cur_tm = _data_[_par['tmIndx']]
			chkDtmOb = dtm.strptime(cur_tm,'%X') ## if cur_tm is not datetime string format, it will error
		except:
			return tm_start, _begin_
		
		err_tm, del_tm = _par['errTm'], _par['delTm']
		cort_tm_now = lambda time: dtm.strftime(time,'%X')
		cort_tm_low = lambda time, err: dtm.strftime(time-err,'%X')
		cort_tm_upr = lambda time, err: dtm.strftime(time+err,'%X')
		tm_logi = lambda _tm, err: ((cur_tm != cort_tm_low(_tm,err))&(cur_tm != cort_tm_now(_tm)))

		## (start time different = 1) test
		if (td1_start&_begin_):
			td1_data = cur_tm == cort_tm_now(tm_start+dtmdt(seconds=1.))
			if td1_data:
				_begin_ = False					## skip the start test
				_fout_.append(nan_dt(tm_start)) ## add nan to current time data
				tm_start += del_tm				## change to the time which has first data
		
		## start test
		if _begin_:
			_begin_ = True if tm_logi(tm_start,_par['stErrTm']) else False
			if _begin_: return tm_start, _begin_

		## wrong current time test
		## check out data time series by current time +/- error time
		if ((cur_tm == cort_tm_low(tm_start,err_tm))|(cur_tm == cort_tm_upr(tm_start,err_tm))):
			cur_tm = dtm.strftime(tm_start,'%X')

		## discountinue data test
		while tm_logi(tm_start,err_tm):
			_fout_.append(nan_dt(tm_start)) ## add nan data to discontinue time
			tm_start += del_tm
			if tm_start>self.final: return None, _begin_

		## data colllect
		_data_[_par['tmIndx']] = tm_start

		_fout_.append(_data_)
		tm_start += del_tm
		if tm_start<=self.final: return tm_start, _begin_
		else: return None, _begin_

	## SMPS
	## change time every 5 min
	## information : 17 lines
	## keys index : 'Date' : 1 (%D)
	## 				'Start Time' : 2 (%X)
	## 				'diameter' : 4~110
	def smps_raw(self):
		path, fout, switch, begin = self.path_smps, [], True, True
		tmStart = self.start

		for file in listdir(path):
			if '.txt' not in file: continue
			with open(pth(path,file),errors='replace') as f:
				## get header and skip instrument information
				if switch==True:
					[ f.readline() for i in range(17) ]
					header = f.readline()[:-1].split('\t')
					switch = False
				else: [ f.readline() for i in range(18) ]

				## collect data from start time, and make nan array for discontinue data
				## [2] is current time index
				par = { 'nanDt'   : lambda _tmStart: n.array([n.nan]*2+[_tmStart]+[n.nan]*len(header[3::])),
						'dtSplt'  : lambda _li: n.array(_li[:-1].split('\t'),dtype=object),
						'delTm'   : dtmdt(minutes=5.),
						'errTm'   : dtmdt(seconds=1.),
						'stErrTm' : dtmdt(seconds=1.),
						'hdrLn'	  : len(header),
						'tmIndx'  : 2 }

				for line in f:
					## time check out and collect data
					if tmStart: tmStart, begin = self.correct_time_data(begin,fout,line,tmStart,**par)
					else: break
			if not tmStart: break
		return dict(zip(header,n.array(fout).T))

	## SMPS_others
	## change time every 5 min
	## information : 15 lines
	## keys index : 'Date' : 1 (%D)
	## 				'Start Time' : 2 (%X)
	## 				'diameter' : 4~110
	def smpsOthers_raw(self):
		path, fout, switch, begin = self.path_smps, [], True, True
		tmStart = self.start

		for file in listdir(path):
			if '.txt' not in file: continue
			with open(pth(path,file),errors='replace') as f:
				## get header and skip instrument information
				if switch==True:
					[ f.readline() for i in range(15) ]
					header = f.readline()[:-1].split('\t')
					switch = False
				else: [ f.readline() for i in range(16) ]

				## collect data from start time, and make nan array for discontinue data
				## [2] is current time index
				par = { 'nanDt'   : lambda _tmStart: n.array([n.nan]*2+[_tmStart]+[n.nan]*len(header[3::])),
						'dtSplt'  : lambda _li: n.array(_li[:-1].split('\t'),dtype=object),
						'delTm'   : dtmdt(minutes=5.),
						'errTm'   : dtmdt(seconds=1.),
						'stErrTm' : dtmdt(seconds=1.),
						'hdrLn'	  : len(header),
						'tmIndx'  : 2 }

				for line in f:
					## time check out and collect data
					if tmStart: tmStart, begin = self.correct_time_data(begin,fout,line,tmStart,**par)
					else: break
			if not tmStart: break
		return dict(zip(header,n.array(fout).T))

	## CCNc
	## name : CCN 100 data %y%m%d%M0000.csv
	## change time every second
	## keys index : '    Time' : 0 (%X)
	## 				' CCN Number Conc' : 45
	def ccn_raw(self):
		path, fout, switch, begin = self.path_ccn, [], True, True
		tmStart = self.start

		for file in listdir(path):
			if '.csv' not in file: continue
			with open(pth(path,file),errors='replace') as f:
				## get header and skip instrument information
				if switch:
					[ f.readline() for i in range(4) ]
					header = f.readline()[:-1].split(',')[:-1]
					f.readline()
					switch = False
				else: [ f.readline() for i in range(6) ]

				## collect data from start time, and make nan array for discontinue data
				## [0] is current time

				par = { 'nanDt'   : lambda _tmStart: n.array([_tmStart]+['nan']*len(header[1::])),
						'dtSplt'  : lambda _li: n.array(_li[:-1].split(','),dtype=object),
						'delTm'   : dtmdt(seconds=1.),
						'errTm'   : dtmdt(seconds=1.),
						'stErrTm' : dtmdt(seconds=0.),
						'hdrLn'	  : len(header),
						'tmIndx'  : 0 }

				for line in f:
					## time check out and collect data
					if tmStart: tmStart, begin = self.correct_time_data(begin,fout,line,tmStart,**par)
					else: break
			if not tmStart: break

		## alarm code
		raw_dt = dict(zip(header,n.array(fout).T))
		alarm_dt = raw_dt[' Alarm Code']!='    0.00'
		raw_dt[' CCN Number Conc'][alarm_dt] = n.nan
		return raw_dt

	## DMA
	## name : %Y%m%d.txt
	## change time every second
	## keys index : 'Date' : 0 (%Y/%m/%d)
	## 				'Time' : 1 (%X)
	## 				'Diameter' : 3
	def dma_raw(self):
		path, fout = self.path_dma, []
		header = ['Date','Time','Diameter','SPD']
		for file in listdir(path):
			if '.txt' not in file: continue
			with open(pth(path,file),errors='replace') as f:
				[ fout.append(n.array(line[:-1].split('\t'))) for line in f ]

		return dict(zip(header,n.array(fout).T))

	## CPC
	## name : %y%m%d_cpc.csv
	## change time every second
	## keys index : 'Time' : 0 (%X)
	## 				'Concentration (#/cm?)' : 1
	def cpc_raw(self):
		path, fout, switch, begin = self.path_cpc, [], True, True
		tmStart = self.start

		for file in listdir(path):
			if '.csv' not in file: continue
			with open(pth(path,file),errors='replace') as f:
				if switch:
					[ f.readline() for i in range(17) ]
					header = f.readline()[:-1].split(',')[:-1]
					switch = False
				else: [ f.readline() for i in range(18) ]

				## collect data from start time, and make nan array for discontinue data
				## [0] is current time
				par = { 'nanDt'   : lambda _tmStart: n.array([_tmStart]+['nan']*len(header[1::])),
						'dtSplt'  : lambda _li: n.array(_li[:-1].split(',')[:-1],dtype=object),
						'delTm'   : dtmdt(seconds=1.),
						'errTm'   : dtmdt(seconds=1.),
						'stErrTm' : dtmdt(seconds=0.),
						'hdrLn'	  : len(header),
						'tmIndx'  : 0 }

				for line in f:
					## check out time and collect data
					if tmStart: tmStart, begin = self.correct_time_data(begin,fout,line,tmStart,td1_start=self.td1DtPrces,**par)
					else: break
			if not tmStart: break

		return dict(zip(header,n.array(fout).T))

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

	def mdfy_data_mesr(self,smps_data=True,cpc_data=True,kappa_data=True,smpsOther_data=False,calib_SS_date=None,outDt=False):
		## smps data
		if smps_data:
			smps = self.smps_raw()
			smpsKey = list(smps.keys())
			smpsBin = n.array([ smps[key] for key in smpsKey[8:115] ],dtype=float) ## [ [ smps[Dp1] ], [ smps[Dp2] ], .... ]
			smpsTm  = smps[smpsKey[2]]
			smpsData = {'time' 	   : smpsTm,
						'bin_data' : smpsBin,
						'bins'	   : n.array(smpsKey[8:115],dtype=float)}
		else: smpsData = None

		## smps other data
		if smpsOther_data:
			smps = self.smps_raw()
			smpsKey = list(smps.keys())
			smpsBin = n.array([ smps[key] for key in smpsKey[4:111] ],dtype=float) ## [ [ smps[Dp1] ], [ smps[Dp2] ], .... ]
			smpsTm  = smps[smpsKey[2]]
			smpsData = {'time' 	   : smpsTm,
						'bin_data' : smpsBin,
						'bins'	   : n.array(smpsKey[4:111],dtype=float)}
		elif not smps_data: smpsData = None

		## cpc data
		if cpc_data:
			cpc = self.cpc_raw()
			cpcData = {'time' : cpc['Time'],
					   'conc' : cpc['Concentration (#/cm�)'].astype(float)}
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
			## import accumulate iterator
			# from itertools import accumulate as accu

			## get ccn, smps, cpc data and take necessary information
			## ccn
			ccn = self.ccn_raw()
			ccnConc = ccn[' CCN Number Conc'].astype(float)
			ccnSS	= ccn[' Current SS'].astype(float)

			## get smps data, if smps_data = False, get raw data
			## first smps bin data could not correspond to activate ratio, neglect it
			
			try:
				smpsBinDp = smpsData['bins']
				smpsBin_perDy = smpsBin.T[1::]
				smpsTm = smpsTm[2::2]
			except:
				if smpsOther_data:
					smps = self.smpsOthers_raw()
					smpsBinDp = n.array(list(smps.keys())[4:111])
				else:
					smps = self.smps_raw()
					smpsBinDp = n.array(list(smps.keys())[8:115])
				smpsBin_perDy = n.array([ smps[key] for key in smpsBinDp ],dtype=float).T[1::]
				smpsTm  = smps[list(smps.keys())[2]][2::2]
				smpsBinDp = smpsBinDp.astype(float)
			
			## get cpc data, if cpc_data = False, get raw data
			try:
				cpcConc = cpcData['conc']
			except:
				cpc = self.cpc_raw()
				cpcConc = cpc['Concentration (#/cm�)'].astype(float)

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
			cpcConc[ccnConc>cpcConc] = n.nan
			actRat 	   = n.nanmean(n.reshape(ccnConc[:-1]/cpcConc[:-1],(-1,300)),axis=1)
			ssPer10min = n.nanmean(ccnSS[:-1].reshape(-1,300),axis=1)[1::2]
			if calib_SS_date is not None:
				with hpFile('output.hdf5','r') as f:
					dset = f['test/Calibration']
					calSS = dset[calib_SS_date]
					for _indSS, _calibSS in zip(calSS['ssSet'],calSS['ssCalib']):
						ssPer10min[n.abs(ssPer10min-float(_indSS))<1e-3] = _calibSS

			## calculate activative data
			## accumulate bins data and normalize to 1, then compare with activative ratio, get the miniumum index
			## normalize bins data with miniumum index divide into activative ratio then multiply 
			## original smps bin diameter with miniumum index to get real activative diameter
			## (ccn/cpc) / (accumulate bin [min]) * (smps bin [min])
			def actDaFunc(_bin_dt,_f_act):
				## nan test
				if (n.isnan(_bin_dt).all()|n.isnan(_f_act)): return n.nan
			
				_accu_bin = n.add.accumulate(_bin_dt[::-1])[::-1]/_bin_dt.sum()
				_ind_da	  = n.abs(_accu_bin-_f_act).argmin()

				try:
					if _accu_bin[_ind_da]>_f_act:
						return n.poly1d(n.polyfit(_accu_bin[_ind_da:_ind_da+2],smpsBinDp[_ind_da:_ind_da+2],1))(_f_act)
					else:
						return n.poly1d(n.polyfit(_accu_bin[_ind_da-1:_ind_da+1],smpsBinDp[_ind_da-1:_ind_da+1],1))(_f_act)
				except:
					return _accu_bin[_ind_da]
			
			actDaDt = n.array([ actDaFunc(bin_dt,f_act) for bin_dt, f_act in 
								zip(smpsBin_perDy[1::2],actRat[1::2]) ])
			actDaDt[(actDaDt<smpsBinDp[0])|(actDaDt>smpsBinDp[-1])] = n.nan
			
			## calculate kappa
			## use simplify if SS < 0.3
			## parameter
			coeA  = 4.*.072/(461.*299.15*997.)*1e9
			limit = ssPer10min<.3 ## use primitive function

			## primitive function
			## calculate critical S
			dEq  = lambda _da : 10**n.arange(n.log10(_da)+.0001,2.64,.0001) ## nm
			criS = lambda _da, _kappa : n.max(((dEq(_da)**3.-_da**3)/(dEq(_da)**3.-_da**3+_kappa*_da**3))*n.exp(coeA/dEq(_da)))

			## kappa approximation
			def kappaApprox(_ss,_da):
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
			kappaPer10min = 4.*coeA**3./(27.*actDaDt**3.)/(n.log(ssPer10min/100.+1.)**2)
			
			## primitive
			kappaPer10min[limit] = n.array([ _k for _k in map(kappaApprox,ssPer10min[limit],actDaDt[limit]) ])
			
			kappaData = {'data_time'  : smpsTm,
						 'time' 	  : cpc['Time'],
						 'act_dia'    : actDaDt,
						 # 'simp_kappa' : simpKappaPer10min,
						 # 'kappa' : simpKappaPer10min,
						 'kappa'	  : kappaPer10min,
						 # 'pri_kappa'  : primKappaPer10min,
						 'SS'		  : ssPer10min}

		else: kappaData = None

		mdfy_data_mesr = {'smps'  : smpsData, 
						  'cpc'   : cpcData, 
						  'kappa' : kappaData}

		## save output data
		if (outDt&smps_data&cpc_data&kappa_data):
			with hpFile('output.hdf5','r+') as f:
				print('it is test program, plz change to other key in real')
				mesr = f['test/Measurement']

				try:
					dset = mesr.require_group(self.start.strftime('%Y%m%d'))
					dset.clear()
				except:
					input('test')
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
	def __init__(self,start,final,data=None,dtDate=None,useOutput=True,**kwarg):
		## set calculating parameter
		default = {'fig_Path' : './'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## func
		def timeIndex(_time):
			return n.where(_time==start)[0][0], n.where(_time==final)[0][0]

		def onHrTicks(_time,split_dt):
			replaceTm = {'day'    : start.day+1, 
						 'hour'   : 0,
						 'minute' : 0,
						 'second' : 0}
			
			_tick1st_indx = n.where(_time==start.replace(**replaceTm))[0][0]
			_tick = n.arange(_tick1st_indx%split_dt,len(_time),int(split_dt))
			_tick_lab = [ _time[indx].strftime('%m/%d%n%X') for indx in _tick ]
			return _tick, _tick_lab
		
		## get data
		if (useOutput&(dtDate is not None)):
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
		self.onHrTicks = onHrTicks
		self.timeIndex = timeIndex
		self.fs = 13.
		self.start = start
		self.final = final
		self.figPath = default['fig_Path']

	## plot smps with date and set log scale data
	def plot_smps2date(self,plotTogether=None,**kwarg):
		from matplotlib.colors import LogNorm
		## set plot parameter
		default = {'splt_hr'  : 6,
				   'cmap'	  : 'jet',
				   'fig_name' : r'mesr_smps2date.png'}
		for key in kwarg:
			if key not in default.keys(): raise TypeError("got an unexpected keyword argument '"+key+"'")
			default.update(kwarg)

		## set plot variable
		if self.smpsData is not None: smps = self.smpsData 
		else: raise ValueError('SMPS data is None !!!')

		start_indx, final_indx = self.timeIndex(smps['time'])
		xTick, xTickLab = self.onHrTicks(smps['time'],default['splt_hr']*12)

		time = smps['time'][start_indx:final_indx+1]
		data = smps['bin_data'][:,start_indx:final_indx+1]
		data[data==0] = 1e-5
		data[n.isnan(data)] = 0.
		fs = self.fs

		## plot
		if plotTogether is None:
			fig, ax = subplots(figsize=(10.,6.),dpi=150.)
		else:
			fig, ax = plotTogether

		## plot data and necessary setting
		pm = ax.pcolormesh(n.arange(len(time)),smps['bins'],data,cmap=default['cmap'],norm=LogNorm(vmin=10**.5,vmax=10**5.3))

		box = ax.get_position()
		ax.set_position([box.x0,box.y0+0.02,box.width,box.height])
		cax = fig.add_axes([.92,box.y0+0.02,.015,box.height])

		cb = fig.colorbar(pm,cax=cax)

		ax.tick_params(which='major',length=6.,labelsize=fs-2.)
		ax.tick_params(which='minor',length=3.5)
		cb.ax.tick_params(which='major',length=5.,labelsize=fs-2.)
		cb.ax.tick_params(which='minor',length=2.5)
		ax.set(xlim=(0.,len(time)),yscale='log',xticks=xTick)

		## single plot
		if plotTogether is None:
			cb.ax.set_title('# conc.\n(#/$cm^3$/$\Delta log D_p$)',fontsize=fs-2.)
			ax.set_ylabel('Electric modify diameter (nm)',fontsize=fs)
			ax.set_xlabel(f"Time({self.start.strftime('%Y')})",fontsize=fs)
			ax.set_xticklabels(xTickLab)

			fig.suptitle(f"SMPS data from ({self.start.strftime('%Y/%m/%d %X')}) to ({self.final.strftime('%Y/%m/%d %X')})",fontsize=fs+2.,style='italic')
			fig.savefig(pth(self.figPath,default['fig_name']))
			close()

		## together plot
		else:
			ax.set_ylabel('diameter (nm)',fontsize=fs-1.)
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

		start_indx, final_indx = self.timeIndex(cpc['time'])
		xTick, xTickLab = self.onHrTicks(cpc['time'],default['splt_hr']*3600)

		time = cpc['time'][start_indx:final_indx+1]
		data = cpc['conc'][start_indx:final_indx+1]/1e3
		fs = self.fs
		ylim = (0.,max(data)+10.) if max(data)<30. else (0.,30.)
		tm_num = len(time)

		## plot
		if plotTogether is None:
			fig, ax = subplots(figsize=(10.,6.),dpi=150.)
		else:
			fig, ax = plotTogether
		box = ax.get_position()

		## plot data and necessary setting
		ax.plot(n.arange(tm_num),data,c='#7396ff',lw=1.5)

		ax.tick_params(which='major',direction='in',length=7,labelright=True,right=True,labelsize=fs-2.5)
		# ax.tick_params(which='minor',direction='in',length=4.5)
		ax.set(ylim=ylim,xlim=(0.,tm_num),xticks=xTick)

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

		start_indx, final_indx = self.timeIndex(kappa['data_time'])
		xTick, xTickLab = self.onHrTicks(kappa['data_time'],default['splt_hr']*6)

		time = kappa['data_time'][start_indx:final_indx+1]
		data = kappa['kappa']
		ss = kappa['SS']
		fs = self.fs
		tm_num = len(time)
		color = ['#7617FF','#FF4924','#0A82FF','#efe934','#17FFA4','#FF9C2F']
		
		## plot
		if plotTogether is None:
			fig, ax = subplots(figsize=(10.,6.),dpi=150.)
		else:
			fig, ax = plotTogether
		collectArt = []
		## plot data and necessary setting
		for _ss, _colr in zip(ss[0:6],color):
			artist, = ax.plot(n.arange(tm_num)[ss==_ss],data[ss==_ss],mec=_colr,ls='',ms=5.,
					mfc='#ffffff',mew=2.,marker='o',label=round(_ss,3))
			collectArt.append(artist)

		ax.tick_params(which='major',direction='in',length=7,right=True,labelsize=fs-2.5)
		# ax.tick_params(which='minor',direction='in',length=4.5)
		ax.set(xlim=(0.,tm_num),ylim=(-.05,1.),xticks=xTick)

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
			cax = fig.add_axes([.935,box.y0,.015,box.height])
			cax.set_axis_off()

			cax.legend(collectArt,[_art.get_label() for _art in collectArt],framealpha=0,fontsize=fs-2.,
					   title='SS (%)',title_fontsize=fs-3.,loc=10)
			ax.set_xticklabels('')
			ax.set_title(r'$\kappa$',fontsize=fs-1.)

			return xTickLab

	def plot_together(self,**kwarg):
		## set plot parameter
		default = {'splt_hr' : 12,
				   'order'   : ['smps','cpc','kappa']}

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

		fig.savefig(pth(self.figPath,'mesr_together2date.png'))


# xtick  ylabel  sep together and non


