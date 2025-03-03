#example syntax: python spec_salt.py 20210501 826 1026
#example syntax if potential PN: python spec_salt.py 20210501 826 1026 PN
from spectroscopic_routines import *
# TED ADDITION: COMMAND-LINE ARGS FOR OBSERVATION DATES
import datetime as dt
import sys
rawdate = sys.argv[1]
data_date = str(sys.argv[1])
print('OBSERVATION DATE: ', data_date)
os.system('mkdir '+data_date)
os.system('mkdir /Users/asishphilipmonai/SALT/reduced/'+data_date)
yr = int(rawdate[0:4])
mn = int(rawdate[4:6])
dy = int(rawdate[6:8])

# TED ADDITION: copy files from salt.data folder to working folder
folder = '/Users/asishphilipmonai/SALT/raw_data/'+data_date+'/product/*'
os.system('cp /Users/asishphilipmonai/SALT/raw_data/'+data_date+'/product/* /Users/asishphilipmonai/SALT/reduction/titus_saurus_rex')
# ADDITION END

d_obs1 = dt.date(yr, mn, dy)
d_obs2 = d_obs1 + dt.timedelta(days=1)

pre_mn = str(d_obs1)
post_mn = str(d_obs2)
time_keyword = 'TIME-OBS'
date_keyword = 'DATE-OBS'
reminder_flag = False
# ADDITION END

run_pipeline = True
apply_flux_cal = False

raw_files_prefix = 'mbxgpP'

grating = 'PG2300'
grating_keyword = 'GRATING'
# TED ADDITION: STUFF FOR DETECTING GRATING ANGLE FROM HEADER
angle_keyword = 'GRTILT'
arc_keyword = 'OBSTYPE'
block_keyword = 'BVISITID'
angle1 = 32.0
angle2 = 30.5
# ADDITION END

identifying_keyword = 'OBSTYPE' # identify science / flat / etc
exp_time_keyword = 'EXPTIME'
science_keyword = 'OBJECT'
object_keyword = 'OBJECT'

gain_keyword = 1.0 #1.507 #'GAINVAL'
readnoise_keyword = 2.225 #'NOISEADU'

spatial_axis = 1
flip_wave_axis = False

#############################
#############################
IMAGES =  raw_files_prefix+'*.fits' # or you can a single image
#ARC = 1	# 1: takes arc after, -1: takes arc before or can take file name eg: 'arc.fits'

#############################
#### Trimming Parameters ####
#############################

Trim = True
xmin = 1
xmax = 3160
ymin = int(sys.argv[2])
ymax = int(sys.argv[3])
Trim_file_prefix = None
Trim_ouput_prefix = 't'

################################
#### Calibration Parameters ####
################################

apply_call = True,
files_prefix = None
cal_file_prefix = 't'
masterbias = None
masterflat = None
cosmicray_removal = True
contrast = 4
cr_threshold = 8
neighbor_threshold = 5
error = None
mask = None
background = None
maxiter = 4
border_mode = 'mirror'


################################
###### Trace Parameters ########
################################

#y0_trace = (ymax-ymin)/2
#yf_trace = (ymax-ymin)/2
trace_prominence = 200
manual_trace = False#True
manual_poly_order = 3
manual_x = [0,1020,1073,2096,2145,3165] # 0 1020 1073 2096 2145 3165
manual_y = [100,100,100,100,100,100] #101.5 100.5 100.5 99.5 99.5 98
tolerance = 1
trace_width = '2.5-fwhm' # or for a factor of 2 of fwhm: '2-fwhm' or just integers for pixels
poly_order = 3
fmask = (1,)
nsteps = 20
recenter = True
prevtrace = (0,)
bigbox = 15
Saxis = spatial_axis
display_trace = False

################################
## Optimal Extract Parameters ##
################################

display_optimal = False

################################
##### Apextract Parameters #####
################################

skysep = 5
skywidth = 10
skydeg = 0
coaddN = 1
display_apextract = False

##########################################
##### Wavelength Solution Parameters #####
##########################################

#reference_spec ='calibrated_arc_cuar32.dat'
line_list = 'calib_cuar.dat'
wave_min = 3800
wave_max = 5300
prominence = 10
order = 2
#parameter_file = 'parameters_new'
view_arc = False
display_wave = False

##########################################
######### Flux Calibration ###############
##########################################

display_flux_cal = False

##########################################
## APM: PN Sky background emmision test ##
##########################################

try:
	if sys.argv[4] == 'PN':
		PN = True
	else:
		PN = False
except IndexError:
	PN = False

##########################################
############## Pipeline ##################
##########################################


if run_pipeline is True:

	if type(IMAGES) is str:
		IMAGES = glob.glob(IMAGES)
		#APM ADDITION: Ignore bp_ag_ff_cr_cg_wr and Flat files; else TypeError
		IMAGES = [x for x in IMAGES if 'wr' not in x]
		IMAGES = [x for x in IMAGES if 'Flat' not in x]
		# print(len(IMAGES))
		# for el in IMAGES:
		# 	print(el)
		# exit()
	trim_images = []

	# TED ADDITION: ONLY SELECT IMAGES FROM CORRECT DATE
	for k in IMAGES:
		frame_date = image_header(k,date_keyword)
		frame_time = image_header(k,time_keyword)
		timedex = int(frame_time[0:2])
		if 12 <= timedex <= 23:
			if frame_date != pre_mn:
				IMAGES.remove(k)
		elif frame_date != post_mn:
			IMAGES.remove(k)
	# IDENTIFY ARC FRAMES
	arc_frames = []

	for k in IMAGES:
		arc_id = image_header(k,arc_keyword)
		ang_no = image_header(k,angle_keyword)
		block_id = image_header(k,block_keyword)
		if arc_id == 'ARC':
			arc_object = image_header(k,object_keyword)
			arc_angle = image_header(k,angle_keyword)
			arc_frames.append([block_id,arc_angle,k])
	#print("arc_frames = ",arc_frames)
	#ADDITION END

	for k in IMAGES:
		gr = image_header(k,grating_keyword)
		if gr == grating:
			trim_images.append(k)

	for k in trim_images:
	    trim(
			apply_trim = Trim,
	    	files_prefix = None,
	    	file = k,
	    	x1 = xmin,
	    	x2 = xmax,
	    	y1 = ymin,
	    	y2 = ymax,
	    	output_prefix = Trim_ouput_prefix)

	for k in IMAGES:
		# TED ADDITION: SWITCH CALIBRATION BASED ON OBJECT AND GRATING ANGLE
		frame_block = image_header(k,block_keyword)
		frame_angle = image_header(k,angle_keyword)
		print("frame_block = ", frame_block)
		print("frame_angle = ", frame_angle)
		print("k = ",k)
		if frame_angle == angle1:
			reference_spec = 'calibrated_arc_cuar32.dat'
			#ARC = arc32
			parameter_file = 'parameters_32'
		elif frame_angle == angle2:
			reference_spec = 'calibrated_arc_cuar305.dat'
			#ARC = arc305
			parameter_file = 'parameters_305'
		for j in arc_frames:
			if j[0] == frame_block:
				if j[1] == frame_angle:
					ARC = j[2]

		# FIND TRACE CENTRE
		if image_header(k, arc_keyword) != 'ARC':
			tf_fits = fits.open(k)
			tf_data = tf_fits['SCI'].data
			centre_col = int(tf_data.shape[1]/2)
			centre_region = tf_data[ymin:ymax,centre_col-100:centre_col+100]
			collapsed = np.empty([centre_region.shape[0],1])
			for i in range(0, len(collapsed)-1):
				collapsed[i] = np.mean(centre_region[i])
			collapsed = np.delete(collapsed, -1) #stops a weird error with certain frames where the final element of collapsed gets set to like 300
			mid_peak = np.argmax(collapsed)
			mid_peak += ymin
			search_lo = mid_peak - 30
			search_hi = mid_peak + 30
			peaks_0 = []
			peaks_f = []
			for i in range(1,11):
				tf_flux_0 = tf_data[search_lo:search_hi, i]
				tf_flux_f = tf_data[search_lo:search_hi, -i]
				peaks_0.append(np.argmax(tf_flux_0))
				peaks_f.append(np.argmax(tf_flux_f))
			peaks_0 = [i for i in peaks_0 if i != 0] # peaks of 0 skew trace downwards
			peaks_f = [i for i in peaks_f if i != 0]
			shift = search_lo-ymin
			y0_trace = (np.median(peaks_0)+shift)
			yf_trace = (np.median(peaks_f)+shift)
			#print('vvvvvvvvvvvvvvvvvvvv')
			#print(k)
			#print('centre_col', centre_col)
			#print('mid_peak', mid_peak)
			#print('search_lo/hi', search_lo, search_hi)
			#print('y0/f_trace', y0_trace, yf_trace)
			#print(peaks_0, peaks_f)
			#print('^^^^^^^^^^^^^^^^^^^^')
			#dummy = np.arange(len(collapsed))
			#plt.plot(dummy, collapsed)
			#plt.title(k)
			#plt.show()
##########################################################
# IF IT WOULD BE EASIER TO JUST FORCE y0, yf, DO IT HERE!

			# y0_trace = 100
			# yf_trace = 100
			# reminder_flag = True
##########################################################
		# ADDITION END
		gr = image_header(k,grating_keyword)
		if gr == grating:
			appy_calibrations(apply_call = apply_call,
	    		identifying_keyword = identifying_keyword,
	    		science_keyword = science_keyword,
	    		files_prefix = None,
	    		file = cal_file_prefix+k,
	    		masterbias = masterbias,
	    		masterflat = masterflat,
	    		cosmicray_removal = cosmicray_removal,
	    		contrast = contrast,
	    		cr_threshold = cr_threshold,
	    		neighbor_threshold = neighbor_threshold,
	    		error = error,
	    		mask = mask,
	    		background = background,
	    		gain_keyword = gain_keyword,
	    		readnoise_keyword = readnoise_keyword,
	    		maxiter = maxiter,
	    		border_mode = border_mode)
	    	#
			exposure_type = image_header(k,identifying_keyword)
			# print(exposure_type,science_keyword)
			if exposure_type == science_keyword:
				# print(k,exposure_type)
				number = int(k.split('.')[0].split(raw_files_prefix)[1])
				if type(ARC) is int:
					arc = 't'+raw_files_prefix+str(number+ARC)+'.fits'
				if type(ARC) is str:
					arc = 't'+ARC
				spec_file = k.split('.fits')[0]
				k = 'ct'+k
				# print (k)
				# input()
				try:
					image_raw, sky_subtracted, sky, xs, ys, nx, ny, yvals = twodsky(k,object_keyword=object_keyword)
					optimal_spec = True
					#
				except ValueError:
					print ('Fit unsuccessful for '+k+' in twodsky, cannot complete optimal extraction')
					optimal_spec = None
				#TED ADDITION: My slightly lazy workaround for the 'index out of range' problem for some frames
				except IndexError:
					print ('Fit unsuccessful for '+k+' in twodsky, cannot complete optimal extraction')
					optimal_spec = None
				#ADDITION END
				trace_output = ap_trace(
					k,
					y0_trace = y0_trace,
					yf_trace = yf_trace,
					trace_prominence = trace_prominence,
					manual_trace = manual_trace,
					manual_poly_order = manual_poly_order,
					manual_x = manual_x,
					manual_y = manual_y,
					tolerance = tolerance,
					trace_width = trace_width,
					poly_order = poly_order,
					object_keyword = object_keyword,
					fmask = fmask,
					nsteps = nsteps,
					recenter = recenter,
					prevtrace = prevtrace,
					bigbox = bigbox,
					Saxis = Saxis,
					display = display_trace)
				my, myfwhm, trace_c = trace_output
				toast = np.arange(200)
				if optimal_spec is True:
					try:
						optimal_spec = optimal(
							image_raw,
							sky = sky,
							xs = xs,
							ys = ys,
							nx = nx,
							ny = ny,
							yvals = yvals,
							trace_c = trace_c,
							display = display_optimal)
					except:
						optimal_spec = None

	    		#
				onedspec, fluxerr, variancespec = ap_extract(
					k,
					trace = my,
					apwidth = round(myfwhm),
					skysep = skysep,
					skywidth = skywidth,
					skydeg = skydeg,
					coaddN = coaddN,
					gain_keyword = gain_keyword,
					readnoise_keyword = readnoise_keyword,
					object_keyword = object_keyword,
					display = display_apextract)
	    		#
				interact,cal_wav = wavelength(						#APM:cal_wav -- calibrated wavelength (np.ndarray)
					onedspec,
	    			onedspec_optimal = optimal_spec,
	    			spec_file_name = spec_file,
	    			arc_file = arc,
	    			reference_spec = reference_spec,
	    		    line_list = line_list,
	    		    trace = my,
	    		    trace_fwhm = myfwhm,
	    		    wave_min = wave_min,
	    		    wave_max = wave_max,
	    		    prominence = prominence,
	    		    order = order,
	    		    parameter_file = parameter_file,
	    		    object_keyword = object_keyword,
	    		    flip_wave_axis = flip_wave_axis,
	    		    view_arc = view_arc,
	    		    display = display_wave)
				view_arc = interact

				if PN:															
					try:
						#APM: sky background with wavlength calibration as fits
						header = fits.getheader(k)
						header.set('FWHM', myfwhm, 'Trace FWHM measured during reduction (pix.)')
						fitsname = 'r'+k[-8:-5]+'_em.fits'												
						fitsdata = [cal_wav,np.mean(sky,axis = 0)]																
						savename1 = os.path.join(data_date, fitsname)					
						savename2 = os.path.join('/Users/asishphilipmonai/SALT/reduced/'+data_date, fitsname)								
						fits.writeto(savename1, fitsdata, header=header, overwrite=True)						
						fits.writeto(savename2, fitsdata, header=header, overwrite=True)
						continue
						#APM: Plot sky background with wavlength calibration						
						# fig, (ax1, ax2) = plt.subplots(2, 1)
						# fig.set_size_inches(18, 7)

						# ax1.plot(cal_wav,np.mean(sky,axis = 0),'-k')
						# ax1.set_title(target_name+' sky background')
						# ax1.axvline(3968,ls = 'dotted',color='r')
						# ax1.axvline(4340,ls = 'dotted',color='r')
						# ax1.axvline(4686,ls = 'dotted',color='r')
						# ax1.axvline(4861,ls = 'dotted',color='r')
						# ax1.axvline(5007,ls = 'dotted',color='blue')
						# ax1.set_ylabel('Counts')
						# ax1.set_xlim(min(cal_wav)-20,max(cal_wav)+20)

						# ax2.plot(cal_wav,onedspec,'-k')
						# ax2.set_title(target_name+' reduced spectrum')
						# ax2.axvline(3968,ls = 'dotted',color='r')
						# ax2.axvline(4340,ls = 'dotted',color='r')
						# ax2.axvline(4686,ls = 'dotted',color='r')
						# ax2.axvline(4861,ls = 'dotted',color='r')
						# ax2.axvline(5007,ls = 'dotted',color='blue')
						# ax2.set_xlabel('Wavelength $\AA$')
						# ax2.set_ylabel('Counts')
						# ax2.set_xlim(min(cal_wav)-20,max(cal_wav)+20)
						# plt.savefig(target_name+'_'+spec_file[-4:]+'_'+'background.png')
						# plt.show()
						# plt.close()
					except NameError:
						print('*************************************************')
						print('\n\nNameError:', 'r'+k[-8:-5]+'_em.fits', '\n\n')
						savename1 = os.path.join(data_date, fitsname)					
						savename2 = os.path.join('/Users/asishphilipmonai/SALT/reduced/'+data_date, fitsname)								
						print(savename1,'\n',savename2,'\n')
						print('*************************************************')
						continue
				#TED ADDITION: OUTPUT AS .FITS
				header = fits.getheader(k)
				header.set('FWHM', myfwhm, 'Trace FWHM measured during reduction (pix.)')
				readname = header[object_keyword]+'_'+spec_file+'_reduced.dat'
				wvl = np.empty(0)
				spc = np.empty(0)
				#sm3 = np.empty(0)
				spc_opt = np.empty(0)
				#sm3_opt = np.empty(0)

				with open(readname, 'r') as file:
					for line in file:
						splitstr = line.split(' ')
						try:
							wvl = np.append(wvl, float(splitstr[0]))
							spc = np.append(spc, float(splitstr[1]))
							#sm3 = np.append(sm3, float(splitstr[2]))
							if hasattr(optimal_spec, '__len__'):
								spc_opt = np.append(spc_opt, float(splitstr[3]))
								#sm3_opt = np.append(sm3_opt, float(splitstr[4]))
						except ValueError:
							pass
				fitsname = 'r'+k[-8:-5]+'.fits'
				if hasattr(optimal_spec, '__len__'):
					#APM: addition to check if optimal fit is too noisy compared to summed fit.
					SN_spc = spc/np.sqrt(variancespec)
					SN_spc_opt = spc_opt/np.sqrt(variancespec)

					p_spc = np.polyval(np.polyfit(wvl, SN_spc, deg=3),wvl)
					p_spc_opt = np.polyval(np.polyfit(wvl, SN_spc_opt, deg=3),wvl)

					if abs(np.std(SN_spc_opt - p_spc_opt) - np.std(SN_spc - p_spc)) > 1.5:
						print(np.std(SN_spc_opt - p_spc_opt) - np.std(SN_spc - p_spc))
						plt.plot(wvl,spc/np.sqrt(variancespec),label='S/N summed')
						plt.plot(wvl,spc_opt/np.sqrt(variancespec),label='S/N optimal')
						plt.plot(wvl,p_spc,label='S/N summed fit')
						plt.plot(wvl,p_spc_opt,label='S/N optimal fit')
						# plt.hlines(np.std(SN_spc - p_spc),min(wvl),max(wvl),ls = 'dotted',color='k')
						# plt.hlines(np.std(SN_spc_opt - p_spc_opt),min(wvl),max(wvl),ls = 'dotted',color='r')
						plt.legend()
						plt.show()
						mask_spc_opt = input('Should optimal trace be ignored? (Y/N)  ')
						if mask_spc_opt in ['y','Y']:
							fitsdata = [wvl, spc]
						else:
							fitsdata = [wvl, spc, spc_opt]
					else:
						fitsdata = [wvl, spc, spc_opt]
				else:
					fitsdata = [wvl, spc]
				savename1 = os.path.join(data_date, fitsname)
				savename2 = os.path.join('/Users/asishphilipmonai/SALT/reduced/'+data_date, fitsname)
				fits.writeto(savename1, fitsdata, header=header, overwrite=True)
				fits.writeto(savename2, fitsdata, header=header, overwrite=True)
				#fitsdata.writeto(savename, overwrite=True)
				#ADDITION END
	    		    #



if apply_flux_cal is True:
	reduced_data_files = glob.glob('*reduced.dat')
	standard = []
	standard_name = []
	for k in reduced_data_files:
		if k.split('_')[0].lower() == 'ltt3218':
			standard.append(k)
			standard_name.append('ltt3218')
		if k.split('_')[0].lower() == 'eg21':
			standard.append(k)
			standard_name.append('eg21')
		if k.split('_')[0].lower() == 'ltt377':
			standard.append(k)
			standard_name.append('ltt377')
		if k.split('_')[0].lower() == 'feige110':
			standard.append(k)
			standard_name.append('feige110')
		if (k.split('_')[0].lower() == 'cd-32-9927') or (k.split('_')[0].lower() == 'cd-32d9927'):
			standard.append(k)
			standard_name.append('cd-32-9927')
		if k.split('_')[0].lower() == 'ltt7379':
			standard.append(k)
			standard_name.append('ltt7379')
		if k.split('_')[0].lower() == 'ltt7987':
			standard.append(k)
			standard_name.append('ltt7987')

	if len(standard) != None:
		for k in reduced_data_files:
			print (k)
			flux_callibration(
		    	standard_reduced_spec = standard[0],
		    	standard_name = standard_name[0],
		    	science_spec = k,
		    	display = display_flux_cal)
	else:
		print('Cannot perform flux callibration')

# TED ADDITION: TIDY FILES INTO REDUCED FOLDER
# APM: if PN only em.fits files moved into folder
if not PN:
	os.system("mv *_traced_arc.png "+data_date+'/')
	os.system("mv *_wavelength_sol.png "+data_date+'/')
	os.system("mv *_reduced.dat "+data_date+'/')
	os.system("mv *_reduced_spectrum.png "+data_date+'/')
	os.system("mv *_fitted_sky.png "+data_date)
	os.system("mv *_Initial_1D_spectrum.png "+data_date+'/')
	os.system("mv *_trace_fit2D.png "+data_date+'/')
	os.system("mv *_trace_fit.png "+data_date+'/')
	os.system("mv *_trace_region.png "+data_date+'/')

	os.system('mv tmbxgpP* '+data_date)
	#os.system('mv cmbxgpP* '+data_date)
	os.system('mv ctmbxgpP* '+data_date)
	os.system('mv sctmbxgpP* '+data_date)

	textpath1 = os.path.join(data_date, 'read_groups.txt')
	textpath2 = os.path.join('/Users/asishphilipmonai/SALT/reduced/'+data_date, 'read_groups.txt')
	PRODUCTS = os.path.join(data_date, 'r*.fits')
	PRODUCTS = glob.glob(PRODUCTS)
	#sort products
	PRODUCTS = sorted(PRODUCTS, key = lambda x: float(x[10:13]))

	#get list of objects
	objects = []
	for k in PRODUCTS:
		obj_id = image_header(k,object_keyword)
		if obj_id not in objects:
			objects.append(obj_id)

	#find which fits matches each object
	with open(textpath1, 'w') as file:
		for j in objects:
			if j[0] == 'J':
				frames_str = "SALT "+j+","+data_date+","
			else:
				frames_str = j+","+data_date+","
			frames_str += "2,2"
			for k in PRODUCTS:
				obj_id = image_header(k,object_keyword)
				if obj_id == j:
					frames_str += (","+k[9:-5])
			frames_str += "\n"
			file.write(frames_str)

	os.system('cp '+textpath1+' '+textpath2)

#Reminder to disable hardcoded y0 and yf
if reminder_flag:
	print('+------------------------------+')
	print('|                              |')
	print('| REMEMBER TO RESET Y0 AND YF! |')
	print('|                              |')
	print('+------------------------------+')
#Reminder to disable manual trace
if manual_trace:
	print('+---------------------------------+')
	print('|                              	 |')
	print('| REMEMBER TO RESET MANUAL TRACE! |')
	print('|                              	 |')
	print('+---------------------------------+')
# ADDITION END

# TED ADDITION: remove previously-copied files from salt.data
os.system('rm mbxgpP*')
# ADDITION END
