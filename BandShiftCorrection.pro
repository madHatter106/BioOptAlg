; docformat:'rst'

;+
; A group of procedures and functions to apply a given set of bandshift corrections to remote sensing reflectances,
; calculating IOPs at the desired wavelengths using a spectral model that starts from -preferably QAA(v5)derived- IOPs at
; reference wavelengths, and feeding the obtained IOPs into QAA(v5) in forward mode to calculate corrections factors.  
;
; :version: 1
;
; :Author:
; 
;   Gert Sclep (European Commission/JRC/IES)
;			
;   Frederic Melin (European Commission/JRC/IES)
;
; :Copyright:
; 
;   All rights reserved
;   
;-   

;+
; :Author:
; 
;   Gert Sclep
;
; :version: 1
;
; :Categories:
; 
;   RRS, bandshift correction, setup
;   
; :Description:
; 
;   Sets up a structure containing information used during bandshift correction of remote sensing reflectances.
;   
;   Contains::
;   
;     input and output wavelengths
;     
;     Bricaud coefficients for those wavelengths (See Bricaud et al. Variability in the chlorophyll-specific absorption 
;     coefficients of natural phytoplankton: Analysis and parameterization. Journal of Geophysical Research. 1995;100(95):321–332)
;     
;     pure water IOPs for those wavelengths  
;     
;     wavelength for which IOPs are fed into the spectral model (deriving IOPs at input and output wavelengths, and using the IOP 
;     ratio to calculate the correction factor, see ATBD), the current implementation uses the blue wavelength (443) as starting point,
;     if one wants to experiment with other starting points (possibily different for each individual bandshift), then changes must 
;     be made within the code. 
;     
; :Pre:
;
;   Each input wavelength (lambda_i) in `correction_context` must be represented once and only once in `rrs_wavelengths`
;   
;   `correction_context` may not contain trivial corrections (output wavelength of lambda_o equals input wavelength of lambda_i) 
;   
;   `correction_context` must at least contain one correction (one elements in lambda_i and lambda_o)
;     
; :Requires:
; 
;   IDL 8.0
;   
; :Uses:
; 
;   readSpectralDataPureWater, getSpectralDataPureWater, get_A_B_Bricaud
;   
; :Params:
; 
;   pw_dir: in, type=string
;			
;      directory containing IOPs for pure water
;			
;   correction_context: out, type=structure
;   
;     Structure containing following fields::
;     
;       sensor: name of the sensor, which determines the input and output wavelengths of the individual corrections
;     
;       lambda_i: array of input wavelengths of to be applied bandshift corrections
;       
;       lambda_o: array of output wavelengths of to be applied bandshift corrections
;       
;       a_i: Bricaud A coefficient (Bricaud et al, 1995, table 2 - see above) for input wavelengths
;       
;       b_i: Bricaud B coefficient (Bricaud et al, 1995, table 2 - see above) for input wavelengths
;       
;       a_o: Bricaud A coefficient (Bricaud et al, 1995, table 2 - see above) for output wavelengths
;       
;       b_o: Bricaud B coefficient (Bricaud et al, 1995, table 2 - see above) for output wavelengths
;       
;       spec_model_start: the wavelength(s) to be used as starting point for spectral model, a scalar if the same is used for each 
;       correction (current implementation), an array with length the number of corrections if different starting wavelengths are used 
;       for the different corrections 
;       
;       a_sms: the Bricaud A coefficient(s) (Bricaud et al, 1995, table 2 - see above) for the chosen spectral model start wavelength(s)
;       
;       b_sms: the Bricaud B coefficient(s) (Bricaud et al, 1995, table 2 - see above) for the chosen spectral model start wavelength(s)
;       
;       aw_i: pure seawater absorption coefficient (m-1) for input wavelengths
;       
;       bbw_i: pure seawater bacskcattering coefficient (m-1) for input wavelengths
;       
;       aw_o: pure seawater absorption coefficient (m-1) for output wavelengths
;       
;       bbw_o: pure seawater bacskcattering coefficient (m-1) for output wavelengths
;              
; :Keywords:
;   sensor: in, optional, type=string, default='MODISA'
;     Name of the satellite sensor for which the bandshift correction structure should be set up. A different sensor name causes different 
;     input and output wavelengths to be selected. 
;   
;-
PRO BandShiftCorrectionSetup, pw_dir, correction_context, sensor=sensor
  
  IF (NOT(KEYWORD_SET(sensor))) THEN sensor='MODISA' ; default sensor
  
  error_status = 11
  
  CASE STRUPCASE(sensor) OF
    'SEAWIFS': BEGIN  
      lambda_i = [412., 490., 510., 555., 555., 555., 670., 670.]
      lambda_o = [413., 488., 531., 531., 547., 560., 665., 667.]
    END
    'MERIS': BEGIN
      lambda_i = [413.,490., 510., 560., 560., 560., 665., 665.]
      lambda_o = [412.,488., 531., 531., 547., 555., 667., 670.]
    END
    'MODISA': BEGIN
      lambda_i = [412., 488., 488., 531., 547., 547., 667., 667.]
      lambda_o = [413., 510., 490., 510., 560., 555., 665., 670.]
    END
  ENDCASE 
  
  subset = INDGEN(N_ELEMENTS(lambda_i))
  readSpectralDataPureWater,pw_dir,sdata ; reads IOPs for pure seawater from directory pw_dir
  
  number_i = N_ELEMENTS(lambda_i)
  number_o = N_ELEMENTS(lambda_o)
  aw_i = FLTARR(number_i)
  bbw_i = FLTARR(number_i)
  a_i =  FLTARR(number_i)
  b_i = FLTARR(number_i)
  aw_o = FLTARR(number_o)
  bbw_o = FLTARR(number_o)
  a_o =  FLTARR(number_o)
  b_o = FLTARR(number_o)
; determine IOPs for pure water and Bricaud parameters at specified wavelengths
  FOR i=0,number_i-1 DO BEGIN 
    status = getSpectralDataPureWater(sdata,lambda_i[i],aw,bw,bandwidth=0)
    IF status EQ 0 THEN BEGIN
      aw_i[i] = aw
      bbw_i[i] = bw/2.0
      get_A_B_Bricaud,lambda_i[i],a_bricaud, b_bricaud
      a_i[i] = reform(a_bricaud)
      b_i[i] = reform(b_bricaud)
    ENDIF ELSE EXIT, STATUS=error_status
  ENDFOR
  FOR i=0,number_o-1 DO BEGIN 
    status = getSpectralDataPureWater(sdata,lambda_o[i],aw,bw,bandwidth=0)
    IF status EQ 0 THEN BEGIN
      aw_o[i] = aw
      bbw_o[i] = bw/2.0
      get_A_B_Bricaud,lambda_o[i],a_bricaud, b_bricaud
      a_o[i] = reform(a_bricaud)
      b_o[i] = reform(b_bricaud)     
    ENDIF ELSE EXIT, STATUS=error_status
  ENDFOR
  
  ;spectral model is applied from blue wavelength onwards, for each of the corrections.
  spec_model_start = 443.
  get_A_B_Bricaud,spec_model_start,a_sms, b_sms
          
  correction_context  = {sensor:sensor, $
                         lambda_i: lambda_i[subset], aw_i: aw_i[subset], bbw_i:bbw_i[subset], a_i:a_i[subset], b_i:b_i[subset],$
                         lambda_o: lambda_o[subset], aw_o: aw_o[subset], bbw_o:bbw_o[subset], a_o:a_o[subset], b_o:b_o[subset],$
                         spec_model_start:spec_model_start,a_sms:a_sms, b_sms:b_sms}
END

;+
; :Author:
; 
;   Gert Sclep
;			
;	Frederic Melin
;
; :version: 1
;   
; :Categories:
; 
;   RRS, bandshift correction
;   
; :Description:
; 
;   Bandshift corrects -fully or partially- a given array or matrix of remote sensing reflectances RRS using IOPs at input
;   and output wavelength derived through a spectral model from values at a reference wavelength, which are then
;   fed in the forward mode of the QAAv5 model. 
;       
; :Pre:
;
;   `correction_context` is created according to `BandShiftCorrectionSetup`
;   
;   Each input wavelength (lambda_i) in `correction_context` must be represented once and only once in `rrs_wavelengths`
;   
;   `correction_context` may not contain trivial corrections (output wavelength of lambda_o equals input wavelength of lambda_i) 
;   
;   `correction_context` must at least contain one correction (one elements in lambda_i and lambda_o)
;
;   `rrs_wavelengths` must contain at least 443nm and the green band (see `sensor`)
;   
;   The order of IOPs in 1st dimension of `qaa_matrix` must be respected: first column: phytoplankton absorption, second column: 
;   detritus-gelbstoff absorption, last column: particle backscattering
;      
;   The resolution of the (QAA(v5) derived) IOP and RRS bin data are supposed to be equal, any resolution conversion need to be performed beforehand
;
; :Requires:
; 
;   IDL 8.1
;   
; :Uses:
; 
;   intersectionSortedUniqueSets, iopSpectralModel
;   
; :Returns:
;   
;   -1 in case an input inconsistency or incompleteness is detected. Currently explicitly checked inconsistencies or incompletenesses::
;   
;     * an input wavelength (lambda_i field) of `correction_context` is is not represented in `rrs_wavelengths`
;     
;     * an input wavelength (lambda_i field) of `correction_context` is represented multiple times in `rrs_wavelengths`
;     
;     * `rrs_wavelengths` does not contain a blue (443) band
;     
;     * `rrs_wavelengths` does not contain a green (SeaWiFS: 555, MODISA: 547, MERIS: 560) band
;     
;   number of bandshift corrected bins, possibly zero, in case no inconsistency or incompleteness is detected. Equals the number of 
;   bins in the intersection of `qaa_bins` and `rrs_bins` that have valid (QAA(v5) derived) IOPs.
;     
; :ToDo:
; 
;   Add debug statements
;   
; :Params:
;
;   rrs_matrix: in, required, type="fltarr(m)/fltarr(m,n)"
;   rrs_wavelengths: in, required, type = fltarr(m)
;   rrs_bins: in, out, required, type="long/lonarr(n)"
;   qaa_matrix: in, required, type="fltarr(3)/fltarr(3,o)/fltarr(3,o,p)"
;   rrs_matrix_corrected: out, type="fltarr(p)/fltarr(p,q)"
;   rrs_wavelengths_corrected: out, type=fltarr(p)
;   correction_context: in, required, type=structure
;
; :Keywords:
;
;   qaa_bins: in, optional, type="long/lonarr(o)", default=:rrs_bins:
;   qaa_min: in, optional, type=float, default=0.
;   qaa_max: in, optional, type=float, default=100.
;   correction_factors: out, type="fltarr(p)/fltarr(p,q)"
;   qaa_valid_index: out, type="lonarr"
;   rrs_valid_index: out, type="lonarr"
;-
FUNCTION BandShiftCorrectionCore, $
rrs_matrix, rrs_wavelengths, rrs_bins, qaa_matrix, correction_context, $
rrs_matrix_corrected, rrs_wavelengths_corrected, $
qaa_bins=qaa_bins, qaa_min=qaa_min, qaa_max = qaa_max, correction_factors=correction_factors, $
qaa_valid_index=qaa_valid_index, rrs_valid_index = rrs_valid_index

; min/max IOP values
IF (NOT(KEYWORD_SET(qaa_min))) THEN qaa_min = 0.
IF (NOT(KEYWORD_SET(qaa_max))) THEN qaa_max = 100.
IF (NOT(KEYWORD_SET(qaa_bins))) THEN qaa_bins = rrs_bins ; IOP bins are assumed those of Rrs if not provided
sensor = correction_context.sensor

; IOP indexes
aph_index  = 0 
acdm_index  = 1
bbp_index = 2

n_qaa_bins = N_ELEMENTS(qaa_bins)
qaa_matrix_dimension = SIZE(qaa_matrix,/dimension)
qaa_matrix_n_dimensions = N_ELEMENTS(qaa_matrix_dimension)
number_correction = N_ELEMENTS(correction_context.lambda_i)
; print,'number_correction = ',number_correction
n_start_bands = N_ELEMENTS(correction_context.spec_model_start)

spec_model_start = correction_context.spec_model_start
a_start = correction_context.a_sms
b_start = correction_context.b_sms
IF qaa_matrix_dimension[0] NE 3 THEN STOP, 'qaa_matrix must have 1st dimension equal to 3'
IF qaa_matrix_n_dimensions EQ 1 THEN BEGIN
  IF n_qaa_bins NE 1 THEN STOP, 'qaa_matrix must have 2nd dimension equal to number of (QAA(v5) derived) IOP bins'
  IF n_start_bands GT 1 THEN STOP, 'qaa_matrix must have 3rd dimension equal to number of spectral model start bands'
  qaa_matrix = REBIN(qaa_matrix,3,1,number_correction)
  spec_model_start = REBIN([spec_model_start],number_correction)
  a_start = REBIN(a_start,number_correction)
  b_start = REBIN(b_start,number_correction)
ENDIF ELSE BEGIN
  IF qaa_matrix_n_dimensions EQ 2 THEN BEGIN
    IF n_start_bands GT 1 THEN STOP, 'qaa_matrix must have 3rd dimension equal to number of spectral model start bands'
    IF ( qaa_matrix_dimension[1] NE n_qaa_bins ) THEN STOP, 'qaa_matrix must have 2nd dimension equal to number of (QAA(v5) derived) IOP bins'
    qaa_matrix = REBIN(qaa_matrix,3,n_qaa_bins,number_correction)
    spec_model_start = REBIN([spec_model_start],number_correction)
    a_start = REBIN(a_start,number_correction)
    b_start = REBIN(b_start,number_correction)
  ENDIF ELSE BEGIN
    IF qaa_matrix_dimension[2] NE number_correction then STOP, 'qaa_matrix must have 3rd dimension equal to number of spectral model start bands'    
  ENDELSE
ENDELSE

; Conversion routine needs RRS at blue and green wavelengths. The blue wavelength
; is always 443, however the green wavelength differs according to the sensor.
blue_index = WHERE(rrs_wavelengths EQ 443,nb_blue)
IF nb_blue NE 1 THEN RETURN, -1 ELSE blue_index = blue_index[0]
CASE STRUPCASE(sensor) OF
  'MODISA': BEGIN
    green_index = WHERE(rrs_wavelengths EQ 547,nb_green)
    IF nb_green NE 1 THEN RETURN, -1 ELSE green_index = green_index[0]
   END
  'MERIS': BEGIN
    green_index = WHERE(rrs_wavelengths EQ 560,nb_green)
    IF nb_green NE 1 THEN RETURN, -1 ELSE green_index = green_index[0]
   END
  'SEAWIFS': BEGIN
    green_index = WHERE(rrs_wavelengths EQ 555,nb_green)
    IF nb_green NE 1 THEN RETURN, -1 ELSE green_index = green_index[0]
   END
ENDCASE

; Determine the indexes of the correction input products and create the correction output product names
input_wavelength_indexes = INTARR(number_correction)
rrs_wavelengths_corrected = FLTARR(number_correction)
FOR li=0,number_correction-1 DO BEGIN
  rrs_prod_position = WHERE(rrs_wavelengths EQ correction_context.lambda_i[li], nb_pos)
  IF nb_pos NE 1 THEN RETURN, -1  
  input_wavelength_indexes[li] = rrs_prod_position
  rrs_wavelengths_corrected[li] = correction_context.lambda_o[li]
; print,input_wavelength_indexes[li],correction_context.lambda_i[li],rrs_wavelengths_corrected[li]
ENDFOR

;constants used in the calculation of the below-water-reflectance
g0 = 0.089
g1 = 0.125

; finds intersection of RRS and QAA bins
intersect_exists = intersectionSortedUniqueSets(rrs_bins, qaa_bins, rrs_intersect_index,qaa_intersect_index)
;in case the QAA bins and RRS bins have no intersection, no corrections can be carried out and 0 will be returned 
nb_corrected = 0

IF intersect_exists THEN BEGIN
  ; Determine which intersection bins have valid IOP values (GT MIN and LT MAX) 
  qaa_intersect_index_valid = MAKE_ARRAY(N_ELEMENTS(qaa_intersect_index),/INTEGER,VALUE=1)

; filter on IOP values
  IF number_correction GT 1 THEN BEGIN
    invalid_aph = WHERE(MIN(qaa_matrix[aph_index,qaa_intersect_index,*],dimension=3) LE qaa_min OR $
                      MAX(qaa_matrix[aph_index,qaa_intersect_index,*],dimension=3) GE qaa_max,nb_invalid_aph)
    invalid_acdm = WHERE(MIN(qaa_matrix[acdm_index,qaa_intersect_index,*],dimension=3) LE qaa_min OR $
                          MAX(qaa_matrix[acdm_index,qaa_intersect_index,*],dimension=3) GE qaa_max,nb_invalid_acdm)  
    invalid_bbp = WHERE(MIN(qaa_matrix[bbp_index,qaa_intersect_index,*],dimension=3) LE qaa_min OR $
                      MAX(qaa_matrix[bbp_index,qaa_intersect_index,*],dimension=3) GE qaa_max,nb_invalid_bbp)
  ENDIF ELSE BEGIN
    invalid_aph = WHERE(qaa_matrix[aph_index,qaa_intersect_index] LE qaa_min OR $
                      qaa_matrix[aph_index,qaa_intersect_index] GE qaa_max,nb_invalid_aph)
    invalid_acdm = WHERE(qaa_matrix[acdm_index,qaa_intersect_index] LE qaa_min OR $
                          qaa_matrix[acdm_index,qaa_intersect_index] GE qaa_max,nb_invalid_acdm)  
    invalid_bbp = WHERE(qaa_matrix[bbp_index,qaa_intersect_index] LE qaa_min OR $
                      qaa_matrix[bbp_index,qaa_intersect_index] GE qaa_max,nb_invalid_bbp)
  ENDELSE
  
  IF nb_invalid_aph GT 0 THEN qaa_intersect_index_valid[invalid_aph] = 0
  IF nb_invalid_acdm GT 0 THEN qaa_intersect_index_valid[invalid_acdm] = 0
  IF nb_invalid_bbp GT 0 THEN qaa_intersect_index_valid[invalid_bbp] = 0
  
  valid_qaa_intersect_index = WHERE(qaa_intersect_index_valid EQ 1, nb_corrected)
  
  ; Only continue if there are intersection bins with a valid value for all of the used IOPs (aph, acdm, bbp)  
  
  IF nb_corrected GT 0 THEN BEGIN    
    qaa_valid_index = qaa_intersect_index[valid_qaa_intersect_index]
    rrs_valid_index = rrs_intersect_index[valid_qaa_intersect_index]  

    Rrs_i = (rrs_matrix[input_wavelength_indexes,*])[*,rrs_valid_index]  ; FLTARR(number of input bands x number of bins)

    ; Below-water reflectance for blue and green wavelengths
    rrs_blue = REFORM((rrs_matrix[blue_index,*])[*,rrs_valid_index])/(0.52+1.7*REFORM((rrs_matrix[blue_index,*])[*,rrs_valid_index]))
    rrs_green = REFORM((rrs_matrix[green_index,*])[*,rrs_valid_index])/(0.52+1.7*REFORM((rrs_matrix[green_index,*])[*,rrs_valid_index]))

    ; Derive the aph, adg and bbp for the correction input wavelengths starting from the blue band
    IOPSpectralModel, correction_context.spec_model_start, correction_context.a_sms, correction_context.b_sms, $
              TRANSPOSE(REFORM(qaa_matrix[0,qaa_valid_index,*],nb_corrected,number_correction)), TRANSPOSE(REFORM(qaa_matrix[1,qaa_valid_index,*],nb_corrected,number_correction)), $
              TRANSPOSE(REFORM(qaa_matrix[2,qaa_valid_index,*],nb_corrected,number_correction)), $
              rrs_blue, rrs_green, correction_context.lambda_i,correction_context.a_i,correction_context.b_i,$
              aph_in, adg_in, bbp_in
    ; Derive the aph, adg and bbp for the correction output wavelengths starting from the blue band          
    IOPSpectralModel, correction_context.spec_model_start, correction_context.a_sms, correction_context.b_sms, $
              TRANSPOSE(REFORM(qaa_matrix[0,qaa_valid_index,*],nb_corrected,number_correction)),TRANSPOSE(REFORM(qaa_matrix[1,qaa_valid_index,*],nb_corrected,number_correction)), $
              TRANSPOSE(REFORM(qaa_matrix[2,qaa_valid_index,*],nb_corrected,number_correction)), $
              rrs_blue, rrs_green, correction_context.lambda_o,correction_context.a_o,correction_context.b_o,$
              aph_out, adg_out, bbp_out 
    ; Calculate the total absorption and backscattering at correction output wavelengths
    a_tot_out = aph_out + adg_out + rebin(correction_context.aw_o,number_correction,nb_corrected)
    bb_tot_out = bbp_out + rebin(correction_context.bbw_o,number_correction,nb_corrected)
    ; Calculate the total absorption and backscattering at correction input wavelengths
    a_tot_in = aph_in + adg_in + rebin(correction_context.aw_i,number_correction,nb_corrected)
    bb_tot_in = bbp_in + rebin(correction_context.bbw_i,number_correction,nb_corrected)
    ; Using the forward QAA mode, calculate the above water RRS for the correction output wavelengths
    QAA_u_out = bb_tot_out/(a_tot_out + bb_tot_out)        
    QAA_rrs_bw_out = (g0 + g1 * QAA_u_out) * QAA_u_out
    QAA_rrs_aw_out = (-1.*0.52*QAA_rrs_bw_out)/((1.7*QAA_rrs_bw_out)-1)
    ; Using the forward QAA model, calculate the above wrater RRS for the correction input wavelengths
    QAA_u_in = bb_tot_in/(a_tot_in + bb_tot_in)        
    QAA_rrs_bw_in = (g0 + g1 * QAA_u_in) * QAA_u_in
    QAA_rrs_aw_in = (-1.*0.52*QAA_rrs_bw_in)/((1.7*QAA_rrs_bw_in)-1)            
    ;Correction factors that, when multiplied with the RRS at correction input wavelengths give RRS at correction output wavelengths 
    correction_factors = QAA_rrs_aw_out/QAA_rrs_aw_in 
    ;Predict RRS at output wavelengths, multiplying with correction factors
    rrs_matrix_corrected = correction_factors * RRS_i ; FLTARR(number of bands x number of bins)
  ENDIF

ENDIF  
; Reset qaa_matrix to original dimensions, if applicable.
IF qaa_matrix_n_dimensions EQ 2 THEN qaa_matrix = REFORM(qaa_matrix[*,*,0])
IF qaa_matrix_n_dimensions EQ 1 THEN qaa_matrix = REFORM(qaa_matrix[*,0,0])

RETURN, nb_corrected 

END


;+
; :Author:
;
;   Gert Sclep
;
; :version: 1
;   
; :Categories:
; 
;   RRS, bandshift correction, weighted average
;   
; :Description:
; 
;   Given an array or matrix of bandshift corrected remote sensing reflectance (RRS) values
;   -wavelengths in columns and bins in rows-, given an array of corresponding wavelengths, 
;   and given a bandshift correction context containing input and output wavelengths of the 
;   bandshift corrections, the RRS columns corresponding to equal wavelengths are united into one single 
;   column calculating a weighted average. The contribution of each of
;   the single output RRS to the weighted average is inversely proportional to the distance of the 
;   corresponding input wavelength to the output wavelength of the correction context.   
;   
;   An example::
;       
;           correction 1: RRS_i_1 from lambda_i_1 to RRS_o_1 of lambda_o
;           ...
;           correction n (in general at most 2 for RRS bs. corr.): RRS_i_n from from lamda_i_n to RRS_o_n of lambda_o 
;           weighted average: ((RRS_o_1 * (1-abs(lambda_o - lambda_i_1)))+ ... +(RRS_o_n * (1-abs(lambda_o - lambda_i_n))))/
;                             ((n-1)*(abs(lambda_o - lambda_i_1)+ ... + abs(lambda_o - lambda_i_n))) 
;           
;           
;     
;   
; :Pre:
; 
;   There must be a 1 one 1 match between the wavelengths in the input `rrs_corrected_wavelengths` and the output
;   wavelengths (field lambda_o) of `correction_context`.
;   
; :Requires:
; 
;   IDL 7.1
;   
; :Params:
;
;   rrs_corrected_matrix: in, out, required, type="fltarr(p)/fltarr(p,q)"; 
;   rrs_corrected_wavelengths: in, out, required, type="fltarr(p)"
;   correction_context: in, required, type=structure 
;-
PRO weightedAverageEqualCorrectionProducts, rrs_corrected_matrix, rrs_corrected_wavelengths, correction_context

  rrs_dimension = SIZE(rrs_corrected_matrix,/dimension)
  IF ( N_ELEMENTS(rrs_dimension)) GT 1 THEN n_bins = (SIZE(rrs_corrected_matrix,/dimension))[1] ELSE n_bins = 1

  n_rrs_corrected_wavelengths = N_ELEMENTS(rrs_corrected_wavelengths)
 ; print,n_rrs_corrected_wavelengths,rrs_corrected_wavelengths ; MODISA:  413 510 490 510 560 555 665 670

  non_doubles = MAKE_ARRAY(n_rrs_corrected_wavelengths,value=1,/INT)
  sorted_rrs_corrected_wavelengths = rrs_corrected_wavelengths(SORT(rrs_corrected_wavelengths))
  unique_indexes = UNIQ(sorted_rrs_corrected_wavelengths)
  n_unique = N_ELEMENTS(unique_indexes)
  ; Create an array with the first index of all unique elements in sorted_rrs_corrected_wavelengths
  compare_indexes = 0
  IF n_unique GT 1 THEN compare_indexes = [compare_indexes,unique_indexes[0:n_unique-2]+1]
  ; Where the last index of the unique elements in sorted_rrs_corrected_wavelengths does not correspond to the first index
  ; of the unique elements in sorted_rrs_corrected_wavelengths, then this means that elements is multiply present. 
  double_index = WHERE(unique_indexes NE compare_indexes,n_doubles) ; n_doubles=1 for MODISA 

  IF ( n_doubles GT 0 ) THEN BEGIN
    doubles = sorted_rrs_corrected_wavelengths[unique_indexes[double_index]] ; 510 for MODISA

    foreach double_wl, doubles DO BEGIN         
          wavelength_index = WHERE(rrs_corrected_wavelengths EQ double_wl)

          non_doubles[wavelength_index] = 0
          non_doubles[wavelength_index[0]] = 1
          wl_index = WHERE(correction_context.lambda_o EQ double_wl, n_wl)
; check on consistency
          IF TOTAL(wl_index - wavelength_index) NE 0 THEN STOP, 'rrs wavelengths must match one on one to output wavelengths of correction context' 
          input_wl = correction_context.lambda_i[wl_index]  
          rel_weight = ( 1. - ( ABS(input_wl - double_wl)/TOTAL(ABS(input_wl-double_wl)) ) ) / (n_wl - 1.) ; 0.488372  0.511628 for MODISA
          rrs_corrected_matrix[wavelength_index[0],*] = TOTAL(REBIN(rel_weight,n_wl,n_bins) * rrs_corrected_matrix[wavelength_index,*],1) ; sum on the rows
     endforeach
     rrs_corrected_matrix = rrs_corrected_matrix[WHERE(non_doubles EQ 1),*]
     rrs_corrected_wavelengths = rrs_corrected_wavelengths[WHERE(non_doubles EQ 1)]
  ENDIF
END


;+
; :Author:
;  Gert Sclep, F. Melin
;
; :version: 1
;  
; :Categories:
;   RRS, bandshift correction, bin data
;   
; :Description:
;   The package contains the code needed to carry out bandshift corrections of a remote sensing reflectance (RRS), 
;   an array of RRS or a matrix of RRS. A spectral model is used to derive IOPs at input and output wavelengths, 
;   starting from the IOP at reference wavelength (preferably as given by QAA(v5)). And these IOPs are fed into 
;   the forward mode of the QAAv5 model to calculate the bandhsift correction factors. 
; 
; :Examples:
;   Required pre-compilation::
;   
;      .compile Intersection.pro
;      .compile IOPSpectralModel.pro
;      .compile PureWater.pro
;      .compile BandShiftCorrection.pro
;      
;   Example syntax::
;     pw_dir = '.' ; local directory for pure water IOP
;     bandShiftCorrectionSetup, pw_dir, correction_context, sensor='MODISA' 
;     ;read RRS data (rrs_matrix, rrs_wavelengths, rrs_bins) from file
;     rrs_matrix = [[0.00719889,0.00567422,0.00468322,0.00254233,0.00196100,0.000256445],$
;                   [0.00709421,0.00560526,0.00464842,0.00256442,0.00196990,0.000251790],$
;                   [0.00711314,0.00559714,0.00459386,0.00249029,0.00189400,0.000241144]]
;     rrs_wavelengths = [412,443,488,531,547,667] 
;     rrs_bins = [284234,284237,284238]
;     ;read QAA IOP data (qaa_matrix (with columns: aph,adg,bbp), qaa_bins) from file
;     qaa_matrix = [[0.0189454,0.00553217,0.0133541],$
;                  [0.0192148,0.00571175,0.0138207],$
;                  [0.0186105,0.00577449,0.0112262]]
;     qaa_bins = [284237,284238,284239]
;     qaa_min = 0.  
;     qaa_max = 5.
;     debug = 1
;     bandshiftCorrection, rrs_matrix, rrs_wavelengths, rrs_bins, qaa_matrix, correction_context, debug=debug, $
;                          qaa_bins=qaa_bins, qaa_min=qaa_min, qaa_max=qaa_max, correction_factors=correction_factors
;	  
;     ;rrs_matrix concatenates input and output values; rrs_wavelengths lists the final wavelength set. 
;       
;   
; :Pre:
;
;   `rrs_wavelengths` must contain at least 443nm and the green band (see `sensor`)
;   
;   The columns of `rrs_wavelengths` and `rrs_matrix` must correspond.
;   
;   The rows of `rrs_bins` and `rrs_matrix` must correspond, as well as the rows of `qaa_bins` and `qaa_matrix`. 
;   
;   The order of IOPs in `qaa_matrix` must be respected: first column: phytoplankton absorption, second column: detritus-gelbstoff 
;   absorption, last column: particle backscattering
;   
;   The IOPs contained in `qaa_matrix` must be of the same wavelength as the spectral model start bands (spec_model_start field) of `conversion_context`
;   
;   The resolution of the (QAA(v5) derived) IOP and RRS bin data are supposed to be equal, any resolution conversion need to be performed beforehand
;   
;   `correction_context` is created according to `bandShiftCorrectionSetup`
;   
;   Each input wavelength (lambda_i) in `correction_context` must be represented once and only once in `rrs_wavelengths`
;
;   `correction_context` may not contain trivial corrections (output wavelength of lambda_o equals input wavelength of lambda_i) 
;   
;   `correction_context` must at least contain one correction (one elements in lambda_i and lambda_o)

;   
; :Requires:
;   IDL 8.1
; 
; :Uses:
;   BandShiftCorrectionCore, weightedAverageEqualCorrectionProducts
;   
; :Params:
;   rrs_matrix: in, out, required, type="fltarr(m)/fltarr(m,n)"
;   
;     In input: remote sensing reflectances that -all or part of them- need to be bandshift corrected, in output: input + bandshift 
;     corrected remote sensing reflectance or reflectances
;          
;     Columns and column dimension m correspond to wavelengths in `rrs_wavelengths`, row(s) and row dimension n 
;     correspond to bin(s) in `rrs_bins`.
;     
;     During procedure execution, dimension m changes: new columns get appended containing bandshift corrected RRS at wavelengths corresponding 
;     to output wavelengths (lambda_o field) of `correction_context` (maintaining the order of the wavelengths in lambda_o, but retaining only 
;     the first appearance of wavelengths multiply present in lambda_o).
;     
;     During procedure execution, also dimension n may change: the rows corresponding to bins that do not have a 
;     valid counterpart in `qaa_bins` are removed. 
;     
;     In case `correction_context` contains different input wavelengths (lambda_i field) corrected to the same output wavelength 
;     (lambda_o field), at first all the individual bandshift corrections having that same output wavelength are carried out, 
;     and then a weighted average is calculated of each of the individually calculated bandshift corrected RRS values. The contribution 
;     of each individual bandshift corrected RRS to the weighted average is inverse proportional to the corresponding distance between 
;     the input and output wavelength of the correction. Only the weighted averaged RRS is retained in the output 
;     (see weightedAverageEqualCorrectionProducts)
;         
;   rrs_wavelengths: in, out, required, type = fltarr(m)
;   
;     Wavelengths corresponding to the columns of `rrs_matrix` (before as well as after procedure execution)
;     
;     m equals m of `rrs_matrix` dimension
;     
;     RRS corresponding to wavelengths appearing in lambda_i are corrected to the matching wavelength in the lambda_o field. RRS corresponding 
;     to wavelengths not appearing in lambda_i are ignored.   
;      
;     During procedure execution, dimension m changes: the unique output wavelengths (lambda_o field) of `correction_context` get appended 
;     (maintaining the order of the wavelengths in lambda_o, but retaining only the first appearance of doubles).
;     
;     In case `correction_context` contains different input wavelengths (lambda_i field) corrected to the same output wavelength (lambda_o field), 
;     the output wavelength multiply present, gets appended once, and the according column position in `rrs_matrix`, contains a weighted averaged 
;     bandshift corrected RRS value (see weightedAverageEqualCorrectionProducts)
;     
;   rrs_bins: in, out, required, type="long/lonarr(n)"
;   
;     Bin numbers (westernmost start bin south pole = 1) that correspond to rows of `rrs_matrix`
;
;     n equals n of `rrs_matrix` dimension, type can be a scalar if also rrs_matrix is 1-dimensional (only 1 bin), during processing however the 
;     scalar will be  converted in lonarr(1)
;         
;     During procedure execution, dimension n may change: in case `qaa_bins` is not equal to rrs_bins, intersection between both bin lists is taken 
;     and rrs_bins is set to this intersection
;     
;   qaa_matrix: in, required, type="fltarr(3)/fltarr(3,o)/fltarr(3,o,p)"
;   
;     IOPs, preferably as calculated by QAA(v5) (See Z. Lee, B. Lubac, J. Werdell, R. Arnone: An update of the quasi-analytical algorithm (QAA_v5): International 
;     Ocean Color Group software report)
;     
;     If `qaa_bins` contains more than 1 bin, there has to be a second dimension o, equalling this number of bins. 
;
;     The third dimension p, if present, should equal the number of spectral model start bands (spec_model_start field) of `correction_context`. If the 
;     number of spectral model start bands of `correction_context` is 1, the third dimension is absent.  
;
;     During procedure execution qaa_matrix is unaltered. 
;
;     The order of IOPs in `qaa_matrix` must be respected: first column: phytoplankton absorption, second column: detritus-gelbstoff 
;     absorption, last column: particle backscattering
;   
;   correction_context: in, required, type=structure
;     
;     Structure containing input and output wavelengths (lambda_i and lambda_o) of desired bandshift corrections, 
;     Bricaud coefficients (a_i, b_i, a_o, b_o) for prediction of phytoplankton absorption at the input (_i) and output (_o) wavelengths, 
;     pure water absorption (aw_i, aw_o) and backscattering (bbw_i, bbw_o) at those input (_i) and output (_o) wavelengths,
;     the start band(s) of the spectral model (spec_model_start), one value applicable to all bandshift corrections or a value for each bandshift 
;     correction, and the Bricaud coefficients (a_sms, b_sms) for these start band(s).
;          
;     correction_context must have the exact same format as the structure returned by `bandShiftCorrectionSetup`
;     
; :Keywords:
;   debug: in, optional, type=int, default=1
;     
;     0: no text output is produced, 1: some text output is produced
;
;   qaa_bins: in, optional, type="long/lonarr(o)", default=:rrs_bins:
;   
;     Bin numbers (westernmost start bin south pole = 1) that correspond to rows of the `qaa_matrix`
;     
;     o equals o of `qaa_matrix` dimension, type can be a scalar if there is only 1 IOP bin, during processing however the 
;     scalar will be converted in lonarr(1)
;     
;     Apart from this conversion, qaa_bins remains unaltered during procedure execution
;     
;   qaa_min: in, optional, type=float, default=0.
;
;     Bins with any (QAA(v5) derived) IOP equal or below qaa_min are discarded 
;   
;   qaa_max: in, optional, type=float, default=100.
;   
;     Bins with any (QAA(v5) derived) IOP equal or above qaa_max are discarded
;   
;   correction_factors: out, type="fltarr(p)/fltarr(p,q)"
;
;     Calculated correction factors for each of the corrections in `correction_context` 
;
;     Dimension p corresponds to the number of corrections in `correction_context` (i.e. number of input wavelengths, lambda_i, 
;     and output wavelengths, lambda_o), and also equals the dimension p of `qaa_matrix` if this dimension is present. 
;     
;     Dimension q corresponds to the number of bins in the intersection of `rrs_bins` and `qaa_bins` that had valid IOP values.
;     
;     Can be 1-dimensional (only p) if only 1 valid bin was found in the intersection of `rrs_bins` and `qaa_bins`. 
;     
;     Factors can also be calculated by dividing output RRS by input RRS, for those output RRS that were not obtained through a 
;     weighted average (see `rrs_matrix` for an explanation on the average weighting mechanism).
;                               
;-
PRO BandShiftCorrection, rrs_matrix, rrs_wavelengths, rrs_bins, qaa_matrix, correction_context, debug=debug, qaa_bins=qaa_bins, qaa_min=qaa_min, qaa_max=qaa_max, correction_factors=correction_factors
  
  ON_ERROR, 1
  nb_corrected = bandShiftCorrectionCore(rrs_matrix, rrs_wavelengths, rrs_bins, qaa_matrix, correction_context, $
    rrs_matrix_corrected, rrs_wavelengths_corrected, $
    qaa_bins = qaa_bins, qaa_min = qaa_min, qaa_max = qaa_max, correction_factors=correction_factors,$
    rrs_valid_index=rrs_valid_index) 
   
  IF KEYWORD_SET(debug) THEN print, 'Corrected ' + STRCOMPRESS(nb_corrected,/REMOVE_ALL) + ' bins.'
  IF nb_corrected LE 0 THEN RETURN  
  
  weightedAverageEqualCorrectionProducts, rrs_matrix_corrected, rrs_wavelengths_corrected, correction_context
  rrs_matrix = [rrs_matrix[*,rrs_valid_index],rrs_matrix_corrected] ; output rrs
  rrs_wavelengths = [rrs_wavelengths,rrs_wavelengths_corrected] ; concatenate wavelengths
  rrs_bins = rrs_bins[rrs_valid_index] ; output bin numbers
  
END




    