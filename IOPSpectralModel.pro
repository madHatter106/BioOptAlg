; docformat:'rst'
 
;+ 
; Procedure and related utility to predict IOPs at given wavelengths starting from known values at other wavelengths, by applying a spectral model.
; The IOPs taken into account are pythoplankton absorption aph, detritus+gelbstoff absorption adg and particle backscattering bbp. 
;
; :version: 3
; 
; :Author:
;
;   Frédéric Mélin (European Commission/JRC/IES/WRES)
;   
;   Gert Sclep (European Commission/JRC/IES/WRES)
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
; :version: 2
; 
; :Categories:
; 
;   IOP, aph, spectral model
;   
; :History:
;   
;   Created by Gert Sclep, 07/2012 
;   
; :Description:
;   
;   Retrieves the A and B Bricaud coefficient for given wavelengths. 
;   
;   See table 2 from Bricaud et al. Variability in the chlorophyll-specific absorption coefficients of natural phytoplankton: Analysis and parameterization. Journal of Geophysical Research. 1995;100(95):321–332.  
;  
; :Pre:
;   Code currently contains only Bricaud coefficients for a selection of wavelengths (412,413,443,488,490,510,531,547,555,560,665,667,670). 
;   Only these wavelengths are allowed as input.   
;
; :Requires:
;   IDL 8.0
;
; :Params:
; 
;    wl: in, required, type="int/float/intarr(n)/fltarr(n)"
;     Wavelength(s) for which the Bricaud A and B coefficients need to be fetched.     
;     
;    a_bricaud: out, required, type="fltarr(n)"
;     Bricaud A coefficients for the wavelengths corresponding to `wl`.
;     
;     If `wl` is a scalar, a_bricaud will have type fltarr(1)     
;    
;    b_bricaud: out, required, type="fltarr(n)"
;     Bricaud B coefficients for the wavelengths corresponding to `wl`.
;     
;     If `wl` is a scalar, b_bricaud will have type fltarr(1)     
;    
;-
PRO get_A_B_Bricaud, wl, a_bricaud, b_bricaud
  ; a and b at 443 derived from values previously in a_bb_prediction
  ; a_443 = 0.0394 and k (= 1/(1-b))= 1.52323  
  a_hash = hash(412,0.0323,443, 0.0394,488,0.0279,510,0.0180,531,0.0115,547,0.00845,667,0.01685,665,0.0152, $
                413,0.032775,490,0.0274, 560, 0.0062, 555, 0.007, 670, 0.0189)
  b_hash = hash(412,0.286,443, 0.3435, 488,0.369,510,0.260,531,0.134,547,0.0625,667,0.140,665,0.134, $
                413,0.28775,490,0.361, 560, 0.016, 555, 0.0315, 670, 0.149)
  n_wl = n_elements(wl)
  a_bricaud = fltarr(n_wl)
  b_bricaud = fltarr(n_wl)
  for i=0,n_wl-1 do begin                         
    a_bricaud[i] = a_hash[wl[i]]
    b_bricaud[i] = b_hash[wl[i]]
  endfor
END

;+
; :Author:
; 
;   Frédéric Mélin 
;   
;   Gert Sclep
;
; :version: 2
;    
; :History:
; 
;   Created by Frédéric Mélin, 01/2011
;   
;   Modified by Gert Sclep, 04/2012::
;     
;     Made procedure 2D matrix-oriented, allowing to process multiple records at once. 
;   
; :Categories:
; 
;   IOP, aph, adg, bbp, spectral model
;   
; :Description:
; 
;   Using a spectral model, IOPs at a given wavelength are converted to a set of desired wavelengths. 
;   This is done for phytoplankton absorption, for detritus-gelbstoff absorption and for particle backscattering.
;   The phytoplankton absorption is converted using the Bricaud formula (see Bricaud et al. Variability in the chlorophyll-specific absorption coefficients 
;   of natural phytoplankton: Analysis and parameterization. Journal of Geophysical Research. 1995;100(95):321–332), and this requires the Bricaud
;   coefficients A and B to be given for both the input as for the output wavelengths. The detritus-gelbstoff and particle backscattering
;   are spectrally evolved using the same spectral slope as used in QAAv5 (See Z. Lee, B. Lubac, J. Werdell, R. Arnone: An update of the quasi-analytical algorithm 
;   (QAA_v5): International Ocean Color Group software report). To calculate these slopes, the below water remote sensing reflectances in the blue and green band
;   need to be given as input. 
;     
; :Requires:
; 
;   IDL 7.1
;   
; :Params:
;   lambda_in: in, required, type="float/fltarr(m)"   
;     start wavelength(s) of the spectral model, in general the start wavelength used is the same for each of the individual conversions
;            
;   a_in: in, required, type="float/fltarr(m)"     
;     Bricaud A coefficient for the wavelengths given in `lambda_in`
;      
;   b_in: in, required, type="float/fltarr(m)"
;     Bricaud B coefficient for the wavelengths given in `lambda_in`     
;           
;   aph_in: in, required, type="float/fltarr(m)/fltarr(m,n)" 
;    phytoplankton absorption for the wavelengths given in `lambda_in` (column dimension m), for 1 or more records (if >1: row dimension n)
;    
;   adg_in: in, required, type="float/fltarr(m)/fltarr(m,n)"
;    CDOM+detritus absorption for the wavelengths given in `lambda_in` (column dimension m), for 1 or more records (if >1: row dimension n)
;     
;   bbp_in: in, required, type="float/fltarr(m)/fltarr(m,n)"
;     particulate back-scattering for the wavelengths given in `lambda_in` (column dimension m), for 1 or more records (if >1: row dimension n)
;   
;   rrs_blue_in: in, required, type="float/fltarr(n)"
;     below water remote sensing reflection in the blue band (443 nm), for 1 or more records (if >1: column dimension n) 
;     
;   rrs_green_in: in, required, type="float/fltarr(n)"
;     below water remote sensing reflection in the green band (547 or 555 or 560 nm, depending on sensor used), for 1 or more records (if >1: column dimension n) 
;     
;   lambda_out: in, required, type="float/fltarr(m)"
;    end wavelengths to which `aph_in`, `adg_in` and `bbp_in` should be evolved towards
;    
;    the dimension m is equal to the dimension m of `lambda_in`
;     
;   a_out: in, required, type="float/fltarr(m)"
;    Bricaud A coefficient for the wavelengths given in `lambda_out`        
;    
;   b_out: in, required, type="float/fltarr(m)"
;    Bricaud B coefficient for the wavelengths given in `lambda_out`
;           
;   aph_out: out, type="float/fltarr(m)/fltarr(m,n)"
;    phytoplankton absorption for the wavelengths given in `lambda_out` (column dimension m), for 1 or more records (if >1: row dimension n), obtained using Bricaud formula starting from values in `aph_in`
;   
;   adg_out: out, type="float/fltarr(m)/fltarr(m,n)"
;    CDOM+detritus absorption for the wavelengths given in `lambda_out` (column dimension m), for 1 or more records (if >1: row dimension n), obtained applying QAAv5 spectral slopes starting from values in `adg_in`.
;     
;   bbp_out: out, type="float/fltarr(m)/fltarr(m,n)"
;     particulate back-scattering for the wavelengths given in `lambda_out` (column dimension m), for 1 or more records (if >1: row dimension n), obtained applying QAAv5 spectral slopes starting from value in `bbp_in`
;
;-			   
PRO IOPSpectralModel, lambda_in, a_in, b_in, aph_in, adg_in, bbp_in, rrs_blue_in, rrs_green_in, $
lambda_out, a_out, b_out, aph_out, adg_out, bbp_out
	
	  dimension_input = size(aph_in,/dimension)
	  nb_conversions = dimension_input[0]
	  IF (N_ELEMENTS(dimension_input) EQ 1) THEN nb_predictions = 1 ELSE nb_predictions = dimension_input[1]
	
    k = 1.52323 ; = 1/(1-b_443)
    a_443 = 0.0394 ; 
    
    S = 0.015
    rat = rrs_blue_in/rrs_green_in 
    sdg = S + 0.002/(0.6+rat)
    
    yy = 2.0 * (1.0 - 1.2 * exp(-0.9*rat))
    
    Chla = (aph_in/REBIN(a_in,nb_conversions,nb_predictions))^(1/(1-REBIN(b_in,nb_conversions,nb_predictions))) ; Chla concentration; Bricaud et al. (1995)
    
    ; IOPs for input wavelengths
    ll = lambda_out - lambda_in
    aph_out = REBIN(a_out,nb_conversions,nb_predictions)*Chla^(1-REBIN(b_out,nb_conversions,nb_predictions)) ; b_out: bricaud coefficient
    adg_out = adg_in * EXP(-TRANSPOSE(REBIN(sdg,nb_predictions,nb_conversions))*REBIN(ll,nb_conversions,nb_predictions))    
    bbp_out = bbp_in * REBIN((lambda_in/lambda_out),nb_conversions,nb_predictions)^TRANSPOSE(REBIN(yy,nb_predictions,nb_conversions))
 
END
