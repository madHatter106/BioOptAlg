import numpy as np

def readSpectralDataPureWater(idir,seadas_mimick = False):
    pass


def bandShiftCorrectionSetup(pw_dir,sensor='MODISA'):
    
    if sensor.upper() == "SEAWIFS":
        lambda_i = array([412.,490.,510.,555.,555.,555.,670.,670.])
        lambda_o = array([413.,488.,531.,531.,547.,560.,665.,667.])
    elif sensor.upper() == "MERIS":
        lambda_i = array([413.,490.,510.,560.,560.,560.,665.,665.])
        lambda_o = array([412.,488.,531.,531.,547.,555.,667.,670.])
    elif sensor.upper() == "MODISA":
        lambda_i = array([412.,488.,488.,531.,547.,547.,667.,667.])
        lambda_o = array([413.,510.,490.,510.,560.,555.,665.,670.])
    
    
    # read IOPs from pure seawater from directory pw_dir
    sdata = readSpectralDataPureWater(pw_dir)
    
    subset = np.arange(lambda_i.size)
    
    aw_i = np.empty_like(lambda_i)
    bbw_i = np.empty_like(lambda_i)
    a_i = np.empty_like(lambda_i)
    b_i = np.empty_like(lambda_i)
    aw_o = np.empty_like(lambda_o)
    bbw_o = np.empty_like(lambda_o)
    a_o = np.empty_like(lambda_o)
    b_o = np.empty_like(lambda_o)
    
    i_cnt = 0
    for wv_i in lambda_i:
        [status,aw,bw] = getSpectralDataPureWater(sdata,wv_i,bandwith=0)
        if status == 0:
            aw_i[i_cnt] = aw
            bbw_i[i_cnt] = bw/2.
            [a_bricaud,b_bricaud] = get_A_B_Bricaud(lambda_i[i_cnt])
            
    