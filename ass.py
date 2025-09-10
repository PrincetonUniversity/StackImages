global version
version = 21
'''
Here are the functions defined herein...
    resize(image, m, n):
    setxybar(*arg):
    color_weights(R, G, B, border=70, help=False):
    color_weights_old(R, G, B, border=20, help=False):
    _zvalue_from_index(arr, ind):
    nan_percentile(arr, q):
    stack(image_stack, outlier_thresh=1.5, help=False):
    createMatcher(method,crossCheck):
    matchKeyPointsBF(featuresA, featuresB, method):
    matchKeyPointsKNN(featuresA, featuresB, ratio, method):
    getHomography(kpsA, kpsB, featuresA, featuresB, matches, reprojThresh):
    detectAndDescribe(image, method=None):
    align_cv(images):
    autoalign3stars(images, thresh=0.1, border=10, max_rotate=30, max_rescale=1):
    autoalign(images, numstars=1, thresh=0.1, border=0, \
    onclick1(event):
    stack_images():
    show_stacked_image():
    ChangeCursor(event):
    UpdateStatusBar(event):
    align(image_stack, method):
    align1star(image_stack):
    onclick(event):
    get_aligned_images():
    get_fwhms():
    onmouse(event):
    slidershow(Img, fignum=1, img_out=[], folder='', fname=''):
    HRdiag_show(fignum=1, img_out=[], folder='', fname=''):
    update(val):
    removeHotPixels(L, thresh=2):
    findstars(L, thresh=0.0006, fignum=2, bkg=False):
    findstars_new(*arg):
    findstars_small_psf(*arg):
    setArcsecPerPixel(app):
    setMaxBoxSize(mbs):
    showAlignStar(TrueFalse, fnames=''):
    setBoxsize(boxsize):
    setExptime(exp_time):
    makeHRdiag(*arg):
    makeHRdiag(RGB, objname='', minmag0=10, shape='rect', shape_params=20, \
    weighted_median(values, weights):
    weighted_percentile(values, weights, p):
    makeHRdiag_new(*arg):
    makeHRdiag_new(RGB, objname='', minmag0=10, shape='rect', shape_params=20, \
    makeHRdiagRRlyrae(*arg):
    makeHRdiagRRlyrae3(RGB, RGB2, RGB3, objname, minmag0, fluxthresh, colorthresh, showit, z):
    makeHRdiagRRlyraeN(RGBs, objname, minmag0, fluxthresh, colorthresh, showit, z):
    gnomonic(ra,dec,ra_mid,dec_mid):
    inv_gnomonic(y,z,ra_mid,dec_mid):
    equatorial_to_galactic(ra,dec):
    fwhm(L,i,j):
    ErosionKF(Image, numpx):
    ErosionKF2(Image, numpx):
    GaussianBlur(Image, r):
    NeatImageCore(Image, StdDev, Blur):
    NeatImage(Image, StdDev, Blur):
    Denoise(im, thresh):
    NeatImageFFTCore(Image, StdDev, Blur):
    NeatImageFFT(Image, StdDev, Blur):
    NeatImageFFT(Image, StdDev, Blur):
    ImageSmooth(*arg):
    ImageSmoothCore(I,thresh_pctl,r,medmean):
    applyFlat(I, flat, bias=0):
    RemoveGradient(I):
    HDR(RGB):
    UnsharpMaskRGB(RGB):
    UnsharpMask(L):
    RLdeconvRGB(RGB, fwhm_target, thresh=0.1, x=0, y=0, std_dev=5, \
    RLdeconvGrid(L, fwhm_target, grid_size=1, thresh=0.1, x=0, y=0, std_dev=5, \
    find_fwhm_star(L, rad=6, thresh=0.1):
    pad(img, rad):
    RLdeconv(L, fwhm_target, psf, thresh=0.1, x=0, y=0, std_dev=5, \
    decimal_to_min_sec(angle):
    min_sec_to_decimal(s, hr, min, sec):
    fits_astrometry(file_location, approx_ra_dec=[0,0], approx_fov=0.5, \
    astrometry(image, approx_ra_dec=[0,0], approx_fov=0.5, \
    radec_to_xyz(radec):
    SiderealTime(month,day,hour,minute,second,longitude):
    getRaDec(obj):
    askopenfilenames(multiple=False, filetypes='fits files (*.fit)|*.fit|(*.fits)|*.fits|(*.fts)|*.fts'):
    get_RGB_image(evt):
    show_fwhm(L):
    get_images(evt):
    merge_paths(p1,p2):
    clear_bias():
    clear_flat():
    get_bias(evt):
    get_flat(evt):
    align_stack_show(evt):
    find_image(evt):
    show_gaia_image(evt):
    show_image(evt):
    min_sec_to_decimal(min_sec_str):
    astrometric_analysis(evt):
    quit(evt):
    deconvolve(evt):
    denoise(evt):
    hdr(evt):
    unsharpmask(evt):
    GUI_HRdiag(evt):
    GUI_makeHRdiag0(evt):
    GUI_makeColorDiag0(evt):
    got_updated():
    set_align_method(event):
    set_showalign(event):
    set_removehotpix(event):
    set_neatimage(event):
    set_flip(event):
    set_showimage(event):
    set_stack_method(event):
    StackImages():
    ShowClusters():    
'''
import math
import numpy as np
from numpy import *
import numpy.ma as ma
from scipy.signal import *
from scipy.stats import *
import matplotlib
from matplotlib.pyplot import *
from matplotlib.colors import LogNorm
#import pyfits   # Tools -> Package Manager -> Available -> pyfits
from astropy.utils.data import download_file
from astropy.io import fits
from matplotlib.widgets import Slider, Button, TextBox
import matplotlib.widgets as wid
#import imageio
import cv2
#from tictoc import tic, toc
#from Tkinter import *
#import png
import csv
from scipy import fftpack
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

import platform

import wx
import wx.lib.agw.aui as aui
import wx.lib.mixins.inspection as wit

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import (
    FigureCanvasWxAgg as FigureCanvas,
    NavigationToolbar2WxAgg as NavigationToolbar)
matplotlib.use('Qt5Agg')


import PIL

global showAlign, FileNames
global ybar0, xbar0
global xmin, xmax, ymin, ymax
global lbl4
global file_location
global tab
global L_idx
global current_fwhm
global nstacked
global plotter, images_frame

plotter = ''
images_frame = ''

current_fwhm = 0
nstacked = 0

tab = ''

global success
success = []

objname = ''

global image_type
image_type = '2D'

global GUI, slash
GUI = False
slash = '/'
print('OS = '+platform.system())
if (platform.system() == 'Windows'): slash = '\\'
else:                                slash = '/'



xmin=0
xmax=0
ymin=0
ymax=0

exptime = 1.
box_size = 80
showAlign = False
FileNames = ''
file_location = ''
MaxBoxSize = 30
ArcsecPerPixel = 1

def resize(image, m, n):
    if (len(image.shape) == 2):
        (m0,n0) = image.shape
        upimage = zeros((m,n), dtype=float32)
    
        xs = (np.array( arange(0,n) )   * n0/n).round()
        ys = (np.array([arange(0,m)]).T * m0/m).round()
        xs = xs.astype(int)
        ys = ys.astype(int)
        upimage[:,:] = image[ys,xs]
        return upimage
    elif (len(image.shape) == 3):
        (m0,n0,N) = image.shape
        upimage = zeros((m,n,N), dtype=float32)
    
        xs = (np.array( arange(0,n) )   * n0/n).round()
        ys = (np.array([arange(0,m)]).T * m0/m).round()
        xs = xs.astype(int)
        ys = ys.astype(int)
        upimage[:,:,:] = image[ys,xs,:]
        return upimage

def setxybar(*arg):
    global ybar0, xbar0
    xbar0 = arg[0]
    ybar0 = arg[1]

######################################################################################
## Code for opening multiple figures as tabs in one window

class Plot(wx.Panel):
    def __init__(self, parent, id=-1, dpi=None, **kwargs):
        super().__init__(parent, id=id, **kwargs)
        self.figure = Figure(dpi=100, figsize=(5, 5))
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)


class PlotNotebook(wx.Panel):
    def __init__(self, parent, id=-1):
        super().__init__(parent, id=id)
        self.nb = aui.AuiNotebook(self)
        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        self.SetSizer(sizer)

    def add(self, name="plot"):
        page = Plot(self.nb)
        self.nb.AddPage(page, name, select=True)
        return page.figure

######################################################################################
## Compute color weights using aligned images

def color_weights(R, G, B, border=70, help=False):

    '''
    r_flux = percentile(R,99.9)-percentile(R,99.0)
    g_flux = percentile(G,99.9)-percentile(G,99.0)
    b_flux = percentile(B,99.9)-percentile(B,99.0)
    '''

    k, LconvL2 = findstars(R,thresh, fignum=0)
    k = k&(R < 0.8*np.max(R));
    flux = LconvL2[k]
    I = argsort(-flux)
    Rflux = flux[I]

    k, LconvL2 = findstars(G,thresh, fignum=0)
    k = k&(R < 0.8*np.max(R));
    flux = LconvL2[k]
    I = argsort(-flux)
    Gflux = flux[I]

    k, LconvL2 = findstars(B,thresh, fignum=0)
    k = k&(R < 0.8*np.max(R));
    flux = LconvL2[k]
    I = argsort(-flux)
    Bflux = flux[I]

    nfluxes = np.min((len(Rflux), len(Gflux), len(Bflux), 20))

    r_flux = sum(Rflux[0:nfluxes])
    g_flux = sum(Gflux[0:nfluxes])
    b_flux = sum(Bflux[0:nfluxes])





    if ( (r_flux > 0) & (g_flux > 0) & (b_flux > 0) ): 
        numerator = (r_flux * g_flux * b_flux)**(1/3)
        r_wt = numerator/r_flux
        g_wt = numerator/g_flux
        b_wt = numerator/b_flux
    else:
        r_wt = 1;
        g_wt = 1;
        b_wt = 1;

    return r_wt, g_wt, b_wt

def color_weights_old(R, G, B, border=20, help=False):

    global xmin, xmax, ymin, ymax

    if (help != False):
        print('\nFunction color_weights parameters: \n', \
                '  R, G, B:  the images for red, green, and blue (Required) \n', \
                'Output:  weight factors for the red, green and blue frames \n', \
                '\n')

    eps = 1e-5
    
    border = np.max([-xmin, xmax, -ymin, ymax])+20


    (m,n) = G.shape

    xs = arange(0,n)
    ys = array([arange(0,m)]).T
    Xs = ones((m,1))*xs
    Ys = ys*ones((1,n))

    thresh = 0.1
    k, GconvL2 = findstars(G,thresh=thresh,fignum=0,bkg=True)
    mask = (Xs > border) & (Xs < n-border) & (Ys > border) & (Ys < m-border) & (G < 0.8*np.max(G))
    k = k&mask;
    Gflux = GconvL2[k]
    '''
    figure(100)
    imshow(G, cmap='gray', norm=LogNorm())
    return 1,1,1
    '''

    while (len(Gflux[Gflux>0]) < 10):
        thresh = thresh/2
        k, GconvL2 = findstars(G,thresh=thresh,fignum=0,bkg=True)
        mask = (G < 0.8*np.max(G)) 
        k = k&mask;
        Gflux = GconvL2[k]

    I = argsort(-Gflux)[0:10]
    Gflux = Gflux[I] + eps

    Laplacian2 = array([[ 0., -1., -1., -1.,  0.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [-1.,  2.,  4.,  2., -1.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [ 0., -1., -1., -1.,  0.]])
    RconvL2 = convolve2d(R, Laplacian2, mode='same')
    BconvL2 = convolve2d(B, Laplacian2, mode='same')

    Rflux = RconvL2[k]
    Bflux = BconvL2[k]

    Rflux = Rflux[I] + eps
    Bflux = Bflux[I] + eps

    kk = (Rflux>0)&(Gflux>0)&(Bflux>0)
    Rflux = Rflux[kk]
    Gflux = Gflux[kk]
    Bflux = Bflux[kk]
    '''
    '''

    numerator = (Rflux*Gflux*Bflux)**(1/3)
    R_wts = numerator/Rflux
    G_wts = numerator/Gflux
    B_wts = numerator/Bflux

    r_wt = np.mean(R_wts)
    g_wt = np.mean(G_wts)
    b_wt = np.mean(B_wts)

    return r_wt, g_wt, b_wt

######################################################################################
## Stack images


# from here:  https://krstn.eu/np.nanpercentile()-there-has-to-be-a-faster-way/
def _zvalue_from_index(arr, ind):
    """private helper function to work around the limitation of np.choose() by employing np.take()
    arr has to be a 3D array
    ind has to be a 2D array containing values for z-indicies to take from arr
    See: http://stackoverflow.com/a/32091712/4169585
    This is faster and more memory efficient than using the ogrid based solution with fancy indexing.
    """
    # get number of columns and rows
    _,nC,nR = arr.shape

    # get linear indices and extract elements with np.take()
#   idx = nC*nR*ind + nC*np.arange(nR)[:,None] + np.arange(nC)
    idx = nC*nR*ind + np.arange(nC*nR).reshape((nC,nR))
    return np.take(arr, idx)

# from here:  https://krstn.eu/np.nanpercentile()-there-has-to-be-a-faster-way/
def nan_percentile(arr, q):
    # valid (non NaN) observations along the first axis
    valid_obs = np.sum(np.isfinite(arr), axis=0)
    # replace NaN with maximum
    max_val = 2*np.nanmax(arr)
    arr_sort = np.zeros(arr.shape)
    arr_sort[:,:] = arr
    arr_sort[np.isnan(arr)] = max_val
    # sort - former NaNs will move to the end
    arr_sort = np.sort(arr_sort, axis=0)

    # loop over requested quantiles
    if type(q) is list:
        qs = []
        qs.extend(q)
    else:
        qs = [q]
    if len(qs) < 2:
        quant_arr = np.zeros(shape=(arr.shape[1], arr.shape[2]))
    else:
        quant_arr = np.zeros(shape=(len(qs), arr.shape[1], arr.shape[2]))

    result = []
    for i in range(len(qs)):
        quant = qs[i]
        # desired position as well as floor and ceiling of it
        k_arr = (valid_obs - 1) * (quant / 100.0)
        f_arr = np.floor(k_arr).astype(np.int32)
        c_arr = np.ceil(k_arr).astype(np.int32)
        fc_equal_k_mask = f_arr == c_arr

        # linear interpolation (like numpy percentile) takes the fractional part of desired position
        floor_val = _zvalue_from_index(arr=arr_sort, ind=f_arr) * (c_arr - k_arr)
        ceil_val = _zvalue_from_index(arr=arr_sort, ind=c_arr) * (k_arr - f_arr)


        quant_arr = floor_val + ceil_val
        quant_arr[fc_equal_k_mask] = _zvalue_from_index(arr=arr_sort, ind=k_arr.astype(np.int32))[fc_equal_k_mask]  # if floor == ceiling take floor value

        result.append(quant_arr)

    return array(result)

def stack(image_stack, outlier_thresh=1.5, help=False):
    global ybar0, xbar0
    global gauge, ntot, nstacked


    if (help != False):
        print('\nFunction stack parameters: \n', \
                '  images:  the set of images to be stacked into a single image (Required) \n', \
                '  outlier_thresh:  ignore pixel values that are further \n'+ \
                '     from the median than this factor times the interquartile range \n'+ \
                '     Default is 1.5 \n', \
                '\n')

    A = image_stack
    A = 1.00000001*A + 0.001
    (N, m, n) = array(A).shape

    if N==0: return array(zeros((m,n)))

    # sadly, nanpercentile is pathetically slow.   nan_percentile is fast!!!
    if (stack_method == 'Average'):
        #####################################################################################
        # Stack images by averaging
        B = (A==0.001)
        A[B] = float('nan')
        A_median = nan_percentile(A, [50])
        A_iqr    = nan_percentile(A, [100]) - nan_percentile(A, [0])
        A[B] = 0
    
        A_median = tile(A_median, (N,1,1))
        A_iqr    = tile(A_iqr,    (N,1,1))
        A = ma.masked_array(A, mask=( (abs(A-A_median)>outlier_thresh*A_iqr) | \
                                      (A<=0) ) )
    elif (stack_method == 'Max'):
        #####################################################################################
        # Stack images by averaging
        B = (A==0.001)
        A[B] = float('nan')
        A_max = nan_percentile(A, [100])
        A[B] = 0
    
        A_max = tile(A_max, (N,1,1))
        A = ma.masked_array(A, mask=( ((0.9*A_max-A)>0) ) )
    else:
        #####################################################################################
        # Stack images by outlier detection (smarter and faster than SD-mask)
        B = (A<=0.001)
        A[B] = float('nan')
        min_bkg = 1e+10
        bkg = zeros((N))
        for j in arange(N):
            bkg[j] = nanpercentile(A[j,:,:], [1])
            if (bkg[j] < min_bkg): min_bkg = bkg[j]
        for j in arange(N):
            A[j,:,:] = A[j,:,:] - bkg[j] + min_bkg
        A_median = nan_percentile(A, [50])
        if (gauge != 0): gauge.SetValue(int(50+50*(nstacked+N/2)/ntot)); wx.Yield()
        A_iqr    = nan_percentile(A, [75]) - nan_percentile(A, [25])
        A[B] = 0

        A_median = tile(A_median, (N,1,1))
        A_iqr    = tile(A_iqr,    (N,1,1))
        mask=( (abs(A-A_median)>outlier_thresh*A_iqr) | (A<=0) ) 
        '''
        print('outlier_thresh = ', outlier_thresh)
        print('hot spot = ', A[:,91,551])
        print('median = ', A_median[0,91,551])
        print('iqr    = ', A_iqr[0,91,551])
        print('min of A = ', np.min(A))
        '''
    
        A = ma.masked_array(A, mask=( (abs(A-A_median)>outlier_thresh*A_iqr) | \
                                      (A<=0) ) )
    A = ma.mean(A, axis=0)

    if (gauge != 0): gauge.SetValue(int(50+50*(nstacked+N)/ntot)); wx.Yield()
    nstacked += N
    
    return A.filled(fill_value=0)

######################################################################################
## Align images on 1 or 2 stars

global fwhms, align_stack

def createMatcher(method,crossCheck):
    "Create and return a Matcher Object"
    
    if method == 'sift' or method == 'surf':
        bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=crossCheck)
    elif method == 'orb' or method == 'brisk':
        bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=crossCheck)
    return bf

def matchKeyPointsBF(featuresA, featuresB, method):
    bf = createMatcher(method, crossCheck=True)
        
    # Match descriptors.
    best_matches = bf.match(featuresA,featuresB)
    
    # Sort the features in order of distance.
    # The points with small distance (more similarity) are ordered first in the vector
    rawMatches = sorted(best_matches, key = lambda x:x.distance)
    print("Raw matches (Brute force):", len(rawMatches))
    return rawMatches

def matchKeyPointsKNN(featuresA, featuresB, ratio, method):
    bf = createMatcher(method, crossCheck=False)
    # compute the raw matches and initialize the list of actual matches
    rawMatches = bf.knnMatch(featuresA, featuresB, 2)
    print("Raw matches (knn):", len(rawMatches))
    matches = []

    # loop over the raw matches
    for m,n in rawMatches:
        # ensure the distance is within a certain ratio of each
        # other (i.e. Lowe's ratio test)
        if m.distance < n.distance * ratio:
            matches.append(m)
    return matches

def getHomography(kpsA, kpsB, featuresA, featuresB, matches, reprojThresh):
    # convert the keypoints to numpy arrays
    kpsA = np.float32([kp.pt for kp in kpsA])
    kpsB = np.float32([kp.pt for kp in kpsB])
    
    if len(matches) > 4:

        # construct the two sets of points
        ptsA = np.float32([kpsA[m.queryIdx] for m in matches])
        ptsB = np.float32([kpsB[m.trainIdx] for m in matches])
        
        # estimate the homography between the sets of points
        (H, status) = cv2.findHomography(ptsA, ptsB, cv2.RANSAC,
            reprojThresh)

        return (matches, H, status)
    else:
        return None

def detectAndDescribe(image, method=None):
    """
    Compute key points and feature descriptors using an specific method
    """
    
    assert method is not None, "You need to define a feature detection method. Values are: 'sift', 'surf'"
    
    # detect and extract features from the image
    if method == 'sift':
        descriptor = cv2.xfeatures2d.SIFT_create()
    elif method == 'surf':
        descriptor = cv2.xfeatures2d.SURF_create()
    elif method == 'brisk':
        descriptor = cv2.BRISK_create()
    elif method == 'orb':
        descriptor = cv2.ORB_create()
        
    # get keypoints and descriptors
    (kps, features) = descriptor.detectAndCompute(image, None)
    
    return (kps, features)

def build_aligned_stack(k,kk,images):
    global i0, i1, j0, j1, m1, n1, Xs, Ys, alphas, thetas, xbar0, ybar0, align_stack
    sin_theta = sin(thetas[k])
    cos_theta = cos(thetas[k])
    img = zeros((m1,n1), dtype=float32)
    img[i0:i1, j0:j1] = images[k,:,:]
    Xold = alphas[k]*(cos_theta*(Xs-xbar0[0])-sin_theta*(Ys-ybar0[0])) + xbar0[k]
    Yold = alphas[k]*(sin_theta*(Xs-xbar0[0])+cos_theta*(Ys-ybar0[0])) + ybar0[k]

    Yold[Yold<0] = 0
    Xold[Xold<0] = 0
    Yold[Yold>=m1-1] = m1-2
    Xold[Xold>=n1-1] = n1-2
    I  = int16(Ys)
    J  = int16(Xs)
    I0 = int16(Yold)
    J0 = int16(Xold)  

    I1 = I0+1
    J1 = J0+1
    px = Xold - J0
    py = Yold - I0
    qx = 1.-px
    qy = 1.-py
    align_stack[kk,I,J] = qy*qx*img[I0,J0] + \
                          qy*px*img[I0,J1] + \
                          py*qx*img[I1,J0] + \
                          py*px*img[I1,J1]

def align_astrometric(images):
    global xbar, ybar, xbar0, ybar0, xbar1, ybar1, stack_size, stack_num, align_stack
    global showAlign, FileNames
    global fwhms
    global progress, gauge
    global Mono_files, RGB_files, NonFit_files
    global success
    global L_idx
    global i0, i1, j0, j1, m1, n1, Xs, Ys, alphas, thetas

    global tab
    print(tab,'Aligning via Astrometry:')
    tab = tab+'    '

    (N,m,n) = images.shape
    ra_ctr  = zeros((N), dtype=float32)
    dec_ctr = zeros((N), dtype=float32)
    dtheta  = zeros((N), dtype=float32)
    scale   = zeros((N), dtype=float32)
    crpix1  = zeros((N), dtype=float32)
    crpix2  = zeros((N), dtype=float32)
    flip_it = zeros((N), dtype=bool)
    success = zeros((N), dtype=bool)
    xbar0   = zeros((N), dtype=float32)
    ybar0   = zeros((N), dtype=float32)

    Files = []
    Files.extend(  Mono_files)
    Files.extend(   RGB_files)
    Files.extend(NonFit_files)
    name = ''
    filetype = Files[0].split('.')[-1]
    if (filetype == 'f'):
      hdul = fits.open(Files[0])   # FYI... HDUL = Header Data Unit List
      for h in hdul[0].header:
        if h.upper() == 'OBJECT':
            name = hdul[0].header['OBJECT']
    if name == '':
        name = objectid.GetValue()
    if name == '':
        name = (((Files[0].split(slash)[-1]).split('.fit')[0]).split('-')[0]).split('_')[0]
    if (name != ''):
        [ra, dec] = getRaDec(name)
        if ((ra == 0) & (dec == 0)):
            ra_str  = objectra.GetValue()
            dec_str = objectdec.GetValue()
            ra  = min_sec_to_decimal(ra_str)
            dec = min_sec_to_decimal(dec_str)
            if ((ra == 0) & (dec == 0)):
                print(tab,'ra/dec approximation failed: ', name)
                return
        print(tab,'approx ra/dec: ', ra, dec)
        approx_ra_dec = [ ra, dec ]
    for j in range(N):
        if (scale[L_idx[j]] == 0):
            print('DOING IT: ', j, Files[j])
            [ra_ctr[j], dec_ctr[j], dtheta[j], scale[j], crpix1[j], crpix2[j], flip_it[j], success[j]] = \
                astrometry(images[L_idx[j],:,:], approx_ra_dec=[ra, dec], \
                                approx_fov=float(objectfov.GetValue())/2, thresh=0.3, showAlign=showAlign)
            ra_ctr[L_idx[j]]     = ra_ctr[j]
            dec_ctr[L_idx[j]]    = dec_ctr[j]
            dtheta[L_idx[j]]     = dtheta[j]
            scale[L_idx[j]]      = scale[j]
            crpix1[L_idx[j]]     = crpix1[j]
            crpix2[L_idx[j]]     = crpix2[j]
            flip_it[L_idx[j]]    = flip_it[j]
            success[L_idx[j]]    = success[j]
        else:
            ra_ctr[j]     = ra_ctr[L_idx[j]]
            dec_ctr[j]    = dec_ctr[L_idx[j]]
            dtheta[j]     = dtheta[L_idx[j]]
            scale[j]      = scale[L_idx[j]]
            crpix1[j]     = crpix1[L_idx[j]]
            crpix2[j]     = crpix2[L_idx[j]]
            flip_it[j]    = flip_it[L_idx[j]]
            success[j]    = success[L_idx[j]]
        if (success[j] == False):
            print(tab,'file ',j,' alignment failed')
            print(tab,'file ',j,': ra = ', ra_ctr[j], \
                               ', dec = ', dec_ctr[j], \
                            ', dtheta = ', dtheta[j], \
                             ', scale = ', scale[j], \
                            ', crpix1 = ', crpix1[j], \
                            ', crpix2 = ', crpix2[j] \
                 )
            print(''); print('');
        else:
            print(tab,'file ',j,': ra = ', ra_ctr[j], \
                           ', dec = ', dec_ctr[j], \
                           ', dtheta = ', dtheta[j], \
                           ', scale = ', scale[j], \
                           ', crpix1 = ', crpix1[j], \
                           ', crpix2 = ', crpix2[j] \
             )
            print(''); print('');
            '''
            if (j==0):
                print('Updating FOV')
                objectfov.SetValue(str(scale[j]))
            '''

    xbar0[0] = crpix1[0]
    ybar0[0] = crpix2[0]
    for j in range(N):
        print('flip_it = ', flip_it[j])

#       dra  =  ra_ctr[j]*cos(dec_ctr[j]*pi/180) -  ra_ctr[0]*cos(dec_ctr[0]*pi/180)
        dra  =  (ra_ctr[j] -  ra_ctr[0])*(cos(dec_ctr[j]*pi/180) + cos(dec_ctr[0]*pi/180))/2
        ddec = dec_ctr[j] - dec_ctr[0]
        dra  *= 15

        if (flip_it[j] == True): dra *= -1

        '''
        print(tab,'dra = ', dra, ', ddec = ', ddec)
        print(tab,'dtheta[0] = ', dtheta[0])
        '''
        dt0 = dtheta[0]*pi/180
        xx =  cos(dt0)*dra + sin(dt0)*ddec
        yy = -sin(dt0)*dra + cos(dt0)*ddec

        if (flip_it[j] == True): 
            xx  *= -1
            yy  *= -1

        xbar0[j] = xbar0[0] + xx*3600/scale[0]
        ybar0[j] = ybar0[0] + yy*3600/scale[0]
        '''
        print(tab,'xbar0[j] - xbar0[0] = ', xbar0[j]-xbar0[0])
        print(tab,'ybar0[j] - ybar0[0] = ', ybar0[j]-ybar0[0])
        '''
    
    xmin = int(np.floor(np.min((xbar0-xbar0[0])[success]))-1)
    xmax = int( np.ceil(np.max((xbar0-xbar0[0])[success]))+1)
    ymin = int(np.floor(np.min((ybar0-ybar0[0])[success]))-1)
    ymax = int( np.ceil(np.max((ybar0-ybar0[0])[success]))+1)
    x_argmin = np.argmin(xbar0)
    x_argmax = np.argmax(xbar0)
    y_argmin = np.argmin(ybar0)
    y_argmax = np.argmax(ybar0)
    print(tab,'xmin, xmax = ', xmin, xmax)
    print(tab,'ymin, ymax = ', ymin, ymax)

#   Show the extreme alignments

    print('')
    print(tab,'   Image       Rescaling    Rotation    Good?')
    print(tab,' Name/Number    percent       degs')
    print(tab,'-------------  ---------    --------    -----')
    progress.SetLabel('Stacking images: ')
    alphas = ones((N), dtype=float32)
    thetas = zeros((N), dtype=float32)
    for k in range(0,N):
        alphas[k]= scale[k]/scale[0]
        thetas[k]= (dtheta[k]-dtheta[0])*pi/180
    alpha_median = percentile(alphas,50)
    alpha_iqr = percentile(alphas, 75) - percentile(alphas, 25)

    m1 = m+ymax-ymin
    n1 = n+xmax-xmin
    align_stack = zeros((N,m1,n1), dtype=float32)
    i0 =  ymax
    i1 =  ymax+m
    j0 =  xmax
    j1 =  xmax+n
    xs = arange(0,n-xmin+xmax)
    ys = array([arange(0,m-ymin+ymax)]).T
    Xs = ones((m1,1))*xs
    Ys = ys*ones((1,n1))
    kk = 0
    for k in range(0,N):
      if (L_idx[k] != -1):
        gauge.SetValue(int(100*k/(N-1)))
        wx.Yield()
        if ((FileNames != '') & (k<len(FileNames))): 
            print(tab,'{fname:13s}    {rescale:9.4f}    {rotate:8.3f}     {good}'.format( \
                            fname = FileNames[k], rescale = abs(alphas[k]-1)*100, \
                            rotate = thetas[k]*180/pi,
                            good = success[k] \
                            ) )
        else:
            print(tab,'   {num:3d}          {rescale:9.4f}    {rotate:8.3f}     {good}'.format( \
                            num = k, rescale = abs(alphas[k]-1)*100, \
                            rotate = thetas[k]*180/pi, \
                            good = success[k] \
                            ) )
#       print(tab,'alphas[k] = ', alphas[k], ', alpha_median= ', alpha_median, ', alpha_iqr = ', alpha_iqr)

        build_aligned_stack(k,kk,images)

        kk += 1
    print('\n')
    print(tab,'Number stacked: ', kk)
#   print(tab,'Number not well aligned: ', N-kk)
    N = kk
    align_stack = align_stack[:, ymax:m+ymin, xmax:n+xmin]
    progress.SetLabel('')
    gauge.SetValue(0)

    tab = tab[0:-4]
    printdone('Done aligning via Astrometry')

def align_cv(images):
    global stack_size, stack_num, align_stack
    global FileNames
    global progress, gauge
    global i0, i1, j0, j1, m1, n1, Xs, Ys, alphas, thetas, xbar0, ybar0

    showAlign = True

    cv2.ocl.setUseOpenCL(False)

    (N,m,n) = images.shape
    stack_num  = 0
    stack_size = N

    # use image 0 as the reference image
    queryImg = images[0,:,:]

    feature_extractor = 'brisk' # one of 'sift', 'surf', 'brisk', 'orb'
    feature_matching = 'bf'

    logquery = np.uint8(np.floor(255*log(queryImg-np.min(queryImg)+1)/np.max(log(queryImg-np.min(queryImg)+1))))
    kpsB, featuresB = detectAndDescribe(logquery, method=feature_extractor)

    progress.SetLabel('Aligning images: ')
    H = zeros((stack_size,3,3), dtype=float)
    j = 0
    for j in range(N):
        gauge.SetValue(int(100*j/(stack_size-1)))
        wx.Yield()
        if ((FileNames != '') & (j<len(FileNames))): 
            print(tab,'working on image ',FileNames[j])
        else:
            print(tab,'working on image ',j)
        trainImg = images[j,:,:]

        logtrain = np.uint8(np.floor(255*log(trainImg-np.min(trainImg)+1)/np.max(log(trainImg-np.min(trainImg)+1))))
        kpsA, featuresA = detectAndDescribe(logtrain, method=feature_extractor)

        if feature_matching == 'bf':
            matches = matchKeyPointsBF(featuresA, featuresB, method=feature_extractor)
        elif feature_matching == 'knn':
            matches = matchKeyPointsKNN(featuresA, featuresB, ratio=0.75, method=feature_extractor)

        M = getHomography(kpsA, kpsB, featuresA, featuresB, matches, reprojThresh=4)
        if M is None:
            print("Error!")
        (matches, H[j,:,:], status) = M

    di = np.int(-np.floor(np.min([np.min(H[:,0,2]),0])))
    dj = np.int(-np.floor(np.min([np.min(H[:,1,2]),0])))
    print('H[:,0,2] = ', H[:,0,2])
    print('H[:,1,2] = ', H[:,1,2])
    print('di, dj = ', di, dj)
    di_max = np.int(np.floor(np.max([np.max(H[:,0,2]),0])))
    dj_max = np.int(np.floor(np.max([np.max(H[:,1,2]),0])))

    width  = queryImg.shape[1] + di + di_max
    height = queryImg.shape[0] + dj + dj_max

    progress.SetLabel('')
    gauge.SetValue(0)

    align_stack = zeros((N,height,width), dtype=float32)

    align_stack[0,dj:queryImg.shape[0]+dj, di:queryImg.shape[1]+di] = queryImg
        
    for j in arange(1,N):
        gauge.SetValue(int(100*j/(N-1)))
        wx.Yield()

        H[j,0,2] += di
        H[j,1,2] += dj
        trainImg = images[j,:,:]
        align_stack[j,:,:] = cv2.warpPerspective(trainImg, H[j,:,:], (width, height))

    if (showAlign):
        rcParams.update({'figure.max_open_warning': 0})
        rcParams['figure.max_open_warning'] = 0
        for j in range(N):
            figure(80+j)
            I = align_stack[j]
            imshow(maximum(minimum(I, percentile(I,99))-percentile(I,50),0), cmap='gray')
            if ((FileNames != '') & (j<len(FileNames))): title(FileNames[j])
            draw()

    progress.SetLabel('')
    gauge.SetValue(0)

def autoalign3stars(images, thresh=0.1, border=10, max_rotate=30, max_rescale=1):
    global xbar, ybar, xbar0, ybar0, xbar1, ybar1, stack_size, stack_num, star, align_stack, img, fig
    global showAlign, FileNames
    global fwhms
    global xmin, xmax, ymin, ymax
    global progress, gauge
    global L_idx
    global i0, i1, j0, j1, m1, n1, Xs, Ys, alphas, thetas

    global tab
    print(tab, 'Autoaligning using 3 stars:')
    tab = tab+'    '

    (N,m,n) = images.shape
    num_stars = 4
    xbar = zeros((num_stars))
    ybar = zeros((num_stars))
    dist = zeros((num_stars,num_stars))
    xbar0 = zeros(N)
    ybar0 = zeros(N)
    xbar1 = zeros(N)
    ybar1 = zeros(N)
    xbar2 = zeros(N)
    ybar2 = zeros(N)
    xbar3 = zeros(N)
    ybar3 = zeros(N)
    quality = zeros(N)
    star = 0
    stack_num  = 0
    stack_size = N

#   print(tab, 'max_rescale set to ', max_rescale)
#   print(tab, 'max_rotate set to ', max_rotate)
    
    # use image 0 as the reference image
    L = images[L_idx[0],:,:]

    # find brightest stars
    k, LconvL2 = findstars(L,thresh=thresh,fignum=0,bkg=True)

    xs = arange(0,n)
    ys = array([arange(0,m)]).T
    Xs = ones((m,1))*xs
    Ys = ys*ones((1,n))

    mask = (Xs > border) & (Xs < n-border) & (Ys > border) & (Ys < m-border) & (L < 0.8*np.max(L)) # \
#          & (LconvL2 > 0.4*np.max(LconvL2)) & (LconvL2 < 0.9*np.max(LconvL2))
    k = k&mask;

    flux = LconvL2[k]
    Xs = Xs[k]
    Ys = Ys[k]
    '''
    print(tab, 'Xs:   ', Xs)
    print(tab, 'Ys:   ', Ys)
    print(tab, 'flux: ', flux)
    '''

    if (len(Xs) < 10):
        thresh = thresh/2
        autoalign3stars(images, thresh=thresh, border=border, \
                        max_rotate=max_rotate, max_rescale=max_rescale)
        return

    I = argsort(-flux)
    flux = flux[I]
    Xs_sort = int16(Xs[I])
    Ys_sort = int16(Ys[I])

    print(tab, 'number of stars found = ', len(Xs))
    num_stars = min([15, len(Xs)])
    if num_stars < 15: 
        thresh = thresh/2
        print(tab, '1: reducing threshold to ', thresh)
        autoalign3stars(images, thresh=thresh, border=border, max_rotate=max_rotate, \
                max_rescale=max_rescale)
        return

    thresh = flux[num_stars-1]/np.max(LconvL2)
    print(tab, 'updated thresh = ', thresh)

    num_stars3 = int(num_stars**3)
    distdata = zeros((num_stars3, 5))
    maxflux = np.max(LconvL2)
    k = 0
    for i0 in arange(0, num_stars):
        if flux[i0] > 0.8*maxflux: continue
        for i1 in arange(i0+1, num_stars):
            if flux[i1] > 0.8*maxflux: continue
            for i2 in arange(i1+1, num_stars):
#               for i3 in arange(i2+1, num_stars):
                    if flux[i2] > 0.8*maxflux: continue
                    distdata[k,0] = i0
                    distdata[k,1] = i1
                    distdata[k,2] = i2
#                   distdata[k,3] = i3
                    distdata[k,4] = sqrt( (Xs_sort[i0]-Xs_sort[i1])**2 + \
                                          (Ys_sort[i0]-Ys_sort[i1])**2 ) \
                                          *flux[i0]*flux[i1] \
                                  + sqrt( (Xs_sort[i1]-Xs_sort[i2])**2 + \
                                          (Ys_sort[i1]-Ys_sort[i2])**2 ) \
                                          *flux[i1]*flux[i2] \
                                  + sqrt( (Xs_sort[i2]-Xs_sort[i0])**2 + \
                                          (Ys_sort[i2]-Ys_sort[i0])**2 ) \
                                          *flux[i2]*flux[i0] 
#                                 + sqrt( (Xs_sort[i2]-Xs_sort[i3])**2 + \
#                                         (Ys_sort[i2]-Ys_sort[i3])**2 ) \
#                                 + sqrt( (Xs_sort[i3]-Xs_sort[i0])**2 + \
#                                         (Ys_sort[i3]-Ys_sort[i0])**2 )

                    equilateral = True
                    if (equilateral == True):
                        theta01 = arctan2(Ys_sort[i1]-Ys_sort[i0],Xs_sort[i1]-Xs_sort[i0]) 
                        theta12 = arctan2(Ys_sort[i2]-Ys_sort[i1],Xs_sort[i2]-Xs_sort[i1]) 
                        theta20 = arctan2(Ys_sort[i0]-Ys_sort[i2],Xs_sort[i0]-Xs_sort[i2]) 
                        theta10 = arctan2(Ys_sort[i0]-Ys_sort[i1],Xs_sort[i0]-Xs_sort[i1]) 
                        theta21 = arctan2(Ys_sort[i1]-Ys_sort[i2],Xs_sort[i1]-Xs_sort[i2]) 
                        theta02 = arctan2(Ys_sort[i2]-Ys_sort[i0],Xs_sort[i2]-Xs_sort[i0]) 
    
                        theta0102 = maximum(theta01,theta02) - minimum(theta01,theta02)
                        theta1210 = maximum(theta12,theta10) - minimum(theta12,theta10)
                        theta2021 = maximum(theta20,theta21) - minimum(theta20,theta21)
    
                        if theta0102 > pi:  theta0102 -= pi
                        if theta1210 > pi:  theta1210 -= pi
                        if theta2021 > pi:  theta2021 -= pi
    
                        distdata[k,4] /= (1 + abs(theta0102-pi/3) \
                                            + abs(theta2021-pi/3) \
                                            + abs(theta1210-pi/3))
                    k = k+1
    I  = argsort(-distdata[0:k,4])
    k  = I[0]
#   k  = argmax(distdata[0:k,2])
    i0 = distdata[k,0]
    i1 = distdata[k,1]
    i2 = distdata[k,2]
#   i3 = distdata[k,3]
    
    xbar0[0] = Xs_sort[int(i0)]
    ybar0[0] = Ys_sort[int(i0)]
    xbar1[0] = Xs_sort[int(i1)]
    ybar1[0] = Ys_sort[int(i1)]
    xbar2[0] = Xs_sort[int(i2)]
    ybar2[0] = Ys_sort[int(i2)]
#   xbar3[0] = Xs_sort[int(i3)]
#   ybar3[0] = Ys_sort[int(i3)]

    if (showAlign):
            rcParams.update({'figure.max_open_warning': 0})
            rcParams['figure.max_open_warning'] = 0
            figure(80)
            I = images[L_idx[0]]
            imshow(maximum(minimum(I, percentile(I,99))-percentile(I,20),0), cmap='gray')
            plot(xbar0[0],ybar0[0],'r.')
            plot(xbar1[0],ybar1[0],'b.')
            plot(xbar2[0],ybar2[0],'g.')
#           plot(xbar3[0],ybar3[0],'m.')
            if (FileNames != ''): title(FileNames[0])
            draw()
#           show()

    dist01 = sqrt((xbar0[0]-xbar1[0])**2 + (ybar0[0]-ybar1[0])**2) 
    dist02 = sqrt((xbar0[0]-xbar2[0])**2 + (ybar0[0]-ybar2[0])**2) 
#   dist03 = sqrt((xbar0[0]-xbar3[0])**2 + (ybar0[0]-ybar3[0])**2) 
    dist12 = sqrt((xbar1[0]-xbar2[0])**2 + (ybar1[0]-ybar2[0])**2) 
#   dist13 = sqrt((xbar1[0]-xbar3[0])**2 + (ybar1[0]-ybar3[0])**2) 
#   dist23 = sqrt((xbar2[0]-xbar3[0])**2 + (ybar2[0]-ybar3[0])**2) 
#   print(tab, 'dists = ', dist01, dist02)

    alpha01 = sqrt((ybar1[0]-ybar0[0])**2 + (xbar1[0]-xbar0[0])**2)
    alpha02 = sqrt((ybar2[0]-ybar0[0])**2 + (xbar2[0]-xbar0[0])**2)
    theta01 = arctan2((ybar1[0]-ybar0[0]),(xbar1[0]-xbar0[0]))
    theta02 = arctan2((ybar2[0]-ybar0[0]),(xbar2[0]-xbar0[0]))
    print(tab, '')
    j = 0
    progress.SetLabel('Aligning images: ')
    while j<stack_size:
      if (xbar0[L_idx[j]] != 0):
        if (L_idx[j] != -1):
          xbar0[j] = xbar0[L_idx[j]] 
          ybar0[j] = ybar0[L_idx[j]] 
          xbar1[j] = xbar1[L_idx[j]] 
          ybar1[j] = ybar1[L_idx[j]] 
          xbar2[j] = xbar2[L_idx[j]] 
          ybar2[j] = ybar2[L_idx[j]] 
      else:
        gauge.SetValue(int(100*j/(stack_size-1)))
        wx.Yield()
#       print(tab, 'j = ', j)
        if ((FileNames != '') & (j<len(FileNames))): 
            print(tab, 'working on image ',FileNames[j])
        else:
            print(tab, 'working on image ',j)
        L = images[L_idx[j],:,:]
#       k, LconvL2 = findstars(L,thresh=thresh/2,fignum=2+j)
        k, LconvL2 = findstars(L,thresh=thresh,fignum=0)

        xs = arange(0,n)
        ys = array([arange(0,m)]).T
        Xs = ones((m,1))*xs
        Ys = ys*ones((1,n))

        mask = (Xs > 11) & (Xs < n-11) & (Ys > 11) & (Ys < m-11) 
        k = k&mask;

        flux = LconvL2[k]
#       print(tab, 'Xs,Ys at [569,280]: ', Xs[569,280], Ys[569,280])
        Xs = Xs[k]
        Ys = Ys[k]
#       print(tab, 'k[569,280] = ', k[569,280])
#       print(tab, 'Xs,Ys at 78: ', Xs[78], Ys[78])

        num_stars = min([120, len(Xs)])
        if num_stars < 3: 
            thresh = thresh/2
            print(tab, '2: reducing threshold to ', thresh)
            continue

        I = argsort(-flux)
        flux = flux[I]
        Xs_sort = int16(Xs[I])
        Ys_sort = int16(Ys[I])
    
#       if (Xs_sort[74] == 280) & (Ys_sort[74] == 569): print(tab, 'it is here: ', I[74])
    
        num_stars2 = int(num_stars**3)
        matchdata = 1000*ones((num_stars2, 9))
        k = 0
#       print(tab, 'number of candidate stars = ',num_stars)
        min_rotate = 1000000
        for i0 in arange(0, num_stars):
            xbar0[j] = Xs_sort[i0]
            ybar0[j] = Ys_sort[i0]
            for i1 in arange(0, num_stars):
                if (i1 == i0): continue
                xbar1[j] = Xs_sort[i1]
                ybar1[j] = Ys_sort[i1]

                alpha01j = sqrt((ybar1[j]-ybar0[j])**2 + (xbar1[j]-xbar0[j])**2)
                if abs(alpha01j-alpha01)/alpha01 > max_rescale/100: continue
#               else: print(tab,  abs(alpha01j-alpha01)/alpha01 )

                for i2 in arange(0, num_stars):
                        if ((i2 == i0) | (i2 == i1)): continue
                        xbar2[j] = Xs_sort[i2]
                        ybar2[j] = Ys_sort[i2]

                        '''
                        if ((abs(xbar0[j]-xbar0[0])<20)
                           &(abs(ybar0[j]-ybar0[0])<20)
                           &(abs(xbar1[j]-xbar1[0])<20)
                           &(abs(ybar1[j]-ybar1[0])<20)
                           &(abs(xbar2[j]-xbar2[0])<20)
                           &(abs(ybar2[j]-ybar2[0])<20)):
                                print(tab, 'yay:  k = ',k,', i0 = ',i0,', i1 = ',i1,', i2 = ',i2)
                        '''

                        theta01j = arctan2((ybar1[j]-ybar0[j]),(xbar1[j]-xbar0[j])) 
                        theta02j = arctan2((ybar2[j]-ybar0[j]),(xbar2[j]-xbar0[j])) 
                        dtheta01 = theta01j - theta01 + pi
                        dtheta02 = theta02j - theta02 + pi
#                       print(tab, 'angle:', abs(dtheta01%(2*pi)-pi)*180/pi)
#                       print(tab, 'angle:', abs(dtheta02%(2*pi)-pi)*180/pi)
                        if abs(dtheta01%(2*pi)-pi)*180/pi < min_rotate:
                            min_rotate = abs(dtheta01%(2*pi)-pi)*180/pi
                        if abs(dtheta01%(2*pi)-pi)*180/pi > max_rotate: continue
                        if abs(dtheta02%(2*pi)-pi)*180/pi > max_rotate: continue
#                       print(tab, '      k = ',k,', i0 = ',i0,', i1 = ',i1,', i2 = ',i2)
#                       print(tab, abs(dtheta01%(2*pi)-pi)*180/pi)
                    
                        matchdata[k,0] = i0
                        matchdata[k,1] = i1
                        matchdata[k,2] = i2

                        dist01j = sqrt((xbar0[j]-xbar1[j])**2 + (ybar0[j]-ybar1[j])**2) 
                        dist02j = sqrt((xbar0[j]-xbar2[j])**2 + (ybar0[j]-ybar2[j])**2) 
                        dist12j = sqrt((xbar1[j]-xbar2[j])**2 + (ybar1[j]-ybar2[j])**2) 

                        matchdata[k,4] = \
                            maximum(
                            abs(dist02j/dist01j-dist02/dist01) , \
                            abs(dist12j/dist01j-dist12/dist01)  
                            )

                        matchdata[k,5] = dist02j/dist01j
                        matchdata[k,6] = dist02 /dist01 
                        matchdata[k,7] = dist12j/dist01j
                        matchdata[k,8] = dist12 /dist01 

                        '''
                        if ((ybar0[j]==472) & (xbar0[j]==679) & (ybar1[j]==222) & (xbar0[j]==1295) ):
                                print(tab, matchdata[k,:])
                        '''

                        k = k+1
#       print(tab, 'number of star trios found = ',k)
#       print(tab, 'minimum rotation = ', min_rotate)
        
        print(tab, 'j = ', j, ', k = ', k);
        if k>0:
            k = argmin(matchdata[0:k,4])
    
            i0 = matchdata[k,0]
            i1 = matchdata[k,1]
            i2 = matchdata[k,2]
#           print(tab, 'i0 = ', int(i0), ', i1 = ', int(i1), ', i2 = ', int(i2))
        
            xbar0[j] = Xs_sort[int(i0)]
            ybar0[j] = Ys_sort[int(i0)]
            xbar1[j] = Xs_sort[int(i1)]
            ybar1[j] = Ys_sort[int(i1)]
            xbar2[j] = Xs_sort[int(i2)]
            ybar2[j] = Ys_sort[int(i2)]

#           print(tab, 'alpha01 = ', matchdata[k,5],', alpha01j = ', matchdata[k,6])
    
            quality[j] = matchdata[k,4]
        else:
            xbar0[j] = xbar0[0]
            ybar0[j] = ybar0[0]
            xbar1[j] = xbar1[0]
            ybar1[j] = ybar1[0]
            xbar2[j] = xbar2[0]
            ybar2[j] = ybar2[0]

        if (showAlign):
            rcParams.update({'figure.max_open_warning': 0})
            rcParams['figure.max_open_warning'] = 0
            figure(80+j)
            I = images[j]
            imshow(maximum(minimum(I, percentile(I,99))-percentile(I,20),0), cmap='gray')
            plot(xbar0[j],ybar0[j],'r.')
            plot(xbar1[j],ybar1[j],'b.')
            plot(xbar2[j],ybar2[j],'g.')
            if ((FileNames != '') & (j<len(FileNames))): title(FileNames[j])
            draw()
#           show()

        xbar0[L_idx[j]] = xbar0[j]
        ybar0[L_idx[j]] = ybar0[j]
        xbar1[L_idx[j]] = xbar1[j]
        ybar1[L_idx[j]] = ybar1[j]
        xbar2[L_idx[j]] = xbar2[j]
        ybar2[L_idx[j]] = ybar2[j]

      j += 1
    progress.SetLabel('')
    gauge.SetValue(0)

    '''
    '''
    r0 = 2
    r1 = 4
    r2 = 8
    box_size0 = 2*r2+1
    fwhms = zeros(stack_size)
    for j in arange(0,stack_size):
        xmin = int(round(xbar0[j]-r2))
        ymin = int(round(ybar0[j]-r2))
        xmax = xmin+box_size0
        ymax = ymin+box_size0
        image = images[j]
        img = float32(array(image))
#
#       Computing the local median corrects errors caused by hot pixels...
        Arect1pile = zeros((box_size0,box_size0,9), dtype=float32)
#       print(tab, 'j = ', j, ', Arect1pile.shape = ', Arect1pile.shape, ', img.shape = ', img.shape)
#       print(tab, 'ymin/ymax = ', ymin, ymax, ',  xmin/max = ', xmin, xmax)
        for ii in arange(-1,2):
            for jj in arange(-1,2):
                Arect1pile[:,:,3*(ii+1)+(jj+1)] = img[ymin+ii:ymax+ii,xmin+jj:xmax+jj]
        Arect1 = percentile(Arect1pile,50,axis=2)

#
#       This is the simple version without the median...
#       Arect1 = img[ymin:ymax,xmin:xmax]
#
        xs =        arange(-r2,r2+1)
        ys = array([arange(-r2,r2+1)]).T
        radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
        disk = (radii < r0)
        annulus = (radii < r2) & (radii > r1)
        vals0 = Arect1[disk]
        vals2 = Arect1[annulus]
    
        maxval = np.max(vals0)
        bkg_avg = np.mean(vals2)
        bkg_std = np.std(vals2)
        bkg_median = percentile(vals2,50)
        bkg_iqr = percentile(vals2, 75) - percentile(vals2, 25)
        flux = np.sum(vals0-bkg_avg)
        Arect0 = Arect1[r2-r0:r2+r0+1,r2-r0:r2+r0+1]
        xs0 =        arange(-r0,r0+1)
        ys0 = array([arange(-r0,r0+1)]).T
        denom = sum(maximum(Arect0-bkg_median,1e-5))
        xbar0[j] = sum(maximum(Arect0-bkg_median,0)*xs0)/denom
        ybar0[j] = sum(maximum(Arect0-bkg_median,0)*ys0)/denom

        if (showAlign):
                figure(200+j)
                imshow(Arect1, cmap='gray', norm=LogNorm())
                plot(r2+xbar0[j], r2+ybar0[j], 'ro')
                if ((FileNames != '') & (j<len(FileNames))): title(FileNames[j])
                draw()
#               show()
        fwhms[j] = fwhm(Arect1,r2,r2)
        xbar0[j] = xmin+r2+xbar0[j]
        ybar0[j] = ymin+r2+ybar0[j]

#   print(tab, '(ybar0[0],xbar0[0]) = (',ybar0[0],',',xbar0[0],')')
    print(tab, '')
    print(tab, 'FWHMs = ',fwhms)
        
    for j in arange(0,stack_size):
        xmin = int(round(xbar1[j]-r2))
        ymin = int(round(ybar1[j]-r2))
        xmax = xmin+box_size0
        ymax = ymin+box_size0
        image = images[j]
        img = float32(array(image))
#
#       Computing the local median corrects errors caused by hot pixels...
        Arect1pile = zeros((box_size0,box_size0,9), dtype=float32)
        for ii in arange(-1,2):
            for jj in arange(-1,2):
                Arect1pile[:,:,3*(ii+1)+(jj+1)] = img[ymin+ii:ymax+ii,xmin+jj:xmax+jj]
        Arect1 = percentile(Arect1pile,50,axis=2)

#
#       This is the simple version without the median...
#       Arect1 = img[ymin:ymax,xmin:xmax]
#
        xs =        arange(-r2,r2+1)
        ys = array([arange(-r2,r2+1)]).T
        radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
        disk = (radii < r0)
        annulus = (radii < r2) & (radii > r1)
        vals0 = Arect1[disk]
        vals2 = Arect1[annulus]
    
        maxval = np.max(vals0)
        bkg_avg = np.mean(vals2)
        bkg_std = np.std(vals2)
        bkg_median = percentile(vals2,50)
        bkg_iqr = percentile(vals2, 75) - percentile(vals2, 25)
        flux = np.sum(vals0-bkg_avg)
        Arect0 = Arect1[r2-r0:r2+r0+1,r2-r0:r2+r0+1]
        xs0 =        arange(-r0,r0+1)
        ys0 = array([arange(-r0,r0+1)]).T
        xbar1[j] = sum(maximum(Arect0-bkg_median,0)*xs0)/sum(maximum(Arect0-bkg_median,0)) 
        ybar1[j] = sum(maximum(Arect0-bkg_median,0)*ys0)/sum(maximum(Arect0-bkg_median,0)) 

        if (showAlign):
                figure(300+j)
                imshow(Arect1, cmap='gray', norm=LogNorm())
                plot(r2+xbar1[j], r2+ybar1[j], 'bo')
                if ((FileNames != '') & (j<len(FileNames))): title(FileNames[j])
                draw()
#               show()
        xbar1[j] = xmin+r2+xbar1[j]
        ybar1[j] = ymin+r2+ybar1[j]
    '''
    '''
    
#   print(tab, xbar0, ybar0)
    '''
    if (showAlign):
        rcParams.update({'figure.max_open_warning': 0})
        rcParams['figure.max_open_warning'] = 0
        for k in range(0,N):
            figure(80+k)
            imshow(images[k], cmap='gray', norm=LogNorm())
            plot(xbar0[k],ybar0[k],'r.')
            plot(xbar1[k],ybar1[k],'b.')
            plot(xbar2[k],ybar2[k],'g.')
#           plot(xbar3[k],ybar3[k],'m.')
            draw()
            show()
    '''

    '''
    print(tab, '(x0,y0) = ', xbar0[0], ybar0[0], ', ', \
          '(x1,y1) = ', xbar0[1], ybar0[1])
    '''
        
    (N,m,n) = images.shape
    xmin = int(np.floor(np.min(xbar0-xbar0[0]))-1)
    xmax = int( np.ceil(np.max(xbar0-xbar0[0]))+1)
    ymin = int(np.floor(np.min(ybar0-ybar0[0]))-1)
    ymax = int( np.ceil(np.max(ybar0-ybar0[0]))+1)
    x_argmin = np.argmin(xbar0)
    x_argmax = np.argmax(xbar0)
    y_argmin = np.argmin(ybar0)
    y_argmax = np.argmax(ybar0)
    '''
    i0 = -ymin
    i1 = -ymin+m
    j0 = -xmin
    j1 = -xmin+n
    '''
    print(tab, 'xmin, xmax = ', xmin, xmax)
    print(tab, 'ymin, ymax = ', ymin, ymax)
    '''
    print(tab, 'x_argmin, x_argmax = ', x_argmin, x_argmax)
    print(tab, 'y_argmin, y_argmax = ', y_argmin, y_argmax)
    '''
    dx_max = int16(ceil(max(abs(xbar0-xbar0[0]))))
    dy_max = int16(ceil(max(abs(ybar0-ybar0[0]))))

#   Show the extreme alignments
    '''
    figure(500)
    I = images[x_argmin]
    imshow(maximum(minimum(I, percentile(I,99))-percentile(I,20),0), cmap='gray')
    plot(xbar0[x_argmin],ybar0[x_argmin],'r.')
    plot(xbar1[x_argmin],ybar1[x_argmin],'b.')
    plot(xbar2[x_argmin],ybar2[x_argmin],'g.')
    if (FileNames != ''): title(FileNames[x_argmin])
    draw()
#   show()
    figure(501)
    I = images[x_argmax]
    imshow(maximum(minimum(I, percentile(I,99))-percentile(I,20),0), cmap='gray')
    plot(xbar0[x_argmax],ybar0[x_argmax],'r.')
    plot(xbar1[x_argmax],ybar1[x_argmax],'b.')
    plot(xbar2[x_argmax],ybar2[x_argmax],'g.')
    if (FileNames != ''): title(FileNames[x_argmax])
    draw()
#   show()
    figure(502)
    I = images[y_argmin]
    imshow(maximum(minimum(I, percentile(I,99))-percentile(I,20),0), cmap='gray')
    plot(xbar0[y_argmin],ybar0[y_argmin],'r.')
    plot(xbar1[y_argmin],ybar1[y_argmin],'b.')
    plot(xbar2[y_argmin],ybar2[y_argmin],'g.')
    if (FileNames != ''): title(FileNames[y_argmin])
    draw()
#   show()
    figure(503)
    I = images[y_argmax]
    imshow(maximum(minimum(I, percentile(I,99))-percentile(I,20),0), cmap='gray')
    plot(xbar0[y_argmax],ybar0[y_argmax],'r.')
    plot(xbar1[y_argmax],ybar1[y_argmax],'b.')
    plot(xbar2[y_argmax],ybar2[y_argmax],'g.')
    if (FileNames != ''): title(FileNames[y_argmax])
    draw()
#   show()
    '''

    dist0 = sqrt((xbar1[0]-xbar0[0])**2+(ybar1[0]-ybar0[0])**2)
    theta0 = arctan2(ybar1[0]-ybar0[0],xbar1[0]-xbar0[0])
    print(tab, '')
    print(tab, '   Image        Quality    Rescaling    Rotation    Good?')
    print(tab, ' Name/Number   0 is best    percent       degs')
    print(tab, '-------------  ---------   ---------    --------    -----')
    progress.SetLabel('Stacking images: ')
    alphas = ones((N), dtype=float32)
    thetas = zeros((N), dtype=float32)
    for k in range(0,N):
        distk = sqrt((xbar1[k]-xbar0[k])**2+(ybar1[k]-ybar0[k])**2)
        alphas[k]= distk/dist0
        thetas[k]= arctan2(ybar1[k]-ybar0[k],xbar1[k]-xbar0[k])-theta0
    alpha_median = percentile(alphas,50)
    alpha_iqr = percentile(alphas, 75) - percentile(alphas, 25)
    '''
    print(tab, 'alphas = ', alphas)
    print(tab, 'alpha_median = ', alpha_median)
    print(tab, 'alpha_iqr = ', alpha_iqr)
    '''
    m1 = m + 2*dy_max
    n1 = n + 2*dx_max
    align_stack = zeros((N,m1,n1), dtype=float32)
    i0 =  dy_max
    i1 =  dy_max+m
    j0 =  dx_max
    j1 =  dx_max+n
    xs = arange(0,n1)
    ys = array([arange(0,m1)]).T
    Xs = ones((m1,1))*xs
    Ys = ys*ones((1,n1))
    kk = 0
    for k in range(0,N):
      if (L_idx[k] != -1):
        gauge.SetValue(int(100*k/(N-1)))
        wx.Yield()
        if ((FileNames != '') & (k<len(FileNames))): 
            print(tab, '{fname:13s}    {quality:7.4f}    {rescale:9.4f}    {rotate:8.3f}     {good}'.format( \
                            fname = FileNames[k], quality = quality[k], rescale = abs(alphas[k]-1)*100, \
                            rotate = thetas[k]*180/pi,
                            good = (abs(alphas[k] - alpha_median) < 3*alpha_iqr) \
                            ) )
        else:
            print(tab, '   {num:3d}          {quality:7.4f}    {rescale:9.4f}    {rotate:8.3f}     {good}'.format( \
                            num = k, quality = quality[k], rescale = abs(alphas[k]-1)*100, \
                            rotate = thetas[k]*180/pi, \
                            good = (abs(alphas[k] - alpha_median) < 3*alpha_iqr) \
                            ) )
        if (abs(alphas[k] - alpha_median) < 3*alpha_iqr):
            build_aligned_stack(k,kk,images)
            kk += 1
    print(tab, '\n')
    print(tab, 'Number stacked: ', kk)
#   print(tab, 'Number not well aligned: ', N-kk)
    N = kk
    align_stack = align_stack[:, dy_max:m+dy_max, dx_max:n+dx_max]
    progress.SetLabel('')
    gauge.SetValue(0)
    tab = tab[0:-4]
    printdone( 'Done autoaligning using 3 stars')

def autoalign(images, numstars=1, thresh=0.1, border=0, \
                max_rotate=30, max_rescale=1, # max_rotate in degrees, max_rescale in percentage \
                help=False, fnames=''):
    global xbar, ybar, xbar0, ybar0, xbar1, ybar1, stack_size, stack_num, star, align_stack, img, fig
    global showAlign, FileNames
    global MaxBoxSize
    global fwhms
    global xmin, xmax, ymin, ymax
    global progress, gauge
    global success
    global L_idx
    global i0, i1, j0, j1, m1, n1, Xs, Ys, alphas, thetas
    twostar = False

    global tab
    print(tab, 'Autoaligning')
    tab = tab+'    '

    if (fnames != []): FileNames = fnames

    name = objectid.GetValue()

    if (help != False):
        print('\nFunction autoalign parameters: \n', \
                '  images:  the stack of images to be aligned (Required) \n', \
                '  numstars:  the number of stars in each image to use to make the alignment '+ \
                '(1, 2 or 3).  \n     Default is 1. \n', \
                '  thresh:  smaller values give more stars from which to identify alignment '+ \
                '\n     Default is 0.1 \n', \
                '  border:  how far away from the edge of the image a star has to be \n'+ \
                '     and how big of a search box to use around a reference star \n'+ \
                '     Default is 30 \n', \
                '  max_rotate:  the maximum amount of rotation in degrees \n'+ \
                '     Default is 30 degrees \n', \
                '  max_rescale:  the maximum amount of image rescaling expressed '+\
                'as a percentage \n     Default is 1% \n', \
                '  fnames: list of file names being aligned \n', \
                '\n' )
    (N,m,n) = images.shape
    if (border == 0):  border = int(minimum(m,n)/5)
    MaxBoxSize = border

    if (numstars==0):
        align_stack = images
        return

    if (numstars==3):
        autoalign3stars(images,thresh=thresh,border=border, \
                max_rotate=max_rotate, max_rescale=max_rescale)
        return

    if (numstars==2):
        twostar = True

    (N,m,n) = images.shape
    num_stars = 4
    xbar = zeros((num_stars))
    ybar = zeros((num_stars))
    dist = zeros((num_stars,num_stars))
    xbar0 = zeros(N)
    ybar0 = zeros(N)
    xbar1 = zeros(N)
    ybar1 = zeros(N)
    star = 0
    stack_num  = 0
    stack_size = N
    
    # use image 0 as the reference image
    L = images[L_idx[0],:,:]

    # find brightest stars
    k, LconvL2 = findstars(L,thresh,fignum=0)
    if sum(k) < 10:
        k, LconvL2 = findstars(L,thresh/10,fignum=0)
        if sum(k) < 10:
            k, LconvL2 = findstars(L,thresh/100,fignum=0)
    xs = arange(0,n)
    ys = array([arange(0,m)]).T
    Xs = ones((m,1))*xs
    Ys = ys*ones((1,n))
    flux = LconvL2[k]
    Xs = Xs[k]
    Ys = Ys[k]

    I = argsort(-flux)
    flux = flux[I]
    Xs = int16(Xs[I])
    Ys = int16(Ys[I])

#   1st star...
    num_stars = len(Xs)
    for j in arange(0, num_stars):
        dist1 = min([Ys[j], m-1-Ys[j], Xs[j], n-1-Xs[j]])
        if (dist1 > MaxBoxSize):
            i0 = Ys[j]
            j0 = Xs[j]
            jj = j
            break

    box_size0 = MaxBoxSize

#   2nd star...
    j1 = -1
    maxdist3=0
    for j in arange(jj+1, num_stars):
        dist2 = min([Ys[j], m-1-Ys[j], Xs[j], n-1-Xs[j]])
        xdiff = 1.*(Xs[j] - Xs[jj])
        ydiff = 1.*(Ys[j] - Ys[jj])
        dist3 = sqrt( xdiff*xdiff + ydiff*ydiff )
        if ((dist2 > MaxBoxSize) & (dist3 > maxdist3)):
            maxdist3 = dist3
            i1 = Ys[j]
            j1 = Xs[j]
    '''
    for j in arange(jj+1, num_stars):
        dist2 = min([Ys[j], m-1-Ys[j], Xs[j], n-1-Xs[j]])
        xdiff = 1.*(Xs[j] - Xs[jj])
        ydiff = 1.*(Ys[j] - Ys[jj])
        dist3 = sqrt( xdiff*xdiff + ydiff*ydiff )

        if ((dist2 > MaxBoxSize) & (dist3 > 0.9*maxdist3)):
            i1 = Ys[j]
            j1 = Xs[j]
            break
    '''
    if (j1 == -1):
        print(tab, 'Only found one star.  Doing one-star alignment...')
        twostar = False

#   box_size1 = 100
    box_size1 = MaxBoxSize

    r0 = 2
    r1 = 4
    r2 = 8
    xmin = int(round(j0-box_size0/2))
    ymin = int(round(i0-box_size0/2))
    xmax = xmin+box_size0
    ymax = ymin+box_size0
    fwhms = zeros(stack_size)
    progress.SetLabel('Aligning images (1st star): ')
    for j in arange(0,stack_size):
      if (xbar0[L_idx[j]] != 0):
        if (L_idx[j] != -1):
          fwhms[j] = fwhms[L_idx[j]] 
          xbar0[j] = xbar0[L_idx[j]] 
          ybar0[j] = ybar0[L_idx[j]] 
      else:
        gauge.SetValue(int(100*j/(stack_size-1)))
        wx.Yield()
        image = images[L_idx[j]]
        img = float32(array(image))
        '''
        '''
#
#       Computing the local median corrects errors caused by hot pixels...
        Arect1pile = zeros((box_size0,box_size0,9), dtype=float32)
        for ii in arange(-1,2):
            for jj in arange(-1,2):
                Arect1pile[:,:,3*(ii+1)+(jj+1)] = img[ymin+ii:ymax+ii,xmin+jj:xmax+jj]
        Arect1 = percentile(Arect1pile,50,axis=2)
#
#       Arect1 = img[ymin:ymax,xmin:xmax]        # This is the simple version without the median...
#
        (m,n) = Arect1.shape
        idx = argmax(Arect1[r2:m-r2,r2:n-r2])
        ii = r2 + int(idx/(n-2*r2))
        jj = r2 + idx%(n-2*r2)


        '''
#       Arect1pile = zeros((box_size0,box_size0,9), dtype=float32)
#       for ii in arange(-1,2):
#           for jj in arange(-1,2):
#               Arect1pile[:,:,3*(ii+1)+(jj+1)] = img[ymin+ii:ymax+ii,xmin+jj:xmax+jj]
#       Arect1 = percentile(Arect1pile,50,axis=2)
#       print('Arect1 shape = ', Arect1.shape)
        Arect1 = img[ymin:ymax,xmin:xmax]
        (m,n) = Arect1.shape
        k, LconvL2 = findstars(Arect1,thresh,fignum=0)
        while sum(k) == 0:
            k, LconvL2 = findstars(Arect1,thresh/10,fignum=0)
            if sum(k) == 0:
                k, LconvL2 = findstars(Arect1,thresh/100,fignum=0)
        xs = arange(0,n)
        ys = array([arange(0,m)]).T
        Xs = ones((m,1))*xs
        Ys = ys*ones((1,n))
        flux = LconvL2[k]
        Xs = Xs[k]
        Ys = Ys[k]

        I = argsort(-flux)
        flux = flux[I]
        Xs = int16(Xs[I])
        Ys = int16(Ys[I])
#       ii = Ys[0]
#       jj = Xs[0]
        print('len(Xs) = ', len(Xs))
        for kk in arange(0, len(Xs)):
            dist1 = min([Ys[kk], m-1-Ys[kk], Xs[kk], n-1-Xs[kk]])
            print('dist1 = ', dist1)
            if (dist1 > r2):
                ii = Ys[kk]
                jj = Xs[kk]
                print('kk = ', kk)
                break
        '''

        Arect = Arect1[ii-r2:ii+r2+1,jj-r2:jj+r2+1]
        xs =        arange(-r2,r2+1)
        ys = array([arange(-r2,r2+1)]).T
        radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
        disk = (radii < r0)
        annulus = (radii < r2) & (radii > r1)
        vals0 = Arect[disk]
        vals2 = Arect[annulus]
    
        maxval = np.max(vals0)
        bkg_avg = np.mean(vals2)
        bkg_std = np.std(vals2)
        bkg_median = percentile(vals2,50)
        bkg_iqr = percentile(vals2, 75) - percentile(vals2, 25)
        flux = np.sum(vals0-bkg_avg)
        Arect0 = Arect1[ii-r0:ii+r0+1,jj-r0:jj+r0+1]
        xs0 =        arange(-r0,r0+1)
        ys0 = array([arange(-r0,r0+1)]).T
        xbar0[j] = sum(maximum(Arect0-bkg_median,0)*xs0)/sum(maximum(Arect0-bkg_median,0)) + jj
        ybar0[j] = sum(maximum(Arect0-bkg_median,0)*ys0)/sum(maximum(Arect0-bkg_median,0)) + ii

        if (showAlign):
                figure(200+j, figsize=(12,9))
                imshow(Arect1, cmap='gray', norm=LogNorm())
                plot(xbar0[j], ybar0[j], 'ro')
                if ((FileNames != '') & (j<len(FileNames))): title(FileNames[j])
                draw()
        fwhms[j] = fwhm(Arect1,int(ybar0[j]),int(xbar0[j]))
        xbar0[j] = xmin+xbar0[j]
        ybar0[j] = ymin+ybar0[j]
        fwhms[L_idx[j]] = fwhms[j]
        xbar0[L_idx[j]] = xbar0[j]
        ybar0[L_idx[j]] = ybar0[j]

    if (showAlign):
        show()

#   print(tab, '(i0,j0) = (',i0,',',j0,'),   (ybar0[0],xbar0[0]) = (',ybar0[0],',',xbar0[0],')')
    print(tab, '')
    print(tab, 'FWHMs = ',fwhms)
        
    dx_max = int16(ceil(max(abs(xbar0-xbar0[0]))))
    dy_max = int16(ceil(max(abs(ybar0-ybar0[0]))))
    
    if twostar:
        print(tab, 'doing two star alignment')
        progress.SetLabel('Aligning images (2nd star): ')
        for j in arange(0,stack_size):
          if (xbar1[L_idx[j]] != 0):
            if (L_idx[j] != -1):
              xbar1[j] = xbar1[L_idx[j]] 
              ybar1[j] = ybar1[L_idx[j]] 
          else:
            gauge.SetValue(int(100*j/(stack_size-1)))
            wx.Yield()
            xmin = int(round(xbar0[j]-j0+j1-box_size1/2))
            ymin = int(round(ybar0[j]-i0+i1-box_size1/2))
            xmax = xmin+box_size1
            ymax = ymin+box_size1
            image = images[L_idx[j]]
            img = float32(array(image))
     
            '''
            '''
            # Computing the local median corrects errors caused by hot pixels...
            Arect1pile = zeros((box_size1,box_size1,9), dtype=float32)
            for ii in arange(-1,2):
                for jj in arange(-1,2):
                    Arect1pile[:,:,3*(ii+1)+(jj+1)] = img[ymin+ii:ymax+ii,xmin+jj:xmax+jj]
            Arect1 = percentile(Arect1pile,50,axis=2)
    
            # This is the simple version without the median...
            # Arect1 = img[ymin:ymax,xmin:xmax]
     
            (m,n) = Arect1.shape
            idx = argmax(Arect1[r2:m-r2,r2:n-r2])
            ii = r2 + int(idx/(n-2*r2))
            jj = r2 + idx%(n-2*r2)
            '''
            Arect1 = img[ymin:ymax,xmin:xmax]
            (m,n) = Arect1.shape
            k, LconvL2 = findstars(Arect1,thresh,fignum=0)
            if sum(k) == 0:
                k, LconvL2 = findstars(Arect1,thresh/10,fignum=0)
                if sum(k) == 0:
                    k, LconvL2 = findstars(Arect1,thresh/100,fignum=0)
            xs = arange(0,n)
            ys = array([arange(0,m)]).T
            Xs = ones((m,1))*xs
            Ys = ys*ones((1,n))
            flux = LconvL2[k]
            Xs = Xs[k]
            Ys = Ys[k]

            I = argsort(-flux)
            flux = flux[I]
            Xs = int16(Xs[I])
            Ys = int16(Ys[I])
#           ii = Ys[0]
#           jj = Xs[0]
            print('len(Xs) = ', len(Xs))
            for kk in arange(0, len(Xs)):
                dist1 = min([Ys[kk], m-1-Ys[kk], Xs[kk], n-1-Xs[kk]])
                print('dist1 = ', dist1)
                if (dist1 > r2):
                    ii = Ys[kk]
                    jj = Xs[kk]
                    print('kk = ', kk)
                    break
            '''

            Arect = Arect1[ii-r2:ii+r2+1,jj-r2:jj+r2+1]
            xs =        arange(-r2,r2+1)
            ys = array([arange(-r2,r2+1)]).T
            radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
            disk = (radii < r0)
            annulus = (radii < r2) & (radii > r1)
            vals0 = Arect[disk]
            vals2 = Arect[annulus]
        
            maxval = np.max(vals0)
            bkg_avg = np.mean(vals2)
            bkg_std = np.std(vals2)
            bkg_median = percentile(vals2,50)
            bkg_iqr = percentile(vals2, 75) - percentile(vals2, 25)
            flux = np.sum(vals0-bkg_avg)
            Arect0 = Arect1[ii-r0:ii+r0+1,jj-r0:jj+r0+1]
            xs0 =        arange(-r0,r0+1)
            ys0 = array([arange(-r0,r0+1)]).T
            xbar1[j] = sum(maximum(Arect0-bkg_median,0)*xs0)/sum(maximum(Arect0-bkg_median,0)) + jj
            ybar1[j] = sum(maximum(Arect0-bkg_median,0)*ys0)/sum(maximum(Arect0-bkg_median,0)) + ii
    
            if (showAlign):
                    figure(300+j, figsize=(12,9))
                    imshow(Arect1, cmap='gray', norm=LogNorm())
                    plot(xbar1[j], ybar1[j], 'bo')
                    if ((FileNames != '') & (j<len(FileNames))): title(FileNames[j])
                    draw()

            xbar1[j] = xmin + xbar1[j]
            ybar1[j] = ymin + ybar1[j]
            xbar1[L_idx[j]] = xbar1[j]
            ybar1[L_idx[j]] = ybar1[j]

#   print(tab, xbar0, ybar0)
    if (showAlign):
        rcParams.update({'figure.max_open_warning': 0})
        rcParams['figure.max_open_warning'] = 0
        for k in range(0,N):
            figure(80+k, figsize=(12,9))
            imshow(1000*images[L_idx[k]]/np.max(images[L_idx[k]])+1, cmap='gray', norm=LogNorm())
            plot(xbar0[k],ybar0[k],'r.')
            plot(xbar1[k],ybar1[k],'b.')
            if ((FileNames != '') & (j<len(FileNames))): title(FileNames[k])
            draw()
            show()

        
#   (xbar0, ybar0) are the (x,y) coordinates of the alignment star in each image.  Image 0 is the reference frame
    (N,m,n) = images.shape
    m1 = m + 2*dy_max
    n1 = n + 2*dx_max
    align_stack = zeros((N,m1,n1), dtype=float32)

    alpha = 1
    theta = 0
    alphas = ones((N), dtype=float32)
    thetas = zeros((N), dtype=float32)
    i0 =  dy_max
    i1 =  dy_max+m
    j0 =  dx_max
    j1 =  dx_max+n
    xs = arange(0,n1)
    ys = array([arange(0,m1)]).T
    Xs = ones((m1,1))*xs
    Ys = ys*ones((1,n1))
    success = ones((N), dtype=bool)
    if twostar:
        dist0 = sqrt((xbar1[0]-xbar0[0])**2+(ybar1[0]-ybar0[0])**2)
        theta0 = arctan2(ybar1[0]-ybar0[0],xbar1[0]-xbar0[0])
        print(tab, '')
        print(tab, '  Image               Rescaling    Rotation    Good?')
        print(tab, 'Name/Number            percent       degs')
        print(tab, '-----------           ---------    --------    -----')
        progress.SetLabel('Stacking images: ')
        for k in range(0,N):
            distk = sqrt((xbar1[k]-xbar0[k])**2+(ybar1[k]-ybar0[k])**2)
            alphas[k]= distk/dist0
            thetas[k]= arctan2(ybar1[k]-ybar0[k],xbar1[k]-xbar0[k])-theta0
        alpha_median = percentile(alphas,50)
        alpha_iqr = percentile(alphas, 75) - percentile(alphas, 25)
        kk = 0
        print(tab, 'len of FileNames = ', len(FileNames), ', N = ', N)
#       print(tab, 'FileNames = ', FileNames)
        for k in range(0,N):
            gauge.SetValue(int(50*k/(N-1)))
            wx.Yield()
            if ((FileNames != '') & (k<len(FileNames))): 
                print(tab, '{fname:17s}    {rescale:9.4f}    {rotate:8.3f}     {good}'.format( \
                            fname = FileNames[k], rescale = abs(alphas[k]-1)*100, rotate = thetas[k]*180/pi, \
                            good = (abs(alphas[k] - alpha_median) < 3*alpha_iqr) \
                            ) )
            else:
                print(tab, '      {num:3d}             {rescale:9.4f}    {rotate:8.3f}     {good}'.format( \
                            num = k, rescale = abs(alphas[k]-1)*100, rotate = thetas[k]*180/pi, \
                            good = (abs(alphas[k] - alpha_median) < 3*alpha_iqr) \
                            ) )
            if (abs(alphas[k] - alpha_median) < 3*alpha_iqr):
                build_aligned_stack(k,kk,images)
                kk += 1
            else:
                success[k] = False

        print(tab, '\n')
        print(tab, 'Number stacked: ', kk)
        print(tab, 'Number not well aligned: ', N-kk)
        N = kk
        align_stack = align_stack[:, dy_max:m+dy_max, dx_max:n+dx_max]
    else:
        progress.SetLabel('Stacking images: ')
#       print(tab, 'bob auto     1-star:  i0,i1 = ', i0, i1, ', j0,j1 = ', j0, j1, ', m1,n1 = ', m1,n1)
        for k in range(0,N):
            gauge.SetValue(int(50*k/(N-1)))
            wx.Yield()
            build_aligned_stack(k,k,images)
    tab = tab[0:-4]
    printdone( 'Done autoaligning')

######################################################################################
## Align images on 1 or 2 stars

def onclick1(event):
    global xbar, ybar, xbar0, ybar0, xbar1, ybar1, stack_size, stack_num, star, align_stack, img, fig
    global fwhms

    global tab
    print(tab, 'Clicking on stars:')
    tab = tab+'    '

    A = img
    (m,n) = A.shape
    if (type(event.xdata) != np.float64) | (type(event.ydata) != np.float64): return
    i = int(floor(event.ydata))
    j = int(floor(event.xdata))

    if (min([i,j])<=box_size) | (i>=m-box_size) | (j>=n-box_size): 
        if (i>0): print(tab, 'Not '+str(box_size)+' pixels from the edge.  Try again.')
        return
    r0 = 10
    r1 = 12
    r2 = 18

    Arect = A[i-r2:i+r2+1,j-r2:j+r2+1]
    xs =        arange(-r2,r2+1)
    ys = array([arange(-r2,r2+1)]).T
    radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
    disk = (radii < r0)
    annulus = (radii < r2) & (radii > r1)
    vals0 = Arect[disk]
    vals2 = Arect[annulus]

    maxval = np.max(vals0)
    bkg_avg = np.mean(vals2)
    bkg_std = np.std(vals2)
    bkg_median = percentile(vals2,50)
    bkg_iqr = percentile(vals2, 75) - percentile(vals2, 25)
    flux = np.sum(vals0-bkg_avg)
    Arect0 = A[i-r0:i+r0+1,j-r0:j+r0+1]
    xs0 =        arange(-r0,r0+1)
    ys0 = array([arange(-r0,r0+1)]).T
    xbar0[0] = sum((Arect0-bkg_median)*xs0)/sum(Arect0-bkg_median) + j
    ybar0[0] = sum((Arect0-bkg_median)*ys0)/sum(Arect0-bkg_median) + i

    xmin = int(round(xbar0[0])-box_size/2)
    ymin = int(round(ybar0[0])-box_size/2)
    xmax = xmin+box_size
    ymax = ymin+box_size
    fwhms = zeros(stack_size)
    for j in arange(0,stack_size):
        image = align_stack[j]
        img = float32(array(image))
        Arect1pile = zeros((box_size,box_size,9), dtype=float32)
#
#       Computing the local median corrects errors caused by hot pixels...
        for ii in arange(-1,2):
            for jj in arange(-1,2):
                Arect1pile[:,:,3*(ii+1)+(jj+1)] = img[ymin+ii:ymax+ii,xmin+jj:xmax+jj]
                
        Arect1 = percentile(Arect1pile,50,axis=2)

#
#       This is the simple version without the median...
#       Arect1 = img[ymin:ymax,xmin:xmax]
#
        (m,n) = Arect1.shape

        idx = argmax(Arect1[r2:m-r2,r2:n-r2])
        ii = r2 + int(idx/(n-2*r2))
        jj = r2 + idx%(n-2*r2)
        Arect = Arect1[ii-r2:ii+r2+1,jj-r2:jj+r2+1]
#        idx = argmax(Arect1)
#        ii = idx/n
#        jj = idx%n
#        Arect = Arect1[ii-r2:ii+r2+1,jj-r2:jj+r2+1]
#        print(tab, j, Arect1.shape, Arect.shape, ii, jj, r2)
        xs =        arange(-r2,r2+1)
        ys = array([arange(-r2,r2+1)]).T
        radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
        disk = (radii < r0)
        annulus = (radii < r2) & (radii > r1)
        vals0 = Arect[disk]
        vals2 = Arect[annulus]
    
        maxval = np.max(vals0)
        bkg_avg = np.mean(vals2)
        bkg_std = np.std(vals2)
        bkg_median = percentile(vals2,50)
        bkg_iqr = percentile(vals2, 75) - percentile(vals2, 25)
        flux = np.sum(vals0-bkg_avg)
        Arect0 = Arect1[ii-r0:ii+r0+1,jj-r0:jj+r0+1]
        xs0 =        arange(-r0,r0+1)
        ys0 = array([arange(-r0,r0+1)]).T
        xbar0[j] = sum((Arect0-bkg_median)*xs0)/sum(Arect0-bkg_median) + jj
        ybar0[j] = sum((Arect0-bkg_median)*ys0)/sum(Arect0-bkg_median) + ii
        if (showAlign):
                figure(30+j)
                imshow(Arect1, cmap='gray', norm=LogNorm())
                plot(xbar0[j],ybar0[j],'ro')
                if ((FileNames != '') & (j<len(FileNames))): title(FileNames[j])
                draw()
        fwhms[j] = fwhm(Arect1,int(ybar0[j]),int(xbar0[j]))
        xbar0[j] = xmin+xbar0[j]
        ybar0[j] = ymin+ybar0[j]

    if (showAlign): show()

    '''
    print(tab, 'dx0',xbar0-xbar0[0])
    print(tab, 'dy0',ybar0-ybar0[0])
    '''
    print(tab, '')
    print(tab, 'FWHMs = ',fwhms)
        
    dx_max = int16(ceil(max(abs(xbar0-xbar0[0]))))
    dy_max = int16(ceil(max(abs(ybar0-ybar0[0]))))
        
    (N,m,n) = align_stack.shape
    for j in range(0,N):
        # pad with zeros
        img = zeros(((m+2*dy_max),(n+2*dx_max)), dtype=float32)
        img[ dy_max:dy_max+m, dx_max:dx_max+n ] = align_stack[j]
    
        dx = xbar0[j]-xbar0[0]
        dy = ybar0[j]-ybar0[0]
        dx0 = int16(dx_max+floor(dx))
        dx1 = int16(dx_max+ceil(dx))
        dy0 = int16(dy_max+floor(dy))
        dy1 = int16(dy_max+ceil(dy))
        px = dx - floor(dx)
        py = dy - floor(dy)
        qx = 1.-px
        qy = 1.-py
        align_stack[j] = qy*qx*img[dy0:dy0+m,dx0:dx0+n] + \
                         qy*px*img[dy0:dy0+m,dx1:dx1+n] + \
                         py*qx*img[dy1:dy1+m,dx0:dx0+n] + \
                         py*px*img[dy1:dy1+m,dx1:dx1+n]

#        figure(40+j)
#        imshow(align_stack[j], cmap='gray', norm=LogNorm())
#        draw()
    tab = tab[0:-4]
    printdone( 'Done clicking on stars')

global suffix
suffix = 'RGB'

def stack_images(): 
    global fnames, name, RGB
    global images, file_location, suffix
    global ntot, nstacked
    global progress, gauge
    global success
    global nL, nR, nG, nB, nHa, nO3, nSpec

    global tab
    print(tab, 'Stacking images:')
    tab = tab+'    '

    if (name == ''): name = objectid.GetValue()

    if align_method != 'no stars': images = get_aligned_images()
    (N,m,n) = images.shape
    if align_method == 'no stars': success = ones((N), dtype=bool)

    if   flip == 'Rotate 180':      images = images[:,::-1,::-1] 
    elif flip == 'Mirror':          images = images[:, :  ,::-1] 
    elif flip == 'Flip vertically': images = images[:,::-1, :  ] 

    '''
    true_success = (success == True)
    print(tab, 'success = ', success)
    print(tab, 'true_success = ', true_success)
    images_L    = images[ true_success[0                  :nL                       ], :, :]
    images_R    = images[ true_success[nL                 :nL+nR                    ], :, :]
    images_G    = images[ true_success[nL+nR              :nL+nR+nG                 ], :, :]
    images_B    = images[ true_success[nL+nR+nG           :nL+nR+nG+nB              ], :, :]
    images_Ha   = images[ true_success[nL+nR+nG+nB        :nL+nR+nG+nB+nHa          ], :, :]
    images_O3   = images[ true_success[nL+nR+nG+nB+nHa    :nL+nR+nG+nB+nHa+nO3      ], :, :]
    images_Spec = images[ true_success[nL+nR+nG+nB+nHa+nO3:nL+nR+nG+nB+nHa+nO3+nSpec], :, :]
    '''

    images_L    = images[0                  :nL                       , :, :][success[0                  :nL                       ], :, :]
    images_R    = images[nL                 :nL+nR                    , :, :][success[nL                 :nL+nR                    ], :, :]
    images_G    = images[nL+nR              :nL+nR+nG                 , :, :][success[nL+nR              :nL+nR+nG                 ], :, :]
    images_B    = images[nL+nR+nG           :nL+nR+nG+nB              , :, :][success[nL+nR+nG           :nL+nR+nG+nB              ], :, :]
    images_Ha   = images[nL+nR+nG+nB        :nL+nR+nG+nB+nHa          , :, :][success[nL+nR+nG+nB        :nL+nR+nG+nB+nHa          ], :, :]
    images_O3   = images[nL+nR+nG+nB+nHa    :nL+nR+nG+nB+nHa+nO3      , :, :][success[nL+nR+nG+nB+nHa    :nL+nR+nG+nB+nHa+nO3      ], :, :]
    images_Spec = images[nL+nR+nG+nB+nHa+nO3:nL+nR+nG+nB+nHa+nO3+nSpec, :, :][success[nL+nR+nG+nB+nHa+nO3:nL+nR+nG+nB+nHa+nO3+nSpec], :, :]

    (nL   ,m,n) = images_L.shape
    (nR   ,m,n) = images_R.shape
    (nG   ,m,n) = images_G.shape
    (nB   ,m,n) = images_B.shape
    (nHa  ,m,n) = images_Ha.shape
    (nO3  ,m,n) = images_O3.shape
    (nSpec,m,n) = images_Spec.shape

    txt = 'Stacking images:  Lum = {nL:0d}, '+\
                            'Red = {nR:0d}, '+\
                          'Green = {nG:0d}, '+\
                           'Blue = {nB:0d}, '+\
                         'Halpha = {nHa:0d}, '+\
                           'OIII = {nO3:0d}, '+\
                        'Spectra = {nSpec:0d} '
    print(tab,txt.format(nL=nL, nR=nR, nG=nG, nB=nB, nHa=nHa, nO3=nO3, nSpec=nSpec))
    tab = tab + '    '

    ntot = nL+nR+nG+nB+nHa+nO3+nSpec
    nstacked = 0
    L    = stack(images_L,    outlier_thresh=1.5)     # stack the images
    R    = stack(images_R,    outlier_thresh=1.5)     # stack the images
    G    = stack(images_G,    outlier_thresh=1.5)     # stack the images
    B    = stack(images_B,    outlier_thresh=1.5)     # stack the images
    Ha   = stack(images_Ha,   outlier_thresh=1.5)     # stack the images
    O3   = stack(images_O3,   outlier_thresh=1.5)     # stack the images
    Spec = stack(images_Spec, outlier_thresh=1.5)     # stack the images
    progress.SetLabel('')
    gauge.SetValue(0)

    #####################################################################################
    # Equalize background

    if align_method != 'no stars':
        L20 = 1e+9; R20 = 1e+9; G20 = 1e+9; B20 = 1e+9; Ha20 = 1e+9; OIII20 = 1e+9; Spec20 = 1e+9;
        if nL    != 0: L20    = percentile(   L[L>0],    1)
        if nR    != 0: R20    = percentile(   R[R>0],    1)
        if nG    != 0: G20    = percentile(   G[G>0],    1)
        if nB    != 0: B20    = percentile(   B[B>0],    1)
        if nHa   != 0: Ha20   = percentile(  Ha[Ha>0],   1)
        if nO3   != 0: OIII20 = percentile(  O3[O3>0],   1)
        if nSpec != 0: Spec20 = percentile(Spec[Spec>0], 1)
        bottom20 = min([L20, R20, G20, B20, Ha20, OIII20, Spec20])
#       print(tab, 'LRGBHaOIII20s = ', L20, R20, G20, B20, Ha20, OIII20)
        if nL    != 0: L    = L    - L20 
        if nR    != 0: R    = R    - R20 
        if nG    != 0: G    = G    - G20 
        if nB    != 0: B    = B    - B20 
        if nHa   != 0: Ha   = Ha   - Ha20 
        if nO3   != 0: O3   = O3   - OIII20 
        if nSpec != 0: Spec = Spec - Spec20 

        '''
        if np.sum(L )   != 0: L    = maximum( L    - percentile(   L[L>0],   20), 0 )
        if np.sum(R )   != 0: R    = maximum( R    - percentile(   R[R>0],   20), 0 )
        if np.sum(G )   != 0: G    = maximum( G    - percentile(   G[G>0],   20), 0 )
        if np.sum(B )   != 0: B    = maximum( B    - percentile(   B[B>0],   20), 0 )
        if np.sum(Ha)   != 0: Ha   = maximum( Ha   - percentile(  Ha[Ha>0],  20), 0 )
        if np.sum(O3)   != 0: O3   = maximum( O3   - percentile(  O3[O3>0],  20), 0 )
        if np.sum(Spec) != 0: Spec = maximum( Spec - percentile(Spec[Spec>0],20), 0 )
        '''

        '''
        L    = RemoveGradient(L)
        R    = RemoveGradient(R)
        G    = RemoveGradient(G)
        B    = RemoveGradient(B)
        Ha   = RemoveGradient(Ha)
        O3   = RemoveGradient(O3)
        Spec = RemoveGradient(Spec)
        '''

    L += 0.002*(R+G+B+Ha+O3+Spec)
    
    #####################################################################################
    # Make color image 

    (m,n) = L.shape

    RGB = zeros((m,n,3), dtype=float32)
    if nL+nR+nHa > 0: RGB[0:m, 0:n, 0] = (nL*L+nR*R+nHa*Ha)/(nL+nR+nHa)
    if nL+nG+nO3 > 0: RGB[0:m, 0:n, 1] = (nL*L+nG*G+nO3*O3)/(nL+nG+nO3)
    if nL+nB+nO3 > 0: RGB[0:m, 0:n, 2] = (nL*L+nB*B+nO3*O3)/(nL+nB+nO3)
    MaxVal = ma.max(RGB, axis=2)
    if nL+nR+nHa == 0: RGB[0:m, 0:n, 0] = MaxVal
    if nL+nG+nO3 == 0: RGB[0:m, 0:n, 1] = MaxVal
    if nL+nB+nO3 == 0: RGB[0:m, 0:n, 2] = MaxVal
    '''
    RGB[0:m, 0:n, 0] = L+R+Ha
    RGB[0:m, 0:n, 1] = L+G+O3
    RGB[0:m, 0:n, 2] = L+B+O3
    '''

    r_wt, g_wt, b_wt = color_weights(RGB[:,:,0],RGB[:,:,1],RGB[:,:,2])

    print('color weights = ', r_wt, g_wt, b_wt)
    '''
    r_wt = 1
    g_wt = 1
    b_wt = 1
    '''
    RGB[:,:,0] = r_wt * RGB[:,:,0]
    RGB[:,:,1] = g_wt * RGB[:,:,1]
    RGB[:,:,2] = b_wt * RGB[:,:,2]

    '''
    maxR  = np.max(abs(RGB[:,:,0]))
    maxG  = np.max(abs(RGB[:,:,1]))
    maxB  = np.max(abs(RGB[:,:,2]))
    if   ((maxG==0) & (maxB==0)): 
        RGB[:,:,1] = RGB[:,:,0]
        RGB[:,:,2] = RGB[:,:,0]
    elif ((maxR==0) & (maxB==0)): 
        RGB[:,:,0] = RGB[:,:,1]
        RGB[:,:,2] = RGB[:,:,1]
    elif ((maxR==0) & (maxG==0)): 
        RGB[:,:,0] = RGB[:,:,2]
        RGB[:,:,1] = RGB[:,:,2]
    '''

    #####################################################################################
    # Save color image as a fits file

    RGBfit = zeros((3,m,n), dtype=float32)
    RGBfit[0, 0:m, 0:n] = RGB[0:m, 0:n, 0]
    RGBfit[1, 0:m, 0:n] = RGB[0:m, 0:n, 1]
    RGBfit[2, 0:m, 0:n] = RGB[0:m, 0:n, 2]

    RGBfit += np.max(RGBfit)/20

    hdu = fits.PrimaryHDU(RGBfit)
    hdulist = fits.HDUList([hdu])
    suffix = ''
    if (nL>0):    suffix += 'L'
    if (nR>0):    suffix += 'R'
    if (nG>0):    suffix += 'G'
    if (nB>0):    suffix += 'B'
    if (nHa>0):   suffix += 'Ha'
    if (nO3>0):   suffix += 'OIII'
    if (nSpec>0): suffix += 'Spectra'
    file_location = folder+name+'-py-'+suffix+'.fit'
    hdulist.writeto(file_location, overwrite=True)
#   hdulist.close()
#   fits.PrimaryHDU(RGBfit).writeto(folder+name+'-py-RGB.fit', overwrite=True)
    
    #####################################################################################
    # Plot image w/ sliderbar
    
    if neatimage == 'Yes':
        L = (RGB[:,:,0]+RGB[:,:,1]+RGB[:,:,2])/3.
        bkd_dev = percentile(L,30)-percentile(L,10)
#       RGB = NeatImage(RGB, bkd_dev, 5)
        RGB = NeatImageFFT(RGB, bkd_dev, 5)
#       RGB = cv2.fastNlMeansDenoisingColored(uint8(np.floor(255.99*RGB/np.max(RGB))),None,10,10,7,21)/256
    
    #RGB = minimum( RGB, percentile(RGB, 99.99) )
    #RGB = RGB/np.max(RGB)

#   btn3.Enable()

    tab = tab[0:-4]
    printdone( 'Done stacking images')
    
def show_stacked_image(): 
    global folder, fnames, name, RGB, suffix

    close(1)
#   print(tab,'calling slidershow')
    fig = slidershow(maximum(65535*RGB/np.max(RGB),0), folder=folder, fname=name+'-py-'+suffix+'.png')


global frame

def ChangeCursor(event):
    global fig_frame
    image = wx.Image(script_dir+'png'+slash+'pointer_cyan.png')
    image.SetOption(wx.IMAGE_OPTION_CUR_HOTSPOT_X, 16)
    image.SetOption(wx.IMAGE_OPTION_CUR_HOTSPOT_Y, 16)
    fig_frame.figure_canvas.SetCursor(wx.Cursor(image))  # THE BEST

def UpdateStatusBar(event):
    global fig_frame
    if event.inaxes:
        fig_frame.statusBar.SetStatusText("x={}  y={}".format(event.xdata, event.ydata))

def align(image_stack, method):
    align1star(image_stack)
#   CanvasFrame()

def align1star(image_stack):
    global xbar, ybar, xbar0, ybar0, xbar1, ybar1, stack_size, stack_num, star, \
        images, align_stack, img, cid
    global i0, i1, j0, j1, m1, n1, Xs, Ys, alphas, thetas
    global xpos

    global fig_frame
    '''
    app = wx.App(redirect=False)
    '''

    (N,m,n) = image_stack.shape
    align_stack = zeros((N,m,n), dtype=float32)
    align_stack[:,:,:] = image_stack
    xbar0 = zeros(N)
    ybar0 = zeros(N)
    xbar1 = zeros(N)
    ybar1 = zeros(N)
    star = 0
    stack_num  = 0
    stack_size = N

    m0=m
    n0=n
    mnfactor = np.min([900/m0, 1200/n0, 1])
    m0 = uint16(mnfactor*m0)
    n0 = uint16(mnfactor*n0)

#   print('m0, n0 = ', m0, n0)
#   fig_frame = wx.Frame(None, title='Align, Stack, Show', size=(m0+100,n0))
    fig_frame = wx.Frame(None, title='Align, Stack, Show', size=(xpos,800))

    fig_frame.figure = Figure()
    fig_frame.axes = fig_frame.figure.add_subplot(111,facecolor='black')
    fig_frame.figure.set_facecolor('black')

    img = float32(array(align_stack[stack_num]))
    img0 = minimum(img,  percentile(img,99.9))
    img0 = maximum(img0, percentile(img0, 0.01))
    fig_frame.axes.imshow(img, cmap='gray', norm=LogNorm())
#   fig_frame.axes.imshow(img, cmap='gray', norm=LogNorm())
#   fig_frame.axes.imshow(arcsinh(img), cmap='gray')
    show()

    fig_frame.figure_canvas = FigureCanvas(fig_frame, -1, fig_frame.figure)

    # Note that event is a MplEvent
    fig_frame.figure_canvas.mpl_connect( 'motion_notify_event', UpdateStatusBar)
    fig_frame.figure_canvas.Bind(wx.EVT_ENTER_WINDOW, ChangeCursor)

    fig_frame.axes.spines['bottom'].set_color('white')
    fig_frame.axes.spines['left'].set_color('white')
    fig_frame.axes.spines['top'].set_color('white')
    fig_frame.axes.spines['right'].set_color('white')
    fig_frame.axes.xaxis.label.set_color('white')
    fig_frame.axes.yaxis.label.set_color('white')
    fig_frame.axes.tick_params(axis='x', colors='white')
    fig_frame.axes.tick_params(axis='y', colors='white')

#   print('fnames = ', fnames)
    thename = fnames[stack_num]
    if (len(thename)> 23):  thename = thename[0:10]+'...'+thename[-10:-1]+thename[-1]
    fig_frame.figure.text(0.2, 0.95, thename, ha='center', \
                    va='bottom', size='medium',color='white')
    fig_frame.figure.text(0.43, 0.95, 'Click on', ha='center', \
                    va='bottom', size='medium',color='white')
    fig_frame.figure.text(0.5, 0.95, '1st', ha='center', \
                    va='bottom', size='medium',color='cyan')
    fig_frame.figure.text(0.55, 0.95, 'star', ha='center', \
                    va='bottom', size='medium',color='white')

    fig_frame.sizer = wx.BoxSizer(wx.VERTICAL)
    fig_frame.sizer.Add(fig_frame.figure_canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
#   fig_frame.SetSizer(fig_frame.sizer)
#   fig_frame.Fit()
    fig_frame.SetSizerAndFit(fig_frame.sizer)

    fig_frame.statusBar = wx.StatusBar(fig_frame, -1)
    fig_frame.SetStatusBar(fig_frame.statusBar)

#   fig_frame.toolbar = NavigationToolbar2Wx(fig_frame.figure_canvas)
#   fig_frame.sizer.Add(fig_frame.toolbar, 0, wx.LEFT | wx.EXPAND)
#   fig_frame.toolbar.Show()

    print('Click on a star at least '+str(box_size)+' pixels from the edge')
    cid = fig_frame.figure.canvas.mpl_connect('button_release_event', onclick)
    fig_frame.Show()

global stacking
stacking = True
global im
def onclick(event):
    global xbar, ybar, xbar0, ybar0, xbar1, ybar1, stack_size, stack_num, star, \
        align_stack, img, cid, btn12, btn13, btn14, btn16, btn20, btn21, lbl7, images
    global fig_frame
    global stacking
    global L_idx
    global i0, i1, j0, j1, m1, n1, Xs, Ys, alphas, thetas
    global tab
    global im
    print(tab,'Click')
    tab = tab + '    '

    [m, n] = img.shape
    ev_x = uint16(event.xdata)
    ev_y = uint16(event.ydata)
    xs = arange(-15,16) + ev_x
    ys = array([arange(-15,16)]).T + ev_y
    Xs = ones((31,1))*xs
    Ys = ys*ones((1,31))
    img0 = img[ ev_y-15:ev_y+16, ev_x-15:ev_x+16 ]

    print(tab,'  1: stack_num = ', stack_num)
    print(tab,'  xbar0', xbar0)
    
    [j, i] = unravel_index(img0.argmax(), img0.shape)
    if (stack_num < stack_size):
        xbar0[      stack_num]  = Xs[j,i]
        ybar0[      stack_num]  = Ys[j,i]
        xbar0[L_idx[stack_num]] = Xs[j,i]
        ybar0[L_idx[stack_num]] = Ys[j,i]
    elif (stack_num < 2*stack_size):
        xbar1[      stack_num-stack_size]  = Xs[j,i]
        ybar1[      stack_num-stack_size]  = Ys[j,i]
        xbar1[L_idx[stack_num-stack_size]] = Xs[j,i]
        ybar1[L_idx[stack_num-stack_size]] = Ys[j,i]





    done = False
    stack_num += 1
    if (stack_num < stack_size):
        print('    2: stack_num = ', stack_num)
        while ((xbar0[L_idx[stack_num]] != 0) & (stack_num < stack_size-1)): 
            if (L_idx[stack_num] != -1):
                xbar0[stack_num] = xbar0[L_idx[stack_num]]
                ybar0[stack_num] = ybar0[L_idx[stack_num]]
                print('  xbar0', xbar0)
            stack_num += 1
            print('      3: stack_num = ', stack_num)
    if (stack_num < stack_size):
            img = float32(array(align_stack[L_idx[stack_num]]))
            img0 = minimum(img,  percentile(img,99.9))
            img0 = maximum(img0, percentile(img0, 0.01))
    
            fig_frame.figure.set_facecolor('black')
            if (stack_num == 1):
                im = fig_frame.axes.imshow(img0, cmap='gray', norm=LogNorm())
#               fig_frame.axes.imshow(img, cmap='gray', norm=LogNorm())
#               fig_frame.axes.imshow(arcsinh(img), cmap='gray')
                show()
    
#           print('stack_num = ', stack_num, ', stack_size = ', stack_size, ', len(fnames) = ', len(fnames))
#           print('(stack_num)%stack_size = ', (stack_num)%stack_size)
    
            thename = fnames[(stack_num-1)%len(fnames)]
            if (len(thename)> 23):  thename = thename[0:10]+'...'+thename[-10:-1]+thename[-1]
            fig_frame.figure.text(0.2, 0.95, thename, ha='center', va='bottom', size='medium',color='black')
            thename = fnames[(stack_num)%len(fnames)]
            if (len(thename)> 23):  thename = thename[0:10]+'...'+thename[-10:-1]+thename[-1]
            fig_frame.figure.text(0.2, 0.95, thename, ha='center', va='bottom', size='medium',color='white')
            fig_frame.figure.text(0.43, 0.95, 'Click on', ha='center', \
                            va='bottom', size='medium',color='white')
            fig_frame.figure.text(0.5, 0.95, '1st', ha='center', \
                            va='bottom', size='medium',color='cyan')
            fig_frame.figure.text(0.55, 0.95, 'star', ha='center', \
                            va='bottom', size='medium',color='white')
            fig_frame.Show()
            if (stack_num == 1):
                fig_frame.figure.canvas.draw()
            else:
                im.set_array(img0)
                fig_frame.figure.canvas.draw()
            image = wx.Image(script_dir+'png'+slash+'pointer_cyan.png')
            image.SetOption(wx.IMAGE_OPTION_CUR_HOTSPOT_X, 16)
            image.SetOption(wx.IMAGE_OPTION_CUR_HOTSPOT_Y, 16)
            fig_frame.figure_canvas.SetCursor(wx.Cursor(image))  # THE BEST
    elif align_method == 'click on 1 star':
            done = True
    elif ((align_method == 'click on 2 stars') & (stack_num < 2*stack_size)):
        print('        4: stack_num = ', stack_num)
        while ((xbar1[L_idx[stack_num-stack_size]] != 0) & (stack_num < 2*stack_size-1)): 
                xbar1[stack_num-stack_size] = xbar1[L_idx[stack_num-stack_size]]
                ybar1[stack_num-stack_size] = ybar1[L_idx[stack_num-stack_size]]
#               print('          5: stack_num = ', stack_num)
                stack_num += 1
        print('        5: stack_num = ', stack_num, '2*stack_size = ', 2*stack_size)
        if (stack_num < 2*stack_size):
            img = float32(array(align_stack[L_idx[stack_num-stack_size]]))
            img0 = minimum(img,  percentile(img,99.9))
            img0 = maximum(img0, percentile(img0, 0.01))

#           fig_frame.figure.set_facecolor('black')
#           fig_frame.axes.imshow(img0, cmap='gray', norm=LogNorm())
#           fig_frame.axes.imshow(img, cmap='gray', norm=LogNorm())
#           fig_frame.axes.imshow(arcsinh(img), cmap='gray')
#           show()

#       print('stack_num = ', stack_num, ', fname = ', fnames[stack_num-stack_size%len(fnames)])
            thename = fnames[(stack_num-stack_size-1)%len(fnames)]
            if (len(thename)> 23):  thename = thename[0:10]+'...'+thename[-10:-1]+thename[-1]
            fig_frame.figure.text(0.2, 0.95, thename, ha='center', va='bottom', size='medium',color='black')
            thename = fnames[(stack_num-stack_size)%len(fnames)]
            if (len(thename)> 23):  thename = thename[0:10]+'...'+thename[-10:-1]+thename[-1]
            fig_frame.figure.text(0.2, 0.95, thename, ha='center', va='bottom', size='medium',color='white')
            fig_frame.figure.text(0.43, 0.95, 'Click on', ha='center', \
                            va='bottom', size='medium',color='white')
            fig_frame.figure.text(0.5, 0.95, '1st', ha='center', \
                            va='bottom', size='medium',color='black')
            fig_frame.figure.text(0.5, 0.95, '2nd', ha='center', \
                            va='bottom', size='medium',color='red')
            fig_frame.figure.text(0.55, 0.95, 'star', ha='center', \
                            va='bottom', size='medium',color='white')
            fig_frame.Show()
            im.set_array(img0)
            fig_frame.figure.canvas.draw()
            image = wx.Image(script_dir+'png'+slash+'pointer_red.png')
            image.SetOption(wx.IMAGE_OPTION_CUR_HOTSPOT_X, 16)
            image.SetOption(wx.IMAGE_OPTION_CUR_HOTSPOT_Y, 16)
            fig_frame.figure_canvas.SetCursor(wx.Cursor(image))  # THE BEST
        else:
            done = True
    else:
        done = True
    print(tab,'done = ', done)
    if (done & stacking):
#   elif (stacking):

        print('xbar0 = ', xbar0)
        print('ybar0 = ', ybar0)
        if (align_method == 'click on 2 stars'):
            print('xbar1 = ', xbar1)
            print('ybar1 = ', ybar1)

        image = wx.Image(script_dir+'png'+slash+'pointer.png')
        image.SetOption(wx.IMAGE_OPTION_CUR_HOTSPOT_X, 16)
        image.SetOption(wx.IMAGE_OPTION_CUR_HOTSPOT_Y, 16)
        fig_frame.figure_canvas.SetCursor(wx.Cursor(image))  # THE BEST
        show()


        stacking = False
        (N,m,n) = align_stack.shape
        xmin = int( -ceil(-np.min(xbar0-xbar0[0])))-1
        xmax = int(  ceil( np.max(xbar0-xbar0[0])))+1
        ymin = int( -ceil(-np.min(ybar0-ybar0[0])))-1
        ymax = int(  ceil( np.max(ybar0-ybar0[0])))+1
#       print('xmin, xmax = ', xmin, xmax)
#       print('ymin, ymax = ', ymin, ymax)
        m1 = m+ymax-ymin
        n1 = n+xmax-xmin
        align_stack = zeros((N,m1,n1), dtype=float32)
        success = ones((N), dtype=bool)
        alpha = 1
        theta = 0
        alphas = ones((N), dtype=float32)
        thetas = zeros((N), dtype=float32)
        i0 =  ymax
        i1 =  ymax+m
        j0 =  xmax
        j1 =  xmax+n
        xs = arange(0,n-xmin+xmax)
        ys = array([arange(0,m-ymin+ymax)]).T
        Xs = ones((m1,1))*xs
        Ys = ys*ones((1,n1))


        if (align_method == 'click on 2 stars'):
            dist0 = sqrt((xbar1[0]-xbar0[0])**2+(ybar1[0]-ybar0[0])**2)
            theta0 = arctan2(ybar1[0]-ybar0[0],xbar1[0]-xbar0[0])
            print('')
            print('  Image               Rescaling    Rotation')
            print('Name/Number            percent       degs')
            print('-----------           ---------    --------')
            for k in range(0,N):
                distk = sqrt((xbar1[k]-xbar0[k])**2+(ybar1[k]-ybar0[k])**2)
                alphas[k] = distk/dist0
                thetas[k] = arctan2(ybar1[k]-ybar0[k],xbar1[k]-xbar0[k])-theta0
#               print('alpha, theta: ',alpha, theta)
#               alpha = 1.0
#               theta = 0.0
                if ((FileNames != '') & (k<len(FileNames))): 
                    print('{fname:17s}    {rescale:9.4f}    {rotate:8.3f}'.format( \
                                fname = FileNames[k], rescale = abs(alphas[k]-1)*100, rotate = thetas[k]*180/pi) )
                else:
                    print('      {num:3d}             {rescale:9.4f}    {rotate:8.3f}'.format( \
                                num = k, rescale = abs(alphas[k]-1)*100, rotate = thetas[k]*180/pi) )
                build_aligned_stack(k,k,images)
            print('\n')
        else:
#           print('bob click on 1-star:  i0,i1 = ', i0, i1, ', j0,j1 = ', j0, j1, ', m1,n1 = ', m1,n1)
            for k in range(0,N):
                build_aligned_stack(k,k,images)
        stack_images()
        show_stacked_image()
        show_fwhm(RGB[:,:,1])

        btn16.Enable()
        btn12.Enable()
        btn13.Enable()
        btn14.Enable()
        btn20.Enable()
#       btn21.Enable()
        lbl7.SetLabel( folder+name+'-py-RGB.fit')
        fig_frame.figure_canvas.mpl_disconnect(cid)
        fig_frame.Destroy()

    tab = tab[0:-4]
    return
    
def get_aligned_images():
    global align_stack, success, FileNames, showAlign
    if (showAlign):
        for j in arange(len(FileNames)):
                img = float32(array(align_stack[j]))
                fig = slidershow(maximum(65535*img/np.max(img),0), folder=folder, fname=name+'-py-'+suffix+'.png')
                if ((FileNames != '') & (j<len(FileNames))): title(FileNames[j])
                draw()
                show()
    return align_stack

def get_fwhms():
    global fwhms
    return fwhms

global i_prev, j_prev

def onmouse(event):
    global Aglobal, xbar, ybar, lbl9, lbl11
    global i_prev, j_prev
    global current_fwhm
    A = Aglobal
    if (len(A.shape) == 2):
        (m,n) = A.shape
        if (type(event.xdata) == np.float64) & (type(event.ydata) == np.float64):
            i = int(floor(event.ydata))
            j = int(floor(event.xdata))
    
            '''
            r0 = 3
            r1 = 5
            r2 = 8
            '''
    
            if (current_fwhm > 0):
                r0 = int(max([round(current_fwhm),2]))
                r1 = r0+2
                r2 = int(1.5*r1)
            else:
                r0 = 10
                r1 = 12
                r2 = 18

            if (min([i,j])>r2) & (i<m-r2) & (j<n-r2):
                Arect = A[i-r2:i+r2+1,j-r2:j+r2+1]
                xs =        arange(-r2,r2+1)
                ys = array([arange(-r2,r2+1)]).T
                radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
                disk = (radii < r0)
                annulus = (radii < r2) & (radii > r1)
                vals0 = Arect[disk]
                vals2 = Arect[annulus]
    
                maxval = np.max(vals0)
                bkg_avg = np.mean(vals2)
                bkg_std = np.std(vals2)
                bkg_median = percentile(vals2,50)
                bkg_iqr = percentile(vals2, 75) - percentile(vals2, 25)
                flux = np.sum(vals0-bkg_avg)
                Arect0 = A[i-r0:i+r0+1,j-r0:j+r0+1]
                xs0 =        arange(-r0,r0+1)
                ys0 = array([arange(-r0,r0+1)]).T
                xbar = sum((Arect0-bkg_median)*xs0)/sum(Arect0-bkg_median+2.7e-18) + j
                ybar = sum((Arect0-bkg_median)*ys0)/sum(Arect0-bkg_median+2.7e-18) + i
                if ((maxval-bkg_median>5*bkg_iqr) & \
                                (abs(i-ybar)<0.6) & \
                                (abs(j-xbar)<0.6) & \
                                ((abs(j-j_prev)>2)|(abs(i-i_prev)>2)) ):
#                   print(' pixel=(%4d,%4d), pix_value=%7.1f' %
#                           (i,j, A[i,j]))
                    print(' centroid: (x,y)=(%7.2f,%7.2f), max_val=%7.1f, flux=%7.1f' %
                            (xbar,ybar, maxval, flux))
                    print(' bkg_avg   =%7.1f, bkg_std=%7.1f ' %
                            (bkg_avg,bkg_std))
                    print(' bkg_median=%7.1f, bkg_iqr=%7.1f ' %
                            (bkg_median,bkg_iqr))
                    print(' A.shape = ', A.shape,', Arect.shape = ',Arect.shape)
                    print(' fwhm = %7.3f \n' %
                            fwhm(A,i,j))
                    if (lbl9 != ''):
                        lbl9.SetLabel('')

                        txt9 = ' Star info:  Centroid: ({x:<7.2f}, {y:<7.2f})  Max Val: {max_val:<7.0f}  Flux: {flux:<7.0f}'
                        txt9 = txt9+ \
                             '\n            Background:  Avg.: {avg:<7.1f}   Std.Dev.: {stddev:<7.1f} '
                        txt9 = txt9+ \
                             '\n                              Median: {median:<7.1f}            IQR: {iqr:<7.1f}'

                        lbl9.SetLabel( txt9.format(x=xbar,y=ybar,max_val=maxval,flux=flux, \
                                                avg=bkg_avg,stddev=bkg_std,median=bkg_median,iqr=bkg_iqr))
                    i_prev = i
                    j_prev = j
    elif (len(A.shape) == 3):
        (m,n,three) = A.shape
        Amax = np.max(A)
        if (type(event.xdata) == np.float64) & (type(event.ydata) == np.float64):
            i = int(floor(event.ydata))
            j = int(floor(event.xdata))
    
            '''
            r0 = 3
            r1 = 5
            r2 = 8
            '''
    
            if (current_fwhm > 0):
                r0 = int(max([round(current_fwhm),2]))
                r1 = r0+2
                r2 = int(1.5*r1)
            else:
                r0 = 10
                r1 = 12
                r2 = 18

            if (min([i,j])>r2) & (i<m-r2) & (j<n-r2):
                R = A[i-r2:i+r2+1,j-r2:j+r2+1,0]
                G = A[i-r2:i+r2+1,j-r2:j+r2+1,1]
                B = A[i-r2:i+r2+1,j-r2:j+r2+1,2]
                L = R+G+B
                xs =        arange(-r2,r2+1)
                ys = array([arange(-r2,r2+1)]).T
                radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
                disk = (radii < r0)
                annulus = (radii < r2) & (radii > r1)
                Rvals0 = R[disk]
                Rvals2 = R[annulus]
                Gvals0 = G[disk]
                Gvals2 = G[annulus]
                Bvals0 = B[disk]
                Bvals2 = B[annulus]
                Lvals0 = L[disk]
                Lvals2 = L[annulus]
    
                maxval = np.max(Rvals0)
                bkg_avg = np.mean(Rvals2)
                bkg_std = np.std(Rvals2)
                bkg_median = percentile(Rvals2,20)
                bkg_iqr = percentile(Rvals2, 75) - percentile(Rvals2, 25)
                Rflux = np.sum(Rvals0-bkg_avg)
                maxval = np.max(Gvals0)
                bkg_avg = np.mean(Gvals2)
                bkg_std = np.std(Gvals2)
                bkg_median = percentile(Gvals2,20)
                bkg_iqr = percentile(Gvals2, 75) - percentile(Gvals2, 25)
                Gflux = np.sum(Gvals0-bkg_avg)
                maxval = np.max(Bvals0)
                bkg_avg = np.mean(Bvals2)
                bkg_std = np.std(Bvals2)
                bkg_median = percentile(Bvals2,20)
                bkg_iqr = percentile(Bvals2, 75) - percentile(Bvals2, 25)
                Bflux = np.sum(Bvals0-bkg_avg)
                maxval = np.max(Lvals0)
                bkg_avg = np.mean(Lvals2)
                bkg_std = np.std(Lvals2)
                bkg_median = percentile(Lvals2,20)
                bkg_iqr = percentile(Lvals2, 75) - percentile(Lvals2, 25)
                Lflux = np.sum(Lvals0-bkg_avg)
                Arect0 = A[i-r0:i+r0+1,j-r0:j+r0+1,0] \
                       + A[i-r0:i+r0+1,j-r0:j+r0+1,1] \
                       + A[i-r0:i+r0+1,j-r0:j+r0+1,2]
                xs0 =        arange(-r0,r0+1)
                ys0 = array([arange(-r0,r0+1)]).T
                denom = sum(maximum(Arect0-bkg_median,1e-5))
                xbar = sum((Arect0-bkg_median)*xs0)/denom + j
                ybar = sum((Arect0-bkg_median)*ys0)/denom + i
                if ((maxval-bkg_median>5*bkg_std) & \
                                (abs(i-ybar)<0.6) & \
                                (abs(j-xbar)<0.6) & \
                                ((abs(j-j_prev)>2)|(abs(i-i_prev)>2)) & \
                                (sum(A[i,j,:])>0.001) ):
                    fwhm_here = fwhm(A[:,:,0]+A[:,:,1]+A[:,:,2],i,j)
                    if (fwhm_here > 0.001):
                        print(' pixel=(%4d,%4d), centroid=(%7.2f,%7.2f), max_val=%7.1f' %
                                (i,j, ybar, xbar, maxval))
                        print(' red_value=%7.1f, green_value=%7.1f, blue_value=%7.1f' %
                                (A[i,j,0], A[i,j,1], A[i,j,2]))
                        print(' red_flux =%7.1f, green_flux =%7.1f, blue_flux =%7.1f, total_flux = %7.1f' %
                                (Rflux, Gflux, Bflux, Rflux+Gflux+Bflux))
                        print(' bkg_avg   =%7.1f, bkg_std=%7.1f ' %
                                (bkg_avg,bkg_std))
                        print(' bkg_median=%7.1f, bkg_iqr=%7.1f \n' %
                                (bkg_median,bkg_iqr))
                        print(' fwhm = %7.3f \n' %
                                (fwhm(A[:,:,0]+A[:,:,1]+A[:,:,2],i,j)))
                        i_prev = i
                        j_prev = j

                        if (lbl11 != ''): lbl11.SetLabel('') 
                        txt11 = ' Star info:  Centroid: ({xbar:<7.2f}, {ybar:<7.2f})  ' \
                               + 'FWHM: {fwhm:<5.2} ' \
                               + 'Fluxes: R={red_flux:<7.2e}, G={green_flux:<7.2e}, B={blue_flux:<7.2e}, Total={total_flux:<7.2e} '
                        lbl11.SetLabel( txt11.format(xbar=xbar,ybar=ybar,fwhm=fwhm_here,red_flux=Rflux,green_flux=Gflux,blue_flux=Bflux,total_flux=Lflux ))
                else:
                    if (lbl9 != ''): lbl9.SetLabel('') 
                    txt9 = ' Pixel info:  (x,y) = ({x:<5.0f}, {y:<5.0f}),  Red: {red_val:<7.5f}   Green: {green_val:<7.5f}   Blue: {blue_val:<7.5f} '
                    lbl9.SetLabel( txt9.format(x=j,y=i,red_val=A[i,j,0]/Amax,green_val=A[i,j,1]/Amax,blue_val=A[i,j,2]/Amax \
                                                ))

global scolsat, sgamma, smini, smaxi, A, B, C, h, maxi, mini, gamma, colsat, axsave, bsave, axname, bname
def slidershow(Img, fignum=1, img_out=[], folder='', fname=''):
    global scolsat, sgamma, smini, smaxi, A, B, C, h, maxi, mini, gamma, colsat
    global axsave, bsave, axname, bname, filename, foldername
    global i_prev, j_prev
    global xpos

    global Aglobal
    global fig

    m = Img.shape[0]
    n = Img.shape[1]
    '''
    print('m,n: ',m,n)
    '''
#   app = wx.App(redirect=False)
    width, height = wx.GetDisplaySize()
    height = int(0.9*height)
    width = int(min([width, n*height/m]))
    dpi = 100

#   A = Img
    A = Img + 0.001*(np.max(Img)-np.min(Img))  # Add a pedestal
    filename = fname
    foldername = folder

    from matplotlib.widgets import Slider
    rcParams.update({'figure.max_open_warning': 0})
    rcParams['figure.max_open_warning'] = 0

    Aglobal = A
#   ion()
    isinteractive()
    fig =figure(fignum, figsize=(width/dpi,height/dpi), dpi=dpi)

    mngr = get_current_fig_manager()
    mngr.window.setGeometry(xpos,26,width,height)

    fig.patch.set_facecolor('gray')

#   fig.canvas.manager.window.move(400,0)    # Doesn't work as hoped

    i_prev = -1
    j_prev = -1
    fig.canvas.mpl_connect('motion_notify_event', onmouse)
    fig.canvas.mpl_connect('button_press_event', onmouse)

#   ion()
    isinteractive()
    subplots_adjust(bottom=0.14)
    subplots_adjust(top=0.95)
    subplots_adjust(left=0.01)
    subplots_adjust(right=0.99)
    subplots_adjust(hspace=0.01)
    
    B = A/np.max(A)
    gamma = 1
    if (len(A.shape) == 2):
        (m,n) = A.shape
        med = percentile(A,50)
        pct_bottom = np.min(A)
        pct_lo = 0.95*percentile(A,10) + 0.05*pct_bottom
        pct_hi = percentile(A,99.9)
#       pct_top = percentile(A,99.9999)
        pct_top = percentile(A,100)

        h = imshow(A, interpolation='nearest', \
                    cmap='gray',
                    origin='upper', \
                    aspect='equal', extent=None)
        h.set_clim(vmin=pct_lo, vmax=pct_hi)
    else:
        (m,n,N) = A.shape

        '''
        Rdim = percentile(B[:,:,0],1)
        Gdim = percentile(B[:,:,1],1)
        Bdim = percentile(B[:,:,2],1)
        Ldim = minimum(Rdim,minimum(Gdim,Bdim))
        Rbright = percentile(B[:,:,0]-Rdim,99.99)
        Gbright = percentile(B[:,:,1]-Gdim,99.99)
        Bbright = percentile(B[:,:,2]-Bdim,99.99)
        B[:,:,0] = (Gbright/Rbright)*(B[:,:,0] - Rdim) + Ldim
        B[:,:,1] =                   (B[:,:,1] - Gdim) + Ldim
        B[:,:,2] = (Gbright/Bbright)*(B[:,:,2] - Bdim) + Ldim
        '''

        Rdim = percentile(B[:,:,0],0.01)
        Gdim = percentile(B[:,:,1],0.01)
        Bdim = percentile(B[:,:,2],0.01)
        Ldim = minimum(Rdim,minimum(Gdim,Bdim))
        B[:,:,0] = B[:,:,0] - Rdim
        B[:,:,1] = B[:,:,1] - Gdim
        B[:,:,2] = B[:,:,2] - Bdim
#       r_wt, g_wt, b_wt = color_weights(B[:,:,0],B[:,:,1],B[:,:,2])
        r_wt = 1; g_wt = 1; b_wt = 1;
        B[:,:,0] = r_wt * B[:,:,0] + Ldim
        B[:,:,1] = g_wt * B[:,:,1] + Ldim
        B[:,:,2] = b_wt * B[:,:,2] + Ldim

        med = percentile(B,50)
#       pct_lo = 0.95*percentile(B,10) + 0.05*np.min(B)
        pct_lo = percentile(B[B>0],0.1)
        pct_hi = percentile(B,99.9)
#       pct_top = percentile(B,99.9999)
        pct_top = 1.0

        h = imshow(B, interpolation='nearest', \
                    # cmap='gray',        ###  ???
                    origin='upper', \
                    aspect='equal', extent=None)
#       h.set_data(minimum(maximum((B-pct_lo+0.2*(pct_hi-pct_lo))/(pct_hi-pct_lo+0.2*(pct_hi-pct_lo)),0),1))
        h.set_data(minimum(maximum((B-pct_lo)/(pct_hi-pct_lo),0),1))

    gca().set_xticks([])
    gca().set_yticks([])
#   print('gama = ', gamma, ', pcts = ',pct_lo,pct_hi,pct_top)

    class Index(object):

        def saveit(self, event):
            global foldername, filename
            colsat = scolsat.val
#           gamma = np.power(2,sgamma.val)
#           gamma = np.power(2,sgamma.val) - 1
            gamma = np.power(3,sgamma.val) - 2
            mini = smini.val**2
            maxi = smaxi.val**2

            if (img_out != []):
                C = img_out
            else:
                C = array(zeros(B.shape, dtype=float32))

            if len(A.shape) == 2:
                if mini<=maxi:
                        C[:,:] = log(   1+minimum(maximum((A-mini)/(maxi-mini),0),1)*gamma)
                else:
                        C[:,:] = log(  2.-minimum(maximum((A-maxi)/(mini-maxi),0),1)*gamma)
            else:
                if mini<=maxi:
                        C = arcsinh(minimum(maximum((B-mini)/(maxi-mini),0),1)*gamma)/arcsinh(gamma)
                        L = (C[:,:,0]+C[:,:,1]+C[:,:,2])/3
                        CL = array(zeros(C.shape, dtype=float32))
                        CL[:,:,0] = C[:,:,0]-L
                        CL[:,:,1] = C[:,:,1]-L
                        CL[:,:,2] = C[:,:,2]-L
                        maxL = np.max(L)
                        C[:,:,0] = L + colsat*(L/maxL)*CL[:,:,0]
                        C[:,:,1] = L + colsat*(L/maxL)*CL[:,:,1]
                        C[:,:,2] = L + colsat*(L/maxL)*CL[:,:,2]
                        C = minimum(maximum(C,0),1)
                else:
                        C[:,:,:] = log( 2.-minimum(maximum((B-maxi)/(mini-maxi),0),1)*gamma)/log(1+gamma)
            if filename != '': 
                if not os.path.isdir( foldername ):
                    if (foldername != ''):
                        print('Creating the folder: ', foldername)
                        os.mkdir( foldername )  # make sure the directory exists
                im = cv2.cvtColor((floor(65535*C)).astype(uint16), cv2.COLOR_RGB2BGR)
                cv2.imwrite(foldername+filename, im)
                print(tab, 'saving it to: ', foldername)
                '''
                figure(100)
                imshow(C)
                draw()
                show()
                '''

        def nameit(self, text):
            global filename
            filename = text
            
    #------ defining sliders ------

    if len(A.shape) == 3:
        axcolsat = axes([0.14, 0.020, 0.52, 0.015], facecolor='darkgray')
    axgamma = axes([0.14, 0.045, 0.52, 0.015], facecolor='darkgray')
    axmini = axes([0.14, 0.070, 0.52, 0.015], facecolor='darkgray')
    axmaxi = axes([0.14, 0.095, 0.52, 0.015], facecolor='darkgray')
    callback = Index()
    axsave = axes([0.94, 0.015, 0.04, 0.04])
    bsave  = wid.Button(axsave, 'Save', color='darkgray')
    bsave.on_clicked(callback.saveit)
    axname = axes([0.78, 0.077, 0.20, 0.025], facecolor='darkgray')
    bname  = TextBox(axname, 'File name: ', initial=fname, color='darkgray')
    bname.on_submit(callback.nameit)
#   bname._rendercursor()
    
    if len(A.shape) == 3:
        scolsat = Slider(axcolsat, 'color saturation',-1., 10, valinit=5)
    sgamma      = Slider(axgamma,  'nonlinear stretch',0., 10, valinit=2)
    smini       = Slider(axmini,   'min',              0.8*sqrt(pct_lo), sqrt(pct_top), valinit=np.power(maximum(pct_lo,0),0.5))
    smaxi       = Slider(axmaxi,   'max',              0.8*sqrt(pct_lo), sqrt(pct_top), valinit=np.power(pct_hi,0.5))
    
    if len(A.shape) == 3:
        scolsat.poly.set_color('dimgray')
    sgamma.poly.set_facecolor('dimgray')
    smini.poly.set_color('dimgray')
    smaxi.poly.set_color('dimgray')
    pause(0.3)
    if len(A.shape) == 3:
        scolsat.on_changed(update)
    sgamma.on_changed(update)
    smini.on_changed(update)
    smaxi.on_changed(update)
    if (slash=='\\'):
        show()   # Without show(), the app crashes on Windows when clicking on a figure window
#                  With show(), the app has to be closed by hitting Ctrl-C in
#                  the Terminal window
#   print('if too bright, comment out update(0) and uncomment out draw()')
#   draw()
    update(0)   # the zero means nothing
    return fig

def HRdiag_show(fignum=1, img_out=[], folder='', fname=''):
    global plotter
    global xpos

    filename = fname
    foldername = folder

    from matplotlib.widgets import Slider
    rcParams.update({'figure.max_open_warning': 0})
    rcParams['figure.max_open_warning'] = 0

    isinteractive()
    if (plotter == ''): 
        images_frame = wx.Frame(None, -1, 'Images', size=[1000,880])
        images_frame.SetPosition(wx.Point(xpos, 26))
        plotter = PlotNotebook(images_frame)
        images_frame.Show(True)
    fig = plotter.add('HR Diagram')
    fig.patch.set_facecolor('gray')

    fig.canvas.mpl_connect('motion_notify_event', onmouse)
    fig.canvas.mpl_connect('button_press_event', onmouse)

    isinteractive()
    subplots_adjust(bottom=0.11)
    subplots_adjust(top=0.96)
    subplots_adjust(left=0.05)
    subplots_adjust(right=0.99)
    subplots_adjust(hspace=0.01)
    
    gca().set_xticks([])
    gca().set_yticks([])

    '''
    class Index(object):

        def saveit(self, event):
            global foldername, filename

            print('clicked on saveit')
            if (img_out != []):
                C = img_out
            else:
                C = array(zeros(B.shape, dtype=float32))
            if filename != '': 
                print(tab, 'saving it to: ', foldername)
                im = cv2.cvtColor((floor(65535*C)).astype(uint16), cv2.COLOR_RGB2BGR)
                cv2.imwrite(foldername+filename, im)

        def nameit(self, text):
            global filename
            filename = text
    '''
            

    '''
    #------ defining sliders ------
    callback = Index()

    axsave = axes([0.94, 0.015, 0.04, 0.025])
    bsave  = wid.Button(axsave, 'Save', color='darkgray')
    bsave.on_clicked(callback.saveit)

    axname = axes([0.58, 0.015, 0.32, 0.025], facecolor='darkgray')
    bname  = TextBox(axname, 'File name: ', initial=fname, color='darkgray')
    bname.on_submit(callback.nameit)
    '''
    
    pause(0.3)
    print('if too bright, comment out update(0) and uncomment out draw()')
#   draw()
    update(0)
    return fig

def update(val):
    global scolsat, sgamma, smini, smaxi, A, B, h, maxi, mini, gamma, colsat
    if len(A.shape) == 3:
        colsat = scolsat.val
    gamma = sgamma.val
    mini = smini.val
    maxi = smaxi.val

    '''
    '''
#   these older versions often makes things too bright
#   gamma = np.power(2,sgamma.val) - 1     # this is good
    gamma = np.power(3,sgamma.val) - 2     # this is better (aka stronger)
    mini = smini.val**2
    maxi = smaxi.val**2
    '''
    '''

    if len(A.shape) == 2:
        if mini<=maxi:
                h.set_data(log(1+minimum(maximum((A-mini)/(maxi-mini),0),1)*gamma)/log(1+gamma))
                h.set_clim(vmin=0, vmax=1)
        else:
                h.set_data(np.power(1.-minimum(maximum((A-maxi)/(mini-maxi),0),1),gamma))
                h.set_clim(vmin=0, vmax=1)
    else:
        if mini<=maxi:
                D = arcsinh(minimum(maximum((B-mini)/(maxi-mini),0),1)*gamma)/arcsinh(gamma)
                L = (D[:,:,0]+D[:,:,1]+D[:,:,2])/3
                DL = array(zeros(D.shape, dtype=float32))
                DL[:,:,0] = D[:,:,0]-L
                DL[:,:,1] = D[:,:,1]-L
                DL[:,:,2] = D[:,:,2]-L
                maxL = np.max(L)
                D[:,:,0] = L + colsat*(L/maxL)*DL[:,:,0]
                D[:,:,1] = L + colsat*(L/maxL)*DL[:,:,1]
                D[:,:,2] = L + colsat*(L/maxL)*DL[:,:,2]
                D = minimum(maximum(D,0),1)
                h.set_data(D)
        else:
#               h.set_data(np.power(1.-minimum(maximum((B-maxi)/(mini-maxi),0),1),gamma))
                h.set_data(log( 2.-minimum(maximum((B-maxi)/(mini-maxi),0),1)*gamma)/log(1+gamma))
    draw()
    '''
    figure(101)
    imshow(D)
    draw()
    show()
    '''

def removeHotPixels(L, thresh=2):
    global progress, gauge
    '''
    (m,n) = L.shape

    Laplacian3 = thresh * array([[ 0.,  1.,  0.], \
                                 [ 1.,  0.,  1.], \
                                 [ 0.,  1.,  0.]]) / 4 
    LconvL3 = convolve2d(L, Laplacian3, mode='same')
    k = (LconvL3 <= L) 
    '''
    
    progress.SetLabel('Removing Hot Pixels: ')
    gauge.SetValue(0)
    wx.Yield()
    if len(L.shape) == 2:
        L_min  = order_filter(1.0*L, ones((3,3)), 0)
        L_median  = medfilt2d(1.0*L, (3,3));
        k = ( thresh*(L_median-L_min) <= L-L_min )
        L[k] = L_median[k]
    else:
        N = L.shape[0]
        for j in arange(0,N):
            gauge.SetValue(int(100*(j+0.5)/N))
            wx.Yield()
            Lj = L[j][:,:]
            L_min  = order_filter(1.0*Lj, ones((3,3)), 0)
            L_median  = medfilt2d(1.0*Lj, (3,3));
            k = ( thresh*(L_median-L_min) <= Lj-L_min )
            L[j,k] = L_median[k]
    progress.SetLabel('')
    gauge.SetValue(0)

def findstars(L, thresh=0.0006, fignum=2, bkg=False, returnbkg=False):
    global tab
    print(tab, 'Finding stars:')
    tab = tab+'    '
    (m,n) = L.shape

    # local max...
    M1  = array([[ 0., 1., -1.]])
    M2  = array([[-1., 1.,  0.]])
    M3  = M1.T
    M4  = M2.T
    M5  = array([[-1.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0.,  0.]])
    M6  = array([[ 0.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0., -1.]])
    M7  = array([[ 0.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [-1.,  0.,  0.]])
    M8  = array([[ 0.,  0., -1.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0.,  0.]])
    '''
    M9  = array([[ 0., 0., 0., 1., -1.]])
    M10 = array([[-1., 1., 0., 0.,  0.]])
    M11 = M9.T
    M12 = M10.T
    '''
    M9  = array([[ 0., 0., 1., 0., -1.]])
    M10 = array([[-1., 0., 1., 0.,  0.]])
    M11 = M9.T
    M12 = M10.T
    M21 = array([[-1., 2., -1.]])
    M22 = M21.T
    Laplacian = array([[ 0., -1.,  0.], \
                       [-1.,  4., -1.], \
                       [ 0., -1.,  0.]])
    '''
    Laplacian2 = array([[-1.,  0., -1.], \
                        [ 0.,  4.,  0.], \
                        [-1.,  0., -1.]])
    '''
    '''
    Laplacian2 = array([[ 0., -1., -1., -1.,  0.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [-1.,  2.,  4.,  2., -1.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [ 0., -1., -1., -1.,  0.]])
    '''
    Laplacian1 = array([[-1., -1., -1., -1., -1., -1., -1., -1., -1.], \
                        [-1., -1., -1.,  0.,  0.,  0., -1., -1., -1.], \
                        [-1., -1.,  0.,  1.,  1.,  1.,  0., -1., -1.], \
                        [-1.,  0.,  1.,  3.,  4.,  3.,  1.,  0., -1.], \
                        [-1.,  0.,  1.,  4.,  4.,  4.,  1.,  0., -1.], \
                        [-1.,  0.,  1.,  3.,  4.,  3.,  1.,  0., -1.], \
                        [-1., -1.,  0.,  1.,  1.,  1.,  0., -1., -1.], \
                        [-1., -1., -1.,  0.,  0.,  0., -1., -1., -1.], \
                        [-1., -1., -1., -1., -1., -1., -1., -1., -1.]
    ])
    Laplacian2 = array([[-1., -1., -1., -1., -1., -1., -1.], \
                        [-1., -1.,  0.,  0.,  0., -1., -1.], \
                        [-1.,  0.,  0.,  5.,  0.,  0., -1.], \
                        [-1.,  0.,  5.,  8.,  5.,  0., -1.], \
                        [-1.,  0.,  0.,  5.,  0.,  0., -1.], \
                        [-1., -1.,  0.,  0.,  0., -1., -1.], \
                        [-1., -1., -1., -1., -1., -1., -1.]])
    Laplacian3 = array([[ 0., -1., -1., -1.,  0.], \
                        [-1.,  1.,  2.,  1., -1.], \
                        [-1.,  2.,  0.,  2., -1.], \
                        [-1.,  1.,  2.,  1., -1.], \
                        [ 0., -1., -1., -1.,  0.]])
    # avoid hot pixels...
    Laplacian4 = array([[ 0.,  3.,  0.], \
                        [ 3., -4.,  3.], \
                        [ 0.,  3.,  0.]])
    # brightest pixel...
    Laplacian5 = array([[ 1.]])

    Lconv1 = convolve2d(L, M1, mode='same')
    Lconv2 = convolve2d(L, M2, mode='same')
    Lconv3 = convolve2d(L, M3, mode='same')
    Lconv4 = convolve2d(L, M4, mode='same')
    Lconv5 = convolve2d(L, M5, mode='same')
    Lconv6 = convolve2d(L, M6, mode='same')
    Lconv7 = convolve2d(L, M7, mode='same')
    Lconv8 = convolve2d(L, M8, mode='same')
    Lconv9 = convolve2d(L, M9, mode='same')
    Lconv10 = convolve2d(L, M10, mode='same')
    Lconv11 = convolve2d(L, M11, mode='same')
    Lconv12 = convolve2d(L, M12, mode='same')
    Lconv21 = convolve2d(L, M21, mode='same')
    Lconv22 = convolve2d(L, M22, mode='same')
    LconvL1 = convolve2d(L, Laplacian, mode='same')
    LconvL2 = convolve2d(L, Laplacian2, mode='same')
    LconvL3 = convolve2d(L, Laplacian3, mode='same')
    LconvL4 = convolve2d(L, Laplacian4, mode='same')

    k = (Lconv1 >= 0) & (Lconv2 >= 0) & (Lconv3 >= 0) & (Lconv4 >= 0) & (Lconv5 >= 0) & \
        (Lconv6 >= 0) & (Lconv7 >= 0) & (Lconv8 >= 0) & \
        (LconvL1 >= thresh*np.max(LconvL1)) & (LconvL2 >= thresh*np.max(LconvL2)) & \
        (LconvL3 >= 0) & (LconvL4 >= 0)

    k[:,0:20] = 0
    k[:,n-19:n] = 0
    k[0:20,:] = 0
    k[m-19:m,:] = 0

    fwhms = array(zeros((m,n), dtype=float32))
    for i in arange(m):
        for j in arange(n):
            if (k[i,j] == 1): fwhms[i,j] = fwhm(L,i,j)
    med_fwhm = percentile(fwhms[k],50)
    k[fwhms<0.7*med_fwhm] = 0
    k[fwhms>1.3*med_fwhm] = 0
    if bkg:
        Mbkg = array([[ 0.,  1.,  1.,  1.,  0.], \
                      [ 1.,  0.,  0.,  0.,  1.], \
                      [ 1.,  0.,  0.,  0.,  1.], \
                      [ 1.,  0.,  0.,  0.,  1.], \
                      [ 0.,  1.,  1.,  1.,  0.]])/12
        LconvBkg = convolve2d(L, Mbkg, mode='same')
        bkg_level = percentile(L,100)
        k = k & (LconvBkg < bkg_level)

    if fignum > 0:
        flux = LconvL2[k]
        logL = sqrt(maximum(L-percentile(L,0),0)+1)
        logL = logL - np.min(logL)
        logL = logL/np.max(logL)
        Lrgb = zeros((m,n,3), dtype=float32)
        Lrgb[:,:,0] = logL - logL*k
        Lrgb[:,:,1] = logL + (1-logL)*k
        Lrgb[:,:,2] = logL - logL*k
        figure(fignum)  # figure(2, figsize=(16,16))
        imshow(Lrgb)
        show()

    tab = tab[0:-4]
    printdone( 'Done finding stars')
    return k, LconvL2

# not as good as the "old" one ...
def findstars_new(*arg):
    
    L = array(arg[0])
    if len(arg) > 1:
        thresh = float32(arg[1])
    else:
        thresh = 0.0006
    (m,n) = L.shape

    r0 = 3
    r1 = 3
    r2 = 4
    xs =        arange(-r2,r2+1)
    ys = array([arange(-r2,r2+1)]).T
    radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
    disk = float32(radii <= r0)
    annulus = float32((radii <= r2) & (radii > r1))
    annulus = sum(disk)*annulus/sum(annulus)
    averager = disk - annulus
    LL = convolve2d(L, averager, mode='same')

    M1  = array([[ 0., 1., -1.]])
    M2  = array([[-1., 1.,  0.]])
    M3  = M1.T
    M4  = M2.T
    M5  = array([[-1.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0.,  0.]])
    M6  = array([[ 0.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0., -1.]])
    M7  = array([[ 0.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [-1.,  0.,  0.]])
    M8  = array([[ 0.,  0., -1.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0.,  0.]])

    Lconv1 = convolve2d(LL, M1, mode='same')
    Lconv2 = convolve2d(LL, M2, mode='same')
    Lconv3 = convolve2d(LL, M3, mode='same')
    Lconv4 = convolve2d(LL, M4, mode='same')
    Lconv5 = convolve2d(LL, M5, mode='same')
    Lconv6 = convolve2d(LL, M6, mode='same')
    Lconv7 = convolve2d(LL, M7, mode='same')
    Lconv8 = convolve2d(LL, M8, mode='same')
    k = (Lconv1 >= 0) & (Lconv2 >= 0) & (Lconv3 >= 0) & (Lconv4 >= 0) & (Lconv5 >= 0) & \
        (Lconv6 >= 0) & (Lconv7 >= 0) & (Lconv8 >= 0) 
    k[:,1:20] = 0
    k[:,n-19:n] = 0
    k[1:20,:] = 0
    k[m-19:m,:] = 0

    flux = LL[k]
#    logL = log(maximum(L,0)+1)
    logL = log(L-np.min(L)+1)
    logL = logL - np.min(logL)
    logL = logL/np.max(logL)
    Lrgb = zeros((m,n,3), dtype=float32)
    Lrgb[:,:,0] = logL - logL*k
    Lrgb[:,:,1] = logL + (1-logL)*k
    Lrgb[:,:,2] = logL - logL*k
#   figure(2)  # figure(2, figsize=(16,16))
#   imshow(Lrgb)
#    figure(3, figsize=(16,16))
#    imshow(logL, cmap='gray')
#   show()
    return k, LL

def findstars_small_psf(*arg):
    
    L = array(arg[0])
    if len(arg) > 1:
        thresh = float32(arg[1])
    else:
        thresh = 0.0006
    (m,n) = L.shape
    M1  = array([[ 0., 1., -1.]])
    M2  = array([[-1., 1.,  0.]])
    M3  = M1.T
    M4  = M2.T
    M5  = array([[-1.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0.,  0.]])
    M6  = array([[ 0.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0., -1.]])
    M7  = array([[ 0.,  0.,  0.], \
                 [ 0.,  1.,  0.], \
                 [-1.,  0.,  0.]])
    M8  = array([[ 0.,  0., -1.], \
                 [ 0.,  1.,  0.], \
                 [ 0.,  0.,  0.]])
    M21 = array([[-1., 2., -1.]])
    M22 = M21.T
    Laplacian = array([[ 0., -1.,  0.], \
                       [-1.,  4., -1.], \
                       [ 0., -1.,  0.]])

    Lconv1 = convolve2d(L, M1, mode='same')
    Lconv2 = convolve2d(L, M2, mode='same')
    Lconv3 = convolve2d(L, M3, mode='same')
    Lconv4 = convolve2d(L, M4, mode='same')
    Lconv5 = convolve2d(L, M5, mode='same')
    Lconv6 = convolve2d(L, M6, mode='same')
    Lconv7 = convolve2d(L, M7, mode='same')
    Lconv8 = convolve2d(L, M8, mode='same')
    Lconv21 = convolve2d(L, M21, mode='same')
    Lconv22 = convolve2d(L, M22, mode='same')
    LconvL1 = convolve2d(L, Laplacian, mode='same')
    k = (Lconv1 >= 0) & (Lconv2 >= 0) & (Lconv3 >= 0) & (Lconv4 >= 0) & (Lconv5 >= 0) & \
        (Lconv6 >= 0) & (Lconv7 >= 0) & (Lconv8 >= 0) & \
        (Lconv21 >= 0.25*Lconv22) & (Lconv22 >= 0.25*Lconv21) & \
        (LconvL1 >= thresh*np.max(LconvL1)) 
    k[:,1:20] = 0
    k[:,n-19:n] = 0
    k[1:20,:] = 0
    k[m-19:m,:] = 0

    flux = LconvL1[k]
#    logL = log(maximum(L,0)+1)
    logL = log(L-np.min(L)+1)
    logL = logL - np.min(logL)
    logL = logL/np.max(logL)
    Lrgb = zeros((m,n,3), dtype=float32)
    Lrgb[:,:,0] = logL - logL*k
    Lrgb[:,:,1] = logL + (1-logL)*k
    Lrgb[:,:,2] = logL - logL*k
#   figure(2)  # figure(2, figsize=(16,16))
#   imshow(Lrgb)
#    figure(3, figsize=(16,16))
#    imshow(logL, cmap='gray')
#   show()
    return k, LconvL1

def setArcsecPerPixel(app):
    global ArcsecPerPixel
    ArcsecPerPixel = app

def setMaxBoxSize(mbs):
    global MaxBoxSize
    MaxBoxSize = mbs

def showAlignStar(TrueFalse, fnames=''):
    global showAlign, FileNames
    showAlign = TrueFalse
    FileNames = fnames
#   print(tab,'FileNames = ', FileNames)

# for the onclick1 function (which doesn't work)
def setBoxsize(boxsize):
    global box_size
    box_size = boxsize

def setExptime(exp_time):
    global exptime
    exptime = exp_time

#def makeHRdiag(RGB, objname, minmag0, xctr, yctr, rad):
#def makeHRdiag(RGB, objname, minmag0, shape, shape_params, fignum, threshold):
    '''
def makeHRdiag(*arg):
    if len(arg) > 6: 
        thresh       = arg[6]; 
    else: 
        thresh       = 0.0001;
    if len(arg) > 5: 
        fignum       = arg[5]; 
    else: 
        fignum       = 1;
    if len(arg) > 4: 
        shape_params = arg[4]; 
    else: 
        shape_params = 20;
    if len(arg) > 3: 
        shape        = arg[3]; 
    else: 
        shape        = 'rect';
    if len(arg) > 2: 
        minmag0      = arg[2]; 
    else: 
        minmag0      = 10;
    if len(arg) > 1: 
        objname      = arg[1]; 
    else:
        objname      = '';
    RGB = arg[0]
    '''
def makeHRdiag(RGB, objname='', minmag0=10, shape='rect', shape_params=20, \
                fignum=1, thresh=0.0001):

    global exptime
    
    RGB = array(RGB)
    R = RGB[:,:,0]/exptime
    G = RGB[:,:,1]/exptime
    B = RGB[:,:,2]/exptime
    L = (R+G+B)/3

    (m,n) = L.shape

    k, LconvL2 = findstars(L,thresh,fignum=0)
    Laplacian2 = array([[-1., -1.,  0., -1., -1.], \
                        [-1.,  1.,  1.,  1., -1.], \
                        [ 0.,  1.,  4.,  1.,  0.], \
                        [-1.,  1.,  1.,  1., -1.], \
                        [-1., -1.,  0., -1., -1.]])
#    Laplacian2 = array([[ 0., -1., -1., -1.,  0.], \
#                        [-1.,  0.,  2.,  0., -1.], \
#                        [-1.,  2.,  4.,  2., -1.], \
#                        [-1.,  0.,  2.,  0., -1.], \
#                        [ 0., -1., -1., -1.,  0.]])
    RconvL2 = convolve2d(R, Laplacian2, mode='same')
    BconvL2 = convolve2d(B, Laplacian2, mode='same')
    k = k & (RconvL2 > thresh*np.max(RconvL2)) & (BconvL2 > thresh*np.max(BconvL2)) 

    if (shape == 'circ'):
        xctr = shape_params[0]
        yctr = shape_params[1]
        rad  = shape_params[2]
        xs = arange(0,n)  -xctr
        ys = array([arange(0,m)]).T-yctr
        r2 = sqrt((ones((m,1))*xs)**2 + (ys*ones((1,n)))**2)
        mask = (r2 < rad)
    elif (shape == 'annulus'):
        xctr = shape_params[0]
        yctr = shape_params[1]
        rad0  = shape_params[2]
        rad1  = shape_params[3]
        xs = arange(0,n)  -xctr
        ys = array([arange(0,m)]).T-yctr
        r2 = sqrt((ones((m,1))*xs)**2 + (ys*ones((1,n)))**2)
        mask = ( (rad0 < r2) & (r2 < rad1) )
    else:
        border = shape_params
        xs = arange(0,n)
        ys = array([arange(0,m)]).T
        Xs = ones((m,1))*xs
        Ys = ys*ones((1,n))
        mask = (Xs > border) & (Xs < n-border) & (Ys > border) & (Ys < m-border)
    k = k&mask;
    flux = LconvL2[k]
    Rflux = RconvL2[k]
    Bflux = BconvL2[k]

    logL = log(maximum(L,0)+1)
    logL = logL - np.min(logL)
    logL = logL/np.max(logL)
    Lrgb = zeros((m,n,3), dtype=float32)
    Lrgb[:,:,0] = logL - logL*k
    Lrgb[:,:,1] = logL + (1-logL)*k
    Lrgb[:,:,2] = logL - logL*k
#   figure(1, figsize=(16,12))
#   imshow((RGB-np.min(RGB)+1)/(np.max(RGB)-np.min(RGB)), norm=LogNorm())
#   figure(30+fignum)  
    figure(30+fignum, figsize=(16,16))
    imshow(Lrgb, cmap='gray')
    draw()
#   show()

    I = argsort(flux)
    fluxsort = flux[I]
    Rfluxsort = Rflux[I]
    Bfluxsort = Bflux[I]
    
    minflux = np.min(fluxsort)
    
    n = len(flux)
    color = 2.5*log10(Rfluxsort/Bfluxsort)
    print('numstars = ', n)
    
    fluxmag = 2.5*log10(fluxsort)
    maxfluxmag = np.max(fluxmag)
    fluxmag = fluxmag-maxfluxmag - minmag0
    fluxmag = -fluxmag
    maxmag = np.max(fluxmag)
    minmag = np.min(fluxmag)
#   fig5 = figure(fignum)
    fig5 = figure(fignum, figsize=(12,9))
    axes = fig5.add_subplot(1,1,1,facecolor='black')

    fig5.set_facecolor('black')
    axes.spines['bottom'].set_color('white')
    axes.spines['left'].set_color('white')
    axes.spines['top'].set_color('white')
    axes.spines['right'].set_color('white')
    axes.xaxis.label.set_color('white')
    axes.yaxis.label.set_color('white')
    axes.tick_params(axis='x', colors='white')
    axes.tick_params(axis='y', colors='white')

#hold(True)
    color    = zeros(n,dtype=float32)
    fluxmag  = zeros(n,dtype=float32)
    marksize = zeros(n,dtype=float32)
    rgb      = zeros((n,4),dtype=float32)
    for k in arange(0,n):
            flxmg = 2.5*log10(fluxsort[k])
            flxmg = flxmg - maxfluxmag - minmag0
            flxmg = -flxmg
            bright = sqrt((flxmg-maxmag)/(minmag-maxmag))
            maxmag1 = max([maxmag, 18])
            diam = ((flxmg-maxmag1)/(12-maxmag1))
            clr = 2.5*log10(Rfluxsort[k]/Bfluxsort[k])
            r = max([min([1 + 0.7*clr,1]),0])
            b = max([min([1 - 0.7*clr,1]),0])
            g = 1 - 0.5*(1-min([r,b]))
            r = bright*r
            g = bright*g
            b = bright*b
            marksize[k] = (0.001+diam*12.);
            color[k] = clr
            fluxmag[k] = flxmg
            rgb[k,0] = r
            rgb[k,1] = g
            rgb[k,2] = b
            rgb[k,3] = 0.8
    scatter(color, fluxmag, facecolor=rgb, s=marksize, edgecolor=rgb/3);
    axis([-2, 2.1, minmag0-0.6, maxmag+0.2])
    xlabel('Color index (B-R)')
    ylabel('Brightness (apparent magnitude)')
    title('Hertzsprung-Russell Diagram for '+objname, color='white')
#hold(False)
    gca().invert_yaxis()
    draw()
#   show()
    return fig5  # in case it's needed

def weighted_median(values, weights):
    ''' compute the weighted median of values list. The 
    weighted median is computed as follows:
    1- sort both lists (values and weights) based on values.
    2- select the 0.5 point from the weights and return the corresponding values as results
    e.g. values = [1, 3, 0] and weights=[0.1, 0.3, 0.6] assuming weights are probabilities.
    sorted values = [0, 1, 3] and corresponding sorted weights = [0.6,     0.1, 0.3] the 0.5 point on
    weight corresponds to the first item which is 0. so the weighted     median is 0.'''

    #convert the weights into probabilities
    sum_weights = sum(weights)
    weights = array([(w*1.0)/sum_weights for w in weights])
    #sort values and weights based on values
    values = array(values)
    sorted_indices = argsort(values)
    values_sorted  = values[sorted_indices]
    weights_sorted = weights[sorted_indices]
    #select the median point
    it = nditer(weights_sorted, flags=['f_index'])
    accumulative_probability = 0
    median_index = -1
    while not it.finished:
        accumulative_probability += it[0]
        if accumulative_probability > 0.5:
            median_index = it.index
            return values_sorted[median_index]
        elif accumulative_probability == 0.5:
            median_index = it.index
            it.iternext()
            next_median_index = it.index
            return mean(values_sorted[[median_index, next_median_index]])
        it.iternext()

    return values_sorted[median_index]

def weighted_percentile(values, weights, p):
    #convert the weights into probabilities
    sum_weights = sum(weights)
    weights = array([(w*1.0)/sum_weights for w in weights])
    #sort values and weights based on values
    values = array(values)
    sorted_indices = argsort(values)
    values_sorted  = values[sorted_indices]
    weights_sorted = weights[sorted_indices]
    #select the percentile point
    it = nditer(weights_sorted, flags=['f_index'])
    accumulative_probability = 0
    pctile_index = -1
    while not it.finished:
        accumulative_probability += it[0]
        if accumulative_probability > p:
            pctile_index = it.index
            return values_sorted[pctile_index]
        elif accumulative_probability == p:
            pctile_index = it.index
            it.iternext()
            next_pctile_index = it.index
            return mean(values_sorted[[pctile_index, next_pctile_index]])
        it.iternext()

    return values_sorted[pctile_index]


'''
def makeHRdiag_new(*arg):
    if len(arg) > 6: 
        thresh       = arg[6]; 
    else: 
        thresh       = 0.0001;
    if len(arg) > 5: 
        fignum       = arg[5]; 
    else: 
        fignum       = 1;
    if len(arg) > 4: 
        shape_params = arg[4]; 
    else: 
        shape_params = 20;
    if len(arg) > 3: 
        shape        = arg[3]; 
    else: 
        shape        = 'rect';
    if len(arg) > 2: 
        minmag0      = arg[2]; 
    else: 
        minmag0      = 10;
    if len(arg) > 1: 
        objname      = arg[1]; 
    else:
        objname      = '';
    RGB = arg[0]
'''

fig5 = 0
fig6 = 0

def makeHRdiag_new(RGB, objname='', minmag0=10, shape='rect', shape_params=20, \
                fignum=1, thresh=0.0001, maxmag0=20):

    global GUI, exptime, name, folder, A
    global xpos
    
    A = RGB
    RGB = array(RGB)
    R = RGB[:,:,0]/exptime
    G = RGB[:,:,1]/exptime
    B = RGB[:,:,2]/exptime
    L = (R+G+B)/3
    L = (R+B)/2

    (m,n) = L.shape

    if (shape == 'circ'):
        xctr = shape_params[0]
        yctr = shape_params[1]
        rad  = shape_params[2]
        xs = arange(0,n)  -xctr
        ys = array([arange(0,m)]).T-yctr
        r2 = sqrt((ones((m,1))*xs)**2 + (ys*ones((1,n)))**2)
        mask = (r2 < rad)
    elif (shape == 'annulus'):
        xctr = shape_params[0]
        yctr = shape_params[1]
        rad0  = shape_params[2]
        rad1  = shape_params[3]
        xs = arange(0,n)  -xctr
        ys = array([arange(0,m)]).T-yctr
        r2 = sqrt((ones((m,1))*xs)**2 + (ys*ones((1,n)))**2)
        mask = ( (rad0 < r2) & (r2 < rad1) )
    elif (shape == 'auto'):
#       p = shape_params

        xs = arange(0,n).T  # values
        ys = arange(0,m).T
        Lmatrix = matrix(L)

        dn = int(round(n/3))  # it would be better to use rad1 instead of n/5
        dm = int(round(m/3))  # it would be better to use rad1 instead of m/5

        Lx = Lmatrix.sum(axis=0).T # weight
        xctr = minimum(maximum(weighted_median(xs,Lx),dm),n-dm-1)

        Ly = Lmatrix.sum(axis=1)
        yctr = weighted_median(ys,Ly)
        yctr = minimum(maximum(weighted_median(ys,Ly),dm),m-dm-1)

#       print('xctr = ', xctr, ', yctr = ',yctr)

        xs = arange(xctr-dn,xctr+dn).T  # values
        ys = arange(yctr-dm,yctr+dm).T
#       print('xctr = ', xctr, ', yctr = ', yctr, ', dm = ', dm)
        Lmatrix = matrix(L[yctr-dm:yctr+dm,xctr-dn:xctr+dn])
#       print('shape of Lmatrix: ', Lmatrix.shape)

        Lx = Lmatrix.sum(axis=0).T # weight
#       print('shape of xs: ', xs.shape)
#       print('shape of Lx: ', Lx.shape)
        xctr = weighted_median(xs,Lx)

        Ly = Lmatrix.sum(axis=1)
#       print('shape of ys: ', ys.shape)
#       print('shape of Ly: ', Ly.shape)
        yctr = weighted_median(ys,Ly)

#       print('xctr = ', xctr, ', yctr = ',yctr)

        xs = arange(0,n)  -xctr
        ys = array([arange(0,m)]).T-yctr
        r2 = sqrt((ones((m,1))*xs)**2 + (ys*ones((1,n)))**2)

        maxr = np.min([xctr, yctr, n-xctr, m-yctr])
        maxr10 = int(round(maxr/10))
        Lr = zeros((maxr10), dtype=float32)
        wt = zeros((maxr10), dtype=float32)
        for i in range(yctr-maxr,yctr+maxr):
            for j in range(xctr-maxr,xctr+maxr):
                k = int(round(r2[i,j]/10))
                if (k<maxr10):
                    Lr[k] += L[i,j]
                    wt[k] += 1
        for k in range(0,maxr10):
            if (wt[k] > 0): Lr[k] = Lr[k]/wt[k]
        Lr_max = np.max(Lr)
        '''
        Lr_min = np.min(Lr)
        '''
        Lr_thresh = percentile(Lr, 10)
        hit_max = False
        rad1 = maxr
        for k in range(0,maxr10):
            if (Lr[k] >= Lr_max): 
                    rad0 = 10*k
                    hit_max = True
            if ((Lr[k] <= Lr_thresh) & hit_max):
                    rad1 = 10*k
                    break

        '''
        r2linear = reshape(r2, m*n)
        Llinear = reshape(L, m*n)
        rad1 = weighted_percentile(r2linear, Llinear, p)
        rad0 = 0.02*rad1
        print('rad0 = ', rad0, ', rad1 = ',rad1)
        '''

        mask = ( (rad0 < r2) & (r2 < rad1) )
    else:
        border = shape_params
        xs = arange(0,n)
        ys = array([arange(0,m)]).T
        Xs = ones((m,1))*xs
        Ys = ys*ones((1,n))
        mask = (Xs > border) & (Xs < n-border) & (Ys > border) & (Ys < m-border)

# k, LconvL2 = findstars_new(L,thresh)
    k, LconvL2 = findstars(L,thresh,fignum=0)


# an alternative method.   works good (2,2,3 is also good)
# an alternative method.   works good (3,3,4 is also good)
    r0 = 2.5
    r1 = 2.5
    r2 = 4
    xs =        arange(-r2,r2+1)
    ys = array([arange(-r2,r2+1)]).T
    radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
    disk = float32(radii <= r0)
    annulus = float32((radii <= r2) & (radii > r1))
#   print('size of disk = ', sum(disk))
#   print('disk = \n', disk)
#   print('size of annulus = ', sum(annulus))
#   print('annulus = \n', annulus)
    annulus = sum(disk)*annulus/sum(annulus)
    averager = disk - annulus
    LconvL2 = convolve2d(L, averager, mode='same')
    RconvL2 = convolve2d(R, averager, mode='same')
    GconvL2 = convolve2d(G, averager, mode='same')
    BconvL2 = convolve2d(B, averager, mode='same')
    '''
    disk = disk/sum(disk)
    LconvL2 = convolve2d(L, disk, mode='same')
    RconvL2 = convolve2d(R, disk, mode='same')
    GconvL2 = convolve2d(G, disk, mode='same')
    BconvL2 = convolve2d(B, disk, mode='same')
    LconvL2 = LconvL2 - medfilt2d(L, (2*r2+1,2*r2+1))
    RconvL2 = RconvL2 - medfilt2d(R, (2*r2+1,2*r2+1))
    GconvL2 = GconvL2 - medfilt2d(G, (2*r2+1,2*r2+1))
    BconvL2 = BconvL2 - medfilt2d(B, (2*r2+1,2*r2+1))
    '''
    
#    Laplacian2 = array([[ 0., -1., -1., -1.,  0.], \
#                        [-1.,  0.,  2.,  0., -1.], \
#                        [-1.,  2.,  4.,  2., -1.], \
#                        [-1.,  0.,  2.,  0., -1.], \
#                        [ 0., -1., -1., -1.,  0.]])
#    RconvL2 = convolve2d(R, Laplacian2, mode='same')
#    BconvL2 = convolve2d(B, Laplacian2, mode='same')

    k = k & (RconvL2 > thresh*np.max(RconvL2)) & \
            (BconvL2 > thresh*np.max(BconvL2)) & \
            (GconvL2 > thresh*np.max(GconvL2)) & \
            (LconvL2 > thresh*np.max(LconvL2)) 

#   k = k&mask;
    flux = LconvL2[k]
    Rflux = RconvL2[k]
    Gflux = GconvL2[k]
    Bflux = BconvL2[k]

    logL = log(maximum(L,0)+1)
    logL = logL - np.min(logL)
    logL = logL/np.max(logL)
    Lrgb = zeros((m,n,3), dtype=float32)
    Lrgb[:,:,0] = logL - logL*k
    Lrgb[:,:,1] = logL + (1-logL)*k
    Lrgb[:,:,2] = logL - logL*k
#   figure(1, figsize=(16,12))
#   imshow((RGB-np.min(RGB)+1)/(np.max(RGB)-np.min(RGB)), norm=LogNorm())
    if (GUI == False):
        figure(30+fignum)   # figure(30+fignum, figsize=(16,16))
        imshow(Lrgb, cmap='gray')
        draw()
#   show()

    y = array([arange(0,m)]).T * ones((1,n))
    x = ones((m,1)) * array([arange(0,n)])
    x = int16(x[k])
    y = int16(y[k])

    I = argsort(flux)
    fluxsort = flux[I]
    Rfluxsort = Rflux[I]
    Gfluxsort = Gflux[I]
    Bfluxsort = Bflux[I]
    xsort = x[I]
    ysort = y[I]
    
    minflux = np.min(fluxsort)
    
    n = len(flux)
    twenty = ones(20)/20.
#   color = 2.5*log10(2.*Rfluxsort/(Gfluxsort+Bfluxsort))
    color = 2.5*log10(Rfluxsort/Bfluxsort)
#color2 = convolve(color, twenty, mode='same')
#color2 = medfilt(color, 21)
    color2 = zeros(n,dtype=float32)
    color2[:] = 1.0*color[:]
#   print('numstars = ', n)
#   print('Lflux = ',L[ysort[n-1],xsort[n-1]])
#   print('fwhm = ',fwhm(L,ysort[n-1],xsort[n-1]))
    
    fluxmag = 2.5*log10(fluxsort)
    maxfluxmag = np.max(fluxmag)
    fluxmag = fluxmag-maxfluxmag - minmag0
    fluxmag = -fluxmag
    maxmag = np.max(fluxmag)
    minmag = np.min(fluxmag)
    print(tab, 'maxmag = ', maxmag, ', minmag = ', minmag)

    clr    = zeros(n,dtype=float32)
    tmp    = zeros(n,dtype=float32)
    fluxmag  = zeros(n,dtype=float32)
    marksize = zeros(n,dtype=float32)
    rgb      = zeros((n,4),dtype=float32)
    text_file = open('txt/color_mag.txt', 'w')
#   print('hello dolly')
    for k in arange(0,n):
            flxmg = 2.5*log10(fluxsort[k])
            flxmg = flxmg - maxfluxmag - minmag0
            flxmg = -flxmg
            dk = int(max([ceil((flxmg-18.0)*0),0]))
            tmp[k] = median(color2[max([k-dk,0]):min([k+dk+1,n+1])])
    median_color = percentile(tmp,70)
#   print('median_color = ', median_color)
    msize = sqrt(640000./n)
#   print('msize = ',msize)
    kk = 0
    maxmag1 = max([maxmag, 20])
    midmag = (maxmag+minmag)/2
    for k in arange(0,n):
            flxmg = 2.5*log10(fluxsort[k])
            flxmg = flxmg - maxfluxmag - minmag0
            flxmg = -flxmg
            if flxmg > maxmag0: continue
            if ((abs(log(Rfluxsort[k]/Gfluxsort[k])) > 2) | (abs(log((Gfluxsort[k]/Bfluxsort[k])) > 2) )): continue
            bright = sqrt(0.00001+(flxmg-maxmag)/(minmag-maxmag))
            diam = ((flxmg-maxmag1)/(midmag-maxmag1))
            dk = int(max(ceil((flxmg-18.0)*0),0))
            clr[kk] = median(color2[max(k-dk,0):min(k+dk+1,n+1)]) - median_color
#           if abs(clr[kk])>1.5: print('color = ',clr[kk],', flxmg = ',flxmg)
#color[k] = mean(color2[max(k-dk,0):min(k+dk+1,n+1)])
            clr_kk = clr[kk]
            r = max([min([1 + 0.7*clr_kk,1]),0])
            b = max([min([1 - 0.7*clr_kk,1]),0])
            g = 1 - 0.5*(1-min([r,b]))
            r = bright*r;  r = min([r,1])
            g = bright*g;  g = min([g,1])
            b = bright*b;  b = min([b,1])
            marksize[kk] = (0.001+diam*msize);
            '''
            if flxmg > 17.5:
                color[k] = median(color2[k-dk:k+dk])
            '''
# color[k] = clr
            fluxmag[kk] = flxmg
            rgb[kk,0] = r
            rgb[kk,1] = g
            rgb[kk,2] = b
            rgb[kk,3] = 0.8
            if flxmg > 17.5:
                text_file.write(''+'{:.4f}'.format(clr_kk)+' '+'{:.4f}'.format(flxmg)+'\n')

            Rflux[kk] = Rfluxsort[k]
            Gflux[kk] = Gfluxsort[k]
            Bflux[kk] = Bfluxsort[k]

            kk += 1
    '''
    averager = zeros(100,dtype=float32) + 1./100.
    color = convolve(color, averager, mode='same')
    '''
    N = kk

    Rflux    = Rflux[0:N]
    Gflux    = Gflux[0:N]
    Bflux    = Bflux[0:N]
    clr      = clr[0:N]
    fluxmag  = fluxmag[0:N]
    rgb      = rgb[0:N, 0:4]
    marksize = marksize[0:N]


    global plotter, fig5, fig6, fig5axes, fig6axes
    global images_frame

    if (fignum==1):
#       fig5 = figure(fignum)   #  
        if (GUI == True):
            fig5 = 0   # delaxes isn't doing what I'd hoped
            if (plotter == ''): 
                images_frame = wx.Frame(None, -1, 'Images', size=[1000,880])
                images_frame.SetPosition(wx.Point(xpos, 26))
                plotter = PlotNotebook(images_frame)
                images_frame.Show(True)
            if (fig5 == 0): 
                fig5 = plotter.add('HR Diagram')
            else:
                print('fig5axes = ', fig5axes)
                fig5.delaxes(fig5axes)
            fig5axes = fig5.add_subplot(1,1,1,facecolor='black')
        else:
            close(fignum)
            fig5 = figure(fignum, figsize=(12,9))

        fig5.set_facecolor('black')
        fig5axes.spines['bottom'].set_color('white')
        fig5axes.spines['left'].set_color('white')
        fig5axes.spines['top'].set_color('white')
        fig5axes.spines['right'].set_color('white')
        fig5axes.xaxis.label.set_color('white')
        fig5axes.yaxis.label.set_color('white')
        fig5axes.tick_params(axis='x', colors='white')
        fig5axes.tick_params(axis='y', colors='white')

        fig5axes.scatter(clr, fluxmag, facecolor=rgb, s=marksize, edgecolor=rgb/3);
#       axis([-2, 2.1, minmag0-0.6, maxmag0])

        colorrange = minimum(maximum(max(clr),-min(clr))+0.4, 2.0)
        fig5axes.axis([-colorrange, colorrange, minmag0-0.6, maxmag0])
        print(tab, 'min/max mag = ', minmag0-0.6, maxmag0)

        fig5axes.set_xlabel('Color index (B-R)')
        fig5axes.set_ylabel('Brightness (apparent magnitude)')
        fig5axes.set_title('Hertzsprung-Russell Diagram for '+objname.upper(), color='white')
        fig5axes.invert_yaxis()
        if (GUI == False): draw()
#       images_frame.Show(True)
#       show()

    if (fignum==6):
        if (GUI == True):
            if (plotter == ''): 
                images_frame = wx.Frame(None, -1, 'Images', size=[1000,880])
                images_frame.SetPosition(wx.Point(xpos, 26))
                plotter = PlotNotebook(images_frame)
                images_frame.Show(True)
            fig6 = 0   # the clf() and cla() aren't doing what I'd hoped
            if (fig6 == 0): 
                fig6 = plotter.add('Color Diag')
            else:
                fig6.clf()
                fig6axes.cla()
        else:
            close(fignum)
            fig6 = figure(fignum, figsize=(12,9))
        fig6axes = fig6.add_subplot(1,1,1,facecolor='black')
    
        fig6.set_facecolor('black')
        fig6axes.spines['bottom'].set_color('white')
        fig6axes.spines['left'].set_color('white')
        fig6axes.spines['top'].set_color('white')
        fig6axes.spines['right'].set_color('white')
        fig6axes.xaxis.label.set_color('white')
        fig6axes.yaxis.label.set_color('white')
        fig6axes.tick_params(axis='x', colors='white')
        fig6axes.tick_params(axis='y', colors='white')
        fig6axes.scatter(log(Gflux/Rflux), log(Bflux/Gflux), facecolor=rgb, s=marksize, edgecolor=rgb/3)
        xmin = percentile(log(Gflux/Rflux),0)
        xmax = percentile(log(Gflux/Rflux),100)
        ymin = percentile(log(Bflux/Gflux),0)
        ymax = percentile(log(Bflux/Gflux),100)
        xmin = xmin - 0.1*(xmax-xmin)
        ymin = ymin - 0.1*(ymax-ymin)
        xmax = xmax + 0.1*(xmax-xmin)
        ymax = ymax + 0.1*(ymax-ymin)
        fig6axes.axis([ xmin, xmax,  ymin, ymax])
        fig6axes.set_xlabel('Log(G/R)')
        fig6axes.set_ylabel('Log(B/G)')
        fig6axes.set_title('Stellar Locus for '+objname, color='white')
#       gca().invert_yaxis()
    
    text_file.close()


    return 


def makeHRdiagRRlyrae(*arg):
    if len(arg) > 8: 
        thresh       = arg[8]; 
    else: 
        thresh       = 0.0001;
    if len(arg) > 7: 
        z            = arg[7]; 
    else: 
        z            = 1.0;
    if len(arg) > 6: 
        showit       = arg[6]; 
    else: 
        showit       = 'merge';
    if len(arg) > 5: 
        colorthresh  = arg[5]; 
    else: 
        colorthresh  = 1.0;
    if len(arg) > 4: 
        fluxthresh   = arg[4]; 
    else: 
        fluxthresh   = 0.5;
    if len(arg) > 3: 
        minmag0      = arg[3]; 
    else: 
        minmag0      = 11.0;
    if len(arg) > 2: 
        objname      = arg[2]; 
    else: 
        objname      = 'M5';

    RGB2         = arg[1]; 
    RGB          = arg[0]; 

    RGB = array(RGB)
    R = RGB[:,:,0]
    G = RGB[:,:,1]
    B = RGB[:,:,2]
    L = (R+G+B)/3

    RGB2 = array(RGB2)
    R2 = RGB2[:,:,0]
    G2 = RGB2[:,:,1]
    B2 = RGB2[:,:,2]
    L2 = (R2+G2+B2)/3

    (m,n) = L.shape

    k, LconvL2 = findstars(L, fignum=0)


# an alternative method.   works okay
    r0 = 3
    r1 = 3
    r2 = 4
    xs =        arange(-r2,r2+1)
    ys = array([arange(-r2,r2+1)]).T
    radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
    disk = float32(radii <= r0)
    annulus = float32((radii <= r2) & (radii > r1))
    annulus = sum(disk)*annulus/sum(annulus)
    averager = disk - annulus
    LconvL2 = convolve2d(L, averager, mode='same')
    RconvL2 = convolve2d(R, averager, mode='same')
    BconvL2 = convolve2d(B, averager, mode='same')
    L2convL2 = convolve2d(L2, averager, mode='same')
    R2convL2 = convolve2d(R2, averager, mode='same')
    B2convL2 = convolve2d(B2, averager, mode='same')



    '''
    Laplacian2 = array([[ 0., -1., -1., -1.,  0.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [-1.,  2.,  4.,  2., -1.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [ 0., -1., -1., -1.,  0.]])
    LconvL2 = convolve2d(L, Laplacian2, mode='same')
    RconvL2 = convolve2d(R, Laplacian2, mode='same')
    BconvL2 = convolve2d(B, Laplacian2, mode='same')
    L2convL2 = convolve2d(L2, Laplacian2, mode='same')
    R2convL2 = convolve2d(R2, Laplacian2, mode='same')
    B2convL2 = convolve2d(B2, Laplacian2, mode='same')
    '''
    k = k & (RconvL2 > thresh*np.max(RconvL2)) & (BconvL2 > thresh*np.max(BconvL2)) \
          & (L2convL2 > thresh*np.max(RconvL2)) \
          & (R2convL2 > thresh*np.max(RconvL2)) & (B2convL2 > thresh*np.max(BconvL2))

    k[:,0:10] = 0
    k[:,n-9:n] = 0
    k[0:10,:] = 0
    k[m-9:m,:] = 0

    y = array([arange(0,m)]).T * ones((1,n))
    x = ones((m,1)) * array([arange(0,n)])
    x = x[k]
    y = y[k]

    Lflux = LconvL2[k]
    Rflux = RconvL2[k]
    Bflux = BconvL2[k]
    L2flux = L2convL2[k]
    R2flux = R2convL2[k]
    B2flux = B2convL2[k]

    print(np.max(k))
    print(np.max(abs(L2flux-Lflux)))

    I = argsort(Lflux)
    Lfluxsort = Lflux[I]
    Rfluxsort = Rflux[I]
    Bfluxsort = Bflux[I]
    L2fluxsort = L2flux[I]
    R2fluxsort = R2flux[I]
    B2fluxsort = B2flux[I]
    xsort = x[I]
    ysort = y[I]
    
    minflux = np.min(Lfluxsort)
    
    n = len(Lflux)
    color = 2.5*log10(Rfluxsort/Bfluxsort)
    
    Lfluxmag = 2.5*log10(Lfluxsort);            L2fluxmag = 2.5*log10(L2fluxsort)                
    maxLfluxmag = np.max(Lfluxmag);             maxL2fluxmag = np.max(L2fluxmag)                
    Lfluxmag = Lfluxmag-maxLfluxmag - minmag0;  L2fluxmag = L2fluxmag-maxL2fluxmag - minmag0
    Lfluxmag = -Lfluxmag;                       L2fluxmag = -L2fluxmag                

    maxmag = np.max(Lfluxmag)
    minmag = np.min(Lfluxmag)

    avgmag = mean(Lfluxmag[n-100:n])
    avgmag2 = mean(L2fluxmag[n-100:n])
    L2fluxmag = L2fluxmag - avgmag2 + avgmag
    
    colormag = 2.5*log10(Rfluxsort/Bfluxsort)
    color2mag = 2.5*log10(R2fluxsort/B2fluxsort)
    avgcol  = percentile(colormag[n-500:n], 40)
    avgcol2 = percentile(color2mag[n-500:n],40)

    A = zeros((n,2), dtype=float32)
    A[:,0] = (Lfluxmag+L2fluxmag)/2
    A[:,1] = array(ones((n)))
    b = log(1e-30+abs(Lfluxmag-L2fluxmag));
    slope, intercept = linalg.lstsq( A, b, rcond=None )[0]
    ab = array([slope, intercept])
    residual = sqrt(mean((b-matmul(A,ab))**2))

    fig9 = figure(9)
    plot(A[:,0], b, 'b.', markersize=0.8)
    Aminmax = array([min(A[:,0]), max(A[:,0])])
    plot(Aminmax, intercept+Aminmax*slope+z*residual, 'r-', linewidth=0.5)

#   fig5 = figure(5)   
    fig5 = figure(5, figsize=(12,12))
    axes = fig5.add_subplot(1,1,1,facecolor='black')

    fig5.set_facecolor('black')
    axes.spines['bottom'].set_color('white')
    axes.spines['left'].set_color('white')
    axes.spines['top'].set_color('white')
    axes.spines['right'].set_color('white')
    axes.xaxis.label.set_color('white')
    axes.yaxis.label.set_color('white')
    axes.tick_params(axis='x', colors='white')
    axes.tick_params(axis='y', colors='white')

    bb = b

#hold(True)
    m = 0
    printflag = 1
    for k in arange(0,n):
        bright = sqrt((Lfluxmag[k]-maxmag)/(minmag-maxmag))
        maxmag = max([maxmag, 18])
        diam = ((Lfluxmag[k]-maxmag)/(12-maxmag))
        color  = 2.5*log10(Rfluxsort[k]/Bfluxsort[k])   - avgcol
        color2 = 2.5*log10(R2fluxsort[k]/B2fluxsort[k]) - avgcol2
        r = max([min([1 + 0.7*color,1]),0])
        b = max([min([1 - 0.7*color,1]),0])
        g = 1 - 0.5*(1-min([r,b]))
        r = bright*r
        g = bright*g
        b = bright*b
        marksize = (0.0001+diam*2.5)    
            
        if showit == 'merge':
            if ( abs(Lfluxmag[k] - L2fluxmag[k]) > fluxthresh ) \
                & ( log(1e-30+abs(Lfluxmag[k] - L2fluxmag[k])) \
                - (ab[1] + ab[0]*(Lfluxmag[k]+L2fluxmag[k])/2) > z*residual ) \
                & ( abs(color - color2) < colorthresh ):
                plot([color, color2], [Lfluxmag[k], L2fluxmag[k]], 'o' \
                    ,color=(r, g, b) \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
                dx = color2-color
                dy = L2fluxmag[k]-Lfluxmag[k]
                arrow(color+diam*0.1*dx, Lfluxmag[k]+diam*0.1*dy, \
                                (1-2*diam*0.1)*dx, (1-2*diam*0.1)*dy \
                    , color=(r, g, b), width=0.005 \
                                )
            else:
                plot((color+color2)/2, (Lfluxmag[k]+L2fluxmag[k])/2, 'o-' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
        elif showit == 'nolines':
            if ( abs(Lfluxmag[k] - L2fluxmag[k]) > fluxthresh ) \
                & ( log(1e-30+abs(Lfluxmag[k] - L2fluxmag[k])) \
                - (ab[1] + ab[0]*(Lfluxmag[k]+L2fluxmag[k])/2) > z*residual ) \
                & ( abs(color - color2) < colorthresh ):
                plot([color, color2], [Lfluxmag[k], L2fluxmag[k]], 'o' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
            else:
                plot((color+color2)/2, (Lfluxmag[k]+L2fluxmag[k])/2, 'o-' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
        else:
            plot([color], [Lfluxmag[k]], 'o-' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
            plot([color2], [L2fluxmag[k]], 'o-' \
                    ,markerfacecolor=(r, g, b/2) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b/2) \
                    )
    axis([-2, 2.1, minmag0-0.2, maxmag+0.2])
    xlabel('Color index (B-R)')
    ylabel('Brightness (apparent magnitude)')
    title('Hertzsprung-Russell Diagram for '+objname, color='white')
#hold(False)
    gca().invert_yaxis()
    draw()
#   show()

    (m,n) = R.shape

    R =R -median(R );  R =maximum(R ,0.)+1;  
    G =G -median(G );  G =maximum(G ,0.)+1;  
    B =B -median(B );  B =maximum(B ,0.)+1;  
    R2=R2-median(R2);  R2=maximum(R2,0.)+1;  
    G2=G2-median(G2);  G2=maximum(G2,0.)+1;  
    B2=B2-median(B2);  B2=maximum(B2,0.)+1;  

    #####################################################################################

    R1 = reshape(R, m*n)-R.min()+1
    G1 = reshape(G, m*n)-G.min()+1
    B1 = reshape(B, m*n)-B.min()+1
    
    R2 = reshape(R2, m*n)-R2.min()+1
    G2 = reshape(G2, m*n)-G2.min()+1
    B2 = reshape(B2, m*n)-B2.min()+1
    
    M = ones(5)
    
    I1 = argsort(R1)
    I2 = argsort(R2)
    #R2[I2] = R1[I1]
    R2[I2] = R2[I2]*convolve(R1[I1]/R2[I2], M, mode='same')
    
    I1 = argsort(G1)
    I2 = argsort(G2)
    #G2[I2] = G1[I1]
    G2[I2] = G2[I2]*convolve(G1[I1]/G2[I2], M, mode='same')
    
    I1 = argsort(B1)
    I2 = argsort(B2)
    #B2[I2] = B1[I1]
    B2[I2] = B2[I2]*convolve(B1[I1]/B2[I2], M, mode='same')
    
    R2 = reshape(R2, (m,n))
    G2 = reshape(B2, (m,n))
    B2 = reshape(G2, (m,n))
    '''
    '''
    
    #####################################################################################

    '''
    R2 = R2*(sum(R2)/sum(R2))*sum(R)/sum(R)
    G2 = G2*(sum(R2)/sum(G2))*sum(G)/sum(R)
    B2 = B2*(sum(R2)/sum(B2))*sum(B)/sum(R)
    '''
    RGBmax  = maximum(maximum( R.max(),  G.max()),  B.max());
    RGB2max = maximum(maximum(R2.max(), G2.max()), B2.max());
    R =R /RGBmax;  logR  = log(255*R +1);              
    G =G /RGBmax;  logG  = log(255*G +1);              
    B =B /RGBmax;  logB  = log(255*B +1);
    R2=R2/RGB2max; logR2 = log(255*R2+1);            
    G2=G2/RGB2max; logG2 = log(255*G2+1);            
    B2=B2/RGB2max; logB2 = log(255*B2+1);

#    logR = logR - median(logR); logG = logG - median(logG); logB = logB - median(logB);
#    logR = logR+1;                logG = logG+1;                logB = logB+1;
    p = 0.9; q = 1.0-p;
    logR = logR/(q*percentile(logR,99)+p*logR.max());   
    logG = logG/(q*percentile(logG,99)+p*logG.max());   
    logB = logB/(q*percentile(logB,99)+p*logB.max());
    logR = minimum(logR,1.0); logG = minimum(logG,1.0); logB = minimum(logB,1.0)
    Lrgb = zeros((m,n,3), dtype=float32)
    Lrgb[:,:,0] = logR;           Lrgb[:,:,1] = logG;           Lrgb[:,:,2] = logB;

#    logR2 = logR2 - median(logR2); logG2 = logG2 - median(logG2); logB2 = logB2 - median(logB2);
#    logR2 = logR2+1;                 logG2 = logG2+1;                 logB2 = logB2+1;
    logR2 = logR2/(q*percentile(logR2,99)+p*logR2.max());
    logG2 = logG2/(q*percentile(logG2,99)+p*logG2.max());
    logB2 = logB2/(q*percentile(logB2,99)+p*logB2.max());
    logR2 = minimum(logR2,1.0); logG2 = minimum(logG2,1.0); logB2 = minimum(logB2,1.0)
    Lrgb2 = zeros((m,n,3), dtype=float32)
    Lrgb2[:,:,0] = logR2;            Lrgb2[:,:,1] = logG2;            Lrgb2[:,:,2] = logB2;
    hFig23 = figure(23, figsize=(14,12))   # hFig23 = figure(23, figsize=(16,12));
    for iter in arange(0,2):
        imshow(Lrgb);  show(); pause(0.5);
        imshow(Lrgb2); show(); pause(0.5);

def makeHRdiagRRlyrae3(RGB, RGB2, RGB3, objname, minmag0, fluxthresh, colorthresh, showit, z):
    
    RGB = array(RGB)
    R = RGB[:,:,0]
    G = RGB[:,:,1]
    B = RGB[:,:,2]
    L = (R+G+B)/3

    RGB2 = array(RGB2)
    R2 = RGB2[:,:,0]
    G2 = RGB2[:,:,1]
    B2 = RGB2[:,:,2]
    L2 = (R2+G2+B2)/3

    RGB3 = array(RGB3)
    R3 = RGB3[:,:,0]
    G3 = RGB3[:,:,1]
    B3 = RGB3[:,:,2]
    L3 = (R3+G3+B3)/3

    (m,n) = L.shape

    k, LconvL2 = findstars(L, fignum=0)
    Laplacian2 = array([[ 0., -1., -1., -1.,  0.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [-1.,  2.,  4.,  2., -1.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [ 0., -1., -1., -1.,  0.]])
    RconvL2 = convolve2d(R, Laplacian2, mode='same')
    GconvL2 = convolve2d(G, Laplacian2, mode='same')
    BconvL2 = convolve2d(B, Laplacian2, mode='same')
    L2convL2 = convolve2d(L2, Laplacian2, mode='same')
    R2convL2 = convolve2d(R2, Laplacian2, mode='same')
    G2convL2 = convolve2d(G2, Laplacian2, mode='same')
    B2convL2 = convolve2d(B2, Laplacian2, mode='same')
    L3convL2 = convolve2d(L3, Laplacian2, mode='same')
    R3convL2 = convolve2d(R3, Laplacian2, mode='same')
    G3convL2 = convolve2d(G3, Laplacian2, mode='same')
    B3convL2 = convolve2d(B3, Laplacian2, mode='same')
    thresh = 0.0006
    k = k & (RconvL2 > thresh*np.max(RconvL2)) & (BconvL2 > thresh*np.max(BconvL2)) \
          & (L2convL2 > thresh*np.max(RconvL2)) \
          & (R2convL2 > thresh*np.max(RconvL2)) & (B2convL2 > thresh*np.max(BconvL2)) \
          & (L3convL2 > thresh*np.max(RconvL2)) \
          & (R3convL2 > thresh*np.max(RconvL2)) & (B3convL2 > thresh*np.max(BconvL2))

    k[:,0:10] = 0
    k[:,n-9:n] = 0
    k[0:10,:] = 0
    k[m-9:m,:] = 0

    y = array([arange(0,m)]).T * ones((1,n))
    x = ones((m,1)) * array([arange(0,n)])
    x = x[k]
    y = y[k]

    Lflux = LconvL2[k]
    Rflux = RconvL2[k]
    Gflux = GconvL2[k]
    Bflux = BconvL2[k]
    L2flux = L2convL2[k]
    R2flux = R2convL2[k]
    G2flux = G2convL2[k]
    B2flux = B2convL2[k]
    L3flux = L3convL2[k]
    R3flux = R3convL2[k]
    G3flux = G3convL2[k]
    B3flux = B3convL2[k]

    print(np.max(abs(L2flux-Lflux)))
    print('R, B',np.mean(Rflux),np.mean(Bflux))

    I = argsort(Lflux)
    Lfluxsort = Lflux[I]
    Rfluxsort = Rflux[I]
    Gfluxsort = Gflux[I]
    Bfluxsort = Bflux[I]
    L2fluxsort = L2flux[I]
    R2fluxsort = R2flux[I]
    G2fluxsort = G2flux[I]
    B2fluxsort = B2flux[I]
    L3fluxsort = L3flux[I]
    R3fluxsort = R3flux[I]
    G3fluxsort = G3flux[I]
    B3fluxsort = B3flux[I]
    xsort = x[I]
    ysort = y[I]
    print('R, B sorted',np.mean(Rfluxsort),np.mean(Bfluxsort))
    
    minflux = np.min(Lfluxsort)
    
    n = len(Lflux)
    color = 2.5*log10(Rfluxsort/Bfluxsort)
    
    Lfluxmag  = 2.5*log10(Lfluxsort);        
    L2fluxmag = 2.5*log10(L2fluxsort);        
    L3fluxmag = 2.5*log10(L3fluxsort);        

    maxLfluxmag  = np.max(Lfluxmag);        
    maxL2fluxmag = np.max(L2fluxmag);        
    maxL3fluxmag = np.max(L3fluxmag);        

    Lfluxmag  = Lfluxmag-maxLfluxmag - minmag0;
    L2fluxmag = L2fluxmag-maxL2fluxmag - minmag0;
    L3fluxmag = L3fluxmag-maxL3fluxmag - minmag0;

    Lfluxmag  = -Lfluxmag;                        
    L2fluxmag = -L2fluxmag;                        
    L3fluxmag = -L3fluxmag;                        

    maxmag = np.max(Lfluxmag)
    minmag = np.min(Lfluxmag)

    avgmag = mean(Lfluxmag[n-100:n])
    avgmag2 = mean(L2fluxmag[n-100:n])
    avgmag3 = mean(L3fluxmag[n-100:n])
    L2fluxmag = L2fluxmag - avgmag2 + avgmag
    L3fluxmag = L3fluxmag - avgmag3 + avgmag
    
    colormag  = 2.5*log10(Rfluxsort/Bfluxsort)
    color2mag = 2.5*log10(R2fluxsort/B2fluxsort)
    color3mag = 2.5*log10(R3fluxsort/B3fluxsort)
    avgcol  = mean(colormag[ n-100:n])
    avgcol2 = mean(color2mag[n-100:n])
    avgcol3 = mean(color3mag[n-100:n])

    A = zeros((n,2), dtype=float32)
    A[:,0] = (Lfluxmag+L2fluxmag)/2
    A[:,1] = array(ones((n)))
    b = log(1e-30+abs(Lfluxmag-L2fluxmag));
    slope, intercept = linalg.lstsq( A, b, rcond=None )[0]
    ab12 = array([slope, intercept])
    residual2 = sqrt(mean((b-matmul(A,ab12))**2))

    A[:,0] = (Lfluxmag+L3fluxmag)/2
    A[:,1] = array(ones((n)))
    b = log(1e-30+abs(Lfluxmag-L3fluxmag));
    slope, intercept = linalg.lstsq( A, b, rcond=None )[0]
    ab13 = array([slope, intercept])
    residual3 = sqrt(mean((b-matmul(A,ab13))**2))

    A[:,0] = (L2fluxmag+L3fluxmag)/2
    A[:,1] = array(ones((n)))
    b = log(1e-30+abs(L2fluxmag-L3fluxmag));
    slope, intercept = linalg.lstsq( A, b, rcond=None )[0]
    ab23 = array([slope, intercept])
    residual23 = sqrt(mean((b-matmul(A,ab23))**2))

#   fig5 = figure(5)   
    fig5 = figure(5, figsize=(12,12))
    axes = fig5.add_subplot(1,1,1,facecolor='black')

    fig5.set_facecolor('black')
    axes.spines['bottom'].set_color('white')
    axes.spines['left'].set_color('white')
    axes.spines['top'].set_color('white')
    axes.spines['right'].set_color('white')
    axes.xaxis.label.set_color('white')
    axes.yaxis.label.set_color('white')
    axes.tick_params(axis='x', colors='white')
    axes.tick_params(axis='y', colors='white')

#hold(True)
    m = 0
    printflag = 1
    rsum = 0
    gsum = 0
    bsum = 0
    print('R, B sorted',np.mean(Rfluxsort),np.mean(Bfluxsort))
    for k in arange(0,n):
        bright = sqrt((Lfluxmag[k]-maxmag)/(minmag-maxmag))
        maxmag = max([maxmag, 18])
        diam = ((Lfluxmag[k]-maxmag)/(12-maxmag))
        color  = 2.5*log10( Rfluxsort[k]/ Bfluxsort[k])
        color2 = 2.5*log10(R2fluxsort[k]/B2fluxsort[k]) - avgcol2 + avgcol
        color3 = 2.5*log10(R3fluxsort[k]/B3fluxsort[k]) - avgcol3 + avgcol
        '''
        r = max([min([1 + 0.7*color,1]),0])
        b = max([min([1 - 0.7*color,1]),0])
        g = 1 - 0.5*(1-min([r,b]))
        '''

        r = Rfluxsort[k]
        g = Gfluxsort[k]
        b = Bfluxsort[k]
        rgb = max([r,g,b])
        r = r/rgb
        g = g/rgb
        b = b/rgb

        r = bright*r
        g = bright*g
        b = bright*b
        r = max([min([r,1]),0])
        g = max([min([g,1]),0])
        b = max([min([b,1]),0])
        marksize = (0.00005+diam*2)    
        if showit == 'merge':
            if  ( (abs( Lfluxmag[k] - L2fluxmag[k]) > fluxthresh ) \
                | (abs( Lfluxmag[k] - L3fluxmag[k]) > fluxthresh ) \
                | (abs(L2fluxmag[k] - L3fluxmag[k]) > fluxthresh ))\
                &(( log(1e-30+abs(Lfluxmag[k] -  L2fluxmag[k])) \
                - ( ab12[1] +  ab12[0]*( Lfluxmag[k]+L2fluxmag[k])/2) > z*residual2 ) \
                | ( log(1e-30+abs(Lfluxmag[k] -  L3fluxmag[k])) \
                - ( ab13[1] +  ab13[0]*( Lfluxmag[k]+L3fluxmag[k])/2) > z*residual3 ) \
                | ( log(1e-30+abs(L3fluxmag[k] - L2fluxmag[k])) \
                - (ab23[1] + ab23[0]*(L2fluxmag[k]+L3fluxmag[k])/2) > z*residual23 ))\
                & ( abs(color  - color2) < colorthresh ) \
                & ( abs(color  - color3) < colorthresh ) \
                & ( abs(color2 - color3) < colorthresh ):
                plot([color, color2, color3, color], [Lfluxmag[k], L2fluxmag[k], L3fluxmag[k], Lfluxmag[k]], 'o-' \
                    ,color=(r, g, b) \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    ,linewidth=1 \
                        )
            else:
                plot((color+color2+color3)/3, (Lfluxmag[k]+L2fluxmag[k]+L3fluxmag[k])/3, 'o-' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
        elif showit == 'nolines':
            if  ( (abs( Lfluxmag[k] - L2fluxmag[k]) > fluxthresh ) \
                | (abs( Lfluxmag[k] - L3fluxmag[k]) > fluxthresh ) \
                | (abs(L2fluxmag[k] - L3fluxmag[k]) > fluxthresh ))\
                &(( log(1e-30+abs(Lfluxmag[k] -  L2fluxmag[k])) \
                - ( ab12[1] +  ab12[0]*( Lfluxmag[k]+L2fluxmag[k])/2) > z*residual2 ) \
                | ( log(1e-30+abs(Lfluxmag[k] -  L3fluxmag[k])) \
                - ( ab13[1] +  ab13[0]*( Lfluxmag[k]+L3fluxmag[k])/2) > z*residual3 ) \
                | ( log(1e-30+abs(L3fluxmag[k] - L2fluxmag[k])) \
                - (ab23[1] + ab23[0]*(L2fluxmag[k]+L3fluxmag[k])/2) > z*residual23 ))\
                & ( abs(color  - color2) < colorthresh ) \
                & ( abs(color  - color3) < colorthresh ) \
                & ( abs(color2 - color3) < colorthresh ):
                plot([color, color2, color3, color], [Lfluxmag[k], L2fluxmag[k], L3fluxmag[k], Lfluxmag[k]], 'o' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
            else:
                plot((color+color2+color3)/3, (Lfluxmag[k]+L2fluxmag[k]+L3fluxmag[k])/3, 'o-' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
        else:
            plot([color], [Lfluxmag[k]], 'o-' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
            plot([color2], [L2fluxmag[k]], 'o-' \
                    ,markerfacecolor=(r, g, b/2) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b/2) \
                    )
            plot([color3], [L3fluxmag[k]], 'o-' \
                    ,markerfacecolor=(r, g/2, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g/2, b) \
                    )
    print(rsum,gsum,bsum)
    axis([-1.5,1.5,10.8,18.5])
    xlabel('Color index (B-R)')
    ylabel('Brightness (apparent magnitude)')
    title('Hertzsprung-Russell Diagram for '+objname, color='white')
#hold(False)
    gca().invert_yaxis()
    show()

    (m,n) = R.shape

    '''
    Annulus = array([   [ 1.,  1., 1.,  1.,  1.,  1., 1., 1., 1., 1., 1.], \
                        [ 1.,  1., 1.,  1.,  1.,  1., 1., 1., 1., 1., 1.], \
                        [ 1.,  1., 1.,  1.,  1.,  1., 1., 1., 1., 1., 1.], \
                        [ 1.,  1., 1.,  0.,  0.,  0., 0., 0., 1., 1., 1.], \
                        [ 1.,  1., 1.,  0.,  0.,  0., 0., 0., 1., 1., 1.], \
                        [ 1.,  1., 1.,  0.,  0.,  0., 0., 0., 1., 1., 1.], \
                        [ 1.,  1., 1.,  0.,  0.,  0., 0., 0., 1., 1., 1.], \
                        [ 1.,  1., 1.,  0.,  0.,  0., 0., 0., 1., 1., 1.], \
                        [ 1.,  1., 1.,  1.,  1.,  1., 1., 1., 1., 1., 1.], \
                        [ 1.,  1., 1.,  1.,  1.,  1., 1., 1., 1., 1., 1.], \
                        [ 1.,  1., 1.,  1.,  1.,  1., 1., 1., 1., 1., 1.]])/96.
    RconvAnnulus  = convolve2d(R , Annulus, mode='same')
    GconvAnnulus  = convolve2d(G , Annulus, mode='same')
    BconvAnnulus  = convolve2d(B , Annulus, mode='same')
    R2convAnnulus = convolve2d(R2, Annulus, mode='same')
    G2convAnnulus = convolve2d(G2, Annulus, mode='same')
    B2convAnnulus = convolve2d(B2, Annulus, mode='same')
    R3convAnnulus = convolve2d(R3, Annulus, mode='same')
    G3convAnnulus = convolve2d(G3, Annulus, mode='same')
    B3convAnnulus = convolve2d(B3, Annulus, mode='same')
    '''

    R_bkg_median  = medfilt2d(R , (11,11));
    G_bkg_median  = medfilt2d(G , (11,11));
    B_bkg_median  = medfilt2d(B , (11,11));
    R2_bkg_median = medfilt2d(R2, (11,11));
    G2_bkg_median = medfilt2d(G2, (11,11));
    B2_bkg_median = medfilt2d(B2, (11,11));
    R3_bkg_median = medfilt2d(R3, (11,11));
    G3_bkg_median = medfilt2d(G3, (11,11));
    B3_bkg_median = medfilt2d(B3, (11,11));

    R2 = R2 - R2_bkg_median + R_bkg_median;
    G2 = G2 - G2_bkg_median + G_bkg_median;
    B2 = B2 - B2_bkg_median + B_bkg_median;
    R3 = R3 - R3_bkg_median + R_bkg_median;
    G3 = G3 - G3_bkg_median + G_bkg_median;
    B3 = B3 - B3_bkg_median + B_bkg_median;

    R =R -median(R );  R =maximum(R ,0.);  R =R /R.max();  logR  = log(255*R +1);              
    G =G -median(G );  G =maximum(G ,0.);  G =G /G.max();  logG  = log(255*G +1);              
    B =B -median(B );  B =maximum(B ,0.);  B =B /B.max();  logB  = log(255*B +1);
    R2=R2-median(R2);  R2=maximum(R2,0.);  R2=R2/R2.max(); logR2 = log(255*R2+1);            
    G2=G2-median(G2);  G2=maximum(G2,0.);  G2=G2/G2.max(); logG2 = log(255*G2+1);            
    B2=B2-median(B2);  B2=maximum(B2,0.);  B2=B2/B2.max(); logB2 = log(255*B2+1);
    R3=R3-median(R3);  R3=maximum(R3,0.);  R3=R3/R3.max(); logR3 = log(255*R3+1);            
    G3=G3-median(G3);  G3=maximum(G3,0.);  G3=G3/G3.max(); logG3 = log(255*G3+1);            
    B3=B3-median(B3);  B3=maximum(B3,0.);  B3=B3/B3.max(); logB3 = log(255*B3+1);

#    logR = logR - median(logR); logG = logG - median(logG); logB = logB - median(logB);
#    logR = logR+1;                logG = logG+1;                logB = logB+1;
    p = 0.9; q = 1.0-p;
    logR = logR/(q*percentile(logR,99)+p*logR.max());   
    logG = logG/(q*percentile(logG,99)+p*logG.max());   
    logB = logB/(q*percentile(logB,99)+p*logB.max());
    logR = minimum(logR,1.0); logG = minimum(logG,1.0); logB = minimum(logB,1.0)
    Lrgb = zeros((m,n,3), dtype=float32)
    Lrgb[:,:,0] = logR;           Lrgb[:,:,1] = logG;           Lrgb[:,:,2] = logB;

#    logR2 = logR2 - median(logR2); logG2 = logG2 - median(logG2); logB2 = logB2 - median(logB2);
#    logR2 = logR2+1;                 logG2 = logG2+1;                 logB2 = logB2+1;
    logR2 = logR2/(q*percentile(logR2,99)+p*logR2.max());
    logG2 = logG2/(q*percentile(logG2,99)+p*logG2.max());
    logB2 = logB2/(q*percentile(logB2,99)+p*logB2.max());
    logR2 = minimum(logR2,1.0); logG2 = minimum(logG2,1.0); logB2 = minimum(logB2,1.0)
    Lrgb2 = zeros((m,n,3), dtype=float32)
    Lrgb2[:,:,0] = logR2;            Lrgb2[:,:,1] = logG2;            Lrgb2[:,:,2] = logB2;

    logR3 = logR3/(q*percentile(logR3,99)+p*logR3.max());
    logG3 = logG3/(q*percentile(logG3,99)+p*logG3.max());
    logB3 = logB3/(q*percentile(logB3,99)+p*logB3.max());
    logR3 = minimum(logR3,1.0); logG3 = minimum(logG3,1.0); logB3 = minimum(logB3,1.0)
    Lrgb3 = zeros((m,n,3), dtype=float32)
    Lrgb3[:,:,0] = logR3;            Lrgb3[:,:,1] = logG3;            Lrgb3[:,:,2] = logB3;

#   hFig23 = figure(23)  
    hFig23 = figure(23, figsize=(14,12));
    for iter in arange(0,2):
        imshow(Lrgb);  show(); pause(0.01);
        imshow(Lrgb2); show(); pause(0.01);
        imshow(Lrgb3); show(); pause(0.01);

    RGB [:,:,0] = R ; RGB [:,:,1] = G ; RGB [:,:,2] = B ;
    RGB2[:,:,0] = R2; RGB2[:,:,1] = G2; RGB2[:,:,2] = B2;
    RGB3[:,:,0] = R3; RGB3[:,:,1] = G3; RGB3[:,:,2] = B3;
    z = (65535*((RGB - RGB.min())/RGB.ptp())).astype(uint16)
    '''
    with open('RGB1.png', 'wb') as f:
        writer = png.Writer(width=z.shape[1], height=z.shape[0], bitdepth=16)
        # Convert z to the Python list of lists expected by
        # the png writer.
        z2list = z.reshape(-1, z.shape[1]*z.shape[2]).tolist()
        writer.write(f, z2list)
    z = (65535*((RGB2 - RGB2.min())/RGB2.ptp())).astype(uint16)
    with open('RGB2.png', 'wb') as f:
        writer = png.Writer(width=z.shape[1], height=z.shape[0], bitdepth=16)
        # Convert z to the Python list of lists expected by
        # the png writer.
        z2list = z.reshape(-1, z.shape[1]*z.shape[2]).tolist()
        writer.write(f, z2list)
    z = (65535*((RGB3 - RGB3.min())/RGB3.ptp())).astype(uint16)
    with open('RGB3.png', 'wb') as f:
        writer = png.Writer(width=z.shape[1], height=z.shape[0], bitdepth=16)
        # Convert z to the Python list of lists expected by
        # the png writer.
        z2list = z.reshape(-1, z.shape[1]*z.shape[2]).tolist()
        writer.write(f, z2list)
    '''

def makeHRdiagRRlyraeN(RGBs, objname, minmag0, fluxthresh, colorthresh, showit, z):
    
    RGBs = array(RGBs)
    Rs = RGBs[:,:,:,0]
    Gs = RGBs[:,:,:,1]
    Bs = RGBs[:,:,:,2]
    Ls = (Rs+Gs+Bs)/3


    (N,m,n) = Ls.shape

    L = median(Ls, axis=0)
    R = median(Rs, axis=0)
    G = median(Gs, axis=0)
    B = median(Bs, axis=0)

    k, LconvL2 = findstars(L, fignum=0)


# an alternative method.   works okay
    r0 = 3
    r1 = 3
    r2 = 4
    xs =        arange(-r2,r2+1)
    ys = array([arange(-r2,r2+1)]).T
    radii = sqrt((ones(((2*r2+1),1))*xs)**2 + (ys*ones((1,(2*r2+1))))**2)
    disk = float32(radii <= r0)
    annulus = float32((radii <= r2) & (radii > r1))
    annulus = sum(disk)*annulus/sum(annulus)
    averager = disk - annulus
    LconvL2 = convolve2d(L, averager, mode='same')
    RconvL2 = convolve2d(R, averager, mode='same')
    GconvL2 = convolve2d(G, averager, mode='same')
    BconvL2 = convolve2d(B, averager, mode='same')
    


    '''
    Laplacian2 = array([[ 0., -1., -1., -1.,  0.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [-1.,  2.,  4.,  2., -1.], \
                        [-1.,  0.,  2.,  0., -1.], \
                        [ 0., -1., -1., -1.,  0.]])
    RconvL2 = convolve2d(R, Laplacian2, mode='same')
    GconvL2 = convolve2d(G, Laplacian2, mode='same')
    BconvL2 = convolve2d(B, Laplacian2, mode='same')
    '''
    LsconvL2 = zeros((N,m,n), dtype=float32)
    RsconvL2 = zeros((N,m,n), dtype=float32)
    BsconvL2 = zeros((N,m,n), dtype=float32)
    for j in arange(0,N): 
        LsconvL2[j,:,:] = convolve2d(Ls[j,:,:], averager, mode='same')
        RsconvL2[j,:,:] = convolve2d(Rs[j,:,:], averager, mode='same')
        BsconvL2[j,:,:] = convolve2d(Bs[j,:,:], averager, mode='same')
        thresh = 0.0006
        k = k & (LsconvL2[j,:,:] > thresh*np.max(LconvL2)) \
              & (RsconvL2[j,:,:] > thresh*np.max(RconvL2)) \
              & (BsconvL2[j,:,:] > thresh*np.max(BconvL2))

    k[:,0:10] = 0
    k[:,n-9:n] = 0
    k[0:10,:] = 0
    k[m-9:m,:] = 0

    y = array([arange(0,m)]).T * ones((1,n))
    x = ones((m,1)) * array([arange(0,n)])
    x = x[k]
    y = y[k]

    Lflux = LconvL2[k]
    Rflux = RconvL2[k]
    Gflux = GconvL2[k]
    Bflux = BconvL2[k]
    '''
    Lsflux = zeros((N,m,n), dtype=float32)
    Rsflux = zeros((N,m,n), dtype=float32)
    Bsflux = zeros((N,m,n), dtype=float32)
    for j in arange(0,N): 
        Lsflux[j,:,:] = LsconvL2[j,k]
        Rsflux[j,:,:] = RsconvL2[j,k]
        Bsflux[j,:,:] = BsconvL2[j,k]
    '''
    Lsflux = LsconvL2[:,k]
    Rsflux = RsconvL2[:,k]
    Bsflux = BsconvL2[:,k]

    I = argsort(Lflux)
    Lfluxsort = Lflux[I]
    Rfluxsort = Rflux[I]
    Gfluxsort = Gflux[I]
    Bfluxsort = Bflux[I]
    Lsfluxsort = Lsflux[:,I]
    Rsfluxsort = Rsflux[:,I]
    Bsfluxsort = Bsflux[:,I]
    xsort = x[I]
    ysort = y[I]
    
    minflux = np.min(Lfluxsort)
    
    n = len(Lflux)
    print('n = ',n)
    color = 2.5*log10(Rfluxsort/Bfluxsort)
    
    Lfluxmag  = 2.5*log10(Lfluxsort);        
    Lsfluxmag = 2.5*log10(Lsfluxsort);        

    maxLfluxmag  = np.max(Lfluxmag);        
    Lfluxmag  = Lfluxmag-maxLfluxmag - minmag0;
    for j in arange(0,N):
        maxLsfluxmag = np.max(Lsfluxmag[j,:]);        
        Lsfluxmag[j] = Lsfluxmag[j]-maxLsfluxmag - minmag0;


    Lfluxmag  = -Lfluxmag;                        
    Lsfluxmag = -Lsfluxmag;                        

    maxmag = np.max(Lfluxmag)
    minmag = np.min(Lfluxmag)

    avgmag = mean(Lfluxmag[n-10:n])
    avgmags = mean(Lsfluxmag[:,n-10:n], axis=1)
    '''
    for j in arange(0,n):
        Lsfluxmag[:,j] = Lsfluxmag[:,j] - avgmags + avgmag
    '''
    
    colormag   = 2.5*log10(Rfluxsort/Bfluxsort)
    colormags  = 2.5*log10(Rsfluxsort/Bsfluxsort)
    avgcol  = mean(colormag[   n-10:n])
    avgcols = mean(colormags[:,n-10:n], axis=1)
    '''
    for j in arange(0,n):
        colormags[:,j] = colormags[:,j] - avgcols + avgcol
    '''

    print('Lsfluxsort = ', Lsfluxsort[:,n-1])
    print('Lsfluxmag = ', Lsfluxmag[:,n-1])
    print('Rsfluxsort = ', Rsfluxsort[:,n-1])
    print('Bsfluxsort = ', Bsfluxsort[:,n-1])
    print('colormags = ', colormags[:,n-1])
    print('x,y = ', xsort[n-1], ysort[n-1])

    '''
    A = zeros((n,2), dtype=float32)
    A[:,0] = (Lfluxmag+L2fluxmag)/2
    A[:,1] = array(ones((n)))
    b = log(1e-30+abs(Lfluxmag-L2fluxmag));
    slope, intercept = linalg.lstsq( A, b, rcond=None )[0]
    ab = array([slope, intercept])
    residual = sqrt(mean((b-matmul(A,ab))**2))
    '''

    A        = zeros((n,2), dtype=float32)
    ab       = zeros((N+1,2), dtype=float32)
    residual = zeros((N+1),   dtype=float32)
    for j1 in arange(0,N):
        if (j1<N-1): 
                j2 = j1+1
        else:
                j2 = 0
        A[:,0] = (Lsfluxmag[j1]+Lsfluxmag[j2])/2
        A[:,1] = array(ones((n)))
        b = log(1e-30+abs(Lsfluxmag[j1]-Lsfluxmag[j2]));
        slope, intercept = linalg.lstsq( A, b, rcond=None )[0]
        ab[j1,:] = array([slope, intercept])
        residual[j1] = sqrt(mean((b-matmul(A,ab[j1,:]))**2))

#        fig9 = figure(9)
#        plot(A[:,0], b, 'b.', markersize=0.8)
#        Aminmax = array([min(A[:,0]), max(A[:,0])])
#        print(Aminmax)
#        print(slope)
#        print(intercept)
#        print(z)
#        print(residual)
#        plot(Aminmax, intercept+Aminmax*slope+z*residual, 'r-', linewidth=0.5)

#   fig5 = figure(5)   # fig5 = figure(5, figsize=(16,12))
    fig5 = figure(5, figsize=(12,12))   
    axes = fig5.add_subplot(1,1,1,facecolor='black')

    fig5.set_facecolor('black')
    axes.spines['bottom'].set_color('white')
    axes.spines['left'].set_color('white')
    axes.spines['top'].set_color('white')
    axes.spines['right'].set_color('white')
    axes.xaxis.label.set_color('white')
    axes.yaxis.label.set_color('white')
    axes.tick_params(axis='x', colors='white')
    axes.tick_params(axis='y', colors='white')

#hold(True)
    m = 0
    printflag = 1
    for k in arange(0,n):
        bright = sqrt((Lfluxmag[k]-maxmag)/(minmag-maxmag))
#       maxmag = max([maxmag, 18])
        diam = ((Lfluxmag[k]-maxmag)/(12-maxmag))
        alph = ((Lfluxmag[k]-maxmag)/(12-maxmag))
        color  = colormag[k]
        colors = colormags[:,k]
        '''
        r = max([min([1 + 0.7*color,1]),0])
        b = max([min([1 - 0.7*color,1]),0])
        g = 1 - 0.5*(1-min([r,b]))
        '''

        r = Rfluxsort[k]
        g = Gfluxsort[k]
        b = Bfluxsort[k]
        rgb = max([r,g,b])
        r = r/rgb
        g = g/rgb
        b = b/rgb

        r = bright*r
        g = bright*g
        b = bright*b
        r = max([min([r,1]),0])
        g = max([min([g,1]),0])
        b = max([min([b,1]),0])
        a = min([alph,1])
        marksize = (0.0001+diam*2)    
        '''
        & ( max( log(1e-30+abs(Lsfluxmag[0:N-1,k] -  Lsfluxmag[1:N,k])) \
        - ( ab[0:N-1,1] +  ab[0:N-1,0]*( Lsfluxmag[0:N-1,k]+Lsfluxmag[1:N,k])/2) ) \
        > z*residual[0:N-1]  )\
        '''
        if showit == 'merge':
            '''
            if (k > n-10):
                print('hello 1: ',mean(abs( Lsfluxmag[0:N-1,k] - Lsfluxmag[1:N,k])))
                print('hello 2: ',mean( abs(colors[0:N-1]  - colors[1:N])))
                print('hello 3: ',min( log(1e-30+abs(Lsfluxmag[0:N-1,k] -  Lsfluxmag[1:N,k])) \
                - ( ab[0:N-1,1] +  ab[0:N-1,0]*( Lsfluxmag[0:N-1,k]+Lsfluxmag[1:N,k])/2) \
                - z*residual[0:N-1] ))
            '''

            if  ( (max(Lsfluxmag[:,k]) - min(Lsfluxmag[:,k])) > fluxthresh ) \
                & ( max( log(1e-30+abs(Lsfluxmag[0:N-1,k] -  Lsfluxmag[1:N,k])) \
                - ( ab[0:N-1,1] +  ab[0:N-1,0]*( Lsfluxmag[0:N-1,k]+Lsfluxmag[1:N,k] )/2 )  \
                - z*residual[0:N-1]*(mean(Lsfluxmag[:,k])/10) ) > 0 ) \
                & (max(colors[:])  - min(colors[:]) < colorthresh ):
                colors2 = []
                Lsfluxmag2 = []
                '''
                for jj in range(N):
                        colors2 = append(colors2, mean(colors))
                        colors2 = append(colors2, colors[jj])
                        Lsfluxmag2 = append(Lsfluxmag2, mean(Lsfluxmag[:,k]))
                        Lsfluxmag2 = append(Lsfluxmag2, Lsfluxmag[jj,k])
                plot(colors2, Lsfluxmag2, 'o-' \
                '''
                print('before: ', colors, Lsfluxmag[:,k])
                angles = arctan2(colors-mean(colors),Lsfluxmag[:,k]-mean(Lsfluxmag[:,k]))
                print('after: ', colors, Lsfluxmag[:,k])
                print('angles: ', angles)
                ii = argsort(angles)
                colors2 = colors[ii]
                Lsfluxmag2 = Lsfluxmag[ii,k]
                print('before: ', Lsfluxmag[:,k])
                print('after: ', Lsfluxmag2[:])
                plot(append(colors2[:], colors2[0]), append(Lsfluxmag2[:], Lsfluxmag2[0]), 'o-' \
#               plot(append(colors[:], colors[0]), append(Lsfluxmag[:,k], Lsfluxmag[0,k]), 'o-' \
                    ,color=(r, g, b, a) \
                    ,markerfacecolor=(r, g, b, a) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b, a) \
                        )
            else:
                plot(mean(colors), mean(Lsfluxmag[:,k]), 'o-' \
                    ,markerfacecolor=(r, g, b, a) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b, a) \
                    )
        elif showit == 'nolines':
            if  ( (mean(abs( Lsfluxmag[0:N-1,k] - Lsfluxmag[1:N,k])) > fluxthresh ) )\
                & ( min( log(1e-30+abs(Lsfluxmag[0:N-1,k] -  Lsfluxmag[1:N,k])) \
                - ( ab[0:N-1,1] +  ab[0:N-1,0]*( Lsfluxmag[0:N-1,k]+Lsfluxmag[1:N,k])/2) \
                - z*residual[0:N-1] ) > 0 )\
                & ( mean( abs(colors[0:N-1]  - colors[1:N]) ) < colorthresh ):
                plot(colors, Lsfluxmag[k], 'o' \
                    ,color=(r, g, b) \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                        )
            else:
                plot(mean(colors), mean(Lsfluxmag[k]), 'o' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
        else:
            for j in arange(0,N):
                plot(colors[j,:], Lsfluxmag[:,k], 'o-' \
                    ,markerfacecolor=(r, g, b) \
                    ,markersize=marksize \
                    ,markeredgecolor=(r, g, b) \
                    )
    axis([-2, 2.1, minmag0, minmag0+8])
    xlabel('Color index (B-R)')
    ylabel('Brightness (apparent magnitude)')
    title('Hertzsprung-Russell Diagram for '+objname, color='white')
#hold(False)
    gca().invert_yaxis()
    show()

    (m,n) = R.shape

    '''
    R_bkg_median  = medfilt2d(R , (11,11));
    G_bkg_median  = medfilt2d(G , (11,11));
    B_bkg_median  = medfilt2d(B , (11,11));
    R2_bkg_median = medfilt2d(R2, (11,11));
    G2_bkg_median = medfilt2d(G2, (11,11));
    B2_bkg_median = medfilt2d(B2, (11,11));
    R3_bkg_median = medfilt2d(R3, (11,11));
    G3_bkg_median = medfilt2d(G3, (11,11));
    B3_bkg_median = medfilt2d(B3, (11,11));

    R2 = R2 - R2_bkg_median + R_bkg_median;
    G2 = G2 - G2_bkg_median + G_bkg_median;
    B2 = B2 - B2_bkg_median + B_bkg_median;
    R3 = R3 - R3_bkg_median + R_bkg_median;
    G3 = G3 - G3_bkg_median + G_bkg_median;
    B3 = B3 - B3_bkg_median + B_bkg_median;

    R =R -median(R );  R =maximum(R ,0.);  R =R /R.max();  logR  = log(255*R +1);              
    G =G -median(G );  G =maximum(G ,0.);  G =G /G.max();  logG  = log(255*G +1);              
    B =B -median(B );  B =maximum(B ,0.);  B =B /B.max();  logB  = log(255*B +1);
    R2=R2-median(R2);  R2=maximum(R2,0.);  R2=R2/R2.max(); logR2 = log(255*R2+1);            
    G2=G2-median(G2);  G2=maximum(G2,0.);  G2=G2/G2.max(); logG2 = log(255*G2+1);            
    B2=B2-median(B2);  B2=maximum(B2,0.);  B2=B2/B2.max(); logB2 = log(255*B2+1);
    R3=R3-median(R3);  R3=maximum(R3,0.);  R3=R3/R3.max(); logR3 = log(255*R3+1);            
    G3=G3-median(G3);  G3=maximum(G3,0.);  G3=G3/G3.max(); logG3 = log(255*G3+1);            
    B3=B3-median(B3);  B3=maximum(B3,0.);  B3=B3/B3.max(); logB3 = log(255*B3+1);

#    logR = logR - median(logR); logG = logG - median(logG); logB = logB - median(logB);
#    logR = logR+1;                logG = logG+1;                logB = logB+1;
    p = 0.9; q = 1.0-p;
    logR = logR/(q*percentile(logR,99)+p*logR.max());   
    logG = logG/(q*percentile(logG,99)+p*logG.max());   
    logB = logB/(q*percentile(logB,99)+p*logB.max());
    logR = minimum(logR,1.0); logG = minimum(logG,1.0); logB = minimum(logB,1.0)
    Lrgb = zeros((m,n,3), dtype=float32)
    Lrgb[:,:,0] = logR;           Lrgb[:,:,1] = logG;           Lrgb[:,:,2] = logB;

#    logR2 = logR2 - median(logR2); logG2 = logG2 - median(logG2); logB2 = logB2 - median(logB2);
#    logR2 = logR2+1;                 logG2 = logG2+1;                 logB2 = logB2+1;
    logR2 = logR2/(q*percentile(logR2,99)+p*logR2.max());
    logG2 = logG2/(q*percentile(logG2,99)+p*logG2.max());
    logB2 = logB2/(q*percentile(logB2,99)+p*logB2.max());
    logR2 = minimum(logR2,1.0); logG2 = minimum(logG2,1.0); logB2 = minimum(logB2,1.0)
    Lrgb2 = zeros((m,n,3), dtype=float32)
    Lrgb2[:,:,0] = logR2;            Lrgb2[:,:,1] = logG2;            Lrgb2[:,:,2] = logB2;

    logR3 = logR3/(q*percentile(logR3,99)+p*logR3.max());
    logG3 = logG3/(q*percentile(logG3,99)+p*logG3.max());
    logB3 = logB3/(q*percentile(logB3,99)+p*logB3.max());
    logR3 = minimum(logR3,1.0); logG3 = minimum(logG3,1.0); logB3 = minimum(logB3,1.0)
    Lrgb3 = zeros((m,n,3), dtype=float32)
    Lrgb3[:,:,0] = logR3;            Lrgb3[:,:,1] = logG3;            Lrgb3[:,:,2] = logB3;

    hFig23 = figure(23, figsize=(16,12));
    for iter in arange(0,5):
        imshow(Lrgb);  show(); pause(0.01);
        imshow(Lrgb2); show(); pause(0.01);
        imshow(Lrgb3); show(); pause(0.01);

    RGB [:,:,0] = R ; RGB [:,:,1] = G ; RGB [:,:,2] = B ;
    RGB2[:,:,0] = R2; RGB2[:,:,1] = G2; RGB2[:,:,2] = B2;
    RGB3[:,:,0] = R3; RGB3[:,:,1] = G3; RGB3[:,:,2] = B3;
    z = (65535*((RGB - RGB.min())/RGB.ptp())).astype(uint16)
    '''

def gnomonic(ra,dec,ra_mid,dec_mid):

    # set angle theta so that 0 is in the middle
    theta = (ra-ra_mid)*pi/12;
    phi = (90-dec)*pi/180;
#    print(theta*180/pi, phi*180/pi)

    # compute cartesian coords
    x0 = sin(phi)*cos(theta);
    y0 = sin(phi)*sin(theta);
    z0 = cos(phi);
#    print(x0, y0, z0)

    # rotate on y-axis so that z=0 at middle
    theta_mid = dec_mid*pi/180;
    x1 =  x0*cos(theta_mid)+z0*sin(theta_mid);
    y1 =  y0;
    z1 = -x0*sin(theta_mid)+z0*cos(theta_mid);
#    print(x1, y1, z1)

    # project onto tangent plane at (1,0,0)
    x = x1/abs(x1);
    y = y1/abs(x1);
    z = z1/abs(x1);
#    print(x, y, z)

#    print(x,y,z)
    return x,y,z;

def inv_gnomonic(y,z,ra_mid,dec_mid):

    x = 1
    r = sqrt(x*x + y*y + z*z)
    x = x/r
    y = y/r
    z = z/r

    theta_mid = dec_mid*pi/180;
    x0 = x*cos(theta_mid) - z*sin(theta_mid)
    y0 = y
    z0 = x*sin(theta_mid) + z*cos(theta_mid)

    phi   = arccos (z0)
    theta = arctan2(y0,x0)

    ra = theta*12/pi + ra_mid
    dec = 90 - phi*180/pi

    return ra, dec;

def equatorial_to_galactic(ra,dec):
 
    ra_G = 192.85948*pi/180
    dec_G = 27.12825*pi/180
    l_NCP = 122.93192*pi/180

    ra_rad  = ra *pi/12
    dec_rad = dec*pi/180

    sin_l_cos_b = sin(dec_rad)*cos(dec_G) - cos(dec_rad)*sin(dec_G)*cos(ra_rad-ra_G)
    cos_l_cos_b =                           cos(dec_rad)*           sin(ra_rad-ra_G)
    sin_b       = sin(dec_rad)*sin(dec_G) + cos(dec_rad)*cos(dec_G)*cos(ra_rad-ra_G)

    cos_b = sqrt( sin_l_cos_b**2 + cos_l_cos_b**2 )
    b = arctan2( sin_b, cos_b )
    l = arctan2( sin_l_cos_b, cos_l_cos_b)
    l = mod( l + l_NCP, 2*pi )

    l = l*180/pi - 90
    b = b*180/pi
    
    return l,b

def fwhm(L,i,j):
    global tab
#   print(tab,'Estimating FWHM')
    tab = tab+'    '
#   print(tab, 'i,j = ', i,j)

    (m,n) = L.shape
    k0 = 0
    k1 = 0
    k2 = 0
    k3 = 0
    width = np.min([i,j,m-1-i,n-1-j,20])
#   print(tab, 'width = ', width)
    bkg = percentile(L[i-width:i+width,j-width:j+width],50)
    v0 = L[i,j]-bkg
    for k in arange(i+1,minimum(i+width,m-1)):
        if (L[k,j]-bkg < v0/2): 
            v1 = L[k-1,j]-bkg
            v2 = L[k  ,j]-bkg
            if (v2 != v1):
                k0 = k-(i+1) + (v0/2-v1)/(v2-v1)
            else:
                k0 = 100
            break
    for k in arange(i-1,maximum(i-width,0),-1):
        if (L[k,j]-bkg < v0/2): 
            v1 = L[k+1,j]-bkg
            v2 = L[k  ,j]-bkg
            if (v2 != v1):
                k1 = (i-1)-k + (v0/2-v1)/(v2-v1)
            else:
                k1 = 100
            break
    for k in arange(j+1,minimum(j+width,n)):
        if (L[i,k]-bkg < v0/2): 
            v1 = L[i,k-1]-bkg
            v2 = L[i,k  ]-bkg
            if (v2 != v1):
                k2 = k-(j+1) + (v0/2-v1)/(v2-v1)
            else:
                k2 = 100
            break
    for k in arange(j-1,maximum(j-width,0),-1):
        if (L[i,k]-bkg < v0/2): 
            v1 = L[i,k+1]-bkg
            v2 = L[i,k  ]-bkg
            if (v2 != v1):
                k3 = (j-1)-k + (v0/2-v1)/(v2-v1)
            else:
                k3 = 100
            break
#   print('i,j = ',i,j,',   fwhm = ', ArcsecPerPixel,'*',(k0+k1+k2+k3)/2)
    tab = tab[0:-4]
    return ArcsecPerPixel*(k0+k1+k2+k3)/2

def ErosionKF(Image, numpx):
    if (len(Image.shape) == 2):
        (m,n) = Image.shape
        return ErosionKF2(Image, numpx)
    elif (len(Image.shape) == 3):
        (m,n,three) = Image.shape
        R = ErosionKF2(Image[:,:,0], numpx)
        G = ErosionKF2(Image[:,:,1], numpx)
        B = ErosionKF2(Image[:,:,2], numpx)
        Image[:,:,0] = R
        Image[:,:,1] = G
        Image[:,:,2] = B
        return Image

def ErosionKF2(Image, numpx):
    (m,n) = Image.shape
    Im = zeros((m,n), dtype=float32)
    d = int(numpx/2)
    I0 = Image[d:m-d, d:n-d]
    for i in arange(-d,d+1):
        for j in arange(-d,d+1):
            I0 = np.minimum(I0, Image[d+i:m-d+i, d+j:n-d+j])
#   Im[d:m-d, d:n-d] = I0
    Im = pad(I0,d)
    return Im

def GaussianBlur(Image, r):
    M = zeros((21,21), dtype=float32)
    for i in arange(0,21):
        for j in arange(0,21):
            M[i,j] = (1/(2*pi*r*r))*exp(-((i-10)**2 + (j-10)**2)/(2*r*r))
    return convolve2d(pad(Image,10), M, mode='valid')

def NeatImageCore(Image, StdDev, Blur):
    background = ErosionKF(Image,7)
    background = GaussianBlur(background,3)
    background = background+5*StdDev
    # figure(40+i); imshow(background)

    diff_pos_part = maximum(Image - background, 0)
    # figure(50+i); imshow(diff_pos_part)

    diff_neg_part = maximum(diff_pos_part + background - Image, 0)
    diff_neg_part = GaussianBlur(diff_neg_part, Blur)
    # figure(60+i); imshow(diff_neg_part)

    result = maximum(diff_pos_part + background - diff_neg_part, 0)
    # figure(70+i); imshow(result)
    return result

def NeatImage(Image, StdDev, Blur):
    if (len(Image.shape) == 2):
        result = NeatImageCore(Image[:,:], StdDev, Blur)
    elif (len(Image.shape) == 3):
        result = zeros(Image.shape, dtype=float32)
        for i in arange(0,3):
            result[:,:,i] = NeatImageCore(Image[:,:,i], StdDev, Blur)
    return result

def Denoise(im, thresh):
    [m, n, three] = shape(im)
    im_denoised = zeros((m,n,three), dtype=float32)
    for j in range(three):
        im_fft = fftpack.fft2(im[:,:,j]/256)
        thresh = percentile(im_fft.real,50)
        k = (im_fft.real**2 + im_fft.imag**2 < thresh**2)
        im_fft[k] = 0
#       im_fft[m-int(round(thresh)):m,                    : ] = 0
#       im_fft[                    : ,n-int(round(thresh)):n] = 0
        z = fftpack.ifft2(im_fft)
        im_denoised[:,:,j] = 256*sqrt(z.real**2 + z.imag**2)
    return im_denoised

def NeatImageFFTCore(Image, StdDev, Blur):
    background = ErosionKF(Image,7)
    background = GaussianBlur(background,3)
    background = background+5*StdDev

    diff_pos_part = maximum(Image - background, 0)

    diff_neg_part = maximum(diff_pos_part + background - Image, 0)
#   diff_neg_part = GaussianBlur(diff_neg_part, Blur)
    im = diff_neg_part
    im_fft = fftpack.fft2(im)
    im_fft2 = im_fft.copy()
    r, c = im_fft2.shape
    xs =        arange(0,c)/c
    ys = array([arange(0,r)]).T/r
    r0 = minimum(r/2,c/2)
    r0 = 1./Blur
    radii00 = sqrt((ones((r,1))*xs)**2 + (ys*ones((1,c)))**2)
    radii01 = radii00[ :  ,::-1]
    radii11 = radii00[::-1,::-1]
    radii10 = radii00[::-1, :  ]

    disk = float32( (radii00 <= r0) | \
                    (radii01 <= r0) | \
                    (radii11 <= r0) | \
                    (radii10 <= r0) )
    im_fft2 = im_fft2*disk
    im_new = fftpack.ifft2(im_fft2).real
    im_new = maximum(im_new, np.min(im))
    diff_neg_part = im_new

    result = maximum(diff_pos_part + background - diff_neg_part, 0)
    return result

def NeatImageFFT(Image, StdDev, Blur):
    if (len(Image.shape) == 2):
        result = NeatImageFFTCore(Image[:,:], StdDev, Blur)
    elif (len(Image.shape) == 3):
        result = zeros(Image.shape, dtype=float32)
        for i in arange(0,3):
            result[:,:,i] = NeatImageFFTCore(Image[:,:,i], StdDev, Blur)
    return result

'''
def NeatImageFFT(Image, StdDev, Blur):
    if (len(Image.shape) == 2):
        (m,n) = Image.shape
        background = ErosionKF(Image[:,:],7)
        background = GaussianBlur(background,3)
        background = background+5*StdDev

        diff_pos_part = maximum(Image[:,:] - background, 0)

        diff_neg_part = maximum(diff_pos_part + background - Image[:,:], 0)
        diff_neg_part = GaussianBlur(diff_neg_part, Blur)

        result = maximum(diff_pos_part + background - diff_neg_part, 0)
    elif (len(Image.shape) == 3):
        (m,n,three) = Image.shape
        result = zeros((m,n,3), dtype=float32)
        for i in arange(0,3):
            background = ErosionKF(Image[:,:,i],7)
            background = GaussianBlur(background,3)
            background = background+5*StdDev
            # figure(40+i); imshow(background)

            diff_pos_part = maximum(Image[:,:,i] - background, 0)
            # figure(50+i); imshow(diff_pos_part)

            diff_neg_part = maximum(diff_pos_part + background - Image[:,:,i], 0)
            diff_neg_part = GaussianBlur(diff_neg_part, Blur)
            # figure(60+i); imshow(diff_neg_part)

            result[:,:,i] = maximum(diff_pos_part + background - diff_neg_part, 0)
            # figure(70+i); imshow(result[:,:,i])
    return result
'''

def ImageSmooth(*arg):
    image = arg[0]
    medmean = arg[3]     if len(arg)>3 else 'median';
    r = int(arg[2])      if len(arg)>2 else 2;
    thresh_pctl = arg[1] if len(arg)>1 else 80;
    
    if   (len(image.shape) == 2):
        smooth_image = ImageSmoothCore(image,thresh_pctl,r,medmean)
    elif (len(image.shape) == 3):
        R = ImageSmoothCore(image[:,:,0],thresh_pctl,r,medmean)
        G = ImageSmoothCore(image[:,:,1],thresh_pctl,r,medmean)
        B = ImageSmoothCore(image[:,:,2],thresh_pctl,r,medmean)
        
        smooth_image = zeros(image.shape, float32)
        smooth_image[:,:,0] = R
        smooth_image[:,:,1] = G
        smooth_image[:,:,2] = B

    return smooth_image

def ImageSmoothCore(I,thresh_pctl,r,medmean):
    d = 2*r+1
    (m,n) = I.shape
    thresh = percentile(I,thresh_pctl)
    Istack = zeros((m-d,n-d,d,d), dtype=float32)
    for i in arange(-r,r+1):
        for j in arange(-r,r+1):
            Istack[:,:,r+i,r+j] = I[r+i:m-r-1+i,r+j:n-r-1+j]
    
    Imin = percentile(I,5)
    Imedian = zeros((m,n), dtype=float32)
    Imean   = zeros((m,n), dtype=float32)
    Imean2  = zeros((m,n), dtype=float32)
    Istd    = zeros((m,n), dtype=float32)
    Imedian[:,:] = I
    Imean[:,:] = I
    Imean2[:,:] = I*I
    Imedian[r:m-r-1,r:n-r-1] = median(Istack,(2,3))
    Imean[r:m-r-1,r:n-r-1]   = mean(Istack,(2,3))
# Imean2[r:m-r-1,r:n-r-1]  = mean(Istack*Istack,(2,3))
# Istd = sqrt(abs(Imean2 - Imean*Imean))
# lo_pxs = (I <= thresh) & (Istd < 40*sqrt(maximum(Imean-Imin,0)))
    lo_pxs = (I <= thresh)
    hi_pxs = ~lo_pxs
    if   (medmean == 'median'):
        I = I*hi_pxs + Imedian*lo_pxs
    elif (medmean ==   'mean'):
        I = I*hi_pxs + Imean*lo_pxs
    return I

def applyFlat(I, flat, bias=0):
    if len(bias.shape) == 0: 
        if bias == 0:
            bias = percentile(I, 1)
    
    flat = maximum(flat, percentile(flat,0.1))  # ignore dead pixels

    scale = median(flat-bias)

    if len(I.shape) == 2:
        I = bias + scale*(I-bias)/(flat-bias+1e-10)
    else:
        N = I.shape[0]
        meanbias = mean(bias)
        print(tab,'Applying flat to ',N,' images')
        for j in arange(0,N):
            I[j,:,:] = meanbias + scale*(I[j,:,:]-bias)/(flat-bias+1e-10)

def RemoveGradient(I):
    global xmin, xmax, ymin, ymax
    (m,n) = I.shape

    I00 = median(I[ -ymin   :50-ymin,  -xmin   :50-xmin])
    I0n = median(I[ -ymin   :50-ymin, n-xmax-50: n-xmax])
    Im0 = median(I[m-ymax-50: m-ymax,  -xmin   :50-xmin])
    Imn = median(I[m-ymax-50: m-ymax, n-xmax-50: n-xmax])
    Iavg = (I00+I0n+Im0+Imn)/4

    xs = arange(0,n)
    ys = array([arange(0,m)]).T

    px = xs/(n-1)
    py = ys/(m-1)
    qx = 1-px
    qy = 1-py

    B = qy*qx*I00 + py*qx*Im0 + qy*px*I0n + py*px*Imn

    I_bkg_adjust = I - B + Iavg
    return I_bkg_adjust
    
'''
def HDR(RGB):
    RGB = RGB-np.min(RGB)+1
    RGB = RGB/np.max(RGB)
    (m,n,three) = RGB.shape
    R = RGB[:,:,0]
    G = RGB[:,:,1]
    B = RGB[:,:,2]
    L = R + G + B

    R1 = np.power(array(arange(0,m*n))/(m*n),3)
    G1 = R1
    B1 = R1
    L1 = R1

    R2 = reshape(R, m*n)
    G2 = reshape(G, m*n)
    B2 = reshape(B, m*n)
    L2 = reshape(L, m*n)

    M = ones(21)/21

    I1 = argsort(R1)
    I2 = argsort(R2)
    R2[I2] = R2[I2]*convolve(L1[I1]/L2[I2], M, mode='same')

    I1 = argsort(G1)
    I2 = argsort(G2)
    G2[I2] = G2[I2]*convolve(L1[I1]/L2[I2], M, mode='same')

    I1 = argsort(B1)
    I2 = argsort(B2)
    B2[I2] = B2[I2]*convolve(L1[I1]/L2[I2], M, mode='same')

    RGB_HDR = array(zeros((m,n,3), dtype=float32))
    RGB_HDR[:,:,0] = reshape(R2, (m,n))
    RGB_HDR[:,:,1] = reshape(G2, (m,n))
    RGB_HDR[:,:,2] = reshape(B2, (m,n))

    return RGB_HDR
'''

def UnsharpMaskRGB(RGB):
    (m,n,three) = RGB.shape
    R = RGB[:,:,0]
    G = RGB[:,:,1]
    B = RGB[:,:,2]

    print('')
    print('Red frame:') 
    R_UM = UnsharpMask(R)
    print('')
    print('Green frame:') 
    G_UM = UnsharpMask(G)
    print('')
    print('Blue frame:') 
    B_UM = UnsharpMask(B)

    RGB_UM = array(zeros((m,n,3), dtype=float32))
    RGB_UM[:,:,0] = R_UM
    RGB_UM[:,:,1] = G_UM
    RGB_UM[:,:,2] = B_UM

    return RGB_UM

def UnsharpMask(L):
    print('Unsharp Masking...')
    rad = 6
    thresh = 0.1
    (m,n) = L.shape
    border = int(minimum(m,n)/5)
    rad2 = 2*rad+1
    k, LconvL2 = findstars(L,thresh=thresh,bkg=True, fignum=0)
    M = ones((rad2,rad2))
    kconv = convolve2d(k, M, mode='same')

    xs = arange(0,n)
    ys = array([arange(0,m)]).T
    Xs = ones((m,1))*xs
    Ys = ys*ones((1,n))

    mask = (Xs > border) & (Xs < n-border) & (Ys > border) & (Ys < m-border) 
    k = k&mask;

    flux = LconvL2[k]
    Xs = Xs[k]
    Ys = Ys[k]

    I = argsort(-flux)
    flux = flux[I]
    Xs_sort = int16(Xs[I])
    Ys_sort = int16(Ys[I])

    for j in arange(0,len(flux)):
        if ( (flux[j] < 0.5*flux[0]) & (kconv[Ys_sort[j],Xs_sort[j]] == 1) ): break
#       if ( (flux[j] < 60000) & (kconv[Ys_sort[j],Xs_sort[j]] == 1) ): break

    x = Xs_sort[j]
    y = Ys_sort[j]
        
    psf = array(zeros((rad2,rad2), dtype=float32))
    psf[0:rad2,0:rad2] = L[y-rad:y+rad+1,x-rad:x+rad+1]
    psf = psf[::-1,::-1]
    bkg = percentile(psf,100/rad2)
    psf = psf - bkg
    psf = psf/np.sum(psf)

    '''
    result      = array(zeros((m,n), dtype=float32))
    result[:,:] = 2*L - convolve2d(L, psf, mode='same')
    '''
    std_dev = np.median(abs(L[1:,1:]-L[:-1,:-1]))
    background = ErosionKF(L, 7)
    background = GaussianBlur(background, 3) + 5*std_dev
    image_high_vals = maximum(L - background,0)
    image_low_vals  = maximum(image_high_vals + background - L,0)
    tmp_result = 2*image_high_vals - convolve2d(image_high_vals, psf, mode='same')
    tmp_result = tmp_result + background - image_low_vals
    result = maximum(tmp_result,0)

    return result

def RLdeconvRGB(RGB, fwhm_target, thresh=0.1, x=0, y=0, std_dev=5, \
                psf_shape='default', rad=8, show_psf='False'):
    global tab

    (m,n,three) = RGB.shape
    R = RGB[:,:,0]
    G = RGB[:,:,1]
    B = RGB[:,:,2]

    L = (R+G+B)/3
    x,y,fwhm_current = find_fwhm_star(L, rad=rad, thresh=thresh)
    fwhm_current = fwhm_current/ArcsecPerPixel
    rad = int(ceil(2*fwhm_current))
    print(tab,'RLdeconvRGB: x, y, fwhm = ',x,y,fwhm_current)
    tab = tab + '    '

    '''
    if psf_shape == 'fat tail':
        psf = ([[ 0.,  1.,  1.,  1.,  0.], \
                     [ 1.,  2.,  8.,  2.,  1.], \
                     [ 1.,  8., 32.,  8.,  1.], \
                     [ 1.,  2.,  8.,  2.,  1.], \
                     [ 0.,  1.,  1.,  1.,  0.]]) / 2
    '''
    rad2 = 2*rad+1
    print(tab,'x = ',x,', y = ',y,', rad = ',rad)

    psf = array(zeros((rad2,rad2), dtype=float32))
    psf[0:rad2,0:rad2] = R[y-rad:y+rad+1,x-rad:x+rad+1]
    psf = psf[::-1,::-1]
    bkg = percentile(psf,100/rad2)
    psf = maximum(psf - bkg,0)
    psf = uint8(255*psf/np.max(psf))
    if (show_psf == True): print(tab,psf)

    print('')
    print(tab,'Red frame:') 
#   R_RL = RLdeconvGrid(R, fwhm_target, thresh=thresh, x=x, y=y, std_dev=std_dev, \
    R_RL = RLdeconv(R, fwhm_target, psf, thresh=thresh, x=x, y=y, std_dev=std_dev, \
                psf_shape=psf_shape, rad=rad, show_psf=show_psf)

    '''
    '''
    psf = array(zeros((rad2,rad2), dtype=float32))
    psf[0:rad2,0:rad2] = G[y-rad:y+rad+1,x-rad:x+rad+1]
    psf = psf[::-1,::-1]
    bkg = percentile(psf,100/rad2)
    psf = maximum(psf - bkg,0)
    psf = uint8(255*psf/np.max(psf))
    '''
    '''
    print('')
    print(tab,'Green frame:') 
    if (np.max(abs(G-R)) == 0):
        print(tab,'  same as red frame')
        G_RL = R_RL
    else:
#       G_RL = RLdeconvGrid(G, fwhm_target, thresh=thresh, x=x, y=y, std_dev=std_dev, \
        G_RL = RLdeconv(G, fwhm_target, psf, thresh=thresh, x=x, y=y, std_dev=std_dev, \
                    psf_shape=psf_shape, rad=rad, show_psf=show_psf)

    '''
    '''
    psf = array(zeros((rad2,rad2), dtype=float32))
    psf[0:rad2,0:rad2] = B[y-rad:y+rad+1,x-rad:x+rad+1]
    psf = psf[::-1,::-1]
    bkg = percentile(psf,100/rad2)
    psf = maximum(psf - bkg,0)
    psf = uint8(255*psf/np.max(psf))
    '''
    '''
    print('')
    print(tab,'Blue frame:') 
    if (np.max(abs(B-R)) == 0):
        print(tab,'  same as red frame')
        B_RL = R_RL
    elif (np.max(abs(B-G)) == 0):
        print(tab,'  same as green frame')
        B_RL = G_RL
    else:
#       B_RL = RLdeconvGrid(B, fwhm_target, thresh=thresh, x=x, y=y, std_dev=std_dev, \
        B_RL = RLdeconv(B, fwhm_target, psf, thresh=thresh, x=x, y=y, std_dev=std_dev, \
                    psf_shape=psf_shape, rad=rad, show_psf=show_psf)

    RGB_RL = array(zeros((m,n,3), dtype=float32))
    RGB_RL[:,:,0] = R_RL
    RGB_RL[:,:,1] = G_RL
    RGB_RL[:,:,2] = B_RL

    #####################################################################################
    # Save color image as a fits file

    RGBfit = zeros((3,m,n), dtype=float32)
    RGBfit[0, 0:m, 0:n] = R_RL
    RGBfit[1, 0:m, 0:n] = G_RL
    RGBfit[2, 0:m, 0:n] = B_RL

    RGBfit += np.max(RGBfit)/20

    hdu = fits.PrimaryHDU(RGBfit)
    hdulist = fits.HDUList([hdu])
    file_location = folder+name+'-py-RL-RGB.fit'
    hdulist.writeto(file_location, overwrite=True)

    #####################################################################################
    # Done

    tab = tab[0:-4]
    printdone('Done deconvolving')

    return RGB_RL

def RLdeconvGrid(L, fwhm_target, grid_size=1, thresh=0.1, x=0, y=0, std_dev=5, \
                psf_shape='default', rad=6, show_psf='False'):
    (m,n) = L.shape
    L0     = array(zeros((m+40,n+40), dtype=float32))
    L0[20:m+20, 20:n+20] = L
    result = array(zeros((m+40,n+40), dtype=float32))
    for i in range(grid_size):
        for j in range(grid_size):
            Lij = L0[int(floor(i*m/grid_size)):int(ceil((i+1)*m/grid_size)+40), \
                     int(floor(j*n/grid_size)):int(ceil((j+1)*n/grid_size)+40)]
            resultij = RLdeconv(Lij, fwhm_target, psf, thresh=thresh/10, x=0, y=0, std_dev=5, \
                                psf_shape='default', rad=6, show_psf='False')
            result[int(floor(i*m/grid_size)):int(ceil((i+1)*m/grid_size)), \
                   int(floor(j*n/grid_size)):int(ceil((j+1)*n/grid_size))] = resultij[20:-20,20:-20]
    return result[20:m+20, 20:n+20]

def find_fwhm_star(L, rad=6, thresh=0.1):
    global tab
    print(tab, 'find_fwhm_star:')
    tab = tab+'    '

    (m,n) = L.shape
    border = int(minimum(m,n)/5)
    rad2 = 2*rad+1
    k, LconvL2 = findstars(L,bkg=True,fignum=0)
    '''
    M = ones((rad2,rad2))
    kconv = convolve2d(k, M, mode='same')
    '''

    xs = arange(0,n)
    ys = array([arange(0,m)]).T
    Xs = ones((m,1))*xs
    Ys = ys*ones((1,n))

    bkgvals = zeros((m,n,28), dtype=float32)
    m0 = border
    m1 = m-border
    n0 = border
    n1 = n-border
    jj = 0
    '''
    for j in arange(-rad,rad+1):
        bkgvals[m0:m1, n0:n1, jj] = L[m0+j:m1+j, n0-rad:n1-rad]; jj += 1
        bkgvals[m0:m1, n0:n1, jj] = L[m0+j:m1+j, n0+rad:n1+rad]; jj += 1
    for j in arange(-rad+1,rad):
        bkgvals[m0:m1, n0:n1, jj] = L[m0-rad:m1-rad, n0+j:n1+j]; jj += 1
        bkgvals[m0:m1, n0:n1, jj] = L[m0+rad:m1+rad, n0+j:n1+j]; jj += 1
    rad1 = rad-1
    bkgvals[m0:m1, n0:n1, jj] = L[m0-rad1:m1-rad1, n0-rad1:n1-rad1]; jj += 1
    bkgvals[m0:m1, n0:n1, jj] = L[m0-rad1:m1-rad1, n0+rad1:n1+rad1]; jj += 1
    bkgvals[m0:m1, n0:n1, jj] = L[m0+rad1:m1+rad1, n0-rad1:n1-rad1]; jj += 1
    bkgvals[m0:m1, n0:n1, jj] = L[m0+rad1:m1+rad1, n0+rad1:n1+rad1]; jj += 1
    '''
    for i in arange(-rad,rad+1):
        for j in arange(-rad,rad+1):
            if ((i*i+j*j>=rad-1) & (i*i+j*j<=rad+1)):
                bkgvals[m0:m1, n0:n1, jj] = L[m0+i:m1+i, n0+j:n1+j]; jj += 1
    noise    = zeros((m,n), dtype=float32)
    noise[k] = percentile(bkgvals[k,:], 90, axis=1) \
             - percentile(bkgvals[k,:], 10, axis=1)

    mask = (Xs > border) & (Xs < n-border) & (Ys > border) & (Ys < m-border) 
    k = k&mask;

    snr = (LconvL2[k]/28)/noise[k]
    Xs = Xs[k]
    Ys = Ys[k]

    I = argsort(-snr)
    snr =  snr[I]
    noise = noise[k][I]
    Xs_sort = int16(Xs[I])
    Ys_sort = int16(Ys[I])
    '''
    print(tab,0, Xs_sort[0], Ys_sort[0], snr[0], noise[0])
    for j in arange(len(I)):
        if ((Xs_sort[j] == 1270) & (Ys_sort[j] == 840)):
            print(tab,j, Xs_sort[j], Ys_sort[j], snr[j], noise[j])
    '''

    x = Xs_sort[0]
    y = Ys_sort[0]
    fwhm_L = fwhm(L,y,x)

    tab = tab[0:-4]
    printdone( 'Done find_fwhm_star')
    return x,y,fwhm_L

def pad(img, rad):
    '''
    (m,n) = img.shape
    padded = array(zeros((m+2*rad,n+2*rad), dtype=float32))
    padded[rad:m+rad, rad:n+rad] = img
    for j in arange(rad):
        padded[rad:m+rad,      j] = img[  :,0]
        padded[rad:m+rad,n+rad+j] = img[  :,n-1]
        padded[      j,rad:n+rad] = img[  0,:]
        padded[m+rad+j,rad:n+rad] = img[m-1,:]

        padded[0:  rad,0:  rad] = img[  0,  0]
        padded[0:  rad,n:n+rad] = img[  0,n-1]
        padded[m:m+rad,0:  rad] = img[m-1,  0]
        padded[m:m+rad,n:n+rad] = img[m-1,n-1]
    return padded
    '''
    return np.pad(img, rad, mode='edge')


def RLdeconv(L, fwhm_target, psf, thresh=0.1, x=0, y=0, std_dev=5, \
                psf_shape='default', rad=6, show_psf='False'):
    global tab
    print(tab,'Richardson-Lucy deconvolution...')
    tab = tab+'    '
    (m,n) = L.shape
    border = int(minimum(m,n)/5)
    rad2 = 2*rad+1

    std_dev = np.median(abs(L[1:,1:]-L[:-1,:-1]))
    print(tab,'std_dev = ', std_dev)
        
    result      = array(zeros((m,n), dtype=float32))
    tmp_result  = array(zeros((m,n), dtype=float32))
    numerator   = array(zeros((m,n), dtype=float32))
    denominator = array(zeros((m,n), dtype=float32))

    result[:,:] = L
    image_num = 0

    '''
    figure(100, figsize=(12,8))
    imshow(L, cmap='gray', norm=LogNorm())
    title('unpadded array')
    show()
    '''
    '''
    figure(101, figsize=(12,8))
    imshow(pad(L,rad), cmap='gray', norm=LogNorm())
    title('padded array')
    show()
    '''

    fwhm_current = fwhm(result,y,x)

    while fwhm_current > fwhm_target:

        background = ErosionKF(result, 7)
        '''
        if (image_num == 0):
            figure(102, figsize=(12,8))
            imshow(background, cmap='gray', norm=LogNorm())
            title('array after erosion')
            show()
        '''
        background = GaussianBlur(background, 3) + 5*std_dev
        '''
        if (image_num == 0):
            figure(103, figsize=(12,8))
            imshow(background, cmap='gray', norm=LogNorm())
            title('array after blur')
            show()
        '''

        image_high_vals = maximum(result - background,0)
        image_low_vals  = maximum(image_high_vals + background - result,0)

        tmp_result[:,:] = image_high_vals

        numerator[:,:]  = image_high_vals

        denominator[:,:] = tmp_result

        fwhm_current = fwhm(result,y,x)

        if psf_shape == 'fat tail':
            nKern = int(fwhm_current**2 / 1.8548)
        else:
            nKern = 1

        if (image_num == 0): print(tab,'x,y = ', x, y)
        print(tab,'  Iteration: ', image_num, ', FWHM: ', fwhm_current)
        image_num += 1

        for i in arange(0,nKern):
            denominator = convolve2d(pad(denominator,rad), psf, mode='valid')

        '''
        figure(100+10*(image_num-1)  , figsize=(16,12))
        imshow(result[y-20:y+20,x-20:x+20]/np.max(result[y-20:y+20,x-20:x+20]), cmap='gray', norm=LogNorm())
        title('result')
        show()
        '''

        numerator[denominator>0] = numerator[denominator>0]/denominator[denominator>0]

#       psf = psf[::-1,::-1]
        for iter in arange(0,nKern):
            numerator   = convolve2d(pad(numerator,rad),   psf, mode='valid')
            '''
            # this is much much slower
            (m0, n0) = psf.shape
            m02 = int(floor(m0/2))
            n02 = int(floor(n0/2))
            tmp = array(zeros((m+m0+1, n+n0+1), dtype=float32))
            tmp[m02:m+m02, n02:n+n02] = numerator
            for i in range(m):
                for j in range(n):
                    for i0 in range(m0):
                        for j0 in range(n0):
                            numerator[i,j] += tmp[i-i0,j-j0]*psf[i0,j0]
            '''

        tmp_result = tmp_result*numerator + background - image_low_vals
        tmp_result = maximum(tmp_result,0)
        result[:,:] = tmp_result

    tab = tab[0:-4]
    printdone('Done doing Richardson-Lucy deconvolution')
    return result

def decimal_to_min_sec(angle):
    if (angle < 0): 
            s = '-'
            angle = -angle
    else:           
            s =  '+'
    angle += 1e-6
    angle_deg =    (floor( angle                                   ))
    angle_min =    (floor((angle - angle_deg                )*  60.))
    angle_sec =    (      (angle - angle_deg - angle_min/60.)*3600. )

    angle_sec = round(100*angle_sec)/100
#   return (s +str(angle_deg)+' : '+str(angle_min)+' : '+str(angle_sec))
    
    txt = '{x1:1s}{x2:02.0f} : {x3:02.0f} : {x4:02.0f}'
    return txt.format(x1=s,x2=angle_deg,x3=angle_min,x4=angle_sec)

def min_sec_to_decimal(s, hr, min, sec):
    if (s == '-'): sgn = -1
    else:          sgn =  1
    return sgn*(hr + (min + sec/60.)/60.)

def fits_astrometry(file_location, approx_ra_dec=[0,0], approx_fov=0.5, \
                thresh=0.02, showAlign=False):
    global slash
    
    image = fits.getdata(file_location).astype(float32)
    image = array(image)
    if (len(shape(image))==3):
        R = image[0,:,:]
        G = image[1,:,:]
        B = image[2,:,:]
        L = R+G+B
    else:
        L = image
    
    #####################################################################################
    # Let's just look at the fits header info...
    
    hdul = fits.open(file_location, mode='update')   # FYI... HDUL = Header Data Unit List
    '''
    print('Info:')
    hdul.info()
    print('\n')
    
    print('Header Cards: ')
    print(repr(hdul[0].header))
    print('\n')
    '''
    
    #####################################################################################
    # Let's do the astrometry
    
    (ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it, success) = \
            astrometry(L, approx_ra_dec=approx_ra_dec, approx_fov=approx_fov, \
                            showAlign=showAlign)
    
    (s,  ra_hr,   ra_min,  ra_sec) = decimal_to_min_sec( ra_ctr)
    (s, dec_deg, dec_min, dec_sec) = decimal_to_min_sec(dec_ctr)

    print(tab,'3:  dec_ctr: ', dec_ctr)

    print('\n')
    print(tab,'Output: ')
    print('  Ra, Dec at center = %3d:%02d:%05.2f, %1s%3d:%02d:%05.2f' %
            (ra_hr, ra_min, ra_sec, s, dec_deg, dec_min, dec_sec))
    print('  Angular rotation  = %7.3f deg' %
            (dtheta))
    print('  Pixel scale    =    %7.3f arcsec/pix' %
            (scale))
    print(' ')
    
    print('4:  dec_ctr: ', dec_deg, dec_min, dec_sec)

    if (success==False): return
    
    #####################################################################################
    # Add astrometry info to fits file
    
    ra = 15*ra_ctr
    dec = dec_ctr
    if (len(image.shape)==2): (m,n) = image.shape
    else:               (three,m,n) = image.shape
    crota2 = -dtheta*pi/180
    cd001001 =  cos(crota2) * -1
    cd001002 =  sin(crota2) * -1
    cd002001 = -sin(crota2) * -1
    cd002002 =  cos(crota2) * -1
    cd1_1 = cd001001*scale/3600
    cd1_2 = cd001002*scale/3600
    cd2_1 = cd002001*scale/3600
    cd2_2 = cd002002*scale/3600

    '''
    crpix1=n/2
    crpix2=m/2
    crval1=ra
    crval2=dec

    print('2: ra,dec = ', ra, dec)
    print('2: ra,dec = ', crval1, crval2)

    crpix1=579.2277
    crpix2=580.2798
    crval1=318.603943485
    crval2=46.2853651935
    '''


    
#   hdul[0].header.clear()
#   hdul[0].header['BITPIX'] = -32  
    hdul[0].header['CTYPE1'] = 'RA---TAN'  # / TAN (gnomic) projection'
    hdul[0].header['CTYPE2'] = 'DEC--TAN'  # / TAN (gnomic) projection'
    hdul[0].header['EQUINOX'] = 2000.0     # / Equatorial coordinates definition (yr)'
    hdul[0].header['CRVAL1'] = ra     # +'         / RA  of reference point'
    hdul[0].header['CRVAL2'] = dec    # +'        / DEC of reference point'
    hdul[0].header['CRPIX1'] = crpix1 # +'        / X reference pixel'
    hdul[0].header['CRPIX2'] = crpix2 # +'        / Y reference pixel'
    hdul[0].header['CUNIT1'] = 'deg     ' # / X pixel scale units'
    hdul[0].header['CUNIT2'] = 'deg     ' # / Y pixel scale units'
    hdul[0].header['CD1_1'] = cd1_1   # +'      / Transformation matrix'
    hdul[0].header['CD1_2'] = cd1_2   # +'      / no comment'
    hdul[0].header['CD2_1'] = cd2_1   # +'      / no comment'
    hdul[0].header['CD2_2'] = cd2_2   # +'      / no comment'
#   hdul[0].header['XPIXSZ'] = scale   # +'      / no comment'
#   hdul[0].header['YPIXSZ'] = scale   # +'      / no comment'
#   hdul[0].header.remove('OBJCTRA',  ignore_missing=True)
#   hdul[0].header.remove('OBJCTDEC', ignore_missing=True)
#   hdul[0].header.remove('OBJCTHA', ignore_missing=True)
    hdul[0].header['IMAGEW'] = n   # +'      / no comment'
    hdul[0].header['IMAGEH'] = m   # +'      / no comment'
    
    image = maximum(image-percentile(image,10),1)
    image = 100*image/np.max(image) + 1
    image = log(image)
    image = 65535*image/np.max(image)
    
    idx = file_location.rfind(slash)
    file_location = file_location[idx+1:]
#   out_file_location = file_location[:-4] + '-astrometry.fit'
#   out_file_location = 'PyWWT'+slash+out_file_location
#   fits.writeto(out_file_location, image, hdul[0].header, overwrite=True)
#   hdul.writeto(out_file_location, image, hdul[0].header, overwrite=True)
    hdul.flush()
#   hdul.close()
    return

global figdumb
figdumb = 400

def astrometry(image, approx_ra_dec=[0,0], approx_fov=0.5, \
                thresh=0.02, showAlign=False):

    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia

#   global xbar, ybar, xbar0, ybar0, xbar1, ybar1, star, align_stack, img, fig
    global tab, lbl6a
    print(tab, 'Astrometry:')
#   tab = tab+'    '

    global figdumb

    (m,n) = image.shape

    # use a shorthand for the image
    L = image

    # find brightest stars
    k, LconvL2 = findstars(L,thresh=thresh,fignum=0,bkg=True)

    xs = arange(0,n)
    ys = array([arange(0,m)]).T
    Xs = ones((m,1))*xs
    Ys = ys*ones((1,n))

    border = 10
    mask = (Xs > border) & (Xs < n-border) & (Ys > border) & (Ys < m-border) 
    k = k&mask;

    flux = LconvL2[k]
    Xs = Xs[k]
    Ys = Ys[k]

    I = argsort(-flux)
    flux = flux[I]
    Xs_sort = int32(Xs[I])
    Ys_sort = int32(Ys[I])

    print(tab,'number of stars found in image = ', len(Xs))
    if (len(Xs)>300):
        thresh=2*thresh
        (ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it, success) = \
            astrometry(image, approx_ra_dec=approx_ra_dec, approx_fov=approx_fov, \
                thresh=thresh, showAlign=showAlign)
        return ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it, success
    num_stars = min([15, len(Xs)])
    if num_stars < 15: 
        thresh = thresh/1.5
        print(tab,'3: reducing threshold to ', thresh)
        (ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it, success) = \
            astrometry(image, approx_ra_dec=approx_ra_dec, approx_fov=approx_fov, \
                        thresh=thresh, showAlign=showAlign)
        return ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it, success

    approx_ra_dec[0] *= 15

    thresh = flux[num_stars-1]/np.max(LconvL2)

    num_stars3 = int(num_stars**3)
    distdata = zeros((num_stars3, 5))
    maxflux = np.max(LconvL2)
    k = 0
    for i0 in arange(0, num_stars):
        if flux[i0] > 0.9*maxflux: continue
        for i1 in arange(i0+1, num_stars):
            if flux[i1] > 0.9*maxflux: continue
            for i2 in arange(i1+1, num_stars):
                    if flux[i2] > 0.9*maxflux: continue
                    distdata[k,0] = i0
                    distdata[k,1] = i1
                    distdata[k,2] = i2
                    distdata[k,4] = sqrt( (Xs_sort[i0]-Xs_sort[i1])**2 + \
                                          (Ys_sort[i0]-Ys_sort[i1])**2 ) \
                                          *flux[i0]*flux[i1] \
                                  + sqrt( (Xs_sort[i1]-Xs_sort[i2])**2 + \
                                          (Ys_sort[i1]-Ys_sort[i2])**2 ) \
                                          *flux[i1]*flux[i2] \
                                  + sqrt( (Xs_sort[i2]-Xs_sort[i0])**2 + \
                                          (Ys_sort[i2]-Ys_sort[i0])**2 ) \
                                          *flux[i2]*flux[i0] 

                    theta01 = arctan2(Ys_sort[i1]-Ys_sort[i0],Xs_sort[i1]-Xs_sort[i0]) 
                    theta12 = arctan2(Ys_sort[i2]-Ys_sort[i1],Xs_sort[i2]-Xs_sort[i1]) 
                    theta20 = arctan2(Ys_sort[i0]-Ys_sort[i2],Xs_sort[i0]-Xs_sort[i2]) 
                    theta10 = arctan2(Ys_sort[i0]-Ys_sort[i1],Xs_sort[i0]-Xs_sort[i1]) 
                    theta21 = arctan2(Ys_sort[i1]-Ys_sort[i2],Xs_sort[i1]-Xs_sort[i2]) 
                    theta02 = arctan2(Ys_sort[i2]-Ys_sort[i0],Xs_sort[i2]-Xs_sort[i0]) 

                    theta0102 = maximum(theta01,theta02) - minimum(theta01,theta02)
                    theta1210 = maximum(theta12,theta10) - minimum(theta12,theta10)
                    theta2021 = maximum(theta20,theta21) - minimum(theta20,theta21)

                    if theta0102 > pi:  theta0102 -= pi
                    if theta1210 > pi:  theta1210 -= pi
                    if theta2021 > pi:  theta2021 -= pi

                    distdata[k,4] /= (1 + abs(theta0102-pi/3) \
                                        + abs(theta2021-pi/3) \
                                        + abs(theta1210-pi/3))
                    k = k+1
    I  = argsort(-distdata[0:k,4])
    k  = I[0]
    i0 = distdata[k,0]
    i1 = distdata[k,1]
    i2 = distdata[k,2]

#   i3, xbar3, ybar3 are to check for "success"
    for i3 in arange(4):
        if ( (i3!=i0) & \
             (i3!=i1) & \
             (i3!=i2) ): break
    
    xbar0 = Xs_sort[int(i0)]
    ybar0 = Ys_sort[int(i0)]
    xbar1 = Xs_sort[int(i1)]
    ybar1 = Ys_sort[int(i1)]
    xbar2 = Xs_sort[int(i2)]
    ybar2 = Ys_sort[int(i2)]

    xbar3 = Xs_sort[i3]
    ybar3 = Ys_sort[i3]

    if showAlign == True:
        figure(figdumb); figdumb += 1
        imshow(minimum(image, percentile(image,99.9)), cmap='gray', norm=LogNorm())
        plot(xbar0, ybar0, 'r.')
        plot(xbar1, ybar1, 'b.')
        plot(xbar2, ybar2, 'g.')
        title('image data')
        draw()
        show()

    dist01 = sqrt((xbar0-xbar1)**2 + (ybar0-ybar1)**2) 
    dist02 = sqrt((xbar0-xbar2)**2 + (ybar0-ybar2)**2) 
    dist12 = sqrt((xbar1-xbar2)**2 + (ybar1-ybar2)**2) 

    ###################################################################
    # Read in Gaia data
    
    print(tab,'Reading Gaia data')

    (xc,yc,zc) = radec_to_xyz(approx_ra_dec)    # c = center (approximate)

    gaia_folder = script_dir+'csv'+slash
    csv_fekete_name  = gaia_folder+'feketeRaDec.csv' # to mag 12
    csv_fekete_name  = ''       # THIS FORCES THE USE OF ONLINE GAIA DATA
    if (os.path.exists(csv_fekete_name)):
        nf = 1000
        fx    = zeros(nf,dtype=float32);
        fy    = zeros(nf,dtype=float32);
        fz    = zeros(nf,dtype=float32);
        dist  = zeros(nf,dtype=float32);
        j = 0
        with open(script_dir+csv_fekete_name, 'rt') as f:
            reader = csv.reader(f)
            for row in reader:
                fekete_ra  = float(row[1])
                fekete_dec = float(row[2])
                fx[j] = cos(fekete_dec*pi/180)*cos(fekete_ra*pi/180)
                fy[j] = cos(fekete_dec*pi/180)*sin(fekete_ra*pi/180)
                fz[j] = sin(fekete_dec*pi/180)
                j = j+1
        for j in range(nf):
            dist[j] = sqrt( (fx[j]-xc)**2 + (fy[j]-yc)**2 + (fz[j]-zc)**2 )
    
        mindist = zeros(3,dtype=float32)
        jmin    = zeros(3,dtype=int)
        for k in range(3):
            mindist[k] = 1000
            for j in range(nf):
                skipit = False
                if ((k==1) & (j+1==jmin[0])): continue;
                if ((k==2) & (j+1==jmin[0])): continue;
                if ((k==2) & (j+1==jmin[1])): continue;
                if (dist[j] < mindist[k]):
                    mindist[k] = dist[j]
                    jmin[k] = j+1
    
        num_gaia_files=1
        if (mindist[0]+approx_fov*pi/180 > mindist[1]): num_gaia_files = 2
        if (mindist[0]+approx_fov*pi/180 > mindist[2]): num_gaia_files = 3
        print(tab,'Number of Gaia files being looked in: ', num_gaia_files)
    
        m = 0
        for k in range(num_gaia_files):
            csv_file_name  = gaia_folder+'data1234short'+slash+'BrighterStars'+str(jmin[k])+'.csv'
            m += len(open(csv_file_name).readlines())
#       print(tab,'Number of bright Gaia stars: ', m)
        ra  = zeros(m,dtype=float32);
        dec = zeros(m,dtype=float32);
        mag = zeros(m,dtype=float32);
        j = 0
        for k in range(num_gaia_files):
            csv_file_name  = gaia_folder+'data1234short'+slash+'BrighterStars'+str(jmin[k])+'.csv'
            with open(csv_file_name, 'rt') as f:
                reader = csv.reader(f)
                for row in reader:
                    ra[j]  = float(row[0])
                    dec[j] = float(row[1])
                    mag[j] = float(row[2])
                    thisdist = abs(ra[j]-approx_ra_dec[0]) + abs(dec[j]-approx_ra_dec[1])
                    if (mag[j] > 13): continue    # not needed
                    (x1,y1,z1) = radec_to_xyz((ra[j],dec[j]))
                    if ( (approx_ra_dec[0] != 0) | (approx_ra_dec[1] != 0) ):
                        dotprod = xc*x1 + yc*y1 + zc*z1
                        crossprod_x = yc*z1-y1*zc
                        crossprod_y = zc*x1-z1*xc
                        crossprod_z = xc*y1-x1*yc
                        crossprod = sqrt(crossprod_x**2 + crossprod_y**2 + crossprod_z**2)
                        angle = 180*arcsin(crossprod)/pi  # angular separation in degrees
                        if ((dotprod<=0) | (angle>approx_fov)): continue  
                    '''
                    ( ra_hr,   ra_min,  ra_sec) = decimal_to_min_sec( ra[j]/15)
                    (dec_deg, dec_min, dec_sec) = decimal_to_min_sec(dec[j])
                    print(tab,'ra,dec = %3d:%02d:%02d, %3d:%02d:%02d,  mag = %7.3f' %
                    (ra_hr, ra_min, ra_sec, dec_deg, dec_min, dec_sec, mag[j]))
                    '''
                    j = j+1
    else:
        the_query = "SELECT  " + \
            "gaia.ra, gaia.dec, gaia.phot_g_mean_mag " + \
            "FROM gaiadr3.gaia_source AS gaia " + \
            "WHERE DISTANCE({ra:14f}, {dec:14f}, gaia.ra, gaia.dec) < {fov:14f} " + \
            "AND gaia.phot_g_mean_mag <= 13 " + \
            "ORDER BY gaia.phot_g_mean_mag " # + \
        print(tab,'The query: ',    the_query.format(ra=approx_ra_dec[0], dec=approx_ra_dec[1], fov=approx_fov))
        print(tab,'Launching query')
        job = Gaia.launch_job_async(the_query.format(ra=approx_ra_dec[0], dec=approx_ra_dec[1], fov=approx_fov))
        r = job.get_results()
        m = shape(r)[0]
        ra       = zeros((m), dtype=float32)
        dec      = zeros((m), dtype=float32)
        mag      = zeros((m), dtype=float32)
        j = 0
        for k in arange(0,m):
            skipit = False
            for ii in arange(0,3): 
                if r[k][ii]=='': skipit = True
                if r[k][ii]=='--': skipit = True
                if type(r[k][ii]) is ma.core.MaskedConstant: skipit = True
            if skipit: continue
            ra[j]  = float(r[k][0])
            dec[j] = float(r[k][1])
            mag[j] = float(r[k][2])
            (x1,y1,z1) = radec_to_xyz((ra[j],dec[j]))
            if ( (approx_ra_dec[0] != 0) | (approx_ra_dec[1] != 0) ):
                dotprod = xc*x1 + yc*y1 + zc*z1
                crossprod_x = yc*z1-y1*zc
                crossprod_y = zc*x1-z1*xc
                crossprod_z = xc*y1-x1*yc
                crossprod = sqrt(crossprod_x**2 + crossprod_y**2 + crossprod_z**2)
                angle = 180*arcsin(crossprod)/pi  # angular separation in degrees
                if ((dotprod<=0) | (angle>approx_fov)): continue  
            j = j+1


    N = j
    NN = N
    ra  =  ra[:N]
    dec = dec[:N]
    mag = mag[:N]
    print(tab,'Number of close-by bright Gaia stars: ', N)
    if (N<4): 
        print(tab,'Fail: too few nearby Gaia stars')
        return 0, 0, 0, 0, 0, 0, False
#       return approx_ra_dec[0], approx_ra_dec[1], 0, 1, m/2, n/2, False
    if (N>60): N = 60

    ra_ctr  = 0
    dec_ctr = 0

    N2 = int(N**3)
    matchdata = zeros((N2,9))
    k = 0
    for i0 in range(N):
        (x0,y0,z0) = radec_to_xyz((ra[i0],dec[i0]))
        for i1 in range(N):
            if (i1==i0): continue
            (x1,y1,z1) = radec_to_xyz((ra[i1],dec[i1]))
            for i2 in range(N):
                if ((i2==i0) | (i2==i1)): continue
                (x2,y2,z2) = radec_to_xyz((ra[i2],dec[i2]))
                matchdata[k,0] = i0
                matchdata[k,1] = i1
                matchdata[k,2] = i2
                gaia_dist01 = sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
                gaia_dist02 = sqrt((x0-x2)**2+(y0-y2)**2+(z0-z2)**2)
                gaia_dist12 = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                '''
                matchdata[k,4] = \
                    maximum(
                    abs(gaia_dist02/gaia_dist01-dist02/dist01) , \
                    abs(gaia_dist12/gaia_dist01-dist12/dist01)  
                    )
                '''
                matchdata[k,4] = \
                    (gaia_dist02/gaia_dist01-dist02/dist01)**2 + \
                    (gaia_dist12/gaia_dist01-dist12/dist01)**2  
                matchdata[k,5] = gaia_dist02/gaia_dist01
                matchdata[k,6] = gaia_dist12/gaia_dist01
                matchdata[k,7] = dist02/dist01
                matchdata[k,8] = dist12/dist01
                k = k+1

    print(tab,'number of star triples = ', k)
    I = argsort(matchdata[0:k,4])
    matchsort = zeros((len(I),9))
    matchsort = matchdata[I,:]
    '''
    print('matchdata = ', matchdata[0:k,4])
    print('matchsort = ', matchsort[0:k,4])
    '''
    '''
    print(tab,'best matchdata: ', matchsort[0:5,4])
    print(tab,'lengths[0]: ', matchsort[0,5:9])
    print(tab,'lengths[1]: ', matchsort[1,5:9])
    '''

    '''
    k = argmin(matchdata[0:k, 4])
    i0 = int(matchdata[k,0])
    i1 = int(matchdata[k,1])
    i2 = int(matchdata[k,2])
    '''
    for iter in range(5):
        i0 = int(matchsort[iter,0])
        i1 = int(matchsort[iter,1])
        i2 = int(matchsort[iter,2])

        ra0  =  ra[i0]
        ra1  =  ra[i1]
        ra2  =  ra[i2]
        dec0 = dec[i0]
        dec1 = dec[i1]
        dec2 = dec[i2]
    
        '''
        print(tab,'star0 ra,dec = %7.3f %7.3f' %
                (ra0,dec0))
        print(tab,'star1 ra,dec = %7.3f %7.3f' %
                (ra1,dec1))
        print(tab,'star2 ra,dec = %7.3f %7.3f' %
                (ra2,dec2))
        
        print(tab,'star0 x,y = %7.3f %7.3f' %
                (xbar0,ybar0))
        print(tab,'star1 x,y = %7.3f %7.3f' %
                (xbar1,ybar1))
        print(tab,'star2 x,y = %7.3f %7.3f' %
                (xbar2,ybar2))
        '''
        
        ###################################################################
        #   Here's the fun part...
        #
        #   Model:  ra = a*x + b*y + c  <-  3 equations in 3 unknowns
        #          dec = d*x + e*y + f  <-  3 equations in 3 unknowns
    
        A = [ [xbar0, ybar0, 1], [xbar1, ybar1, 1], [xbar2, ybar2, 1] ]
        b = [ ra0, ra1, ra2 ]
        abc = np.linalg.solve(A,b)
    
        b = [ dec0, dec1, dec2 ]
        efg = np.linalg.solve(A,b)
    
    
        (m,n) = image.shape
        x_ctr = n/2
        y_ctr = m/2
        ra_ctr  = abc[0]*x_ctr + abc[1]*y_ctr + abc[2]
        dec_ctr = efg[0]*x_ctr + efg[1]*y_ctr + efg[2]

        ra_ctr /= 15   # convert from degrees back to hours

        flip_it = False
    
        if (showAlign):
                rcParams.update({'figure.max_open_warning': 0})
                rcParams['figure.max_open_warning'] = 0
                fig = figure(1)
                image[0,0] = 0
    #           imshow(minimum(image, percentile(image,99.9)), cmap='gray')
                slidershow(image, fignum=1)
                plot(xbar0,ybar0,'r.')
                plot(xbar1,ybar1,'b.')
                plot(xbar2,ybar2,'g.')
                draw()
                show()
    
    
                ###################################################################
                # Load Psf File
        
                folder = 'fit'+slash
                PSF = fits.getdata(folder+'sxv_psf4.fit').astype(float32)
                w = 199
                PSF = 1+PSF[250-w:250+w,250-w:250+w]
                PSF = PSF - percentile(PSF,80)
                PSF = maximum(PSF,0.0)
                PSF = PSF/PSF.max()
    
                ###################################################################
                # Figure 3         --  Gaia stars aligned
        
                fig = figure(3)
                ax = fig.add_subplot(1,1,1)
                ax.patch.set_facecolor('black')
                
                img = array(zeros((m+2*w,n+2*w), dtype=float32)) + 100
                img[0,0] = 0
                maxmag = np.max(mag)
                fluxfactor = 3
    
                for k in range(NN):
                        A = [ [abc[0], abc[1]], [efg[0], efg[1]] ]
                        b = [ ra[k] - abc[2], dec[k] - efg[2] ]
                        xy = np.linalg.solve(A,b)
                        flux = 10**((maxmag-mag[k]+fluxfactor)/2.5)
                        j = int(xy[0]) + w
                        i = int(xy[1]) + w
                        if ((i>=w) & (i<m+w) & (j>=w) & (j<n+w)):
                            img[i-w:i+w,j-w:j+w] += flux*PSF
                img = img[w:-w,w:-w]
                imshow(minimum(img, percentile(img,99.9)), cmap='gray', norm=LogNorm())
                for k in range(NN):
                        A = [ [abc[0], abc[1]], [efg[0], efg[1]] ]
                        b = [ ra[k] - abc[2], dec[k] - efg[2] ]
                        xy = np.linalg.solve(A,b)
                        if   (k==i0):
                            plot(xy[0], xy[1], 'r.')
                        elif (k==i1):
                            plot(xy[0], xy[1], 'b.')
                        elif (k==i2):
                            plot(xy[0], xy[1], 'g.')
                title('Alignment stars in Gaia data')
                draw()
    
                ###################################################################
                # Figure 4         --  Gaia stars ra/dec plot
        
                fig = figure(4)
                ax = fig.add_subplot(1,1,1)
                ax.patch.set_facecolor('black')
                
                ra_tmp = array(zeros(NN), dtype=float32)
                for k in range(NN):
                    ra_tmp[k] = ra[k]
                    if (ra[k]+180 < 15*ra_ctr): ra_tmp[k] += 360
                    if (ra[k]-180 > 15*ra_ctr): ra_tmp[k] -= 360
                ra_min  = np.min(ra_tmp)
                ra_max  = np.max(ra_tmp)
                dec_min = np.min(dec)
                dec_max = np.max(dec)
    
                maxmag = np.max(mag)
                fluxfactor = 3
    
                '''
                mn = minimum(m,n)
                img = array(zeros((m+2*w,n+2*w), dtype=float32)) + 100
                img[0,0] = 0
                denom = maximum( ra_max-ra_min, dec_max-dec_min )
                for k in range(N):
                        flux = 10**((maxmag-mag[k]+fluxfactor)/2.5)
                        j = int(w+mn*( ra_max- ra_tmp[k])/denom)
                        i = int(w+mn*(dec_max-dec[k])/denom)
                        if ((i>=w) & (i<m+w) & (j>=w) & (j<n+w)):
                            img[i-w:i+w,j-w:j+w] += flux*PSF
                imshow(minimum(img, percentile(img,99.9)), cmap='gray', norm=LogNorm())
                for k in range(N):
                        j = int(w+mn*( ra_max- ra_tmp[k])/denom)
                        i = int(w+mn*(dec_max-dec[k])/denom)
                        if   (k==i0): plot(j, i, 'r.')
                        elif (k==i1): plot(j, i, 'b.')
                        elif (k==i2): plot(j, i, 'g.')
                '''
                x,y,zmax = gnomonic(ra_ctr,    dec_max, ra_ctr,    dec_ctr)
                x,y,zmin = gnomonic(ra_ctr,    dec_min, ra_ctr,    dec_ctr)
                x,ymax,z = gnomonic(ra_max/15, dec_ctr, ra_ctr,    dec_ctr)
                x,ymin,z = gnomonic(ra_min/15, dec_ctr, ra_ctr,    dec_ctr)
                m = int(n*(zmax-zmin)/(ymax-ymin))
                img = array(zeros((m+2*w,n+2*w), dtype=float32)) + 100
                img[0,0] = 0
                for k in range(NN):
                        flux = 10**((maxmag-mag[k]+fluxfactor)/2.5)
                        x,y,z = gnomonic(ra_tmp[k]/15,dec[k],ra_ctr,   dec_ctr)
                        j = w+n-int(round(n*(y-ymin)/(ymax-ymin)))
                        i = w+m-int(round(m*(z-zmin)/(zmax-zmin)))
                        if ((i>=w) & (i<m+w) & (j>=w) & (j<n+w)):
                            img[i-w:i+w,j-w:j+w] += flux*PSF
                imshow(minimum(img, percentile(img,99.9)), cmap='gray', norm=LogNorm())
                for k in range(NN):
                        x,y,z = gnomonic(ra_tmp[k]/15,dec[k],ra_ctr,   dec_ctr)
                        j = w+n-int(round(n*(y-ymin)/(ymax-ymin)))
                        i = w+m-int(round(m*(z-zmin)/(zmax-zmin)))
                        if   (k==i0): plot(j, i, 'r.')
                        elif (k==i1): plot(j, i, 'b.')
                        elif (k==i2): plot(j, i, 'g.')
                title('Alignment stars in Gaia data')
                draw()
                show()
    
        '''
        '''
    #   crpix1=579.2277
    #   crpix2=580.2798
    #   crval1=318.603943485
    #   crval2=46.2853651935
        '''
        '''
    
        (m,n) = image.shape
        crpix1=n/2
        crpix2=m/2
        ra_ctr  = abc[0]*crpix1 + abc[1]*crpix2 + abc[2]
        dec_ctr = efg[0]*crpix1 + efg[1]*crpix2 + efg[2]
    #   print(tab,'ra,dec = ', ra_ctr, dec_ctr)
        ra_ctr /= 15
    
        '''
        (s,  ra_hr,   ra_min,  ra_sec) = decimal_to_min_sec( ra_ctr)
        (s, dec_deg, dec_min, dec_sec) = decimal_to_min_sec(dec_ctr)
        '''
    
        theta0 = arctan2(-ybar1+ybar0,xbar1-xbar0)          # flip vertical axis coords
        theta1 = arctan2(dec[i1]-dec[i0], 
                        (-ra[i1]+ ra[i0])/cos(dec[i0]*pi/180)) # flip horizontal coords
        dtheta = (theta0-theta1)*180/pi
    
        pixel_dist = sqrt((xbar1-xbar0)**2+(ybar1-ybar0)**2) # in pixels
        x0,y0,z0 = gnomonic(ra[i0]/15, dec[i0], ra_ctr, dec_ctr)
        x1,y1,z1 = gnomonic(ra[i1]/15, dec[i1], ra_ctr, dec_ctr)
        y0 = -y0; y1 = -y1;                                 # flip vertical axis coords
        theta1 = arctan2(z1-z0,y1-y0)
        dtheta = (theta0-theta1)*180/pi
        angle_dist = sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)  # in radians
        angle_dist *= 180/pi                                 # in degrees
        angle_dist *= 3600                                   # in arcsecs
        scale = angle_dist/pixel_dist

        ###################################################################
        #   Let's check the solution
    
        ra_3  = abc[0]*xbar3 + abc[1]*ybar3 + abc[2]
        dec_3 = efg[0]*xbar3 + efg[1]*ybar3 + efg[2]
        stardist = 100000
        for k in range(NN):
                distk = (ra[k] - ra_3)**2 + (dec[k] - dec_3)**2
                if (distk < stardist): stardist = distk
        stardist = 3600*sqrt(stardist)
        print(tab,'')
        if (stardist < 6*scale): 
                success=True;  print(tab,'Success') 
                det = linalg.det([[abc[0], abc[1]], [efg[0], efg[1]]])
#               print(tab,'determinant = ', det)
                if (det < 0): flip_it = True
                else:         flip_it = False

                break
        else:                    
                success=False; print(tab,'Fail.  4th-star distance : ', \
                                stardist, 'arcsecs.  Matchsort : ', matchsort[iter,4], \
                                ', iter = ',iter)

    ###################################################################
    #   Let's find the brightest star

    xbar0 = Xs_sort[0]
    ybar0 = Ys_sort[0]
    ra_0  = abc[0]*xbar0 + abc[1]*ybar0 + abc[2]
    dec_0 = efg[0]*xbar0 + efg[1]*ybar0 + efg[2]
    stardist = 100000
    for k in range(NN):
            distk = (ra[k] - ra_0)**2 + (dec[k] - dec_0)**2
            if (distk < stardist): 
                stardist = distk
                k0 = k
    mag0.SetValue(str(mag[k0]))
    
#   tab = tab[0:-4]
#   print('')

    return ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it, success

def radec_to_xyz(radec):
    ra  = radec[0]*pi/180
    dec = radec[1]*pi/180
    x = cos(dec)*cos(ra)
    y = cos(dec)*sin(ra)
    z = sin(dec)
    return x,y,z

#####################################################################################
# Sidereal time = Right Ascension of Zenith at time of observation
    
def SiderealTime(month,day,hour,minute,second,longitude):
    days = 31*(month>1) + \
           28*(month>2) + \
           31*(month>3) + \
           30*(month>4) + \
           31*(month>5) + \
           30*(month>6) + \
           31*(month>7) + \
           31*(month>8) + \
           30*(month>9) + \
           31*(month>10) + \
           30*(month>11) + \
           day-1 + (hour+(minute+(second/60.))/60.)/24.
    
#   ra  = 7+(10+1/60.)/60. # at midnight on Jan. 1 at longitude=0.0
    ra  = 7+(41+7/60.)/60. # at midnight on Jan. 1 at longitude=0.0
    ra  = ra + (days/365)*24 + hour
    ra  = ra + longitude/15
    ra  = (24+ra)%24

    return ra

#####################################################################################
# Ra and Dec of catalog objects.  Ra in hrs, Dec in degs.  Digital values

def getRaDec(obj):
    global objname

    objname = obj.upper()

    print(tab,'obj = ', obj)
    if   (obj.lower()[0]=='m'):     num = int(obj[1:]); fname='NI2021_MessierObjects.csv'; j0=1;
    elif (obj.lower()[0:3]=='ngc'): num = int(obj[3:]); fname=    'NI2021_NGCObjects.csv'; j0=3;
    elif (obj.lower()[0:2]=='ic'):  num = int(obj[2:]); fname=                'ICObj.csv'; j0=2;
    else: print(tab,'Unrecognized object: ', obj); return 0,0

#   print(tab,'num = ', num)
    
    with open(script_dir+'csv'+slash+fname, 'rt') as f:
        reader = csv.reader(f)
        next(reader,None)  # skip first line
        for row in reader:
            j = int(row[0][j0:])
            if (j==num):
                row1 = float(row[1])
                row2 = float(row[2])
                row3 = float(row[3])
                row4 = float(row[4])
                row5 = float(row[5])
                row6 = float(row[6])
#               print(tab,row1, row2, row3, row4, row5, row6)
                ra  = row1+(row2+row3/60)/60
                if (row[4][0] != '-'):
                    dec = row4+(row5+row6/60)/60
                else:
                    dec = row4-(row5+row6/60)/60
                return ra,dec
    print(tab,'Not found: ', obj); return 0,0



#####################################################################################
# Functions that used to be part of StackImages.py

import os.path
from os import path

script_dir = os.path.dirname(os.path.abspath(__file__))+slash

global files, L_files, R_files, G_files, B_files, Ha_files, O3_files, Spectra_files, RGB_files, NonFit_files
global ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it
global fullname, fnames, folder
global flat, flat_file
global bias, bias_file
global panel

folder=''
name=''
flip_it = False

#def askopenfilenames(multiple=False, filetypes='fits files (*.fit)|*.fit|(*.fits)|*.fits|(*.fts)|*.fts'):
#def askopenfilenames(multiple=False, filetypes='JPG and PNG files (*.jpg;*.png)|*.jpg;*.png|FIT files (*.fit)|*.fit|(*.fits)|*.fits|(*.fts)|*.fts'):
#def askopenfilenames(multiple=False, filetypes='JPG and PNG files (*.jpg;*.png;*.fit;*.fits;*.fts;)'):
def askopenfilenames(multiple=False, filetypes='image files (;*.jpeg;*.jpg;*.png;*.fit;*.fits;*.fts;)'):
    global panel
    if (multiple==True):
        with wx.FileDialog(panel, 'Open image files', wildcard=filetypes,
                 style=wx.FD_MULTIPLE | wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return      # the user changed his/her mind
            else:
                return fileDialog.GetPaths()
    else:
        with wx.FileDialog(panel, 'Open image file', wildcard=filetypes,
                        style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:
            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return      # the user changed his/her mind
            else:
                return fileDialog.GetPath()

fullname = ''
fnames=[]
fnames.append('None')

L_files=[]; R_files=[]; G_files=[]; B_files=[]; Ha_files=[]; O3_files=[]; Spectra_files=[]; RGB_files=[]; NonFit_files=[];

global dist

# This function will be used to select an already stacked RGB image
# Only Fits files 
def get_RGB_image(evt):
    global RGB, folder, name, dist
    global objectra, objectdec
    global slash
    global file_location
    global lbl6a
    global btn18
    global dtheta
    global ra_ctr, dec_ctr, crpix1, crpix2, scale

    global tab

    print(tab+'Openning an RGB fits file:')
    tab = tab+'    '

    file_location = askopenfilenames(multiple=False)
    print(tab,'file_location = ', file_location)
    k = 0
    happy = False
    while(not happy):
        image = array(fits.open(file_location, mode='update')[k].data)
        if (len(image.shape)==3):
             (three,m,n) = image.shape
             R = image[0,:,:]
             G = image[1,:,:]
             B = image[2,:,:]
             RGB = zeros((m,n,3), dtype=float32)
             RGB[:,:,0] = R
             RGB[:,:,1] = G
             RGB[:,:,2] = B
             happy = True
        elif (len(image.shape)==2):
             (m,n) = image.shape
             RGB = zeros((m,n,3), dtype=float32)
             RGB[:,:,0] = image
             RGB[:,:,1] = image
             RGB[:,:,2] = image
             happy = True
        k = k+1

    crpix1 = n/2
    crpix2 = m/2

    '''
    Mac_splits = file_location.split('/')
    Win_splits = file_location.split('\\')
    if (len(Mac_splits) >= len(Win_splits)): slash = "/"
    else:                                    slash = "\\"
    '''

    hdul = fits.open(file_location)   # FYI... HDUL = Header Data Unit List
    name = ''
    ra = 0
    dec = 0
    cd1_1 = 0
    cd1_2 = 0
    scale = 0
    dtheta = 0
    for h in hdul[0].header:
        if h.upper() == 'OBJECT':
            name = hdul[0].header['OBJECT']
            name = name.split('-')[0].split('_')[0]
        elif h.upper() == 'CRVAL1':
            ra  = float(hdul[0].header['CRVAL1'])/15
            ra_ctr = ra
        elif h.upper() == 'CRVAL2':
            dec = float(hdul[0].header['CRVAL2'])
            dec_ctr = dec
        elif h.upper() == 'CD1_1':
            cd1_1 = float(hdul[0].header['CD1_1'])
        elif h.upper() == 'CD1_2':
            cd1_2 = float(hdul[0].header['CD1_2'])
    if ((cd1_1 != 0) & (cd1_2 != 0)):
        scale = 3600*sqrt(cd1_1*cd1_1 + cd1_2*cd1_2)
        tancrota2 = cd1_2/cd1_1
        crota2 = arctan2(cd1_2,cd1_1)
        dtheta = -180*crota2/pi

    '''
    if name == '':
        splits = file_location.split(slash)
        NN = len(splits)
        folder = ''
        for j in range(0,NN-1): folder += splits[j]+slash
        name = (((file_location.split(slash)[-1]).split('.fit')[0]).split('-')[0]).split('_')[0]
    '''
    if name == '':
        splits = file_location.split(slash)
        NN = len(splits)
        folder = ''
        for j in range(0,NN-1): folder += splits[j]+slash
        name = (((file_location.split(slash)[-1]).split('.fit')[0]).split('-')[0]).split('_')[0]

    lbl7.SetLabel(merge_paths(folder,file_location))
    if (name != ''): 
        lbl6a.SetLabel(' Default FOV:')
        objectid.SetValue(name)
    if ((ra == 0) & (dec == 0)):
        [ra, dec] = getRaDec(name)

    with open(script_dir+'csv'+slash+'Best500DeepSkyObjects.csv', 'rt') as f:
        reader = csv.reader(f)
        firstrow=True;
        for row in reader:
            if (firstrow == True):
                firstrow = False
                continue
#           print(name, row[0])
            if   ((row[0].split('/')[0] == name.upper()) | (row[1].split(' ')[0] == name.upper()) ):
                dist = 0.306601*float(row[4])   # distance in parsecs
                print(tab,'distance to ', name, ' is ', dist, ' pc')
                break

    objectra.SetValue( decimal_to_min_sec( ra))
    objectdec.SetValue(decimal_to_min_sec(dec))
    lbl_ra.SetLabel( ' Ra  (hh:mm:ss):')
    lbl_dec.SetLabel(' Dec (dd:mm:ss):')
    if (scale != 0):
        objectfov.SetValue(str(round(100*n*scale/3600)/100))
        lbl6a.SetLabel(' FOV (degrees):')
        btn18.Enable()
    '''
    [s,  ra_hr,   ra_min,  ra_sec] = decimal_to_min_sec( ra)
    [s, dec_deg, dec_min, dec_sec] = decimal_to_min_sec(dec)
    ra_sec = round(100*ra_sec)/100
    dec_sec = round(10*dec_sec)/10
    objectra.set(' '+str( ra_hr )+':'+str( ra_min)+':'+str( ra_sec))
    objectdec.set(s+str(dec_deg)+':'+str(dec_min)+':'+str(dec_sec))
    '''
    '''
    objectra.set(round(100000*ra)/100000)
    objectdec.set(round(10000*dec)/10000)
    '''

    lbl4.SetLabel('')
    lbl5.SetLabel('')
    lbl6.SetLabel('')
    lbl6.SetLabel('')
    pixres.SetLabel('')

    show_fwhm(RGB[:,:,1])

    btn16.Enable()
    btn12.Enable()
    btn13.Enable()
    btn14.Enable()
    btn20.Enable()
#   btn21.Enable()

    '''
    hdu = fits.PrimaryHDU(image)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(folder+name+'-py-RGB.fit', overwrite=True)
#   hdulist.close()
#   fits.PrimaryHDU(image).writeto(folder+name+'-py-RGB.fit', overwrite=True)
    '''
    show_stacked_image()

    tab = tab[0:-4]
    printdone('Done openning an RGB fits file')

#   astrometric_analysis(evt)
#   make3D(evt)

def show_fwhm(L): 
    global current_fwhm
    global tab
    print(tab, 'Estimating FWHM:')
    tab = tab+'    '

    global lbl8, lbl8a, pixres, deconvFWHMa, deconvFWHM
    [x,y,current_fwhm] = find_fwhm_star(L)
    if (pixres.GetValue() != ''): 
        lbl8a.SetLabel(' Current FWHM (arcsecs):')
        deconvFWHMa.SetLabel(' Target FWHM (arcsecs): ')
    else:
        lbl8a.SetLabel(' Current FWHM (pixels): ')
        deconvFWHMa.SetLabel(' Target FWHM (pixels):  ')
    lbl8.SetLabel('')
    txt8 = ' {fwhm:<7.3f}'
    lbl8.SetLabel( txt8.format(fwhm=current_fwhm))
    deconvFWHM.SetValue(str(round(10*current_fwhm*3/4)/10))

    tab = tab[0:-4]
    printdone( 'Done estimating FWHM')

# This function will be used to select the files to be stacked
# Only Fits files 
def get_images(evt): 
    global files, L_files, R_files, G_files, B_files, Ha_files, O3_files, Spectra_files, RGB_files, NonFit_files
    global nL, nR, nG, nB, nHa, nO3, nSpec, nRGB, nNonFit, ntot, nstacked
    global images, images1, images2, name, fullname, bias_file, flat_file, fnames, folder
    global objectid, objectra, objectdec
    global showimage
    global slash
    global Mono_files
    global success
    global L_idx
    global lbl6a

    global tab
    print(tab, 'Getting images:')
    tab = tab+'    '

    newfiles = askopenfilenames(multiple=True)
    if shape(newfiles)[0] > 0:
        newfiles = array(newfiles)
        for j in arange(0,len(newfiles)):
          filetype = newfiles[j].split('.')[-1]
          if (filetype[0] == 'f'):
            hdul = fits.open(newfiles[j])   # FYI... HDUL = Header Data Unit List
            k = len(newfiles[j])-4
            filter  = newfiles[j].split('.')[-2].split('-')[-1].split('_')[-1].split(' ')[-1];
            foundit = False
            for h in hdul[0].header:
                if h.upper() == 'FILTER':
                    print(tab,'filter = '+ hdul[0].header['FILTER'])
                    foundit = True
#                   print(tab,'filter = |'+hdul[0].header['FILTER']+'|')
                    if   hdul[0].header['FILTER'][0:3].lower() == 'lum':       L_files.append(newfiles[j])
                    elif hdul[0].header['FILTER'][0:3].lower() == 'red':       R_files.append(newfiles[j])
                    elif hdul[0].header['FILTER'][0:5].lower() == 'green':     G_files.append(newfiles[j])
                    elif hdul[0].header['FILTER'][0:4].lower() == 'blue':      B_files.append(newfiles[j])
                    elif hdul[0].header['FILTER'][0:2].lower() == 'ha':       Ha_files.append(newfiles[j])
                    elif hdul[0].header['FILTER'][0:2].lower() == 'o3':       O3_files.append(newfiles[j])
                    elif hdul[0].header['FILTER'][0:4].lower() == 'oiii':     O3_files.append(newfiles[j])
                    else:                                             foundit = False
                elif ((h.upper() == 'NAXIS') & (hdul[0].header['NAXIS'] == 3)):
                    foundit = True
                    RGB_files.append(newfiles[j])
                    break
            if foundit == False:
                if   filter.upper() == 'L' :        L_files.append(newfiles[j])
                elif filter.upper() == 'R' :        R_files.append(newfiles[j])
                elif filter.upper() == 'G' :        G_files.append(newfiles[j])
                elif filter.upper() == 'B' :        B_files.append(newfiles[j])
                elif filter.upper() == 'HA':       Ha_files.append(newfiles[j])
                elif filter.upper() == 'O3':       O3_files.append(newfiles[j])
                elif filter.upper() == 'OIII':     O3_files.append(newfiles[j])
                elif filter.upper() == 'SPECTRUM': Spectra_files.append(newfiles[j])
                elif newfiles[j][k-1:k].lower() == 'l' :    L_files.append(newfiles[j])
                elif newfiles[j][k-1:k].lower() == 'r' :    R_files.append(newfiles[j])
                elif newfiles[j][k-1:k].lower() == 'g' :    G_files.append(newfiles[j])
                elif newfiles[j][k-1:k].lower() == 'b' :    B_files.append(newfiles[j])
                elif newfiles[j][k-2:k].lower() == 'ha':   Ha_files.append(newfiles[j])
                elif newfiles[j][k-2:k].lower() == 'o3':   O3_files.append(newfiles[j])
                elif newfiles[j][k-4:k].lower() == 'oiii': O3_files.append(newfiles[j])
                elif newfiles[j][k-4:k].lower() == 'trum': Spectra_files.append(newfiles[j])
                elif newfiles[j].find('Light')   > 0:  L_files.append(newfiles[j])
                elif newfiles[j].find('Lum')   > 0:    L_files.append(newfiles[j])
                elif newfiles[j].find('Red')   > 0:    R_files.append(newfiles[j])
                elif newfiles[j].find('Green') > 0:    G_files.append(newfiles[j])
                elif newfiles[j].find('Blue')  > 0:    B_files.append(newfiles[j])
                else: L_files.append(newfiles[j])
            hdul.close()
          else:
            NonFit_files.append(newfiles[j])

        nL      = len( L_files)
        nR      = len( R_files)
        nG      = len( G_files)
        nB      = len( B_files)
        nHa     = len(Ha_files)
        nO3     = len(O3_files)
        nSpec   = len(Spectra_files)
        nRGB    = len(RGB_files)
        nNonFit = len(NonFit_files)

        Mono_files = []
        Mono_files.extend(      L_files)
        Mono_files.extend(      R_files)
        Mono_files.extend(      G_files)
        Mono_files.extend(      B_files)
        Mono_files.extend(     Ha_files)
        Mono_files.extend(     O3_files)
        Mono_files.extend(Spectra_files)
        nMono = len(Mono_files)
        print(tab,'Mono_files = ', Mono_files)

        images0 = [ fits.open(n)[0].data.copy()       for n in Mono_files ]
        images1 = [ fits.open(n)[0].data.copy()       for n in RGB_files ]
        images2 = [ matplotlib.image.imread(n) for n in NonFit_files ]
        minm = 1000000
        minn = 1000000
        print(tab,'len of Mono_files = ',   len(  Mono_files))
        print(tab,'len of RGB_files = ',    len(   RGB_files))
        print(tab,'len of NonFit_files = ', len(NonFit_files))
        for j in range(len(Mono_files)): 
                (m, n) = images0[j].shape
                if (m < minm): minm = m
                if (n < minn): minn = n
        for j in range(len(RGB_files)): 
                (three, m, n) = images1[j].shape
                if (m < minm): minm = m
                if (n < minn): minn = n
        for j in range(len(NonFit_files)): 
                print('j = ', j, ', shape = ', images2[j].shape)
                if (len(images2[j].shape) == 3):
                    (m, n, three) = images2[j].shape
                else:
                    (m, n) = images2[j].shape
                if (m < minm): minm = m
                if (n < minn): minn = n

        N  = len(Mono_files)+3*len(RGB_files)+3*len(NonFit_files)
        M  = len(RGB_files)+len(NonFit_files)
        images  = zeros((nMono+4*nRGB+4*nNonFit,minm,minn))
        L_idx   = zeros( nMono+4*nRGB+4*nNonFit           , dtype='int16')
        jj = 0
        kk = N
        for j in range(0,nL+nR):                     images[jj]  = images0[j][  0:minm,0:minn  ]; L_idx[jj] = jj; jj=jj+1
        for j in range(0,nRGB):                      images[jj]  = images1[j][0,0:minm,0:minn  ]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
        if (nNonFit > 0):
            if (len(shape(images2[0])) == 3):
                for j in range(0,nNonFit):           images[jj]  = images2[j][  0:minm,0:minn,0]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
            else:
                for j in range(0,nNonFit):           images[jj]  = images2[j][  0:minm,0:minn  ]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
        kk = N
        for j in range(nL+nR,nL+nR+nG):              images[jj]  = images0[j][  0:minm,0:minn  ]; L_idx[jj] = jj; jj=jj+1
        for j in range(0,nRGB):                      images[jj]  = images1[j][1,0:minm,0:minn  ]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
        if (nNonFit > 0):
            if (len(shape(images2[0])) == 3):
                for j in range(0,nNonFit):           images[jj]  = images2[j][  0:minm,0:minn,1]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
            else:
                for j in range(0,nNonFit):           images[jj]  = images2[j][  0:minm,0:minn  ]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
        kk = N
        for j in range(nL+nR+nG,nL+nR+nG+nB):        images[jj]  = images0[j][  0:minm,0:minn  ]; L_idx[jj] = jj; jj=jj+1
        for j in range(0,nRGB):                      images[jj]  = images1[j][2,0:minm,0:minn  ]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
        if (nNonFit > 0):
            if (len(shape(images2[0])) == 3):
                for j in range(0,nNonFit):           images[jj]  = images2[j][  0:minm,0:minn,2]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
            else:
                for j in range(0,nNonFit):           images[jj]  = images2[j][  0:minm,0:minn  ]; L_idx[jj] = kk; jj=jj+1; kk=kk+1
        kk = N

        for j in range(nL+nR+nG+nB,len(Mono_files)): images[jj]  = images0[j][  0:minm,0:minn  ]; L_idx[jj] = jj; jj=jj+1

        for j in range(0,nRGB):                      images[jj]  = images1[j][:,0:minm,0:minn  ].sum(axis=0); L_idx[jj] = -1; jj=jj+1
        if (nNonFit > 0):
            if (len(shape(images2[0])) == 3):
                for j in range(0,nNonFit):           images[jj]  = images2[j][  0:minm,0:minn,:].sum(axis=2); L_idx[jj] = -1; jj=jj+1

        images  = array(images)


        nR += nRGB + nNonFit
        nG += nRGB + nNonFit
        nB += nRGB + nNonFit
        files = []
        files.extend(      L_files)
        files.extend(      R_files)
        files.extend(      G_files)
        files.extend(      B_files)
        files.extend(     Ha_files)
        files.extend(     O3_files)
        files.extend(Spectra_files)
        files.extend(    RGB_files)
        files.extend( NonFit_files)

        success = ones((N), dtype=bool)

        '''
        Mac_splits = files[0].split('/')
        Win_splits = files[0].split('\\')
        if (len(Mac_splits) >= len(Win_splits)): slash = "/"
        else:                                    slash = "\\"
        '''

        ra = 0; dec = 0
        name = objectid.GetValue()
        filetype = newfiles[j].split('.')[-1]
        if ((files != [])&(filetype[0] == 'f')): 
            hdul = fits.open(files[0])   # FYI... HDUL = Header Data Unit List
            for h in hdul[0].header:
                if h.upper() == 'OBJECT':
                    name = hdul[0].header['OBJECT']
                    name = (name.split('-')[0]).split('_')[0]
                    if (name != ''): 
                        objectid.SetValue(name)
                        lbl6a.SetLabel(' Default FOV:')
                    break
            fullname = files[0]
            gotName = True
            if (name == ''):
                gotName = False
            elif ((name[0] != 'm') & (name[0] != 'M') & (name[0] != 'n') & (name[0] != 'N') & (name[0] != 'i') & (name[0] != 'I')):
                gotName = False
            if (gotName == False):
                name = (((fullname.split(slash)[-1]).split('.fit')[0]).split('-')[0]).split('_')[0]
                if (name != ''): 
                    objectid.SetValue(name)
                    lbl6a.SetLabel(' Default FOV:')
            [ra, dec] = getRaDec(name)

        objectra.SetValue( decimal_to_min_sec( ra))
        objectdec.SetValue(decimal_to_min_sec(dec))
        '''
        [s,  ra_hr,   ra_min,  ra_sec] = decimal_to_min_sec( ra)
        [s, dec_deg, dec_min, dec_sec] = decimal_to_min_sec(dec)
        ra_sec = round(100*ra_sec)/100
        dec_sec = round(10*dec_sec)/10
        objectra.set(' '+str( ra_hr )+':'+str( ra_min)+':'+str( ra_sec))
        objectdec.set( s+str(dec_deg)+':'+str(dec_min)+':'+str(dec_sec))
        '''
        '''
        objectra.set(round(100000*ra)/100000)
        objectdec.set(round(10000*dec)/10000)
        '''

        splits = fullname.split(slash)
        NN = len(splits)
        folder = ''
        for j in range(0,NN-1): folder += splits[j]+slash

        if len(bias_file) > 0:
          if (fullname != ''):
            lbl2.SetLabel(merge_paths(fullname,bias_file))
          else:
            lbl2.SetLabel(bias_file.split(slash)[-1]+' ')
        if len(flat_file) > 0:
          if (fullname != ''):
            lbl3.SetLabel(merge_paths(fullname,flat_file))
          else:
            lbl3.SetLabel(flat_file.split(slash)[-1]+' ')

        fnames=[]
        for k in range(len(files)):
                fnames.append((files[k].split(slash)[-1]).split('.fit')[0])
    
        if (fnames != []):
            option5.Clear()
            showimage = fnames[0]
            for string in fnames: option5.Append(string)
            option5.SetSelection(0)

        btn8.Enable()

        name = objectid.GetValue()
    tab = tab[0:-4]
    printdone( 'Done getting images')

bias = zeros((0,0))
flat = zeros((0,0))

def merge_paths(p1,p2):
    global slash
    split1 = p1.split(slash)
    split2 = p2.split(slash)
    N1 = len(split1)
    N2 = len(split2)
    for k in range(minimum(N1, N2)):
        if split1[k] != split2[k]: break
    tmp1 = ''
    for j in range(k,N1-1): tmp1 += '..'+slash
    tmp2 = ''
    for j in range(k,N2-1): tmp2 += split2[j]+slash
    tmp2 += split2[N2-1]
    return(tmp1+tmp2+' ')
    
def clear_bias(): 
    global bias, bias_file, fullname
    bias_file = ''
    bias = zeros((0,0))
    lbl2.SetLabel('')
#   btn10.Disable()

def clear_flat(): 
    global flat, flat_file, fullname
    flat_file = ''
    flat = zeros((0,0))
    lbl3.SetLabel('')
#   btn11.Disable()

def get_bias(evt): 
    if btn6.GetLabel() == 'Select Bias Frame':
        global bias, bias_file, fullname
        bias_file = askopenfilenames(multiple=False)
#       if len(bias_file) == 0:
        if bias_file == None: 
            clear_bias()
            return
    
        bias = array( fits.open(bias_file)[0].data )
        bias_text = open('bias_file.txt','w')
        bias_text.write(bias_file)
        bias_text.close()
    
        if (fullname != ''):
            lbl2.SetLabel(merge_paths(fullname,bias_file))
        else:
            lbl2.SetLabel(bias_file.split(slash)[-1]+' ')
        btn6.SetLabel( 'DeSelect' )
#   btn10.Enable()
    else:
        clear_bias()
        btn6.SetLabel( 'Select Bias Frame' )

def get_flat(evt): 
    if btn7.GetLabel() == 'Select Flat Frame':
        global flat, flat_file, fullname
        flat_file = askopenfilenames(multiple=False)
        if len(flat_file) == 0: 
            clear_flat()
            return
    
        flat = array( fits.open(flat_file)[0].data )
        flat_text = open('flat_file.txt','w')
        flat_text.write(flat_file)
        flat_text.close()
    
        if (fullname != ''):
            lbl3.SetLabel(merge_paths(fullname,flat_file))
        else:
            lbl3.SetLabel(flat_file.split(slash)[-1]+' ')
        btn7.SetLabel( 'DeSelect' )
#   btn11.Enable()
    else:
        clear_flat()
        btn7.SetLabel( 'Select Flat Frame' )

def align_stack_show(evt): 
    global files, L_files, R_files, G_files, B_files, Ha_files, O3_files, Spectra_files
    global bias, flat
    global nL, nR, nG, nB, nHa, nO3, nSpec, nRGB, ntot, nstacked
    global images, fnames, name, option4, option5, option6, showimage
    global folder, name
    global thresh, pixres
    global progress
    global tab
    global success

    print(tab, 'Align, Stack, Show:')
    tab = tab+'    '

    #####################################################################################
    # Calibrate 

    '''
    flat_image = array( fits.open('../calib_2020/flat_2x2-L.fit')[0].data )
    bias_image = array( fits.open('../calib_2020/bias_2x2.fit')[0].data )
    '''
    '''
    print(tab,'shape(bias) = ', shape(bias))
    print(tab,'shape(flat) = ', shape(flat))
    print(tab,'shape(image0) = ', shape(images[0,:,:]))
    '''

    if ((shape(bias)[0] == 0)):
        bias = percentile(images[0,:,:],0.1)* ones(shape(images[0,:,:]))
    if ((shape(bias)[0] > 0) & (shape(flat)[0] > 0)):
        if ( (shape(bias) == shape(images[0,:,:])) & 
             (shape(flat) == shape(images[0,:,:])) ):
            print(tab,'Calibrating with bias and flat images')
            applyFlat(images, flat, bias=bias)
        else:
            print(tab,'Bias and/or flat file has wrong size')
            print(tab,'Bias = ', shape(bias), 
                ', Flat = ', shape(flat), 
                ', Image = ', shape(images))
            calib_ok = 0
            [m0, n0] = shape(bias)
            [m1, n1] = shape(flat)
            [m2, n2] = shape(images[0,:,:])
            if ( (m0 == 2*m2) & (n0 == 2*n2) ):
                calib_ok += 1
                bias = (bias[1:,1:]+bias[:-1,1:]+bias[1:,:-1]+bias[:-1,:-1])/4
            elif ( (2*m0 == m2) & (2*n0 == n2) ):
                calib_ok += 1
                tmp = array(zeros((m2,n2)))
                tmp[0::2,0::2] = bias
                tmp[1::2,0::2] = bias
                tmp[0::2,1::2] = bias
                tmp[1::2,1::2] = bias
                bias = tmp
            if ( (m1 == 2*m2) & (n1 == 2*n2) ):
                calib_ok += 1
                flat = (flat[1:,1:]+flat[:-1,1:]+flat[1:,:-1]+flat[:-1,:-1])/4
            elif ( (2*m1 == m2) & (2*n1 == n2) ):
                calib_ok += 1
                tmp = array(zeros((m2,n2)))
                tmp[0::2,0::2] = flat
                tmp[1::2,0::2] = flat
                tmp[0::2,1::2] = flat
                tmp[1::2,1::2] = flat
                flat = tmp
            if (calib_ok == 2):
                print(tab,'Calibrating with resized bias and flat files')
                applyFlat(images, flat, bias=bias)

    #####################################################################################
    # Remove hot pixels

#   if (stack_method == 'Average'):
    if removehotpix == 'Yes':
        print(tab,'Removing hot pixels')
        removeHotPixels(images, thresh=2)
    else:
        print(tab,'Not Removing hot pixels')

    #####################################################################################
    # Align the images

    setMaxBoxSize(150)
    px = pixres.GetValue()
    if (px != ''): setArcsecPerPixel(float(px))
    else:          setArcsecPerPixel(1.0)

    if showalign == 'Yes':
        showAlignStar(True,fnames)
    else:
        showAlignStar(False,fnames)

    print(tab,'Aligning images based on: ', align_method)

    char0 = align_method[0]
    if (len(align_method)>10): char9 = align_method[9]
    else:                      char9 = ''

    [m2, n2] = shape(images[0,:,:])
    if   ((char0=='1') | (char0=='2') | (char0=='3')):
        autoalign(images,int(char0),thresh=float(thresh),border=0,fnames=fnames)
    elif ((char9=='1') | (char9=='2')):
        align(images,int(char9))
    elif (align_method == 'OpenCV'):
        align_cv(images)
    elif (align_method == 'astrometric'):
        align_astrometric(images)

    btn16.Enable()
    btn12.Enable()
    btn13.Enable()
    btn14.Enable()
    btn20.Enable()
#   btn21.Enable()

    if ((char9!='1') & (char9!='2')):
        stack_images()
        show_stacked_image()
        show_fwhm(RGB[:,:,1])

    lbl7.SetLabel(merge_paths(fullname, folder+name+'-py-LRGB.fit'))
    tab = tab[0:-4]
    printdone( 'Done aligning, stacking, showing')

def find_image(evt): 
    global objectra, objectdec, objname
    keycode = evt.GetKeyCode()
    if keycode == wx.WXK_RETURN or keycode == wx.WXK_NUMPAD_ENTER or keycode == wx.WXK_TAB:
        find_image_coords()
    evt.Skip()

def find_image_coords(): 
    global objectra, objectdec, objname
    newname = objectid.GetValue()
    oldname = objname
    if (newname != ''):
        [ra, dec] = getRaDec(newname)
        if ((ra != 0) | (dec != 0)):
            objectra.SetValue( decimal_to_min_sec( ra))
            objectdec.SetValue(decimal_to_min_sec(dec))
            if newname.upper() != oldname:  
                if oldname != '': objname = oldname + ', ' + newname.upper()
                else:             objname = newname.upper()
        else:
            objectra.SetValue( decimal_to_min_sec(0.0))
            objectdec.SetValue(decimal_to_min_sec(0.0))


def get_star_info(ra_mid,dec_mid,fov,maxmag):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia

    global coloridx, mag, absmag, pm_ra, pm_dec
    global image_type
    global ra, dec, parallax, gflux, rflux, bflux

    '''
    the_query = "SELECT  " + \
        "gaia_source.ra, gaia_source.dec, gaia_source.parallax, gaia_source.phot_g_mean_mag, gaia_source.bp_rp " + \
        "FROM gaiadr3.gaia_source " + \
        "WHERE " + \
        "CONTAINS(" + \
                "POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec), " + \
                "CIRCLE('ICRS',{ra:14f},{dec:14f}, {fov:14f})" + \
        ")=1 " + \
        "AND gaia_source.phot_g_mean_mag < {maxmag:5f} " # + \
#       "AND gaia_source.bp_rp < 3 AND gaia_source.bp_rp > -3 " + \
#       "AND gaia_source.parallax >= -1000"
    print(the_query.format(ra=ra_mid*15,dec=dec_mid,fov=fov/sqrt(2),maxmag=maxmag))
    print(tab,'launching query')
    job = Gaia.launch_job_async(the_query.format(ra=ra_mid*15,dec=dec_mid,fov=fov/sqrt(2),maxmag=maxmag))
    '''

    '''
    the_query = "SELECT  " + \
        "ra, dec, parallax, phot_g_mean_mag, bp_rp " + \
        "FROM gaiadr3.gaia_source " + \
        "WHERE " + \
        "CONTAINS(" + \
                "POINT({ra:14f},{dec:14f}, " + \
                "CIRCLE(ra,dec,{fov:14f})" + \
        ")=1 " + \
        "AND phot_g_mean_mag < {maxmag:5f} " + \
        "AND bp_rp < 3 AND bp_rp > -3 " + \
        "AND parallax >= -1000"
    print(the_query.format(ra=ra_mid*15,dec=dec_mid,fov=fov/sqrt(2),maxmag=maxmag))
    print(tab,'launching query')
    job = Gaia.launch_job_async(the_query.format(ra=ra_mid*15,dec=dec_mid,fov=fov/sqrt(2),maxmag=maxmag))
    '''

    '''
    the_query = "SELECT  " + \
        "gaia_source.ra, gaia_source.dec, gaia_source.parallax, gaia_source.phot_g_mean_mag, gaia_source.bp_rp, gaia_source.pmra, gaia_source.pmdec, gaia_source.phot_g_mean_flux " + \
        "FROM gaiadr3.gaia_source " + \
        "WHERE " + \
        "CONTAINS(" + \
                "POINT('ICRS',{ra:14f},{dec:14f}), " + \
                "CIRCLE('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec, {fov:14f})" + \
        ")=1 " + \
        "AND gaia_source.phot_g_mean_mag < {maxmag:5f} " # + \
    print(the_query.format(ra=ra_mid*15,dec=dec_mid,fov=fov/sqrt(2),maxmag=maxmag))
    print(tab,'launching query')
    job = Gaia.launch_job_async(the_query.format(ra=ra_mid*15,dec=dec_mid,fov=fov/sqrt(2),maxmag=maxmag))
    '''

    '''
    fov2 = fov/2
    ramin=(ra_mid*15-fov2/cos(dec_mid*pi/180))
    ramax=(ra_mid*15+fov2/cos(dec_mid*pi/180))
    decmin=dec_mid-fov2
    decmax=dec_mid+fov2
    the_query = "SELECT  " + \
        "gaia.ra, gaia.dec, gaia.parallax, gaia.phot_g_mean_mag, gaia.bp_rp, gaia.pmra, gaia.pmdec, gaia.phot_g_mean_flux " + \
        "FROM gaiadr3.gaia_source AS gaia " + \
        "WHERE " + \
        "gaia.ra  BETWEEN  {ramin:14f} AND  {ramax:14f} AND " + \
        "gaia.dec BETWEEN {decmin:14f} AND {decmax:14f} " + \
        "AND gaia.phot_g_mean_mag < {maxmag:5f} " # + \
    print(                      the_query.format(ramin=ramin, ramax=ramax, decmin=decmin, decmax=decmax, maxmag=maxmag))
    print(tab,'launching query')
    job = Gaia.launch_job_async(the_query.format(ramin=ramin, ramax=ramax, decmin=decmin, decmax=decmax, maxmag=maxmag))
    '''

    '''
    the_query = "SELECT  " + \
        "gaia.ra, gaia.dec, gaia.parallax, gaia.phot_g_mean_mag, gaia.bp_rp, gaia.pmra, gaia.pmdec, gaia.phot_g_mean_flux " + \
        "FROM gaiadr3.gaia_source AS gaia " + \
        "WHERE DISTANCE({ra:14f}, {dec:14f}, gaia.ra, gaia.dec) < {fov:14f} " + \
        "AND gaia.phot_g_mean_mag < {maxmag:5f} " # + \
    '''
    the_query = "SELECT  " + \
        "gaia.ra, gaia.dec, gaia.parallax, gaia.phot_g_mean_mag, gaia.bp_rp, gaia.pmra, gaia.pmdec, gaia.phot_g_mean_flux, " + \
        "gaia.phot_rp_mean_flux, gaia.phot_bp_mean_flux " + \
        "FROM gaiadr3.gaia_source AS gaia " + \
        "WHERE DISTANCE({ra:14f}, {dec:14f}, gaia.ra, gaia.dec) < {fov:14f} " + \
        "AND gaia.phot_g_mean_mag < {maxmag:5f} " # + \
    print(tab,'The query: ',    the_query.format(ra=ra_mid*15, dec=dec_mid, fov=fov/sqrt(2), maxmag=maxmag))
    print(tab,'Launching query')
    if ((ra_mid == 0) & (dec_mid == 0)):  print("WARNING:  RA = 0 and Dec = 0.  You probably didn't hit 'Enter' after entering the Object Catalog Name")
    job = Gaia.launch_job_async(the_query.format(ra=ra_mid*15, dec=dec_mid, fov=fov/sqrt(2), maxmag=maxmag))

#   printdone('almost done')
    r = job.get_results()
#   printdone('done')
    m = shape(r)[0]

    ra       = zeros((m+119615), dtype=float32)
    dec      = zeros((m+119615), dtype=float32)
    parallax = zeros((m+119615), dtype=float32)
    mag      = zeros((m+119615), dtype=float32)
    absmag   = zeros((m+119615), dtype=float32)
    coloridx = zeros((m+119615), dtype=float32)
    pm_ra    = zeros((m+119615), dtype=float32)
    pm_dec   = zeros((m+119615), dtype=float32)
    gflux    = zeros((m+119615), dtype=float32)
    rflux    = zeros((m+119615), dtype=float32)
    bflux    = zeros((m+119615), dtype=float32)
    '''
    ra       = {}
    dec      = {}
    mag      = {}
    coloridx = {}
    '''

    '''
    plx = ones((m), dtype=float32)
    j = 0
    for k in arange(0,m):
        plx[j] = float(r[k][2])
        if (isnan(plx[j])): 
            j += 0
        else:
            j += 1
    print(tab,'plx = ', percentile(plx,0), percentile(plx,0.1), percentile(plx,1))
    '''

    j = 0
    for k in arange(0,m):
#       print(tab,'k = ',k, r[k][2], r[k][5], r[k][6])
        skipit = False
        for ii in arange(0,7): 
            if r[k][ii]=='': skipit = True
            if r[k][ii]=='--': skipit = True
            if type(r[k][ii]) is ma.core.MaskedConstant: skipit = True
#           if math.isnan(r[k][ii]): skipit = True
        if skipit: continue
        ra[j]       = float(r[k][0])/15
        dec[j]      = float(r[k][1])
        parallax[j] = abs(float(r[k][2]))
        mag[j]      = float(r[k][3])
        gflux[j]    = float(r[k][7])
        rflux[j]    = float(r[k][8])
        bflux[j]    = float(r[k][9])
        absmag[j]   = mag[j] - 2.5*log10(100**2/parallax[j]**2)
#       absmag[j]   = 16 - 2.5*log10(flux[j]/parallax[j]**2)
        coloridx[j] = float(r[k][4])
        pm_ra[j]    = float(r[k][5])/15
        pm_dec[j]   = float(r[k][6])
        if (math.isnan(ra[j]) | math.isnan(dec[j]) | math.isnan(mag[j]) | math.isnan(coloridx[j])): # | (parallax[j]<0.002)):
#           print(tab,'Ouch: ',k, ', ', r[k][0], ', ', r[k][1], ', ', r[k][2], ', ', r[k][3], ', ', r[k][4], ', the end')
            j += 0
        elif ((mag[j] >= 4.5) & (coloridx[j]<3.0)):      # Gaia data goes far into the IR
            j += 1

    '''
    k = argmin(absmag)
    if (absmag[k] < -7):
            print(tab,k,len(absmag),'bright star. ra = ', ra[k], 'dec = ', dec[k], 'mag = ', mag[k], 'absmag = ', absmag[k], 'parallax = ', parallax[k])
    '''

    #####################################################################################
    # Read in Hipparcos data
    
    (xc,yc,zc) = radec_to_xyz((ra_mid,dec_mid))    # c = center (approximate)
    '''
    print(tab,'reading Hipparcos data')
    '''
    with open(script_dir+'csv'+slash+'hygdata_v3.csv', 'rt') as f:
        reader = csv.reader(f)
        firstrow=True;
        for row in reader:
            if (firstrow == True):
                firstrow = False
                continue
            if   (row[7] != '') \
               & (row[8] != '') \
               & (row[13] != '') \
               & (row[16] != ''):
                ra[j]       = float(row[7])
                dec[j]      = float(row[8])
                mag[j]      = float(row[13])  # apparent mag
                absmag[j]   = float(row[14])  # absolute mag
                coloridx[j] = float(row[16])
                pm_ra[j]    = float(row[10])/15
                pm_dec[j]   = float(row[11])
                (x1,y1,z1) = radec_to_xyz((ra[j],dec[j]))
                dotprod = xc*x1 + yc*y1 + zc*z1
                crossprod_x = yc*z1-y1*zc
                crossprod_y = zc*x1-z1*xc
                crossprod_z = xc*y1-x1*yc
                crossprod = sqrt(crossprod_x**2 + crossprod_y**2 + crossprod_z**2)
                angle = 180*arcsin(crossprod)/pi  # angular separation in degrees
                if ((dotprod<=0) | (angle>fov*sqrt(2))): continue  
                if ((mag[j] < 5.5) & (mag[j] > -5)): 
                    j += 1

    num_gaia_stars = j

    ra       = ra[0:j]
    dec      = dec[0:j]
    mag      = mag[0:j]
    absmag   = absmag[0:j]
    coloridx = coloridx[0:j]
    parallax = parallax[0:j]
    pm_ra    = pm_ra[0:j]
    pm_dec   = pm_dec[0:j]
    gflux    = gflux[0:j]
    rflux    = rflux[0:j]
    bflux    = bflux[0:j]

    return ra, dec, mag, absmag, coloridx, parallax, pm_ra, pm_dec 


def show_gaia_image(evt): 
#   import astroquery.gaia
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia

    find_image_coords()

    global minmag0, maxmag0, RGB, coloridx, mag, absmag, btn14, objname, abs_vs_app, pm_ra, pm_dec
    global image_type
    global ra, dec, parallax
    global gflux, rflux, bflux

    w = 250
    folder = 'fit'+slash
    sxv_psf = fits.getdata(folder+'sxv_psf4-RGB.fit').astype(float32)
    PSF = zeros((500,500,3), dtype='float32')
    PSF[:,:,0] = 1+sxv_psf[0,250-w:250+w,250-w:250+w]
    PSF[:,:,1] = 1+sxv_psf[1,250-w:250+w,250-w:250+w]
    PSF[:,:,2] = 1+sxv_psf[2,250-w:250+w,250-w:250+w]
    PSF = maximum(PSF - percentile(PSF,80),0)
    PSF = PSF/PSF.max()
    PSF_R = PSF
    '''
    folder = 'png'+slash
    PSF = PIL.Image.open(folder+'stellar_spectrum.png')
    '''
    '''
    folder = 'fit'+slash
    stellar_spectrum = fits.getdata(folder+'stellar_spectrum-RGB.fit').astype(float32)
    PSF = zeros((500,500,3), dtype='float32')
    PSF[:,:,0] = 1+stellar_spectrum[0,250-w:250+w,250-w:250+w]
    PSF[:,:,1] = 1+stellar_spectrum[1,250-w:250+w,250-w:250+w]
    PSF[:,:,2] = 1+stellar_spectrum[2,250-w:250+w,250-w:250+w]
    PSF = maximum(PSF - percentile(PSF,80),0)
    PSF = PSF/PSF.max()

    stellar_spectrum = fits.getdata(folder+'stellar_spectrum-RGB-right-eye.fit').astype(float32)
    PSF_R = zeros((500,500,3), dtype='float32')
    PSF_R[:,:,0] = 1+stellar_spectrum[0,250-w:250+w,250-w:250+w]
    PSF_R[:,:,1] = 1+stellar_spectrum[1,250-w:250+w,250-w:250+w]
    PSF_R[:,:,2] = 1+stellar_spectrum[2,250-w:250+w,250-w:250+w]
    PSF_R = maximum(PSF_R - percentile(PSF_R,80),0)
    PSF_R = PSF_R/PSF_R.max()
    '''

    ra_str  = objectra.GetValue()
    dec_str = objectdec.GetValue()

    ra_mid  = min_sec_to_decimal( ra_str)
    dec_mid = min_sec_to_decimal(dec_str)

    fov = float(objectfov.GetValue())
    maxmag = float(objectmaxmag.GetValue())
    if (fov <= 5):
        abs_vs_app = 'Apparent'
    else:
        abs_vs_app = 'Abs'

    (ra, dec, mag, absmag, coloridx, parallax, pm_ra, pm_dec) \
        = get_star_info(ra_mid=ra_mid,dec_mid=dec_mid,fov=fov,maxmag=maxmag)

    num_gaia_stars = len(ra)

    maxmag  = max([np.max(mag), 18])

    '''
    print('shape(dec)', shape(dec))
    print('ra_min,  ra_max: ',  np.min(ra),  np.max(ra))
    print('dec_min, dec_max: ', np.min(dec), np.max(dec))
    print('done with this stuff')

    print('min negative parallax = ', np.min(parallax[parallax<0]))
    print('max negative parallax = ', np.max(parallax[parallax<0]))
    print('min positive parallax = ', np.min(parallax[parallax>0]))
    print('max positive parallax = ', np.max(parallax[parallax>0]))
    print('number of negative parallaxes = ', shape(parallax[parallax<0]))
    print('number of positive parallaxes = ', shape(parallax[parallax>0]))
    print('number of zero parallaxes = ', shape(parallax[parallax==0]))
    '''

    n = 4000
    m = n
    n2w = n+2*w
    m2w = m+2*w

    if ((image_type == '2D') | (image_type == '3D')): A1s = zeros((m,n,3), dtype=float32)
    if (image_type == '3D'):                          A2s = zeros((m,n,3), dtype=float32)
    if (image_type == 'future'):                      A3s = zeros((m,n,3), dtype=float32)

    vertshift = 0.
    push = 0.1
    maglimit = 22
    fluxfactor = -6

    #####################################################################################
    ## Make star chart
    #
    
    if ((image_type == '2D') | (image_type == '3D')): A10 = zeros((m2w,n2w,3), dtype=float32)
    if (image_type == '3D'):                          A20 = zeros((m2w,n2w,3), dtype=float32)
    if (image_type == 'future'):                      A30 = zeros((m2w,n2w,3), dtype=float32)

    dec_max = dec_mid+fov/2
    dec_min = dec_mid-fov/2

    x,y,zmax = gnomonic(ra_mid,dec_max,ra_mid,dec_mid)
    x,y,zmin = gnomonic(ra_mid,dec_min,ra_mid,dec_mid)
    zmid = (zmax+zmin)/2
    size = zmax-zmid
    ymax =  size
    ymin = -size
    zmax = zmid + size
    zmin = zmid - size

    font = {'family': 'Trade Gothic LT Std',
        'color':  'white',
        'weight': 'normal',
        'fontsize': 4
    }

    num_stars = 0
    maxmag0 = -30
    minmag0 =  30
    maxdec = dec_mid
    mindec = dec_mid
    keepit = zeros((num_gaia_stars), dtype=bool)
    maxpar = max(parallax)
    medianpar = np.median(parallax)
    par90pct = percentile(parallax,90)
    print('min, med, 90pct, max parallax = ', min(parallax), medianpar, par90pct, maxpar)
    rflux = rflux/median(rflux[rflux>0])
    gflux = gflux/median(gflux[gflux>0])
    bflux = bflux/median(bflux[bflux>0])
    for k in arange(num_gaia_stars):
        if ( mag[k] < maglimit ):
            x,y,z    = gnomonic(ra[k],         dec[k],          ra_mid,dec_mid)
        
            keepit[k] = False
            if (x>0):
#           if ((x>0) & (rflux[k]>0) & (gflux[k]>0) & (bflux[k]>0)):
                j  = w+m-int(round(m*(z-zmin)/(zmax-zmin)))
                i  = w+n-int(round(n*( y  - ymin)/( ymax-ymin)))
                if ((i>=w) & (i<n+w) & (j>=w) & (j<m+w)):
                    color = (coloridx[k]-0.7)
                    r = 10**(0.5*color/2.5)  # red
                    b = 1/r
                    g = 1
                    flux = 10**((maxmag-mag[k]+fluxfactor-2)/2.5)
                    '''  
                    # Hipparcos data doesn't include fluxes
                    r = rflux[k]
                    g = gflux[k]
                    b = bflux[k]
                    flux = 1
                    '''
                    if ((image_type == '2D') | (image_type == '3D')):
                        A10[j-w:j+w,i-w:i+w,0] += r*flux*PSF[:,:,0]
                        A10[j-w:j+w,i-w:i+w,1] += g*flux*PSF[:,:,1]
                        A10[j-w:j+w,i-w:i+w,2] += b*flux*PSF[:,:,2]
                    if (mag[k] > maxmag0): maxmag0 = mag[k]
                    if (mag[k] < minmag0): minmag0 = mag[k]
                    if (dec[k] < mindec): mindec = dec[k]
                    if (dec[k] > maxdec): maxdec = dec[k]
                    keepit[k] = True
                    num_stars += 1

                j = w+m-int(round(m*(z-zmin)/(zmax-zmin)))
                i = w+n-int(round(n*(y-ymin)/(ymax-ymin)) + 20*(parallax[k]-2))
#               i = w+n-int(round(n*(y-ymin)/(ymax-ymin)) + 100*parallax[k])
#               i = w+n-int(round(n*(y-ymin)/(ymax-ymin)) + 100*(parallax[k]-par90pct))
#               i = w+n-int(round(n*(y-ymin)/(ymax-ymin)) + 100*abs(parallax[k]-par90pct))
#               i = w+n-int(round(n*(y-ymin)/(ymax-ymin)) + 10*arcsinh(10*(parallax[k]-par90pct)))
                if ((i>=w) & (i<n+w) & (j>=w) & (j<m+w)):
                    '''
                    j = max([j2,w]); j = min([j,m2w-w-1])
                    i = max([i2,w]); i = min([i,n2w-w-1])
                    '''
                    if (image_type == '3D'):
                        A20[j-w:j+w,i-w:i+w,0] += r*flux*PSF_R[:,:,0]
                        A20[j-w:j+w,i-w:i+w,1] += g*flux*PSF_R[:,:,1]
                        A20[j-w:j+w,i-w:i+w,2] += b*flux*PSF_R[:,:,2]

#               1000 years in the future
                j  = w+m-int(round(m*(z-zmin)/(zmax-zmin)))
                i  = w+n-int(round(n*( y  - ymin)/( ymax-ymin)))
                dx,dy,dz = gnomonic(ra[k]+pm_ra[k]/3600,dec[k]+pm_dec[k]/3600,ra_mid,dec_mid)
                j2 = w+m-int(round(m*(dz-zmin)/(zmax-zmin)))
                i2 = w+n-int(round(n*( dy  - ymin)/( ymax-ymin)))
                if ((i>=w) & (i<n+w) & (j>=w) & (j<m+w)):
                    j = max([j2,w]); j = min([j,m2w-w-1])
                    i = max([i2,w]); i = min([i,n2w-w-1])
                    if (image_type == 'future'):
                        A30[j-w:j+w,i-w:i+w,0] += r*flux*PSF[:,:,0]
                        A30[j-w:j+w,i-w:i+w,1] += g*flux*PSF[:,:,1]
                        A30[j-w:j+w,i-w:i+w,2] += b*flux*PSF[:,:,2]
    '''
                else:
                    print(tab,'3: i = ', i, ', j = ', j)
            else:
                print(tab,'2: x = ', x)
        else:
            print(tab,'1: mag[k]', mag[k])
    '''

    ra       =       ra[keepit]
    dec      =      dec[keepit]
    parallax = parallax[keepit]
    mag      =      mag[keepit]
    absmag   =   absmag[keepit]
    coloridx = coloridx[keepit]
    pm_ra    =    pm_ra[keepit]
    pm_dec   =   pm_dec[keepit]
    gflux    =    gflux[keepit]
    rflux    =    rflux[keepit]
    bflux    =    bflux[keepit]

    print(tab,'fov = ', maxdec-mindec)


    print(tab,'number of plotted stars = ', num_stars)
    txt4 = '{space:32s}Number of stars: {space:2s}{ns:d} '
    lbl4.SetLabel( txt4.format(space=' ', ns= num_stars) )

    txt10 = '{space:32s}Brightest star: {space:4s}{minmag:7.2f} '
    txt10fmt = txt10.format(space=' ', minmag= minmag0)
    lbl10.SetLabel(txt10fmt)

    btn14.Enable()
    btn20.Enable()

    if ((image_type == '2D') | (image_type == '3D')):
        A1r = log(A10[w:m+w,w:n+w,0] + 20) - log(20)
        A1g = log(A10[w:m+w,w:n+w,1] + 20) - log(20)
        A1b = log(A10[w:m+w,w:n+w,2] + 20) - log(20)
    
        A1s[:,:,0] = A1r
        A1s[:,:,1] = A1g
        A1s[:,:,2] = A1b

    if (image_type == '3D'):
        A2r = log(A20[w:m+w,w:n+w,0] + 20) - log(20)
        A2g = log(A20[w:m+w,w:n+w,1] + 20) - log(20)
        A2b = log(A20[w:m+w,w:n+w,2] + 20) - log(20)
    
        A2s[:,:,0] = A2r
        A2s[:,:,1] = A2g
        A2s[:,:,2] = A2b

    if (image_type == 'future'):
        A3r = log(A30[w:m+w,w:n+w,0] + 20) - log(20)
        A3g = log(A30[w:m+w,w:n+w,1] + 20) - log(20)
        A3b = log(A30[w:m+w,w:n+w,2] + 20) - log(20)
    
        A3s[:,:,0] = A3r
        A3s[:,:,1] = A3g
        A3s[:,:,2] = A3b

    deltaj = int(round(vertshift*m/1600))

    z = array([arange(0,m)]).T * ones((1,n))
    y = ones((m,1)) * array([arange(0,n)])
    z = zmin + (m-z+deltaj)*(zmax-zmin)/m
    y = ymin + (n-y       )*(ymax-ymin)/n
    raij,decij = inv_gnomonic(y,z, ra_mid, dec_mid)

    if ((image_type == '2D') | (image_type == '3D')):
        A1 = A1s
    if (image_type == 'future'):
        A3 = A3s

    if (image_type == '3D'):
        (m,n,three) = shape(A2s)
        n2 = int(2.05*n)
        A2 = ones((m,n2,3), dtype='float32') * (np.max(A1s) * 0.01)
        A2[:,   0:n ,:] = A1s
        A2[:,n2-n:n2,:] = A2s

#   if (image_type == '3D'): print(tab,'mid value of A2 = ', A2[int(m/2),int(n2/2),:])

    if (image_type == '2D'):
        RGB = A1

    close(2)
    if (image_type == 'future'):
        slidershow(A3, fignum=2, folder=script_dir+'pngs'+slash, fname=objname+'.png')
        suptitle(objname+' -- 1000 Years in the Future', color='black')
#       draw()
    elif (image_type == '3D'):
        slidershow(A2, fignum=2, folder=script_dir+'pngs'+slash, fname=objname+'_splayeyed.png')
        suptitle(objname+' -- 3D Splay-Eyed', color='black')
#       draw()
    elif (image_type == '2D'):
        slidershow(A1, fignum=2, folder=script_dir+'pngs'+slash, fname=objname+'.png')
        suptitle(objname, color='black')
#       draw()
#   show()

def show_image(evt): 
    global files, L_files, R_files, G_files, B_files, Ha_files, O3_files, Spectra_files, RGB_files, NonFit_files
    global images1, images2
    global bias, flat
    global nL, nR, nG, nB, nHa, nO3, nSpec, nRGB, ntot, nstacked
    global images, fnames, name, option4, option5, option6, showimage, folder

    imagename = showimage
    for k in range(len(RGB_files)):
        print(tab,imagename,', ', (RGB_files[k].split(slash)[-1]).split('.fit')[0])
        if (imagename == (RGB_files[k].split(slash)[-1]).split('.fit')[0]):
            close(2)
            (three,m,n) = shape(images1[k])
            image1 = zeros((m,n,3), dtype=float32)
            image1[:,:,0] = images1[k][0,:,:]
            image1[:,:,1] = images1[k][1,:,:]
            image1[:,:,2] = images1[k][2,:,:]
            print('Bob1')
            fig = slidershow(image1, folder=folder, fignum=2)
            return
    for k in range(len(NonFit_files)):
        print(tab,imagename,', ', (NonFit_files[k].split(slash)[-1]).split('.')[0])
        if (imagename == NonFit_files[k].split(slash)[-1]):
            close(2)
            (m,n,three) = shape(images2[k])
            image2 = zeros((m,n,3), dtype=float32)
            image2[:,:,:] = images2[k][:,:,:]
            print('Bob2')
            fig = slidershow(image2, folder=folder, fignum=2)
            return
    for k in range(len(fnames)):
        if imagename == fnames[k]: 
            close(2)
            fig = slidershow(images[k], folder=folder, fignum=2)
            return
#   figure(10)
#   imshow(minimum(images[k], 0.1*np.max(images[k])), cmap='gray')
#   show()
    
def min_sec_to_decimal(min_sec_str):
    splits = min_sec_str.split(':')
    sgn    = splits[0][0]
    hour   = int(splits[0][1:])
    minute = int(splits[1])
    second = float(splits[2])
    if (sgn != '-'): 
        decimal =   hour + (minute + second/60)/60
    else:
        decimal = -(hour + (minute + second/60)/60)
    return decimal



def astrometric_analysis(evt):
    global RGB, folder, name
    global objectra, objectdec
    global file_location
    global btn18
    global lbl6a
    global tab
    global ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it

    print(tab, 'Astrometric Analysis:')
    tab = tab+'    '

    objname = objectid.GetValue()
    if (objname != ''):
        [ra, dec] = getRaDec(objname)
        if ((ra != 0) | (dec != 0)):
            objectra.SetValue( decimal_to_min_sec( ra))
            objectdec.SetValue(decimal_to_min_sec(dec))

    ra_str  = objectra.GetValue()
    dec_str = objectdec.GetValue()

    ra = min_sec_to_decimal(ra_str)
    '''
    splits = ra_str.split(':')
    ra_sgn =       splits[0][0]
    ra_hr  =   int(splits[0][1:])
    ra_min =   int(splits[1])
    ra_sec = float(splits[2])
    if (ra_sgn != '-'): 
        ra =   ra_hr + (ra_min + ra_sec/60)/60
    else:
        ra = -(ra_hr + (ra_min + ra_sec/60)/60)
    '''

    dec = min_sec_to_decimal(dec_str)
    '''
    splits = dec_str.split(':')
    dec_sgn =       splits[0][0]
    dec_hr  =   int(splits[0][1:])
    dec_min =   int(splits[1])
    dec_sec = float(splits[2])
    if (dec_sgn != '-'):
        dec =   dec_hr + (dec_min + dec_sec/60)/60
    else:
        dec = -(dec_hr + (dec_min + dec_sec/60)/60)
    '''

    fov = float(objectfov.GetValue())

    L = RGB[:,:,0]+RGB[:,:,1]+RGB[:,:,2]
    (ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it, success) = \
            astrometry(L, approx_ra_dec=[ra, dec], approx_fov=fov, \
                            showAlign=False)
    
    '''
    (s,  ra_hr,   ra_min,  ra_sec) = decimal_to_min_sec( ra_ctr)
    (s, dec_deg, dec_min, dec_sec) = decimal_to_min_sec(dec_ctr)
    '''

    '''
    objectra.set(round(100000*ra_ctr)/100000)
    objectdec.set(round(10000*dec_ctr)/10000)
    '''

    print(tab,'')
    print(tab,'Output: ')
    print(tab,'  RA, Dec at center = %s, %s' %
            (decimal_to_min_sec( ra_ctr), decimal_to_min_sec(dec_ctr)))
    print(tab,'  Angular rotation  = %7.3f deg' %
            (dtheta))
    print(tab,'  Pixel scale    =    %7.3f arcsec/pix' %
            (scale))
    print(tab,' ')
    
    txt4 = '  Success:  RA = {ra:14s},   '+\
                     'Dec = {dec:14s}'
    txt5 = '                   Angular rotation = {rot:7.3f} deg, '
    txt6 = '                   Pixel scale = {pix:7.3f} arcsec/pix'
    if (success==False): 
        lbl4.SetLabel('Fail (maybe try a different FOV)     ')
        lbl5.SetLabel('')
        lbl6.SetLabel('')
        btn18.Disable()
        return
    else:                
        lbl4.SetLabel( \
           txt4.format(ra= decimal_to_min_sec( ra_ctr), \
                      dec= decimal_to_min_sec(dec_ctr) ))
        lbl5.SetLabel( txt5.format(rot=dtheta))
        lbl6.SetLabel( txt6.format(pix=scale))
        pixres.SetValue(str(round(1000*scale)/1000))
        setArcsecPerPixel(float(pixres.GetValue()))
        objectfov.SetValue(str(round(10*float(pixres.GetValue())*2*sqrt(crpix1**2+crpix2**2)/3600)/10))
        lbl6a.SetLabel(' Current FOV (degrees): ')
        show_fwhm(RGB[:,:,1])

        #####################################################################################
        # Add astrometry info to fits file
        
        image = fits.getdata(file_location).astype(float32)
        image = array(image)
        hdul = fits.open(file_location, mode='update')   # FYI... HDUL = Header Data Unit List
    
        '''
        ra = 15*(ra_hr + (ra_min + ra_sec/60)/60)
        if (s != '-'): dec =   dec_deg + (dec_min + dec_sec/60)/60
        else:          dec = -(dec_deg + (dec_min + dec_sec/60)/60)
        '''
        ra  = ra_ctr * 15
        dec = dec_ctr
        if (len(image.shape)==2): (m,n) = image.shape
        else:               (three,m,n) = image.shape
        crota2 = -dtheta*pi/180
        cd001001 =  cos(crota2) * -1
        cd001002 =  sin(crota2) * -1
        cd002001 = -sin(crota2) * -1
        cd002002 =  cos(crota2) * -1
        cd1_1 = cd001001*scale/3600
        cd1_2 = cd001002*scale/3600
        cd2_1 = cd002001*scale/3600
        cd2_2 = cd002002*scale/3600
    
        hdul0header = hdul[0].header
        hdul[0].header['CTYPE1'] = 'RA---TAN'  # / TAN (gnomic) projection'
        hdul[0].header['CTYPE2'] = 'DEC--TAN'  # / TAN (gnomic) projection'
        hdul[0].header['EQUINOX'] = 2000.0     # / Equatorial coordinates definition (yr)'
        hdul[0].header['CRVAL1'] = ra     # +'         / RA  of reference point'
        hdul[0].header['CRVAL2'] = dec    # +'        / DEC of reference point'
        hdul[0].header['CRPIX1'] = crpix1 # +'        / X reference pixel'
        hdul[0].header['CRPIX2'] = crpix2 # +'        / Y reference pixel'
        hdul[0].header['CUNIT1'] = 'deg     ' # / X pixel scale units'
        hdul[0].header['CUNIT2'] = 'deg     ' # / Y pixel scale units'
        hdul[0].header['CD1_1'] = cd1_1   # +'      / Transformation matrix'
        hdul[0].header['CD1_2'] = cd1_2   # +'      / no comment'
        hdul[0].header['CD2_1'] = cd2_1   # +'      / no comment'
        hdul[0].header['CD2_2'] = cd2_2   # +'      / no comment'
#       hdul[0].header['XPIXSZ'] = scale   # +'      / no comment'
#       hdul[0].header['YPIXSZ'] = scale   # +'      / no comment'
        hdul[0].header['IMAGEW'] = n   # +'      / no comment'
        hdul[0].header['IMAGEH'] = m   # +'      / no comment'
        
        '''
        image = maximum(image-percentile(image,10),1)
        image = 100*image/np.max(image) + 1
        image = log(image)
        image = 65535*image/np.max(image)
        '''
    
#       out_file_location = folder+name+'-py-RGB.fit'
#       fits.writeto(out_file_location, image, hdul[0].header, overwrite=True)
#       hdul.writeto(out_file_location, image, hdul[0].header, overwrite=True)
#       hdul.close()
        hdul.flush()

        btn18.Enable()
#       make3D(evt)

    tab = tab[0:-4]
    printdone( 'Done doing an Astrometric Analysis')

'''
def quit(evt):
    global app
    close('all')
    app.Destroy()
'''

def make3D(evt):
    global RGB, file_location
    global ra_ctr, dec_ctr, dtheta, scale, crpix1, crpix2, flip_it
    global dist

    folder_file = file_location.split('.')[0]

    close(1)
    
#   ###################################################################################
#    Combine the Red, Green, and Blue frames into a single Luminance frame
    
    (m,n) = shape(RGB[:,:,1])
    R = array(zeros((m,n), dtype=float32))
    G = array(zeros((m,n), dtype=float32))
    B = array(zeros((m,n), dtype=float32))
    R[:,:] = RGB[:,:,0]
    G[:,:] = RGB[:,:,1]
    B[:,:] = RGB[:,:,2]
    L = R+G+B       
    
    (m,n) = shape(L)
    image = array(zeros((m,n,3), dtype=float32))
    image[:,:,0] = R
    image[:,:,1] = G
    image[:,:,2] = B
    
    rightimage = array(zeros((m,n,3), dtype=float32))   # right-eye image
    rightimage = image[:,:,:]
    
#   --------------------------------------------------------------
#   ------------------ Find and Remove Stars ---------------------
    
    thresh = 0.01
    '''
    unhappy = True
    while (unhappy):
        K, LconvL2 = findstars(L,thresh, fignum=0)
        close(2)
        num_stars = np.sum(K)
        print(tab,'num stars = ', num_stars)
        if (num_stars > 500):
            thresh *= 1.5
        else:
            unhappy = False
    '''
    K, LconvL2 = findstars(L,thresh, fignum=23)
    num_stars = np.sum(K)
#   close(22)


    '''
    K[0:30,:] = 0
    K[m-30:m,:] = 0
    K[0,0:30] = 0
    K[0,n-30:n] = 0
    K[LconvL2==0] = 0
    num_stars = np.sum(K)
    print(tab,'num stars = ', num_stars)
    fwhms = []
    y = array([arange(0,m)]).T * ones((1,n))
    x = ones((m,1)) * array([arange(0,n)])
    x = x[K].astype(int)
    y = y[K].astype(int)
    L2flux = LconvL2[K]
    I = argsort(L2flux)
    L2fluxsort = L2flux[I]
    xsort = x[I]
    ysort = y[I]
    figure(23); 
    for j in arange(len(L2fluxsort)):
        fwhms.append(fwhm(L,ysort[j],xsort[j]))
        plot(j,fwhms[j],'r.')
    figure(24); 
    for j in arange(len(L2fluxsort)):
        plot(j,L2fluxsort[j],'r.')
    show();
    med_fwhm = np.median(fwhms)
    '''

    
    I  = (ones((1,n)).T*arange(m)).T
    J  =  ones((1,m)).T*arange(n)
    
    Kflat = K.flat
    Iflat = array((I.flat)).astype(int)
    Jflat = array((J.flat)).astype(int)
    
    rmax  = 20
    rmax2 = rmax+1
    rmax3 = rmax+2
    stars = array(zeros((num_stars, 2*rmax2+1,2*rmax2+1,3), dtype=float32))
    r2s = array(zeros(num_stars, dtype=int))
    ii  = array(zeros(num_stars, dtype=int))
    jj  = array(zeros(num_stars, dtype=int))
    ra  = array(zeros(num_stars, dtype=float32))
    dec = array(zeros(num_stars, dtype=float32))
    star = 0
    
    maxG = np.max(G)
    Gorig = array(zeros((m,n), dtype=float32))
    Gorig[:,:] = G
    maxflux = np.max(LconvL2)
    print(tab,'Making starless image')
    '''
    figure(23); imshow(G[600-rmax2:600+rmax3, 707-rmax2:707+rmax3], cmap='gray'); title('More Before')
    '''
    fwhms = zeros(num_stars)
    for iter in arange(len(Kflat)):                 # for each pixel in the image...
        k = Kflat[iter]
        i = Iflat[iter]
        j = Jflat[iter]
        if k:                                       # if there' a star on this pixel...
            fwhms[star] = fwhm(L,i,j)
            ii[star] = i
            jj[star] = j
            star += 1;
    med_fwhm = np.median(fwhms)

    min_flux = 1e+10
    for star in arange(num_stars):                 # for each pixel in the image...
            i = ii[star]
            j = jj[star]
            flux = LconvL2[i,j]
            if (flux < min_flux): min_flux = flux

    for star in arange(num_stars):                 # for each pixel in the image...
            i = ii[star]
            j = jj[star]
#           print(tab,star, i, j)
            '''
#           if (star == 263):
            if (star == 191):
                figure(24); imshow(G[i-rmax2:i+rmax3, j-rmax2:j+rmax3], cmap='gray'); title('Before')
                print(tab,'191 flux = ', LconvL2[i,j])
                print(tab,'191 fwhm = ', fwhms[star])
            '''

            flux = LconvL2[i,j]
            r0 = min([max([int(2.5*arcsinh(flux/min_flux)*med_fwhm),4]),rmax])
            r1 = r0 + 0
            r2 = r0 + 1
            r2s[star] = r2
#           print(tab,star, ', i,j = ', i,j,': r2 = ', r2)
            xs =        arange(-r2,r2+1)
            ys = array([arange(-r2,r2+1)]).T
            Xs = ones(((2*r2+1),1))*xs
            Ys = ys*ones((1,(2*r2+1)))
            radii = sqrt( Xs**2 + Ys**2 )
            disk = (radii <= r0)
#           annulus = float32((radii <= r2) & (radii > r1))
            annulus = ((radii <= r2) & (radii > r1))
    
            Dflat = disk.flat
            Aflat = annulus.flat
            Xflat = Xs.flat
            Yflat = Ys.flat
    
            for px0 in arange(len(Dflat)):          # for each pixel around the star...
                dd = Dflat[px0]
                if dd:                              # if the pixel is on the disk...
                    x0 = int(Xflat[px0])
                    y0 = int(Yflat[px0])
                    if ( (i+y0>=0) & (i+y0<m) & (j+x0>=0) & (j+x0<n) ):
                        numR = 0
                        numG = 0
                        numB = 0
                        den  = 0
                        '''
                        for px1 in arange(len(Aflat)):  # for another pixel around the star...
                            aa = Aflat[px1]
                            if aa:                      # if the second pixel in in the annulus...
                                x1 = int(Xflat[px1])
                                y1 = int(Yflat[px1])
                                numR += R[i+y1,j+x1]/sqrt( (x1-x0)**2 + (y1-y0)**2 )
                                numG += G[i+y1,j+x1]/sqrt( (x1-x0)**2 + (y1-y0)**2 )
                                numB += B[i+y1,j+x1]/sqrt( (x1-x0)**2 + (y1-y0)**2 )
                                den  +=      1      /sqrt( (x1-x0)**2 + (y1-y0)**2 )
                        '''
                        tmpR = R[i-r2:i+r2+1,j-r2:j+r2+1]
                        numR = np.sum(tmpR[annulus]/sqrt( (Xs[annulus]-x0)**2 + (Ys[annulus]-y0)**2 ))
                        tmpG = G[i-r2:i+r2+1,j-r2:j+r2+1]
                        numG = np.sum(tmpG[annulus]/sqrt( (Xs[annulus]-x0)**2 + (Ys[annulus]-y0)**2 ))
                        tmpB = B[i-r2:i+r2+1,j-r2:j+r2+1]
                        numB = np.sum(tmpB[annulus]/sqrt( (Xs[annulus]-x0)**2 + (Ys[annulus]-y0)**2 ))
                        den  = np.sum(     1       /sqrt( (Xs[annulus]-x0)**2 + (Ys[annulus]-y0)**2 ))

                        bkg = numR/den
                        stars[star, rmax2+y0, rmax2+x0, 0] = R[i+y0,j+x0] - bkg
                        R[i+y0,j+x0] = bkg
        
                        bkg = numG/den
                        stars[star, rmax2+y0, rmax2+x0, 1] = G[i+y0,j+x0] - bkg
                        G[i+y0,j+x0] = bkg
        
                        bkg = numB/den
                        stars[star, rmax2+y0, rmax2+x0, 2] = B[i+y0,j+x0] - bkg
                        B[i+y0,j+x0] = bkg

            '''
#           if (star == 263):
            if (star == 191):
                figure(25); imshow(G[i-rmax2:i+rmax3, j-rmax2:j+rmax3], cmap='gray'); title('After')
            '''

            cos_dt = cos(dtheta*pi/180)
            sin_dt = sin(dtheta*pi/180)

            if (flip_it == False):
                ra[star]  =  ra_ctr - ( cos_dt*(j-crpix1) - sin_dt*(i-crpix2))*scale/(3600*15)
                dec[star] = dec_ctr - ( sin_dt*(j-crpix1) + cos_dt*(i-crpix2))*scale/ 3600
            else:
                ra[star]  =  ra_ctr + ( cos_dt*(j-crpix1) - sin_dt*(i-crpix2))*scale/(3600*15)
                dec[star] = dec_ctr - ( sin_dt*(j-crpix1) + cos_dt*(i-crpix2))*scale/ 3600

            '''
            if (Gorig[i,j] > 0.5*maxG):
                print(tab,'Gorig[i,j] = ',Gorig[i,j], ', maxG = ', maxG)
                print(tab,'i,j = ', i, j, \
                      '   ra = ', decimal_to_min_sec( ra[star]), \
                      ', dec = ', decimal_to_min_sec(dec[star]))
            '''

    
    starlessimage = array(zeros((m,n,3), dtype=float32))
    starlessimage[:,:,0] = R
    starlessimage[:,:,1] = G
    starlessimage[:,:,2] = B
    
    tmp = starlessimage - np.min(starlessimage)
    tmp = tmp/np.max(tmp)
    tmp[tmp<0] = 0
    matplotlib.image.imsave(folder_file+'_nostars.png', tmp)
    
    starlessfitimage = array(zeros((3,m,n), dtype=float32))
    starlessfitimage[0,:,:] = R
    starlessfitimage[1,:,:] = G
    starlessfitimage[2,:,:] = B
    hdu = fits.PrimaryHDU(starlessfitimage)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(folder_file+'_nostars.fit', overwrite=True)
    
    '''
    figure(3, figsize=(18, 14), dpi=80)
    imshow(starlessimage/np.max(starlessimage))
    '''
#   show()

    if (path.exists(folder_file+'_depth_PS.png')):
        Lum = PIL.Image.open(folder_file+'_depth_PS.png')
        Lum = Lum - np.min(Lum)
        Lum = Lum/np.max(Lum)
    else:
        L = (R+G+B)/3
        Lum = array(zeros((m,n,3), dtype=float32))
        Lum[:,:,0] = L
        Lum[:,:,1] = L
        Lum[:,:,2] = L
    
        Lum = Lum - np.min(Lum)
        Lum = Lum/np.max(Lum)
        matplotlib.image.imsave(folder_file+'_depth.png', Lum)
    
    
#   --------------------------------------------------------------
#   -------------------- Make 3D nebula  -------------------------
    
    leftimage = array(zeros((m,n,3), dtype=float32))        # left-eye image
    max_shift = 30
    for j in arange(max_shift):
        k = ((Lum >= (j/max_shift)**4) & (Lum <= ((j+1)/max_shift)**4))
        print(tab,'j = ', j, ', sum(k) = ', np.sum(k))
#       leftimage[:,0:n-j,:][k[:,0:n-j,:]] = starlessimage[:,j:n,:][k[:,0:n-j,:]]
        leftimage[:,j:n,:][k[:,j:n,:]] = starlessimage[:,0:n-j,:][k[:,j:n,:]]
    leftimage[:,0:n-10,:] = leftimage[:,10:n,:]
    leftimage[:,n-10:n,:] = np.min(starlessimage[10:-10,n-10,:])

    ra_max  =  np.max(ra)
    ra_min  =  np.min(ra)
    dec_max =  np.max(dec)
    dec_min =  np.min(dec)
    ra_max  =   ra_ctr + 1.2*( ra_max- ra_ctr)
    ra_min  =   ra_ctr + 1.2*( ra_min- ra_ctr)
    dec_max =  dec_ctr + 1.2*(dec_max-dec_ctr)
    dec_min =  dec_ctr + 1.2*(dec_min-dec_ctr)

    (hyg_ra, hyg_dec, hyg_mag, hyg_absmag, hyg_coloridx, hyg_parallax, hyg_pm_ra, hyg_pm_dec) \
        = get_star_info(ra_mid=ra_ctr,dec_mid=dec_ctr,fov=1.5*(dec_max-dec_min),maxmag=13)

    num_hyg_stars = len(hyg_ra)


    parallax_max    = np.percentile(hyg_parallax, 99)
    parallax_min    = np.percentile(hyg_parallax,  1)
    parallax_median = np.percentile(hyg_parallax, 50)
    parallax_iqr    = np.percentile(hyg_parallax, 75) - np.percentile(hyg_parallax, 25)

    print('dist = ', dist)
    obj_parallax = 1000/dist
    print(tab,'obj parallax = ', obj_parallax, ', median parallax = ', parallax_median)



    pos_diff = zeros(num_stars, dtype='float32')
    for star in arange(num_stars):
        i = ii[star]
        j = jj[star]  # + int16(floor(random.uniform(0,30))) - 10
        hyg_star = argmin(15*abs(hyg_ra           - ra[star]) + abs(hyg_dec           - dec[star]))
        pos_diff[star]  = 15*abs(hyg_ra[hyg_star] - ra[star]) + abs(hyg_dec[hyg_star] - dec[star])
        dj = int(10* (hyg_parallax[hyg_star] - obj_parallax) / parallax_median) - 20
        j += dj
#       print(tab,'i,j = ', i, j, ', leftimage shape = ', shape(leftimage[i-rmax2:i+rmax3,j-rmax2:j+rmax3,:]), \
#                       ', stars shape = ', shape(stars[star,:,:,:]))
        if ((i-rmax2 > 0) & (i+rmax3 < n) & \
            (j-rmax2 > 0) & (j+rmax3 < n)):
            leftimage[i-rmax2:i+rmax3,j-rmax2:j+rmax3,:] += stars[star,:,:,:]

#   Check things out (good for debugging)
        i = ii[star]
        j = jj[star]
        if (Gorig[i,j] > 0.5*maxG):
            print(tab,'star = ',star,': i,j = ', i, j, ',  Gorig[i,j]/maxG = ',Gorig[i,j]/maxG)
            print(tab,'           ra = ', decimal_to_min_sec( ra[star]), \
                  ',         dec = ', decimal_to_min_sec(dec[star]))
            print(tab,'       hyg_ra = ', decimal_to_min_sec( hyg_ra[hyg_star]), \
                  ',     hyg_dec = ', decimal_to_min_sec(hyg_dec[hyg_star]))

    
#   --------------------------------------------------------------
#   -------------------- Normalize and Show ----------------------
    
    '''
#   leftimage = leftimage - np.min(leftimage[leftimage>0])
    leftimage = leftimage/np.percentile(leftimage,99.9)
    leftimage[leftimage>1] = 1
    leftimage[leftimage<0] = 0
    
#   rightimage = rightimage - np.min(rightimage[rightimage>0])
    rightimage = rightimage/np.percentile(rightimage,99.9)
    rightimage[rightimage>1] = 1
    rightimage[rightimage<0] = 0
    '''
    
    n2 = 2*n + 30
    m2 =   m + 20
    Img3Dsplay = 0.8*np.percentile(rightimage,99.9)*ones((m2,n2,3), dtype='float32')
#   Img3Dsplay[10:-10,     10:n +10,:] =  leftimage
#   Img3Dsplay[10:-10,n2-n-10:n2-10,:] = rightimage
    Img3Dsplay[10:-10,     10:n +10,:] = rightimage
    Img3Dsplay[10:-10,n2-n-10:n2-10,:] =  leftimage
    
    '''
    figure(6, figsize=(18, 14), dpi=80)
    title('M27  --   View Splayeyed')
    imshow(Img3Dsplay)
    '''
    fig = slidershow(Img3Dsplay, folder=folder, fname=name+'-3D.png')
#   fig.suptitle(name.upper()+'   -   3D View Splay-Eyed', fontsize=12)
    fig.suptitle(name.upper()+'   -   3D View Cross-Eyed', fontsize=12)
    
    '''
    Img3Dcross = ones((m,n2,3), dtype='float32')
    Img3Dcross[:,   0:n ,:] = rightimage
    Img3Dcross[:,n2-n:n2,:] =  leftimage
    
    figure(7, figsize=(18, 14), dpi=80)
    title('M27  --   View Crosseyed')
    imshow(Img3Dcross)
    
    fig = figure(10, figsize=(5.6, (7.225/14.3)*5.6))
    fig.suptitle('M27   --   3D View Splay-eyed', fontsize=12)
    
    axes1 = fig.add_subplot(1,2,2)  # last param 1 for splay 2 for cross
    axes1.set_xlim([0,n])
    axes1.set_ylim([0,m])
    gca().invert_yaxis()
    axes1.imshow(rightimage)
    axes1.axis('off')
    
    axes2 = fig.add_subplot(1,2,1)  # last param 2 for splay 1 for cross
    axes2.set_xlim([0,n])
    axes2.set_ylim([0,m])
    gca().invert_yaxis()
    axes2.imshow(leftimage)
    axes2.axis('off')
    
    tight_layout()
    
    subplots_adjust(top=1.0, bottom=0.0, left=0.15/14, right=1.0-0.15/14, hspace= 0.0, wspace= 0.15/7)
    '''
    
    
    '''
#   --------------------------------------------------------------
#   -------------------- put the stars back in -------------------
    
#   Good for debugging
    for iter in arange(num_stars):
            i = ii[iter]
            j = jj[iter]
            R[i-rmax2:i+rmax3,j-rmax2:j+rmax3] += stars[iter,:,:,0]
            G[i-rmax2:i+rmax3,j-rmax2:j+rmax3] += stars[iter,:,:,1]
            B[i-rmax2:i+rmax3,j-rmax2:j+rmax3] += stars[iter,:,:,2]
    
    starryimage = array(zeros((m,n,3), dtype=float32))
    starryimage[:,:,0] = R
    starryimage[:,:,1] = G
    starryimage[:,:,2] = B
    starryimage = starryimage - np.min(starryimage)
    starryimage = starryimage/np.max(starryimage)
    starryimage[starryimage<0] = 0
    
    figure(4, figsize=(18, 14), dpi=80)
    imshow(starryimage)
    '''
    
    '''
#   --------------------------------------------------------------
#   --------------------- if using Gaia data ---------------------
    
    hdul = fits.open(folder+file+'.fit')
#   print('hdul = ', hdul[0].header)
    if ('EQUINOX' not in hdul[0].header):
        print('send the fits file to astrometry.net')
    else:
        ctype1  = hdul[0].header['CTYPE1']  # TAN (gnomic) projection
        ctype2  = hdul[0].header['CTYPE2']  # TAN (gnomic) projection
        equinox = hdul[0].header['EQUINOX'] # Equatorial coordinates definition (yr)
        crval1  = hdul[0].header['CRVAL1']  # RA  of reference point
        crval2  = hdul[0].header['CRVAL2']  # DEC of reference point
        crpix1  = hdul[0].header['CRPIX1']  # X reference pixel
        crpix2  = hdul[0].header['CRPIX2']  # Y reference pixel
        cunit1  = hdul[0].header['CUNIT1']  # X pixel scale units
        cunit2  = hdul[0].header['CUNIT2']  # Y pixel scale units
        cd1_1   = hdul[0].header['CD1_1']   # Transformation matrix
        cd1_2   = hdul[0].header['CD1_2']   
        cd2_1   = hdul[0].header['CD2_1']   
        cd2_2   = hdul[0].header['CD2_2']   
        imagew  = hdul[0].header['IMAGEW']  
        imageh  = hdul[0].header['IMAGEH']  
        
        print('CTYPE1  = ', hdul[0].header['CTYPE1'])  
        print('CTYPE2  = ', hdul[0].header['CTYPE2'])  
        print('EQUINOX = ', hdul[0].header['EQUINOX']) 
        print('CRVAL1  = ', hdul[0].header['CRVAL1'])  
        print('CRVAL2  = ', hdul[0].header['CRVAL2'])  
        print('CRPIX1  = ', hdul[0].header['CRPIX1'])  
        print('CRPIX2  = ', hdul[0].header['CRPIX2'])  
        print('CUNIT1  = ', hdul[0].header['CUNIT1'])  
        print('CUNIT2  = ', hdul[0].header['CUNIT2'])  
        print('CD1_1   = ', hdul[0].header['CD1_1'])   
        print('CD1_2   = ', hdul[0].header['CD1_2'])   
        print('CD2_1   = ', hdul[0].header['CD2_1'])   
        print('CD2_2   = ', hdul[0].header['CD2_2'])   
        print('IMAGEW  = ', hdul[0].header['IMAGEW'])  
        print('IMAGEH  = ', hdul[0].header['IMAGEH'])  
        
        cd1_1 *= -3600
        cd1_2 *= -3600
        cd2_1 *= -3600
        cd2_2 *= -3600
    
        w = 250
        folder = 'fit'+slash
        sxv_psf = fits.getdata(folder+'sxv_psf4-RGB.fit').astype(float32)
        PSF = zeros((500,500,3), dtype='float32')
        PSF[:,:,0] = 1+sxv_psf[0,250-w:250+w,250-w:250+w]
        PSF[:,:,1] = 1+sxv_psf[1,250-w:250+w,250-w:250+w]
        PSF[:,:,2] = 1+sxv_psf[2,250-w:250+w,250-w:250+w]
        PSF = maximum(PSF - percentile(PSF,80),0)
        PSF = PSF/PSF.max()
        PSF_R = PSF
    
#   --------------------------------------------------------------
    '''
    
    '''
    depth = PIL.Image.open(folder+file+'_depth.png')
    figure(5, figsize=(18, 14), dpi=80)
    imshow(depth)
    '''
    
    
    
    '''
    close(2)
    close(3)
    close(4)
    close(5)
    
    show()
    '''

def deconvolve(evt):
    global folder, name, RGB
    global tab

    print(tab,'Deconvolving:')
    tab = tab + '    '

    RGB_RL = RLdeconvRGB(RGB, float(deconvFWHM.GetValue()))
    close(1)
    fig = slidershow(RGB_RL, folder=folder, fname=name+'-py-RL-RGB.png')
    show_fwhm(RGB_RL[:,:,1])

    tab = tab[0:-4]
    printdone('Done deconvolving')

def denoise(evt):
    global folder, name, RGB

    RGB_Denoise = Denoise(RGB, float(deconvFWHM.GetValue()))
    close(1)
    fig = slidershow(RGB_Denoise, folder=folder, fname=name+'-py-Denoise-RGB.png')
    show_fwhm(RGB_Denoise[:,:,1])
    
def hdr(evt):
    global folder, name, RGB

    RGB_HDR = HDR(RGB)
    close(1)
    fig = slidershow(RGB_HDR, folder=folder, fname=name+'-py-HDR-RGB.png')
#   show_fwhm(RGB_HDR[:,:,1])

def unsharpmask(evt):
    global folder, name, RGB

    RGB_UM = UnsharpMaskRGB(RGB)
    close(1)
    fig = slidershow(RGB_UM, folder=folder, fname=name+'-py-UM-RGB.png')
    show_fwhm(RGB_UM[:,:,1])
    
def GUI_HRdiag(evt):
    global name, RGB, mag0, mag1
    if (mag0.GetValue() != ''): minmag0 = float(mag0.GetValue())
    else:                       minmag0 = 10
    if (mag1.GetValue() != ''): maxmag0 = float(mag1.GetValue())
    else:                       maxmag0 = 20
    makeHRdiag_new(RGB, objname=name, minmag0=minmag0, shape='auto', fignum=1, thresh=0.00001, maxmag0=maxmag0)

def GUI_ColorDiag(evt):
    global name, RGB, mag0, mag1
    if (mag0.GetValue() != ''): minmag0 = float(mag0.GetValue())
    else:                       minmag0 = 10
    if (mag1.GetValue() != ''): maxmag0 = float(mag1.GetValue())
    else:                       maxmag0 = 20
    makeHRdiag_new(RGB, objname=name, minmag0=minmag0, shape='auto', fignum=6, thresh=0.00001, maxmag0=maxmag0)

def GUI_makeHRdiag0(evt):
#   global RGB, minmag0, maxmag0, objname
#   makeHRdiag_new(RGB, objname=objname, minmag0=minmag0, shape='auto', fignum=100, thresh=0.00001, maxmag0=maxmag0)
        
    global objname, coloridx, mag, absmag, abs_vs_app, parallax
    global plotter, fig5, fig5axes
    global images_frame
    global xpos

    if (abs_vs_app == 'Apparent'): mag0 = mag
    else:                          mag0 = absmag

    N = len(mag0)

    if (abs_vs_app == 'Absolute'):
            k = argmin(mag0)
            mag0[k] = mag0[N-1]
            coloridx[k] = coloridx[N-1]
            print('k = ',k, mag0[k], coloridx[k])
#           N = N-1
#           mag0     = mag0[0:N]
#           coloridx = coloridx[0:N]

    maxmag  = max(mag0)
    minmag  = min(mag0)
    maxmag1 = maxmag

    medcol = median(coloridx)
    clridx = 1.25*(coloridx-1.0)
    color = zeros((N,4),dtype=float32)
    diam  = zeros( N,   dtype=float32)
    alph  = zeros( N,   dtype=float32)
    for k in arange(1,N):
            diam[k] = 12*sqrt(min([mag0[k]-maxmag,0])/(minmag-maxmag))
            alph[k] = ((mag0[k]-maxmag)/(minmag-maxmag))
            col = 0.6*clridx[k]
            r = maximum(minimum(1 + col,1),0)
            b = maximum(minimum(1 - col,1),0)
            g = (r+b)/2
            color[k,0] = r
            color[k,1] = g
            color[k,2] = b
            color[k,3] = minimum(alph[k],1)

#   close(6)
#   figure(6, figsize=(10,10))
#   figure(6)
    '''
    if 'fig6' not in globals(): 
        fig6 = figure(6, figsize=(10, 10))
        fig6axes = fig6.add_subplot(1,1,1,facecolor='black')
    '''
    if (GUI == True):
        fig5 = 0   # delaxes isn't doing what I'd hoped
        if (plotter == ''): 
            images_frame = wx.Frame(None, -1, 'Images', size=[1000,880])
            images_frame.SetPosition(wx.Point(xpos, 26))
            plotter = PlotNotebook(images_frame)
            images_frame.Show(True)
        if (fig5 == 0): 
            fig5 = plotter.add('HR Diagram')
        else:
            print('fig5axes = ', fig5axes)
            fig5.delaxes(fig5axes)
        fig5axes = fig5.add_subplot(1,1,1,facecolor='black')
    else:
        close(fignum)
        fig5 = figure(fignum, figsize=(12,9))
#   fig6 = figure(6, figsize=(10, 10))
#   fig6axes = fig6.add_subplot(1,1,1,facecolor='black')

    fig5.set_facecolor('black')
    fig5axes.spines['bottom'].set_color('white')
    fig5axes.spines['left'].set_color('white')
    fig5axes.spines['top'].set_color('white')
    fig5axes.spines['right'].set_color('white')
    fig5axes.xaxis.label.set_color('white')
    fig5axes.yaxis.label.set_color('white')
    fig5axes.tick_params(axis='x', colors='white')
    fig5axes.tick_params(axis='y', colors='white')
    fig5axes.scatter(clridx, mag0, alpha=0.5, facecolor=color, s=diam, edgecolors='none')
    colorrange = minimum(maximum(max(clridx),-min(clridx)), 2.0)
    fig5axes.axis([-colorrange, colorrange, minmag-0.2, maxmag])
    fig5axes.set_xlabel('Color index (B-R)')
    if (abs_vs_app == 'Apparent'):
        fig5axes.set_ylabel('Brightness (apparent magnitude)')
    else:
        fig5axes.set_ylabel('Brightness (absolute magnitude)')
    fig5axes.set_title(objname+'   Hertzsprung-Russell Diagram using Gaia Data', color='white')
    fig5axes.invert_yaxis()
    if (GUI == False): draw()
#   show()
    
def GUI_makeColorDiag0(evt):
    global objname, coloridx, mag, absmag, abs_vs_app, parallax, gflux, rflux, bflux, ra, dec
    global plotter, fig6, fig6axes
    global images_frame
    global xpos

    mag0 = mag
    maxmag  = max(mag0)
    minmag  = min(mag0)

    rflux = rflux/median(rflux)
    gflux = gflux/median(gflux)
    bflux = bflux/median(bflux)

    N = len(mag)

    totflux = rflux+gflux+bflux

    color = zeros((N,4),dtype=float32)
    color[:,0] = rflux/totflux
    color[:,1] = gflux/totflux
    color[:,2] = bflux/totflux
    color[:,3] = 1   

    diam  = zeros( N,   dtype=float32)
    for k in arange(1,N):
            diam[k] = 12*sqrt(min([mag0[k]-maxmag,0])/(minmag-maxmag))

    if (GUI == True):
        fig6 = 0   # delaxes isn't doing what I'd hoped
        if (plotter == ''): 
            images_frame = wx.Frame(None, -1, 'Images', size=[1000,880])
            images_frame.SetPosition(wx.Point(xpos, 26))
            plotter = PlotNotebook(images_frame)
            images_frame.Show(True)
        if (fig6 == 0): 
            fig6 = plotter.add('Color Diagram')
        else:
            print('fig6axes = ', fig6axes)
            fig6.delaxes(fig6axes)
        fig6axes = fig6.add_subplot(1,1,1,facecolor='black')
    else:
        close(fignum)
        fig6 = figure(fignum, figsize=(12,9))
#   fig6 = figure(6, figsize=(10, 10))
#   fig6axes = fig6.add_subplot(1,1,1,facecolor='black')

    fig6.set_facecolor('black')
    fig6axes.spines['bottom'].set_color('white')
    fig6axes.spines['left'].set_color('white')
    fig6axes.spines['top'].set_color('white')
    fig6axes.spines['right'].set_color('white')
    fig6axes.xaxis.label.set_color('white')
    fig6axes.yaxis.label.set_color('white')
    fig6axes.tick_params(axis='x', colors='white')
    fig6axes.tick_params(axis='y', colors='white')
    dotsize = 4*sqrt(400/N)
#   scatter((gflux-rflux)/totflux, (bflux-gflux)/totflux, alpha=0.5, facecolor=color, s=4, edgecolors='none')
    fig6axes.scatter(log(gflux/rflux), log(bflux/gflux), alpha=1.0, facecolor=color, s=diam, edgecolors='none')

    xmin = min(log(gflux/rflux))
    xmax = max(log(gflux/rflux))
    ymin = min(log(bflux/gflux))
    ymax = max(log(bflux/gflux))
    xmin = xmin - 0.1*(xmax-xmin)
    ymin = ymin - 0.1*(ymax-ymin)
    xmax = xmax + 0.1*(xmax-xmin)
    ymax = ymax + 0.1*(ymax-ymin)
    '''
    c00 = array(( exp(-xmin), 1, exp(ymin) ))
    c01 = array(( exp(-xmin), 1, exp(ymax) ))
    c10 = array(( exp(-xmax), 1, exp(ymin) ))
    c11 = array(( exp(-xmax), 1, exp(ymax) ))
    print('colors: ', c00, c01, c10, c11)
    c00 = c00/max(c00)
    c01 = c01/max(c01)
    c10 = c10/max(c10)
    c11 = c11/max(c11)
    print('colors: ', c00, c01, c10, c11)
    '''
    fig6axes.axis([ xmin, xmax,  ymin, ymax])
    fig6axes.set_xlabel('Log(G/R)')
    fig6axes.set_ylabel('Log(B/G)')
    fig6axes.set_title(objname+'   Stellar Locus', color='white')
#   gca().invert_yaxis()
#   draw()
#   show()
    
def got_updated():
    show_fwhm()
    btn16.Enable()




#=======================================================

global btn0, btn4, btn6, btn7, btn8, btn10, btn11, btn12, btn13, btn14, btn15, btn16, btn17, btn18, btn20, btn21
global lbl1, lbl2, lbl3, lbl4, lbl5, lbl6, lbl6a, lbl7, lbl8, lbl8a, lbl9, lbl11, lbl_ra, lbl_dec
global progress, gauge
global option1, option2, option3, option4, option5, option6, option7, option8, option9
global objectid, objectra, objectdec, objectfov
global pixres, showalign, align_method, stack_method, flip, removehotpix, neatimage, deconvFWHM, deconvFWHMa, mag0, mag1
global thresh

gauge = 0

stack_method = 'Outliers Removed'

lbl9  = ''
lbl11 = ''

def set_align_method(event):
    global align_method
    align_method = option1.GetString(option1.GetSelection())

def set_showalign(event):
    global showalign 
    showalign = option2.GetString(option2.GetSelection())

def set_removehotpix(event):
    global removehotpix
    removehotpix = option7.GetString(option7.GetSelection())

def set_neatimage(event):
    global neatimage
    neatimage = option3.GetString(option3.GetSelection())

def set_flip(event):
    global flip
    flip = option4.GetString(option4.GetSelection())

def set_showimage(event):
    global showimage 
    showimage = option5.GetString(option5.GetSelection())

def set_stack_method(event):
    global stack_method 
    stack_method = option6.GetString(option6.GetSelection())

def set_image_type(event):
    global image_type
    image_type = option8.GetString(option8.GetSelection())

def set_mag_type(event):
    global abs_vs_app
    abs_vs_app = option9.GetString(option9.GetSelection())

def StackImages():
    global version
    global btn0, btn4, btn6, btn7, btn8, btn10, btn11, btn12, btn13, btn14, btn15, btn16, btn17, btn18, btn20, btn21
    global progress, gauge
    global lbl1, lbl2, lbl3, lbl4, lbl5, lbl6, lbl6a, lbl7, lbl8, lbl8a, lbl9, lbl11, lbl_ra, lbl_dec
    global option1, option2, option3, option4, option5, option6, option7
    global objectid, objectra, objectdec, objectfov
    global flat, flat_file
    global bias, bias_file
    global pixres, showalign, align_method, stack_method, flip, removehotpix, neatimage, deconvFWHM, deconvFWHMa, mag0, mag1
    global thresh
    global frame, panel
    global objectid, objectra, objectdec
    global sizer
    global GUI
    global app
    global tab
    global plotter
    global xpos
#   global images_frame

    print('ass.py version: ', version)
    print(tab, 'Stacking Images:')
    tab = tab+'    '

    GUI = True
    xpos = 640

    app = wx.App(redirect=False)
    appsize = minimum(array(wx.DisplaySize())-50,(xpos,880))
    print(tab,'appsize = ', appsize)

    frame = wx.Frame(None, title='Align, Stack, Show', size=appsize)
    frame.SetPosition(wx.Point(0, 26))
#   panel = wx.Panel(frame)
    panel = wx.ScrolledWindow(frame)
    sizer = wx.GridBagSizer(0,0)
        
    row = 0
    
    lbl1 = wx.StaticText(panel, label='')
    sizer.Add(lbl1, pos=(row,0), flag=wx.ALL, border=2)
    lbl1 = wx.StaticText(panel, label='')
    sizer.Add(lbl1, pos=(row,2), flag=wx.ALL, border=2)
    
    row += 1
    lbl1 = wx.StaticText(panel, label='')
    sizer.Add(lbl1, pos=(row,0), flag=wx.ALL, border=2)
    
    btn0 = wx.Button(panel, label='Select Images')  #, command = lambda:get_images())
    panel.Bind(wx.EVT_BUTTON, get_images, btn0)
    sizer.Add(btn0, pos=(row, 1), flag=wx.ALL, border=2)
    
    row += 1
    '''
    showimage = tk.StringVar()
    showimage.set('None')
    wx.StaticText(panel, label='Select one to show: '), pos=(row,0), flag=wx.ALL, border=2)
    option5 = tk.OptionMenu(frame, showimage, *fnames)
    option5.grid(row=row, column=3, sticky='nsew')
    '''
    sizer.Add(wx.StaticText(panel, label=' Select one to show: '), pos=(row, 0), flag=wx.ALL, border=2)
    option5 = wx.Choice(panel, choices = ['None                          '])
    sizer.Add(option5, pos=(row,1), flag=wx.ALL, border=2)
    option5.Bind(wx.EVT_CHOICE, set_showimage)
    showimage='None'
    option5.SetSelection(0)
    
    btn15 = wx.Button(panel, label='Show')  #, command = lambda:show_image())
    panel.Bind(wx.EVT_BUTTON, show_image, btn15)
    sizer.Add(btn15, pos=(row, 2), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row,0), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' Calibration...  Bias:  '), pos=(row,0), flag=wx.ALL, border=2)
    
    btn6 = wx.Button(panel, label='Select Bias Frame')  #, command = lambda:get_bias())
    panel.Bind(wx.EVT_BUTTON, get_bias, btn6)
    sizer.Add(btn6, pos=(row, 1), flag=wx.ALL, border=2)
    
    if path.exists('bias_file.txt'):
        bias_text = open('bias_file.txt','r')
        bias_file = bias_text.readline()
        bias_text.close()
        if (os.path.exists(bias_file)):
            bias = array( fits.open(bias_file)[0].data )
            btn6.SetLabel( '     DeSelect        ' )
        else:
            bias_file = ''
            btn6.SetLabel( 'Select Bias Frame' )
    else:
        bias_file = ''
        btn6.SetLabel( 'Select Bias Frame' )
#   lbl2 = wx.StaticText(panel, label=bias_file+' ')
#   lbl2 = wx.StaticText(panel, label=bias_file[0:5]+'...'+bias_file[-min([38,len(bias_file)]):]+' ')
#   lbl2 = wx.StaticText(panel, label='...'+bias_file[-min([43,len(bias_file)]):]+' ')
    lbl2 = wx.StaticText(panel, label=bias_file.split(slash)[-1]+' ')
    sizer.Add(lbl2, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label='                         Flat:  '), pos=(row,0), flag=wx.ALL, border=2)
    
    btn7 = wx.Button(panel, label='Select Flat Frame')  #, command = lambda:get_flat())
    panel.Bind(wx.EVT_BUTTON, get_flat, btn7)
    sizer.Add(btn7, pos=(row, 1), flag=wx.ALL, border=2)
    
    if path.exists('flat_file.txt'):
        flat_text = open('flat_file.txt','r')
        flat_file = flat_text.readline()
        flat_text.close()
        if (os.path.exists(flat_file)):
            flat = array( fits.open(flat_file)[0].data )
            btn7.SetLabel( '     DeSelect        ' )
        else:
            flat_file = ''
            btn7.SetLabel( 'Select Flat Frame' )
    else:
        flat_file = ''
        btn7.SetLabel( 'Select Flat Frame' )
#   lbl3 = wx.StaticText(panel, label=flat_file+' ')
    lbl3 = wx.StaticText(panel, label=flat_file.split(slash)[-1]+' ')
    sizer.Add(lbl3, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    #frame.grid_rowconfigure(0, weight=1)
    #frame.grid_rowconfigure(1, weight=1)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row,0), flag=wx.ALL, border=2)
    
    row += 1
    '''
    align_method = tk.StringVar()
    align_method.set('3 stars')
    wx.StaticText(panel, label='Align images using: ').grid(row=row, column=0)
    option1 = tk.OptionMenu(frame, align_method, \
            'no stars', '1 star', '2 stars', '3 stars', 'click on 1 star', 'click on 2 stars')
    option1.grid(row=row, column=1, sticky='nw')
    
    thresh = tk.StringVar()
    thresh.set('0.1')
    '''
    thresh = 0.1
    sizer.Add(wx.StaticText(panel, label=' Align images using: '), pos=(row, 0), flag=wx.ALL, border=2)
    option1 = wx.Choice(panel, choices = \
            ['no stars', '1 star', '2 stars', '3 stars', 'click on 1 star',
            'click on 2 stars' \
            , 'astrometric' \
#           , 'OpenCV' \
            ])
    sizer.Add(option1, pos=(row,1), flag=wx.ALL, border=2)
    option1.Bind(wx.EVT_CHOICE, set_align_method)
    align_method = '1 star'
    option1.SetSelection(1)

    
    row += 1
    '''
    showalign = tk.StringVar()
    showalign.set('No')
    wx.StaticText(panel, label='Show alignment frames? ').grid(row=row, column=0)
    option2 = tk.OptionMenu(frame, showalign, 'No', 'Yes')
    option2.grid(row=row, column=1, sticky='w')
    '''
    sizer.Add(wx.StaticText(panel, label=' Show alignment frames? '), pos=(row, 0), flag=wx.ALL, border=2)
    option2 = wx.Choice(panel, choices = ['No', 'Yes'])
    sizer.Add(option2, pos=(row,1), flag=wx.ALL, border=2)
    option2.Bind(wx.EVT_CHOICE, set_showalign)
    showalign = 'No'
    option2.SetSelection(0)
    
    row += 1
    '''
    neatimage = tk.StringVar()
    neatimage.set('No')
    wx.StaticText(panel, label='Apply NeatImage filter? ').grid(row=row, column=0)
    option3 = tk.OptionMenu(frame, neatimage, 'No', 'Yes')
    option3.grid(row=row, column=1, sticky='w')
    '''
    sizer.Add(wx.StaticText(panel, label=' Apply Neatimage filter? '), pos=(row, 0), flag=wx.ALL, border=2)
    option3 = wx.Choice(panel, choices = ['No', 'Yes'])
    sizer.Add(option3, pos=(row,1), flag=wx.ALL, border=2)
    option3.Bind(wx.EVT_CHOICE, set_neatimage)
    neatimage = 'No'
    option3.SetSelection(0)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' Remove hot pixels? '), pos=(row, 0), flag=wx.ALL, border=2)
    option7 = wx.Choice(panel, choices = ['No', 'Yes'])
    sizer.Add(option7, pos=(row,1), flag=wx.ALL, border=2)
    option7.Bind(wx.EVT_CHOICE, set_removehotpix)
    removehotpix = 'No'
    option7.SetSelection(0)
    
    row += 1
    '''
    flip = tk.StringVar()
    flip.set('Neither')
    wx.StaticText(panel, label='Flip/rotate? ').grid(row=row, column=0)
    option3 = tk.OptionMenu(frame, flip, 'Neither', 'Rotate 180', \
                    'Flip horizontally', 'Flip vertically')
    option3.grid(row=row, column=1, sticky='w')
    '''
    sizer.Add(wx.StaticText(panel, label=' Flip/mirror/rotate? '), pos=(row, 0), flag=wx.ALL, border=2)
    option4 = wx.Choice(panel, choices = \
            ['No', 'Flip vertically' , 'Mirror', 'Rotate 180'])
    sizer.Add(option4, pos=(row,1), flag=wx.ALL, border=2)
    option4.Bind(wx.EVT_CHOICE, set_flip)
    flip = 'No'
    option4.SetSelection(0)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' Stacking Method '), pos=(row, 0), flag=wx.ALL, border=2)
    option6 = wx.Choice(panel, choices = \
            ['Outliers Removed', 'Average', 'Max'])
    sizer.Add(option6, pos=(row,1), flag=wx.ALL, border=2)
    option6.Bind(wx.EVT_CHOICE, set_stack_method)
    stack_method = 'Outliers Removed'
    option6.SetSelection(0)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)
    
    row += 1
    btn8 = wx.Button(panel, label='Align/Stack/Show')  #, command = lambda:align_stack_show())
    panel.Bind(wx.EVT_BUTTON, align_stack_show, btn8)
    sizer.Add(btn8, pos=(row, 1), flag=wx.ALL, border=2)
    btn8.Disable()
    progress = wx.StaticText(panel, label='Progress: ')
    sizer.Add(progress, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
    gauge = wx.Gauge(panel, range=100, size=(200,15), style=wx.GA_HORIZONTAL)
    sizer.Add(gauge, pos=(row, 2), span=(1,3), flag=wx.ALL, border=2)
#   gauge.Hide()
    
    row += 1
    dashes  = '________________________________'
    dashes += '________________________________'
    dashes += '________________________________'
    dashes += '________________________________'
    sizer.Add(wx.StaticText(panel, label=dashes), pos=(row, 0), span=(1,5), flag=wx.ALL, border=2)
    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' Open an RGB fits file: '), pos=(row, 0), flag=wx.ALL, border=2)
    btn17 = wx.Button(panel, label='Select')  #, command = lambda:get_RGB_image())
    panel.Bind(wx.EVT_BUTTON, get_RGB_image, btn17)
    sizer.Add(btn17, pos=(row, 1), flag=wx.ALL, border=2)
    lbl7 = wx.StaticText(panel, label='')
    sizer.Add(lbl7, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
#   objectid = tk.StringVar()
    sizer.Add(wx.StaticText(panel, label=' Object Catalog Name: '), pos=(row, 0), flag=wx.ALL, border=2)
    objectid = wx.TextCtrl(panel, value='', size=(108,22))
    sizer.Add(objectid, pos=(row, 1), flag=wx.ALL, border=2)
    #objectid.bind("<Return>", lambda x: got_updated())
    btn16 = wx.Button(panel, label='Astrometric Analysis')  #, command = lambda:astrometric_analysis())
    panel.Bind(wx.EVT_BUTTON, astrometric_analysis, btn16)
    sizer.Add(btn16, pos=(row, 2), flag=wx.ALL, border=2)
    btn16.Disable()
    btn18 = wx.Button(panel, label='Show 3D')  #, command = lambda:unsharpmask())
    panel.Bind(wx.EVT_BUTTON, make3D, btn18)
    sizer.Add(btn18, pos=(row, 3), flag=wx.ALL, border=2)
    btn18.Disable()

    '''
    btn19 = wx.Button(panel, label='Sky Chart')  #, command = lambda:unsharpmask())
    panel.Bind(wx.EVT_BUTTON, GUI_SkyChart, btn19)
    sizer.Add(btn19, pos=(row, 3), flag=wx.ALL, border=2)
    btn19.Disable()
    '''
    
    row += 1
#   objectra = tk.StringVar()
    lbl_ra = wx.StaticText(panel, label=' Approximate RA: ')
    sizer.Add(lbl_ra, pos=(row, 0), flag=wx.ALL, border=2)
#   objectra = wx.TextCtrl(panel, value='')
    objectra = wx.TextCtrl(panel, value='', size=(108,22))
    sizer.Add(objectra, pos=(row, 1), flag=wx.ALL, border=2)
    #objectra.bind("<Return>", lambda x: got_updated())
    lbl4 = wx.StaticText(panel, label='')
    sizer.Add(lbl4, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
#   objectdec = tk.StringVar()
    lbl_dec = wx.StaticText(panel, label=' Approximate Dec: ')
    sizer.Add(lbl_dec, pos=(row, 0), flag=wx.ALL, border=2)
    objectdec = wx.TextCtrl(panel, value='', size=(108,22))
    sizer.Add(objectdec, pos=(row, 1), flag=wx.ALL, border=2)
    #objectdec.bind("<Return>", lambda x: got_updated())
    lbl5 = wx.StaticText(panel, label='')
    sizer.Add(lbl5, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
#   objectfov = tk.StringVar()
    lbl6a = wx.StaticText(panel, label=' Approximate FOV: ')
    sizer.Add(lbl6a, pos=(row, 0), flag=wx.ALL, border=2)
    objectfov = wx.TextCtrl(panel, value='0.5', size=(108,22))
    sizer.Add(objectfov, pos=(row, 1), flag=wx.ALL, border=2)
    lbl6 = wx.StaticText(panel, label='')
    sizer.Add(lbl6, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=dashes), pos=(row, 0), span=(1,5), flag=wx.ALL, border=2)
    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)
    
    row += 1
#   pixres = tk.StringVar()
#   sizer.Add(wx.StaticText(panel, label=' Arcseconds per pixel: '), pos=(row, 0), flag=wx.ALL, border=2)
    sizer.Add(wx.StaticText(panel, label=' Pixel size (arcsec/pix): '), pos=(row, 0), flag=wx.ALL, border=2)
    pixres = wx.TextCtrl(panel, value='', size=(108,22))
    sizer.Add(pixres, pos=(row, 1), flag=wx.ALL, border=2)
    btn12 = wx.Button(panel, label='Deconvolve')  #, command = lambda:deconvolve())
    panel.Bind(wx.EVT_BUTTON, deconvolve, btn12)
    sizer.Add(btn12, pos=(row, 2), flag=wx.ALL, border=2)
    btn12.Disable()
    '''
    '''
    btn13 = wx.Button(panel, label='Unsharp Mask')  #, command = lambda:unsharpmask())
    panel.Bind(wx.EVT_BUTTON, unsharpmask, btn13)
    sizer.Add(btn13, pos=(row, 3), flag=wx.ALL, border=2)
    btn13.Disable()

    '''
    btn21 = wx.Button(panel, label='HDR Stretch')  
    panel.Bind(wx.EVT_BUTTON, hdr, btn21)
    sizer.Add(btn21, pos=(row+1, 3), flag=wx.ALL, border=2)
    btn21.Disable()
    '''
    
    row += 1
    lbl8a = wx.StaticText(panel, label=' Current FWHM (pixels): ')
    sizer.Add(lbl8a, pos=(row, 0), flag=wx.ALL, border=2)
    lbl8 = wx.StaticText(panel, label='')
    sizer.Add(lbl8, pos=(row,1), flag=wx.ALL, border=2)
    
    row += 1
#   deconvFWHM = tk.StringVar()
    deconvFWHMa = wx.StaticText(panel, label=' Target FWHM (pixels):  ')
    sizer.Add(deconvFWHMa, pos=(row, 0), flag=wx.ALL, border=2)
    deconvFWHM = wx.TextCtrl(panel, value='1.6', size=(108,22))
    sizer.Add(deconvFWHM, pos=(row, 1), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=dashes), pos=(row, 0), span=(1,5), flag=wx.ALL, border=2)
    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' Brightest star (mag): '), pos=(row, 0), flag=wx.ALL, border=2)
    mag0 = wx.TextCtrl(panel, value='', size=(108,22))
    sizer.Add(mag0, pos=(row, 1), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' Faintest star (mag): '), pos=(row, 0), flag=wx.ALL, border=2)
    mag1 = wx.TextCtrl(panel, value='20', size=(108,22))
    sizer.Add(mag1, pos=(row, 1), flag=wx.ALL, border=2)

    btn14 = wx.Button(panel, label='HR diagram')  #, command = lambda:unsharpmask())
    panel.Bind(wx.EVT_BUTTON, GUI_HRdiag, btn14)
    sizer.Add(btn14, pos=(row, 2), flag=wx.ALL, border=2)
    btn14.Disable()
    '''
    btnQuit = wx.Button(panel, label='Quit')  #, command = lambda:unsharpmask())
    panel.Bind(wx.EVT_BUTTON, quit, btnQuit)
    sizer.Add(btnQuit, pos=(row, 4), flag=wx.ALL, border=2)
    '''
    
    btn20 = wx.Button(panel, label='Color Diagram')  #, command = lambda:unsharpmask())
    panel.Bind(wx.EVT_BUTTON, GUI_ColorDiag, btn20)
    sizer.Add(btn20, pos=(row, 3), flag=wx.ALL, border=2)
    btn20.Disable()

    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)
    
    row += 1
    lbl9 = wx.StaticText(panel, label='')
    sizer.Add(lbl9, pos=(row,0), span=(1,5), flag=wx.ALL, border=2)
    
    row += 1
    lbl11 = wx.StaticText(panel, label=' Star Info: ')
    sizer.Add(lbl11, pos=(row,0), span=(3,5), flag=wx.ALL, border=2)
    
    row += 3
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)

    panel.SetSizerAndFit(sizer)
#   panel.SetScrollbars(20,20,55,40)
    panel.SetScrollbars(1,1,1,1)

    frame.Show(True)
    '''
    frame.Raise() # doesn't seem to help
    panel.Raise() # doesn't seem to help
    '''
#   panel.Center()
#   panel.Show()
#   panel.Fit()




    '''
    images_frame = wx.Frame(None, -1, 'Images', size=[1000,880])
#   images_frame.SetPosition(wx.Point(xpos, 26))
    plotter = PlotNotebook(images_frame)
    images_frame.Show(True)

    class MainFrame(wx.Frame):
        def __init__(self, parent, id, title):
            wx.Frame.__init__(self, None, -1, title, size=(1700,800))

            splitter = wx.SplitterWindow(self, -1)
            panel1 = wx.Panel(splitter, -1)
            b1 = wx.BoxSizer(wx.HORIZONTAL)
            b1.Add(frame, 1, wx.EXPAND)
            panel1.SetSizerAndFit(b1)

            panel2 = wx.Panel(splitter, -1)
            b2 = wx.BoxSizer(wx.HORIZONTAL)
            b2.Add(images_frame ,1)
            panel2.SetSizer(b2)
                
            splitter.SplitVertically(panel1,panel2)
            self.Centre()
            self.Show()
            panel1.Show()
            panel2.Show()

    mainframe = MainFrame(None, -1, 'wx.SplitterWindow')
    mainframe.Show(True)
    '''


    app.MainLoop()

    tab = tab[0:-4]
    printdone( 'Done stacking images')

global abs_vs_app
abs_vs_app = 'Apparent'

def printdone(itsdone):
    global tab
    printing = True
    if (printing): print(tab,itsdone)

def ShowClusters():
    global version
    global btn0, btn4, btn6, btn7, btn8, btn10, btn11, btn12, btn13, btn14, btn15, btn16, btn17, btn20
    global progress, gauge
    global lbl1, lbl2, lbl3, lbl4, lbl5, lbl6, lbl7, lbl8, lbl8a, lbl9, lbl10, lbl11
    global option8, option9
    global objectid, objectra, objectdec, objectfov, objectmaxmag
    global flat, flat_file
    global bias, bias_file
    global pixres, showalign, align_method, flip, removehotpix, neatimage, deconvFWHM, deconvFWHMa, mag0, mag1
    global thresh
    global frame, panel
    global objectid, objectra, objectdec
    global sizer
    global GUI
    global app
    global xpos

    print('ass.py version: ', version)

    GUI = True
    xpos = 360

    app = wx.App(redirect=False)
    appsize = minimum(array(wx.DisplaySize())-50,(xpos,450))
    frame = wx.Frame(None, title='Show Clusters w/ Gaia Data', size=appsize)
    frame.SetPosition(wx.Point(0, 26))
#   panel = wx.Panel(frame)
    panel = wx.ScrolledWindow(frame)
    sizer = wx.GridBagSizer(0,0)
        
    row = 0
    
    lbl1 = wx.StaticText(panel, label='')
    sizer.Add(lbl1, pos=(row,0), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label='Object Catalog Name: '), pos=(row, 0), flag=wx.EXPAND|wx.ALL, border=2)
    objectid = wx.TextCtrl(panel, style=wx.TE_PROCESS_ENTER, value='', size=(108,22))
#   objectid.Bind(wx.EVT_TEXT_ENTER, find_image)
#   objectid.Bind(wx.EVT_TEXT, find_image)
    objectid.Bind(wx.EVT_KEY_DOWN, find_image)
    objectid.SetFocus()
    sizer.Add(objectid, pos=(row, 1), flag=wx.ALL, border=2)
    '''
    btn16 = wx.Button(panel, label='Find')  #, command = lambda:find_image())
    panel.Bind(wx.EVT_BUTTON, find_image, btn16)
    sizer.Add(btn16, pos=(row, 2), flag=wx.ALL, border=2)
    '''
    
    row += 1
    sizer.Add(wx.StaticText(panel, label='Right Ascension (RA): '), pos=(row, 0), flag=wx.ALL, border=2)
    objectra = wx.TextCtrl(panel, value='', size=(108,22))
    objectra.SetValue('+00 : 00 : 00')
    objectra.Bind(wx.EVT_KEY_UP, find_image)
    sizer.Add(objectra, pos=(row, 1), flag=wx.ALL, border=2)
    #objectra.bind("<Return>", lambda x: got_updated())
    lbl2 = wx.StaticText(panel, label='h:m:s ')
    sizer.Add(lbl2, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label='Declination (Dec): '), pos=(row, 0), flag=wx.ALL, border=2)
    objectdec = wx.TextCtrl(panel, value='', size=(108,22))
    objectdec.SetValue('+00 : 00 : 00')
    sizer.Add(objectdec, pos=(row, 1), flag=wx.ALL, border=2)
    #objectdec.bind("<Return>", lambda x: got_updated())
    lbl5 = wx.StaticText(panel, label='d:m:s ')
    sizer.Add(lbl5, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label='Field of View (FOV): '), pos=(row, 0), flag=wx.ALL, border=2)
    objectfov = wx.TextCtrl(panel, value='0.3', size=(108,22))
    sizer.Add(objectfov, pos=(row, 1), flag=wx.ALL, border=2)
    lbl6 = wx.StaticText(panel, label='degrees ')
    sizer.Add(lbl6, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label='Maximum Magnitude: '), pos=(row, 0), flag=wx.ALL, border=2)
    objectmaxmag = wx.TextCtrl(panel, value='16', size=(108,22))
    sizer.Add(objectmaxmag, pos=(row, 1), flag=wx.ALL, border=2)
    lbl7 = wx.StaticText(panel, label='')
    sizer.Add(lbl7, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)

    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label='Select Image Type: '), pos=(row, 0), flag=wx.ALL, border=2)
    option8 = wx.Choice(panel, choices = ['2D', '3D', 'future'])
    sizer.Add(option8, pos=(row,1), flag=wx.ALL, border=2)
    option8.Bind(wx.EVT_CHOICE, set_image_type)
    image_type = '2D'
    option8.SetSelection(0)

    row += 1
    btn15 = wx.Button(panel, label='Show Stars')  #, command = lambda:show_gaia_image())
    panel.Bind(wx.EVT_BUTTON, show_gaia_image, btn15)
    sizer.Add(btn15, pos=(row, 1), flag=wx.ALL, border=2)
    lbl8 = wx.StaticText(panel, label='')
    sizer.Add(lbl8, pos=(row,2), span=(1,3), flag=wx.ALL, border=2)

    row += 1
    lbl4 = wx.StaticText(panel, label='', size=(108,22))
    sizer.Add(lbl4, pos=(row,0), span=(1,4), flag=wx.ALL, border=2)
    
    row += 1
    lbl10 = wx.StaticText(panel, label='', size=(108,22))
    sizer.Add(lbl10, pos=(row,0), span=(1,4), flag=wx.ALL, border=2)

    '''
    row += 1
    sizer.Add(wx.StaticText(panel, label='Brightest star (mag): '), pos=(row, 0), flag=wx.ALL, border=2)
    mag0 = wx.TextCtrl(panel, value='')
    sizer.Add(mag0, pos=(row, 1), flag=wx.ALL, border=2)
    '''
    
    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)
    
    row += 1
    sizer.Add(wx.StaticText(panel, label='Select Magnitude Type: '), pos=(row, 0), flag=wx.ALL, border=2)
    option9 = wx.Choice(panel, choices = ['Apparent', 'Absolute'])
    sizer.Add(option9, pos=(row,1), flag=wx.ALL, border=2)
    option9.Bind(wx.EVT_CHOICE, set_mag_type)
    abs_vs_app = 'Apparent'
    option9.SetSelection(0)

    row += 1
    btn14 = wx.Button(panel, label='HR Diagram')  #, command = lambda:unsharpmask())
    panel.Bind(wx.EVT_BUTTON, GUI_makeHRdiag0, btn14)
    sizer.Add(btn14, pos=(row, 1), flag=wx.ALL, border=2)
    btn14.Disable()

    row += 1
    sizer.Add(wx.StaticText(panel, label=' '), pos=(row, 0), flag=wx.ALL, border=2)
    
    row += 1
    btn20 = wx.Button(panel, label='Color Diagram')  #, command = lambda:unsharpmask())
    panel.Bind(wx.EVT_BUTTON, GUI_makeColorDiag0, btn20)
    sizer.Add(btn20, pos=(row, 1), flag=wx.ALL, border=2)
    btn20.Disable()

    panel.SetSizerAndFit(sizer)
    panel.SetScrollbars(1,1,1,1)

    frame.Show(True)
#   panel.Center()
#   panel.Show()
#   panel.Fit()
    
    app.MainLoop()

#StackImages()
