# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 14:35:49 2024

@author: staechelin
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit




import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.linewidth'] = 2.0
mpl.rcParams['grid.linestyle'] = "--"
mpl.rcParams['grid.alpha'] = 0.75
mpl.rcParams['font.size'] = '16'
mpl.rcParams['legend.framealpha'] = '1'
mpl.rcParams['figure.figsize'] = (4.8, 4.8)




class Raman():
    
    
    
    def __init__(self, path, file):
        """ 
        initialize object by loading data
        """
        self.path = path
        self.file = file
        self.sample = ""
        self.data = np.genfromtxt(f"{path}\\{file}", delimiter = None)
        self.shift = self.data[::-1,0]
        self.spec = self.data[::-1,1]
        self.load_Si_spec()
        
        
        
        
        self.smooth()
        
    def smooth(self):
        
        kernel_size = 3
        kernel = np.ones(kernel_size)/kernel_size
        data_smoothed = np.convolve(self.spec, kernel, mode="same")
        
        
        
        
        # plt.plot(self.shift, self.spec)
        # plt.plot(self.shift, data_smoothed)
        
        self.spec = data_smoothed
        
    def plot(self, xlim=[]):
        """
        plot Raman spectrum
        """
        
        fig, ax = plt.subplots()
        ax.plot(self.shift, self.spec)
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_xlabel(r"Raman Shift (cm$^{-1}$)")
        
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        plt.grid()
        plt.show()
        
        
        
    def load_Si_spec(self):
        
        Si_data = np.genfromtxt(r"L:\Yannic\Messdaten\Raman\20240228\Si chip\Si_spectrum.txt", 
                                delimiter = None)
        
        self.Si_spec = Si_data[::-1,1]
        
        kernel_size = 3
        kernel = np.ones(kernel_size)/kernel_size
        data_smoothed = np.convolve(self.Si_spec, kernel, mode="same")
        
        self.Si_spec = data_smoothed
        
        self.Si_spec = self.Si_spec/np.max(self.Si_spec)
        
    def plot_Si_spec(self, xlim=[]):
        
        fig, ax = plt.subplots()
        ax.plot(self.shift, self.Si_spec)
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_xlabel(r"Raman Shift (cm$^{-1}$)")
        
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        plt.grid()
        plt.show()
        
        
        
        
        
    def crop(self, xlim=[]):
        
        idx1 = next(i for i,v in enumerate(self.shift) if v > xlim[0])
        idx2 = next(i for i,v in enumerate(self.shift) if v > xlim[1])
        
        self.spec = self.spec[idx1:idx2]
        self.shift = self.shift[idx1:idx2]
        
        
        
        
        
        
        
        
    def subtract_Si(self, plot=False):
        
        value_Si = self.spec[962]*0.95
        
        Si_spec_scaled = self.Si_spec*value_Si
        
        self.spec = self.spec - Si_spec_scaled
        
        if plot == True:
                
            plt.figure()
            plt.plot(self.shift, self.spec)
            plt.ylabel("Intensity (a.u.)")
            plt.xlabel("Raman Shift (cm$^{-1}$)")
            plt.grid()
            plt.show()
            
        
    def fit(self, guesses, xlim=[], plot_bulk=[], export=True):
        
        idx1 = next(i for i,v in enumerate(self.shift) if v > xlim[0])
        idx2 = next(i for i,v in enumerate(self.shift) if v > xlim[1])
        
        
        shift_crop = self.shift[idx1:idx2]
        spec_crop = self.spec[idx1:idx2]    
    
        
    
    
        # define guesses
        cguess=[0]
        guezz = []
        for guess in guesses:
            idx = next(i for i,v in enumerate(shift_crop) if v > guess)
            guez = [spec_crop[idx], shift_crop[idx], 50]
            guezz = guezz + guez
        guezz = guezz + cguess 
    
    
    
        # define boundaries
        upper_bounds = []
        lower_bounds = []
        for guess in guesses:
            lower = [0,guess-30,0]
            lower_bounds = lower_bounds + lower
            # upper = [np.inf, np.inf, np.inf]
            upper = [np.max(spec_crop), guess+30, np.inf]
            upper_bounds = upper_bounds + upper
        
        lower_bounds = lower_bounds + [0]    
        upper_bounds = upper_bounds + [np.inf]
            
        boundz = (lower_bounds, upper_bounds)
    
        # print(boundz)
    
        # define which fit function to use
        if len(guesses) == 1:
            def lorentzian(x, amp1,cen1,wid1, c):
                return (amp1*wid1**2/((x-cen1)**2+wid1**2)) + c
    
        elif len(guesses) == 2:
            def lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, c):
                return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                    (amp2*wid2**2/((x-cen2)**2+wid2**2)) + c 
        
        elif len(guesses) == 3:
            def lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, amp3,cen3,wid3, c):
                return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                    (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                    (amp3*wid3**2/((x-cen3)**2+wid3**2)) + c
        
        elif len(guesses) == 4:
            def lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, amp3,cen3,wid3, amp4,cen4,wid4, c):
                return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                    (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                    (amp3*wid3**2/((x-cen3)**2+wid3**2)) +\
                    (amp4*wid4**2/((x-cen4)**2+wid4**2)) + c   

        elif len(guesses) == 5:
            def lorentzian(x, amp1,cen1,wid1, amp2,cen2,wid2, amp3,cen3,wid3, amp4,cen4,wid4, amp5,cen5,wid5, c):
                return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
                    (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                    (amp3*wid3**2/((x-cen3)**2+wid3**2)) +\
                    (amp4*wid4**2/((x-cen4)**2+wid4**2)) +\
                    (amp5*wid5**2/((x-cen5)**2+wid5**2)) + c   



        # print(guezz)
        popt, pcov = curve_fit(lorentzian, shift_crop, spec_crop, p0 = guezz, bounds = boundz)
        fit = lorentzian(shift_crop, *popt)
        
        # plt.figure()
        # plt.plot(shift_crop, spec_crop, alpha=0.5)
        # plt.plot(shift_crop, fit, color="tab:blue")
        # plt.show()
        
        
        
        
        def single_lorentzian(x, amp,cen,wid):
            return (amp*wid**2/((x-cen)**2+wid**2))
        
        
        plt.figure(figsize=[4.8,2.4])
        plt.plot(shift_crop, spec_crop, alpha=1)
        plt.plot(shift_crop, fit, color="tab:blue", alpha=1)
        plt.fill_between(shift_crop, fit, alpha=0.2)
        # plt.plot(shift_crop, np.full((len(shift_crop)),popt[-1]), color="black")
        plt.plot(shift_crop, single_lorentzian(shift_crop, popt[0], popt[1], popt[2]), 
                 color="tab:blue", alpha=0.5)
        if len(guesses) == 2 :
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[3], popt[4], popt[5]), 
                     color="tab:blue", alpha=0.5)
        if len(guesses) == 3:
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[3], popt[4], popt[5]), 
                     color="tab:blue", alpha=0.5)
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[6], popt[7], popt[8]), 
                     color="tab:blue", alpha=0.5)
        if len(guesses) == 4:
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[3], popt[4], popt[5]), 
                     color="tab:blue", alpha=0.5)
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[6], popt[7], popt[8]), 
                     color="tab:blue", alpha=0.5)
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[9], popt[10], popt[11]), 
                     color="tab:blue", alpha=0.5)
        if len(guesses) == 5:
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[3], popt[4], popt[5]), 
                     color="tab:blue", alpha=0.5)
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[6], popt[7], popt[8]), 
                     color="tab:blue", alpha=0.5)
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[9], popt[10], popt[11]), 
                     color="tab:blue", alpha=0.5)
            plt.plot(shift_crop, single_lorentzian(shift_crop, popt[12], popt[13], popt[14]), 
                     color="tab:blue", alpha=0.5)  
        # plt.grid(True)
        plt.xlabel(r"Raman Shift (cm$^{-1}$)")
        plt.yticks([])
        # plt.yaxis("off")
        if plot_bulk != []:
            for value in plot_bulk:
                plt.vlines(value, 0, np.max(spec_crop), color="tab:grey")
        plt.show()
        
        
        
        fit_export = np.zeros((len(shift_crop), len(guesses)+3))
        fit_export[:,0] = shift_crop
        fit_export[:,1] = spec_crop
        fit_export[:,2] = fit
        
        print(fit_export)
        
        if len(guesses) == 1:
            print(f"Peak No 1: center: {np.round(popt[1],1)} cm-1, width: {np.round(2*popt[2],1)} cm-1, amplitude: {np.round(popt[0],1)}")
            fit_export[:,3] = single_lorentzian(shift_crop, popt[0], popt[1], popt[2])
        elif len(guesses) == 2:
            print(f"Peak No 1: center: {np.round(popt[1],1)} cm-1, width: {np.round(2*popt[2],1)} cm-1, amplitude: {np.round(popt[0],1)}")
            print(f"Peak No 2: center: {np.round(popt[4],1)} cm-1, width: {np.round(2*popt[5],1)} cm-1, amplitude: {np.round(popt[3],1)}")
            fit_export[:,3] = single_lorentzian(shift_crop, popt[0], popt[1], popt[2])
            fit_export[:,4] = single_lorentzian(shift_crop, popt[3], popt[4], popt[5])
        elif len(guesses) == 3:
            print(f"Peak No 1: center: {np.round(popt[1],1)} cm-1, width: {np.round(2*popt[2],1)} cm-1, amplitude: {np.round(popt[0],1)}")
            print(f"Peak No 2: center: {np.round(popt[4],1)} cm-1, width: {np.round(2*popt[5],1)} cm-1, amplitude: {np.round(popt[3],1)}")
            print(f"Peak No 3: center: {np.round(popt[7],1)} cm-1, width: {np.round(2*popt[8],1)} cm-1, amplitude: {np.round(popt[6],1)}")
            fit_export[:,3] = single_lorentzian(shift_crop, popt[0], popt[1], popt[2])
            fit_export[:,4] = single_lorentzian(shift_crop, popt[3], popt[4], popt[5])
            fit_export[:,5] = single_lorentzian(shift_crop, popt[6], popt[7], popt[8])
        elif len(guesses) == 4:
            print(f"Peak No 1: center: {np.round(popt[1],1)} cm-1, width: {np.round(2*popt[2],1)} cm-1, amplitude: {np.round(popt[0],1)}")
            print(f"Peak No 2: center: {np.round(popt[4],1)} cm-1, width: {np.round(2*popt[5],1)} cm-1, amplitude: {np.round(popt[3],1)}")
            print(f"Peak No 3: center: {np.round(popt[7],1)} cm-1, width: {np.round(2*popt[8],1)} cm-1, amplitude: {np.round(popt[6],1)}")
            print(f"Peak No 4: center: {np.round(popt[10],1)} cm-1, width: {np.round(2*popt[11],1)} cm-1, amplitude: {np.round(popt[9],1)}")
            fit_export[:,3] = single_lorentzian(shift_crop, popt[0], popt[1], popt[2])
            fit_export[:,4] = single_lorentzian(shift_crop, popt[3], popt[4], popt[5])
            fit_export[:,5] = single_lorentzian(shift_crop, popt[6], popt[7], popt[8])
            fit_export[:,6] = single_lorentzian(shift_crop, popt[9], popt[10], popt[11])
        elif len(guesses) == 5:
            print(f"Peak No 1: center: {np.round(popt[1],1)} cm-1, width: {np.round(2*popt[2],1)} cm-1, amplitude: {np.round(popt[0],1)}")
            print(f"Peak No 2: center: {np.round(popt[4],1)} cm-1, width: {np.round(2*popt[5],1)} cm-1, amplitude: {np.round(popt[3],1)}")
            print(f"Peak No 3: center: {np.round(popt[7],1)} cm-1, width: {np.round(2*popt[8],1)} cm-1, amplitude: {np.round(popt[6],1)}")
            print(f"Peak No 4: center: {np.round(popt[10],1)} cm-1, width: {np.round(2*popt[11],1)} cm-1, amplitude: {np.round(popt[9],1)}")
            print(f"Peak No 5: center: {np.round(popt[13],1)} cm-1, width: {np.round(2*popt[14],1)} cm-1, amplitude: {np.round(popt[12],1)}")
            fit_export[:,3] = single_lorentzian(shift_crop, popt[0], popt[1], popt[2])
            fit_export[:,4] = single_lorentzian(shift_crop, popt[3], popt[4], popt[5])
            fit_export[:,5] = single_lorentzian(shift_crop, popt[6], popt[7], popt[8])
            fit_export[:,6] = single_lorentzian(shift_crop, popt[9], popt[10], popt[11])
            fit_export[:,7] = single_lorentzian(shift_crop, popt[12], popt[13], popt[14])
        print(f"c: {np.round(popt[-1],1)}")
        
        
        if export == True:
           np.savetxt(f"{self.path}\\{self.file[:-4]}_fit.txt", fit_export, delimiter="\t")     
        
        
        
class RamanRow():
    
    
    def __init__(self):
        self.measurements = []
        
    def add_measurements(self, measurements):
        """
        add previously imported measurement to measurement row
        import multiple measurements at a time
        """
        for measurement in measurements:
            self.measurements.append(measurement)    
        
        
        
    def compare_spectra(self, xlim=[]):
        
        fig, ax = plt.subplots(figsize=(4.8,len(self.measurements)))
        
        num = len(self.measurements)
        
        
        r = np.linspace(31, 214, num)
        g = np.linspace(119, 39, num)
        b = np.linspace(180, 40, num)
        
        for idx, meas in enumerate(self.measurements):
            
            
            ax.hlines(y=num-idx*1, xmin=meas.shift[0], xmax=meas.shift[-1], 
                      linewidth=0.5, linestyle="--", alpha=0.65, color="grey")
            ax.plot(meas.shift, num-idx*1+(meas.spec/np.max(meas.spec)), label = meas.sample,
                    color=(r[idx]/255, g[idx]/255, b[idx]/255))
            
            
        ax.set_ylabel("Intensity (a.u.)")
        ax.set_xlabel("Raman Shift (cm$^{-1}$)")
        ax.set_yticks([])
        
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        
        ax.legend( prop = { "size": 13 }, bbox_to_anchor=(1, 0.5))
        
        plt.grid()
        plt.show()
        
        # mpl.rcParams['axes.linewidth'] = 2.0
        # mpl.rcParams['grid.linestyle'] = "--"
        # mpl.rcParams['grid.alpha'] = 0.75
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        