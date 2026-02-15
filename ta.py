
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 19:33:40 2021

@author: Yannic
"""

"""

-------------------------------------------------------------------
-                                                                 -
-                                                                 -
-      Data analysis of TA data                                   -
-        YS 2022                                                  -
-                                                                 -
-                                                                 -
-------------------------------------------------------------------

"""




import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date
from scipy.optimize import curve_fit
import scipy.constants as sc
from scipy.interpolate import interp1d
today = date.today()
today = today.strftime("%y%m%d")
import matplotlib as mpl

mpl.rcParams['grid.linestyle'] = "--"
mpl.rcParams['grid.alpha'] = 0.75
mpl.rcParams['font.size'] = '16'
mpl.rcParams['legend.framealpha'] = '1'

class TAData:
    
    def __init__(self, path, measurement):
        """ 
        initialize object by loading data
        """
        self.path = path
        self.measurement = measurement
        self.data = self.load_data(path, measurement)
        self.metadata, self.data = self.sep_metadata(self.data)
        self.data = np.nan_to_num(self.data)
        self.wl, self.tau = self.sep_wl_tau(self.data)
        self.data = self.data[1:,1:]
        self.break_metadata()
    
    
    def load_data(self, path, measurement):
        """
        load data from measurement file

        """
        load_data_from = path + "\\" + measurement + ".csv"
        self.data = pd.read_csv(load_data_from, sep = ",", header = None, encoding='latin-1')
        self.data = self.data.to_numpy()
        return self.data
        
    
    
    def sep_metadata(self, data):
        """ 
        separate metadata from actual data
        
        to do:
            nichts tun, falls "file" nicht in der Datei vorkommt
        """
        for i in range(len(self.data[:,0])):
            if self.data[i, 0][0:4] == "file":
                self.metadata = self.data[i:,:]
                self.data = self.data[:i, :]
                self.data = np.array(self.data, dtype = np.float)
                break
        return self.metadata, self.data        

    
    def break_metadata(self):
        """ 
        break metadata up into own self.variable for sample, date etc
        """
        if self.metadata[1,0][:7] == "Sample:":
            self.sample = self.metadata[1,0][8:]
        if self.metadata[3,0][:11] == "Pump energy":
            self.pumpfluence = self.metadata[3,0][18:]
        if self.metadata[4,0][:15] == "Pump wavelength":
            self.pumpwavelength = self.metadata[4,0][22:]
        self.find_timeunit()
    
    def edit_metadata(self, new_sample_name):
        """
        edit metadata ... make it possible to edit every single entry... save edited metadata to original measurement .csv file
        
        """
        self.sample = new_sample_name
        
        
    def find_timeunit(self):
        """
        find time unit of measurement in metadata to distinguish EOS and HELIOS measurements 
        """

        for entry in self.metadata[:,0]:
            if entry[:11] == "Time units:":
                self.timeunit = entry[-2:]
        if self.timeunit == "us":
            self.timeunit = "$\mathrm{\mu}$s"
            self.ta_type = "EOS"
        if self.timeunit == "ps":
            self.ta_type = "HELIOS"
    
    def convert_timeunit_to_ns(self):
        """
        convert timeunit of EOS measurement from us to ns
        """
        if self.timeunit == "us" or self.timeunit == "$\mathrm{\mu}$s":
            self.tau = self.tau*10**3
            self.timeunit = "ns"
            
        elif self.timeunit == "ps":
            print("Only possible for EOS measurements.")
            
        elif self.timeunit == "ns":
            print("Already converted to ns.")

    
    def sep_wl_tau(self, data):
        """
        get wavelengths and delay times from data
        """
        self.wl = self.data[1:,0]
        self.tau = self.data[0, 1:]
        return self.wl, self.tau
    
    def get_sample_name(self):
        return self.sample
    
    def get_pumpfluence(self):
        return self.pumpfluence
    
    def get_pumpwavelength(self):
        return self.pumpwavelength
    
    def get_wl(self):
        return self.wl
    
    def get_tau(self):
        return self.tau
    
    def get_data(self):
        return self.data
    
    def get_metadata(self):
        return self.metadata
    
    def get_one_dynamic(self, wl_to_get):
        """
        export the dynamic of one specific probe wavelength
        """
        index  = next(i for i,v in enumerate(self.wl) if v > wl_to_get)
        return index, self.data[index,:]
    
    def get_spectra(self, taus_to_get):
        """ 
        export spectra at specific delay times
        """
        if type(taus_to_get) == int:
            taus_to_get = [taus_to_get]
        
        indices = []
        spectra = []
        for tau in taus_to_get:
            index = next(i for i,v in enumerate(self.tau) if v > tau)
            indices.append(index)
            spectra.append(self.data[:, index])
        return spectra
    
    
    def get_one_spectrum(self, tau_to_get):
        """ 
        export spectrum at specific delay time
        """
        
        index = next(i for i,v in enumerate(self.tau) if v > tau_to_get)
        spectrum = self.data[:, index]
        return spectrum        
        
    
    
    def get_deltaA_value(self, delaytime, wavelength):
        """
        return delta A value at specified delaytime and wavelength
        """
        index1 = next(i for i,v in enumerate(self.tau) if v > delaytime)
        index2 = next(i for i,v in enumerate(self.wl) if v > wavelength)
        
        deltaA_value = self.data[index2, index1]
        
        return deltaA_value
    
    
    def get_deltaA_value_av(self, delaytime_range, wavelength_range):
        """
        return delta A value averaged over a specified delay time range and wavelength range

        """
    
        index1 = next(i for i,v in enumerate(self.tau) if v > delaytime_range[0])
        index11 = next(i for i,v in enumerate(self.tau) if v > delaytime_range[1])
        
        index2 = next(i for i,v in enumerate(self.wl) if v > wavelength_range[0])
        index22 = next(i for i,v in enumerate(self.wl) if v > wavelength_range[1])
        
        deltaA_value_av = np.mean(self.data[index2:index22, index1:index11])
        deltaA_value_std = np.std(self.data[index2:index22, index1:index11])
        
        return deltaA_value_av, deltaA_value_std
    
    
    def load_corresponding_uvvis(self, path_uvvis, filename_uvvis):
        """
        load uvvis of sample, input path and filename
        Uv-vis and TA data come with different wavelength values
        --> interpolate uvvis to wavelength array of TA
        """        
        wl_uvvis, data_uvvis = self.import_uvvis(path_uvvis, filename_uvvis)
        
        # interpolate uvvis to TA wl array
        g = interp1d(wl_uvvis, data_uvvis, kind = "cubic")
        self.data_uvvis_int = g(self.wl)
        # self.data_uvvis_int = self.data_uvvis_int-0.015
                
        
    def import_uvvis(self, path_uvvis, filename_uvvis):
        """
        import UV-Vis data and average it if multiple measurements were taken

        """
        uvvis = np.genfromtxt(path_uvvis + '\\' + filename_uvvis, 
                         delimiter = ',', skip_header = 2)
        wl_uvvis = uvvis[:,0]
        
        uvvis = uvvis[:,1::2] # only use uneven columns
        
        data_uvvis = np.zeros(len(uvvis))
        
        for column in range(len(uvvis[0,:])): # build average of measured spectra
            data_uvvis = data_uvvis + uvvis[:,column]
        data_uvvis = data_uvvis/len(uvvis[0,:])   
        
        return wl_uvvis, data_uvvis
        
    
    
    def plot_uvvis(self, xlim=[]):
        """
        make simple plot of uvvis
        """
        
        fig, ax = plt.subplots()
        
        # ax.figure()
        ax.plot(self.wl, self.data_uvvis_int, label = f"{self.sample}")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Absorption (OD)")
        ax.set_ylim(0, np.max(self.data_uvvis_int))
        ax.legend()
        ax.grid(True)
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        fig.show()
    
        return fig, ax
    
    
    
    def calc_nonlinear_abs(self):
        """
        calculate the nonlinear absorption (pump induced absorption / absorption spectrum of
        the excited state) of the sample. 
        
        DeltaA + A0 = A_exc
        
        if A_exc < 0: stimulated emission
        
        not yet finished.. why is there sometimes an offset in absorption spectrum?
        
        """
        
        self.data_nonlinear = np.zeros((len(self.wl), len(self.tau)))
        
        for i in range(len(self.tau)):    
            self.data_nonlinear[:,i] = self.data[:, i] + self.data_uvvis_int
        
        # self.plot_map_nonlinear()
            
        
    def plot_map_nonlinear(self, ylim=[-1,200], save=False, log=False):
        """
        plot map on nonlinear absorption
        
        not yet finished
        
        """
        
        fig, ax = plt.subplots()
        plt.set_cmap("seismic")
        plt.pcolormesh(self.wl, self.tau, self.data_nonlinear.T, 
                       vmin=-0.01, vmax = 0.01,
                       rasterized=True, shading = "auto")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel(f"Delay Time ({self.timeunit})")
       
        
        # plt.plot(self.wl, self.data_uvvis_int)
        cbar = plt.colorbar()
        cbar.set_label(r"$\Delta$ A (mOD)")
        
        if ylim != []:
            plt.ylim(ylim[0], ylim[1])
    
        if log == True:
            plt.yscale("log")
    
        if save == True:
            plt.savefig(f"{today}_{self.name}_map_nonlinear.svg", bbox_inches = "tight")
        plt.show()
        
        
        
        
    def plot_spectra_nonlinear(self, tau_to_be_plotted=[1,10,100], xlim=[], ylim=[]):
        
        
        indices = []
        for tau in tau_to_be_plotted:
            index = next(i for i,v in enumerate(self.tau) if v > tau)
            indices.append(index)
            
        fig, ax = plt.subplots()
        for n in range(len(indices)):
            ax.plot(self.wl, 1000*self.data_nonlinear[:,indices[n]], 
                     label = str(round(self.tau[indices[n]])) + f" {self.timeunit}")
        ax.legend(title = "Delay Times", framealpha = 1)
        ax.grid(True)
        ax.set_ylabel("A$_{nonlinear}$ (mOD)")
        ax.set_xlabel("Wavelength (nm)")
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])   
        ax.axhline(y = 0, color = "k", linewidth = 0.5)  
        
          # plot grey area if data was deleted due to scattered pump light
        spectrum = self.get_spectra(tau_to_be_plotted[0])
        isnan = np.argwhere(np.isnan(spectrum))
        if isnan.size > 0:
            ax.axvspan(self.wl[isnan[0,1]], self.wl[isnan[-1,1]], 
                        facecolor = "0.2", alpha = 0.25)
        
        fig.show()
        return fig, ax
        
        
        
        
        
    
    def plot_dynamics(self, wl_to_be_plotted, xlim=[], ylim=[], 
                      normalize = "", log = False, save=False, invert_sign=False):
        """
        plot dynamics at certain wavelength
        provide wavelengths as list, eg. self.plot_dynamics([400,500,600])
        """ 
        
        if invert_sign == False:
            invert = 1
        elif invert_sign == True:
            invert = -1
            
            
            
        indices = []
        for wl in wl_to_be_plotted:
            index = next(i for i,v in enumerate(self.wl) if v > wl)
            indices.append(index)
   
        fig, ax = plt.subplots()
        
        if normalize == "max":
            for n in range(len(indices)):
                if np.abs(np.min(self.data[indices[n],:])) > np.max(self.data[indices[n],:]):
                    ax.plot(self.tau, invert*self.data[indices[n],:]/np.min(self.data[indices[n],:]), 
                             label = str(round(self.wl[indices[n]], 1))[:-2] + " nm")
                
                else:
                    ax.plot(self.tau, invert*self.data[indices[n],:]/np.max(self.data[indices[n],:]), 
                             label = str(round(self.wl[indices[n]], 1))[:-2] + " nm")
        
            ax.set_ylabel("$\Delta$A (OD) (normalized)")
        
        elif normalize == "last":
            for n in range(len(indices)):
                ax.plot(self.tau, invert*self.data[indices[n],:]/self.data[indices[n],-1],
                         label = str(round(self.wl[indices[n]], 1))[:-2] + " nm")
            ax.set_ylabel(r"$\Delta$A (OD) (normalized to long dynamic)")
                
        
        else:
            for n in range(len(indices)):
                ax.plot(self.tau, invert*1000*self.data[indices[n],:], 
                         label = str(round(self.wl[indices[n]], 1))[:-2] + " nm")
                ax.set_ylabel("$\Delta$A (mOD)")
                
                
        ax.legend(title = "Probe Wavelengths", framealpha = 1)
        ax.grid(True)
        ax.set_title(f"{self.sample}")
        ax.set_xlabel(f"Delay Time ({self.timeunit})")
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])
        
        if log == True:
            ax.set_xscale("log")
        
        # plot grey area if data was deleted due to artifact
        _, dynamic = self.get_one_dynamic(wl_to_be_plotted[0])
        isnan = np.argwhere(np.isnan(dynamic))
        if isnan.size > 0:
            ax.axvspan(self.tau[isnan[0]], self.tau[isnan[-1]], 
                        facecolor = "0.2", alpha = 0.25)    
            
            
        if save == True:
            fig.savefig(f"{today}_{self.sample}_dynamics.svg", bbox_inches = "tight")
        # fig.show()
    
        return fig, ax
    
    
    
    
    def plot_dynamics_av(self, spectral_range, xlim=[], ylim=[]):
        """
        plot dynamics average over a spectral range
        """
        dynamics_av = self.get_dynamics_averaged(spectral_range)
        
        fig, ax = plt.subplots()
        ax.plot(self.tau, 1000*dynamics_av)
        ax.set_title(f"averaged dynamics from {spectral_range[0]} to {spectral_range[1]} nm")
        ax.grid(True)
        ax.set_xlabel(f"Delay Time ({self.timeunit})")
        ax.set_ylabel(r"$\Delta$A (mOD)")
        # plt.xscale("log")
        # plt.yscale("log")
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])
        fig.show()
        
        return fig, ax
    
    

    def get_dynamics_averaged(self, spectral_range):
        """
        get dynamics averaged over a certain spectral range
        """
        indices = []
        for wl in spectral_range:
            index = next(i for i,v in enumerate(self.wl) if v > wl)
            indices.append(index)

        # print(indices)

        dynamics_av = np.zeros(len(self.tau))
        for n in range(indices[0], indices[1]):
            dynamics_av += self.data[n,:]
        dynamics_av = dynamics_av/(indices[1]-indices[0])

        return dynamics_av


    def get_spectrum_averaged(self, timerange):
        """
        get spectrum averaged over a certain timerange
        """
        indices = []
        for tau in timerange:
            index = next(i for i,v in enumerate(self.tau) if v > tau)
            indices.append(index)
        
        
        spectrum_av = np.zeros(len(self.wl))
        for n in range(indices[0], indices[1]):
            spectrum_av += self.data[:,n]
        spectrum_av = spectrum_av/(indices[1]-indices[0])
        
        return spectrum_av
        

        
    def plot_spectra(self, 
                     tau_to_be_plotted=[1,10,100], 
                     xlim=[], 
                     ylim=[],
                     normalize=False,
                     save=False):
        """ 
        plot spectra at certain delay times
        provide delay times as list, eg. self.plot_spectra([1,10,100])
        """
        
        indices = []
        for tau in tau_to_be_plotted:
            index = next(i for i,v in enumerate(self.tau) if v > tau)
            indices.append(index)
            
        fig, ax = plt.subplots()
        
        if normalize == False:
            for n in range(len(indices)):
                ax.plot(self.wl, 1000*self.data[:,indices[n]], 
                         label = str(round(self.tau[indices[n]])) + f" {self.timeunit}")
            ax.set_ylabel("$\Delta$A (mOD)")
            
        elif normalize == True:
            for n in range(len(indices)):
                ax.plot(self.wl, -self.data[:, indices[n]]/np.nanmin(self.data[:, indices[n]]),
                         label = str(round(self.tau[indices[n]])) + f" {self.timeunit}")
            ax.set_ylabel("$\Delta$A (OD) (normalized)")
            
        ax.legend(title = "Delay Times", framealpha = 1)
        ax.grid(True)
        ax.set_xlabel("Wavelength (nm)")
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])   
        ax.axhline(y = 0, color = "k", linewidth = 0.5)  
        
        # plot grey area if data was deleted due to scattered pump light
        spectra = self.get_spectra(tau_to_be_plotted)
        spectrum = spectra[0]
        isnan = np.argwhere(np.isnan(spectrum))
        
        if isnan.size > 0:
            ax.axvspan(self.wl[isnan[0]], self.wl[isnan[-1]], 
                        facecolor = "0.2", alpha = 0.25)
        
        
        
        if save == True:
            fig.savefig(f"{today}_{self.name}_spectra.svg", bbox_inches = "tight")
        plt.title(f"{self.sample}")
        fig.show()
        
        return fig, ax
        
            
   
    
    
    def plot_map(self, ylim=[], save=False):
        """
        plot TA map of measurement
        """
        
    
    
        fig, ax = plt.subplots()
        plt.title(f"{self.sample}")
        # plt.set_cmap("seismic")
        im = ax.pcolormesh(self.wl, self.tau, self.data.T*1000, 
                       vmin= -np.max(np.abs(self.data*1000)), 
                       vmax = np.max(np.abs(self.data*1000)),
                       rasterized=True, shading = "auto", cmap = "seismic")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel(f"Delay Time ({self.timeunit})")
        cbar = plt.colorbar(im)
        cbar.set_label(r"$\Delta$A (mOD)")
        
        # ax.set_yscale("log")
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])
    
        if save == True:
            fig.savefig(f"{today}_{self.name}_map.svg", bbox_inches = "tight")
        fig.show()
        
        return fig, ax
    
    
    
    def plot_map_linlog(self, linlog=2, xlim=[], save=False):
        """
        plot TA map of measurement with broken axis. first few ps linear, rest logarithmic
        linlog: delay time of transition between linear and logarithmic scale
        """
        
        # initialize figure
        import matplotlib as mpl
        fig, (ax, ax2) = plt.subplots(2,1, sharex=True, figsize = (4,4), 
                                      gridspec_kw={'height_ratios': [3, 2]})
        plt.set_cmap("seismic")
        
        # plot logarithmic part
        ax.pcolormesh(self.wl, self.tau, self.data.T, 
                       vmin= -np.max(np.abs(self.data)), 
                       vmax = np.max(np.abs(self.data)),
                       rasterized=True, shading = "auto")
        ax.set_ylim(linlog, np.max(self.tau))
        ax.set_yscale("log")
     
        
        # plot linear part
        ax2.pcolormesh(self.wl, self.tau, self.data.T, 
                       vmin= -np.max(np.abs(self.data)), 
                       vmax = np.max(np.abs(self.data)),
                       rasterized=True, shading = "auto")
        ax2.set_ylim(-0.5, linlog)
        
        
        # set axis
        plt.xlabel("Wavelength (nm)")
        plt.ylabel(f"Delay Time ({self.timeunit})")
        if xlim != []:
            plt.xlim(xlim[0],xlim[1])
        
        
        # colorbar
        cmap = mpl.cm.seismic
        norm = mpl.colors.Normalize(vmin= -np.max(np.abs(self.data*10**3)), 
                       vmax = np.max(np.abs(self.data*10**3)))
        cax = plt.axes([0.95, 0.15, 0.04, 0.7])
        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
        cbar.set_label(r"$\Delta$A (mOD)")
        
        # set distance between subplots to zero
        plt.subplots_adjust(hspace=0)
    
        if save == True:
            plt.savefig(f"{today}_{self.sample}_map.svg", bbox_inches = "tight")
        plt.show()
    
        return fig, ax
    
    
    
    def plot_map_dyn_spec(self, wl_to_be_plotted, tau_to_be_plotted):
        """
        plot TA map of measurement with dynamics on the left and spectra down below
        not yet finished, just for the one example suitable which i wrote it for
        """
        
    
        fig, ax = plt.subplots(2, 2, 
                       gridspec_kw={
                           'width_ratios': [3, 1.5],
                           'height_ratios': [3, 2]})
        
        
        # plot map
        plt.set_cmap("seismic")
        # ax[0][0].set_aspect("equal")
        ax[0][0].pcolormesh(self.wl, self.tau, self.data.T, 
                       vmin= -np.max(np.abs(self.data)), 
                       vmax = np.max(np.abs(self.data)), 
                       rasterized=True, shading = "auto")
        ax[0][0].set_xlabel("Wavelength (nm)")
        ax[0][0].set_ylabel(f"Delay Time ({self.timeunit})")
        ax[0][0].set_ylim(-0.25,1)
        ax[0][0].set_xlim(450,600)
        # cbar = ax[0][0].colorbar()
        # cbar.set_label(r"$\Delta$ A (mOD)")
        
        
        #plot dynamics
        indices = []
        for wl in wl_to_be_plotted:
            index = next(i for i,v in enumerate(self.wl) if v > wl)
            indices.append(index)
            
        for n in range(len(indices)):
            ax[0][1].plot(self.data[indices[n],:], self.tau,  label = str(round(self.wl[indices[n]], 0))[:-2] + " nm")
            
        ax[0][1].set_xlim(0.005, -0.085)
        ax[0][1].set_ylim(-3, 300)
        ax[0][1].set_xlabel(r"$\Delta$A")
        ax[0][1].legend()
        
        
        
        
        # plot spectra
        indices = []
        for tau in tau_to_be_plotted:
            index = next(i for i,v in enumerate(self.tau) if v > tau)
            indices.append(index)
            
   
        for n in range(len(indices)):
            ax[1][0].plot(self.wl, self.data[:,indices[n]], 
                     label = str(round(self.tau[indices[n]])) + f" {self.timeunit}")
        
        ax[1][0].axhline(y = 0, color = "k", linewidth = 0.3)
        ax[1][0].set_xlim(420,700)
        ax[1][0].set_ylabel(r"$\Delta$A")
        ax[1][0].set_xlabel(r"Wavelength (nm)")
        ax[1][0].legend()
        
        
        wl_to_be_plotted2 = [505]
        
        indices = []
        for wl in wl_to_be_plotted2:
            index = next(i for i,v in enumerate(self.wl) if v > wl)
            indices.append(index)
            
        for n in range(len(indices)):
            ax[1][1].plot( self.tau, self.data[indices[n],:],  label = str(round(self.wl[indices[n]], 0))[:-2] + " nm")
   
        
        
        ax[1][1].axhline(y = 0, color = "k", linewidth = 0.3)
        ax[1][1].set_xlim(-0.25,5)
        ax[1][1].set_ylim(-0.02,0.01)
        ax[1][1].set_ylabel(r"$\Delta$A")
        ax[1][1].set_xlabel(r"Delay Time (ps)")
        ax[1][1].legend()
       
    

        plt.show()
    
    
    
    def delete_data_scattered_pump(self, limits_in_nm):
        """
        delete data in the region of the pump pulse, if scattered pump pulse distorts data
        
        input: give to two wavelengths limiting the region you want to delete as list
        """
        indices = []
        for wl in limits_in_nm:
            index = next(i for i,v in enumerate(self.wl) if v > wl)
            indices.append(index)
            
        self.data[indices[0]:indices[1], :] = np.nan    
    
    
    def delete_data_artifact(self, limits_in_ps=[0,1]):
        """
        delete data in the temporal domain of the artifact, if it distorts the data
        
        input: give time points limiting the artifact
        default is from 0 to 1 ps
        """
        indices = []
        for delaytime in limits_in_ps:
            index = next(i for i,v in enumerate(self.tau) if v > delaytime)
            indices.append(index)
        
        self.data[:, indices[0]:indices[1]] = np.nan
        
        
        
    def correct_zeropoint(self, zero_point):
        """
        correct zero point, input zero-point in appropriate time unit
        """
    
        self.tau = self.tau - zero_point
    
  
    def smooth_dynamics(self, degree = 3):
        """
        smooth dynamics
        degree: 3 not so much
                21 very much
        """
        from scipy.signal import savgol_filter
        
        smoothed = np.zeros(self.data.shape)
        for idx, row in enumerate(self.data):
            dat = savgol_filter(row[:], degree, 1)
            smoothed[idx,:] = dat
        
        self.data = smoothed
        
      
  
    
    def analyse_bleach_shift(self, limits, wl_to_compare=[], plot_fit=False):
        """
        follow shift of bleach by fitting gaussian functions to data 
        e.g. for TA measurements of AuNP (bleach maximum shifts blue with delay time)
        
        inputs: 
            limits: give spectral area of bleach as list
         
        
        """
        indices = []
        for wl in limits:
            index = next(i for i,v in enumerate(self.wl) if v > wl)
            indices.append(index)
            
        data_cropped = - self.data[indices[0]:indices[1], :]    
        wl_cropped = self.wl[indices[0]:indices[1]]
    
    
        self.wl_bleach = []
        self.deltaA_bleach = []
    
        self.wl_bleach_fit = [] 
        self.deltaA_bleach_fit = [] 
        
          
        # def gaussian(x, amp, cen, wid):
        #     return amp * np.exp(-((x-cen)/wid)**2)
        
        # def gaussian2(x, amp1, cen1, wid1, amp2, cen2, wid2):
        #     return amp1 * np.exp(-((x-cen1)/wid1)**2) + amp2 * np.exp(-((x-cen2)/wid2)**2)
        
        # def gaussian3(x, amp1, cen1, wid1, amp2, cen2, wid2, amp3, cen3, wid3):
        #     return amp1 * np.exp(-((x-cen1)/wid1)**2) + amp2 * np.exp(-((x-cen2)/wid2)**2) + amp3 * np.exp(-((x-cen3)/wid3)**2)
        
        def poly(x, A0, A1, A2, A3, A4, A5):
            return A0 + A1*x + A2*x**2 + A3*x**3 + A4*x**4 + A5*x**5
        
        for i, spectrum in enumerate(data_cropped.T):
        
            try:   
                if np.max(spectrum)*1000000 > 1: 
                        
                    # guess = [np.max(spectrum), wl_cropped[np.argmax(spectrum)], 15,np.max(spectrum), wl_cropped[np.argmax(spectrum)], 15]
                    popt, pcov = curve_fit(poly, wl_cropped, spectrum)
                    
                    fit = poly(wl_cropped, *popt)
                    
                    self.deltaA_bleach.append(np.max(spectrum))
                    self.wl_bleach.append(wl_cropped[np.argmax(spectrum)])
                    
                    
                    self.deltaA_bleach_fit.append(np.max(fit))
                    self.wl_bleach_fit.append(wl_cropped[np.argmax(fit)])
                    
                    if plot_fit == True:
                        if i == 90 or i == 190:
                            plt.figure()
                            plt.title(f"{self.tau[i]} ps")
                            plt.plot(wl_cropped, spectrum)
                            plt.plot(wl_cropped, fit)
                            plt.show()
                else:
                    self.deltaA_bleach_fit.append(np.nan)
                    self.wl_bleach_fit.append(np.nan)
                    self.deltaA_bleach.append(np.max(spectrum))
                    self.wl_bleach.append(wl_cropped[np.argmax(spectrum)])
                    
                 
            except:
                 
                self.deltaA_bleach_fit.append(np.nan)
                self.wl_bleach_fit.append(np.nan)
                self.deltaA_bleach.append(np.max(spectrum))
                self.wl_bleach.append(wl_cropped[np.argmax(spectrum)])
            
        plt.figure()
        plt.plot(self.tau, self.wl_bleach, label = "max")
        plt.plot(self.tau, self.wl_bleach_fit, label = "fit")
        plt.legend()
        plt.ylabel("wavelength of bleach max")
        plt.xlabel("Delay Time (ps)")
        # plt.xlim(-1,10)
        plt.grid(True)
        plt.show()
        
        plt.figure()
        plt.plot(self.tau, self.deltaA_bleach, label = "max")
        plt.plot(self.tau, self.deltaA_bleach_fit, label = "fit")
        plt.legend()
        plt.ylabel(r"$\Delta$ A of bleach max")
        plt.xlabel("Delay Time (ps)")
        # plt.xlim(-10,1000)
        plt.grid(True)
        plt.show()
        
        
        
        # compare dynamics of bleach maximum with dynamics at one wavelength
        if wl_to_compare != []:
            plt.figure()
            plt.plot(self.tau, self.deltaA_bleach, label = "max")
            plt.plot(self.tau, self.deltaA_bleach_fit, label = "fit")
            plt.plot(self.tau, -self.get_one_dynamic(wl_to_compare)[1], label = f"{wl_to_compare} nm")
            plt.legend()
            plt.ylabel(r"$\Delta$ A of bleach max")
            plt.xlabel("Delay Time (ps)")
            # plt.xlim(-10,1000)
            plt.grid(True)
            plt.show()
    

    
    def analyse_max_bleach(self, limits):
        """
        follow shift of bleach by just checking maximum bleach value
        e.g. for TA measurements of AuNP (bleach maximum shifts blue with delay time)
        
        inputs: 
            limits: give spectral area of bleach as list
         
        
        """
        indices = []
        for wl in limits:
            index = next(i for i,v in enumerate(self.wl) if v > wl)
            indices.append(index)
            
        data_cropped = - self.data[indices[0]:indices[1], :]    
        wl_cropped = self.wl[indices[0]:indices[1]]
    
    
        self.wl_bleach = []
        self.deltaA_bleach = []
    
        
        
        for i, spectrum in enumerate(data_cropped.T):
            
                    self.deltaA_bleach.append(np.max(spectrum))
                    self.wl_bleach.append(wl_cropped[np.argmax(spectrum)])
            
        
        
        
        
            
        # plt.figure()
        # plt.plot(self.tau, self.wl_bleach, label = "max")
        # plt.legend()
        # plt.ylabel("wavelength of bleach max")
        # plt.xlabel("Delay Time (ps)")
        # plt.xlim(-1,10)
        # plt.grid(True)
        # plt.show()
        
        plt.figure()
        plt.title(self.pumpfluence)
        plt.plot(self.tau, self.deltaA_bleach, label = "max")
        plt.legend()
        plt.ylabel(r"$\Delta$ A of bleach max")
        plt.xlabel("Delay Time (ps)")
        # plt.xlim(-10,1000)
        plt.grid(True)
        plt.show()
        
        
        
        
        
        
        
        
        
        
        
    def fit_exp(self, wl_to_be_fitted=[], min_delay_time=[], max_delay_time = None, 
                log = False, av=False, spectral_range=[], xlim=[], ylim=[], bounds=[]):
        """
        fit exponential decay to dynamics at one probe wavelength
        or to dynamics averaged over a certain spectral range
        give delay time range you want to fit 
        
        set wl_to_be_fitted to "bleach_max" to fit the max of bleach as evaluated by self.analyse_bleach_shift
        
        """
        
        
        if av == True:
        
            indices = []
            for wl in spectral_range:
                index = next(i for i,v in enumerate(self.wl) if v > wl)
                indices.append(index)
    
            print(indices)
    
            dynamics_av = np.zeros(len(self.tau))
            for n in range(indices[0], indices[1]):
                dynamics_av += self.data[n,:]
            dynamic = dynamics_av/(indices[1]-indices[0])
        
        elif wl_to_be_fitted == "bleach_max":
            dynamic = self.deltaA_bleach_fit
        
        else:
            dynamic = self.get_one_dynamic(wl_to_be_fitted)[1]
        
        
        
        if np.abs(np.min(dynamic)) > np.max(dynamic):
            dynamic = -dynamic
        
        
        if max_delay_time == None:
            max_delay_time = self.tau[-2]
        
        if min_delay_time != None and max_delay_time != None:
            index1 = next(i for i,v in enumerate(self.tau) if v > min_delay_time)
            index2 = next(i for i,v in enumerate(self.tau) if v > max_delay_time)
            tau = self.tau[index1:index2]
            dynamic = dynamic[index1:index2]
            if np.abs(np.min(dynamic)) > np.max(dynamic):
                dynamic = -dynamic
        
            
        else: 
            tau = self.tau[np.nanargmax(dynamic):]
            dynamic = dynamic[np.nanargmax(dynamic):]   
        
        

        def exp_decay(t, A, tau1):
            return A * np.exp(-t/tau1)
            
        
        boundz  = ((0,0), (np.inf, np.inf))
        
        popt, pcov = curve_fit(exp_decay, tau, dynamic, bounds = boundz)    
        
        tau1 = popt[1]
        perr = np.sqrt(np.diag(pcov))
        tau1sd = perr[1] 
        
        tau1_str = f"({round(tau1,3)} $\pm$ {round(tau1sd, 3)}) {self.timeunit}"
        
        fit = exp_decay(tau, *popt)
        
        
        # plot figure
        plt.figure()
        plt.plot(tau, dynamic, color = "tab:blue", alpha = 0.5, label = "data")
        plt.plot(tau, fit, color = "tab:orange", label = "fit")
        plt.xlabel(f"Delay Time ({self.timeunit})")
        plt.ylabel("$\Delta$A (mOD)")
        plt.legend(title = f"{self.sample}, " 
                    + "$\lambda_{probe}$: "  
                    + f"{wl_to_be_fitted} nm \n" + r"$\tau$ = " + tau1_str, 
                    framealpha = 1)
        plt.grid()
        if log == True:
            plt.xscale("log")
        if xlim != []:
            plt.xlim(xlim[0], xlim[1])
        if ylim != []:
            plt.ylim(ylim[0], ylim[1])   
        plt.show()
        
        return tau1, tau1sd
        
        
    def fit_2exp(self, wl_to_be_fitted=[], min_delay_time=None, max_delay_time=None,
                 log=False, av=False, spectral_range=[], xlim=[], ylim=[]):
        
        """
        fit exponential decay to dynamics at one probe wavelength
        or to dynamics averaged over a certain spectral range
        give delay time range you want to fit using min_delay_time and max_delay_time
        
        """
        
        
        if av == True:
        
            indices = []
            for wl in spectral_range:
                index = next(i for i,v in enumerate(self.wl) if v > wl)
                indices.append(index)
    
            print(indices)
    
            dynamics_av = np.zeros(len(self.tau))
            for n in range(indices[0], indices[1]):
                dynamics_av += self.data[n,:]
            dynamic = dynamics_av/(indices[1]-indices[0])
        
        elif wl_to_be_fitted == "bleach_max":
            dynamic = self.deltaA_bleach_fit
        
        else:
            dynamic = self.get_one_dynamic(wl_to_be_fitted)[1]
        
        if np.abs(np.min(dynamic)) > np.max(dynamic):
            dynamic = -dynamic
        
        
        
        
        if min_delay_time != None:
            index = next(i for i,v in enumerate(self.tau) if v > min_delay_time)
            tau = self.tau[index:]
            dynamic = dynamic[index:]
            if np.abs(np.min(dynamic)) > np.max(dynamic):
                dynamic = -dynamic
        
        
        elif max_delay_time != None and min_delay_time == None:
            index = next(i for i,v in enumerate(self.tau) if v > max_delay_time)
            tau = self.tau[np.nanargmax(dynamic):index]
            dynamic = dynamic[np.nanargmax(dynamic):index]
            
        else: 
            tau = self.tau[np.nanargmax(dynamic):]
            dynamic = dynamic[np.nanargmax(dynamic):]   
        




        def exp_decay2(t, a1, tau1, a2, tau2):
            return a1 * np.exp(-t/tau1) + a2 * np.exp(-t/tau2)
            
        
        boundz  = ((0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf))
        
        popt, pcov = curve_fit(exp_decay2, tau, dynamic, bounds = boundz) 
        
        a1 = popt[0]
        tau1 = popt[1]
        a2 = popt[2]
        tau2 = popt[3]
        
        perr = np.sqrt(np.diag(pcov))
        
        tau1_sd = perr[1]
        tau2_sd = perr[3]
        
        if tau2 < tau1:
            tau1, tau2 = tau2, tau1
            a1, a2 = a2, a1
            tau1_sd, tau2_sd = tau2_sd, tau1_sd
        
        tau1_str = f"{round(tau1,1)} $\pm$ {round(tau1_sd, 1)} {self.timeunit}"
        tau2_str = f"{round(tau2)} $\pm$ {round(tau2_sd)} {self.timeunit}"
        
        fit = exp_decay2(tau, *popt)
        
        
        # plot figure
        plt.figure()
        plt.plot(tau, dynamic, color = "tab:blue", alpha = 0.5, label = "data")
        plt.plot(tau, fit, color = "tab:orange", 
                 label = "fit")
        plt.xlabel(f"Delay Time ({self.timeunit})")
        plt.ylabel("$\Delta$A (mOD)")
        plt.legend(title =  f"{self.sample}, " 
                   + "$\lambda_{probe}$: "  
                   + f"{wl_to_be_fitted} nm \n" + r"$\tau_{1}$ = " + tau1_str + "\n" + 
                 r"$\tau_{2} = $" + tau2_str,
                   framealpha = 1)
        plt.grid()
        if log == True:
            plt.xscale("log")
            # plt.yscale("log")
        if xlim != []:
            plt.xlim(xlim[0], xlim[1])
        if ylim != []:
            plt.ylim(ylim[0], ylim[1])       
        plt.show()
            
        print(tau1, tau2)
        
        # def bi_exp():
        return tau1, tau1_sd    
        
    def average_spectra(self, timerange, xlim=[], ylim=[]):
        """
        average spectra in a certain timerange
        """    
        
        spectrum_av = self.get_spectrum_averaged(timerange)
        
        fig, ax  = plt.subplots()
        ax.set_title(f"time-averaged TA spectrum of {self.measurement} \n averaged from {timerange[0]} to {timerange[1]} {self.timeunit}")
        ax.plot(self.wl, 1000*spectrum_av)
        ax.grid(True)
        ax.set_ylabel("$\Delta$A (mOD)")
        ax.set_xlabel("Wavelength (nm)")
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])   
        ax.axhline(y = 0, color = "k", linewidth = 0.5)  
        fig.show()
        return fig, ax    
        
    def nm_to_ev(self):
        """
        Convert nm scale to eV scale using Jacobian transformation
        
        see dx.doi.org/10.1021/jz401508t
        """
        
        
    
        self.e = (sc.h * sc.c )/(sc.e * self.wl*10**-9)
        self.data_e = np.zeros(self.data.shape)
        for i in range(len(self.tau)):
            spec_e = self.data[:,i] * ((sc.h * sc.c)/self.e**2)
            scale = np.min(spec_e) / np.min(self.data[:,i])
            self.data_e[:,i] = spec_e/scale
        
        
        
        
        
    def calc_carrier_temp(self, xlim=[], tau_range=[0.1,10]):
        """
        calculate carrier temperature from high-energy tail of TA spectra
        by fitting Boltzmann distribution
        
        input:
            xlim in nm
            delay time range to be evaluated
            
        to do:
            clean up code
            take into account excitonic features (how?)
            xlim depending on spectral shape at each delay time
            
        """

        # convert nm to eV
        self.nm_to_ev()
        
        
        # normalize spectra
        data_normalized = np.zeros(self.data_e.shape)
        for i in range(len(self.tau)):
            dat = self.data_e[:,i] / np.abs(np.min(self.data_e[:,i]))
            data_normalized[:,i] = -dat
        
        
        
        
        # crop to xlim
        xlim = np.array(xlim)
        xlim_e = (sc.h * sc.c )/(sc.e * xlim*10**-9) 
       
        
        indices_x = []
        for x in xlim_e:
            index = next(i for i,v in enumerate(self.e[::-1]) if v > x)
            index = len(self.e)-index
            indices_x.append(index)
        
        
        
        e_fit = self.e[indices_x[0]: indices_x[1]]
        data_fit = data_normalized[indices_x[0]:indices_x[1],:]
        
        
        
        
        # plt.figure()
        # plt.plot(self.e, data_normalized[:,600])
        # plt.plot(e_fit, data_fit[:,600], 'o')
        # plt.show()
        
        
        
        
        
        # crop to tau_range
        indices_tau = []
        for tau in tau_range:
            index = next(i for i,v in enumerate(self.tau) if v > tau)
            indices_tau.append(index)
            
        print(indices_tau)
        
        tau_fit = self.tau[indices_tau[0]:indices_tau[1]]
        data_fit = data_fit[:, indices_tau[0]:indices_tau[1]]        
        data_normalized = data_normalized[:, indices_tau[0]:indices_tau[1]]
        
        
        
        
        k = 8.617E-5 # Boltzmann constant in eV/K
        
        def boltzmann(e, A_0, T_e, c):
            return A_0 * np.exp(-(e-c)/(k*T_e)) 
        
        
        
        
        T_e = []
        fits = []
        for i in range(len(tau_fit)):
            dat = data_fit[:, i]
            dat_normalized = data_normalized[:, i]
            
            if i == 1:
                plt.figure()
                plt.plot(self.e, dat_normalized)
                plt.plot(e_fit, dat, 'o')
                plt.show()
            
            boundz = ((-1, 0, -5), (1, 1E6, 5))
            guess = (np.max(dat), 1000, 2)
            
            
            popt, pcov = curve_fit(boltzmann, e_fit, dat, bounds = boundz, p0 = guess)
            fit = boltzmann(e_fit, *popt)
            fits.append(fit)
            T_e.append(popt[1])
            
            if i == 0 or i == 10 or i == 100:
                plt.figure()
                plt.plot(e_fit, dat, label = "data")
                plt.plot(e_fit, fit, label = "fit")
                plt.show()
          
        fits = np.array(fits)  
        print(fits.shape)  
        
        tau_to_be_plotted = [0.1, 0.2, 0.5, 0.8, 1, 1.5, 2, 3,200,3000]
        indices = []
        for tau in tau_to_be_plotted:
            index = next(i for i,v in enumerate(tau_fit) if v > tau)
            indices.append(index)
            
        plt.figure()
        for index in indices:
            plt.plot(self.e, data_normalized[:,index], color = "tab:blue")
            plt.plot(e_fit, fits[index,:], color = "tab:orange")
        plt.show()    
            
            
        
        plt.figure()
        plt.plot(tau_fit, T_e, 'o')
        plt.ylabel("Carrier Temperature (K)")
        plt.xlabel("Delay Time (ps)")
        plt.xscale("log")
        plt.grid(True)
        plt.show()
        
        print(T_e[-1], "K: final temperature")
        
        
        
        # fit Boltzmann-distribution from Cao et al, Nanoscale, 2017, 9, 17133–17142
        # plot time dependent carrier temperature
        
        return
        
        
        
    def calc_nonthermal_TA(self, tau=3.2, wavelength=490):
        """
        calculate TA of thermal and nonthermal carriers for AuNPs
        cf. O'Keeffe, 2024, ACS Photonics, https://doi.org/10.1021/acsphotonics.4c00591
        """
        
        spec_thermal = self.get_one_spectrum(tau)
        dyn_thermal = self.get_one_dynamic(wavelength)[1]
        
        # map_thermal = np.zeros((len(spec_thermal), len(dyn_thermal)))
        
        # print(np.array(map_thermal))
        
        map_thermal = np.outer(spec_thermal, dyn_thermal)
        map_thermal = (map_thermal/np.max(np.abs(map_thermal)))*np.max(np.abs(self.data))
        
        print(map_thermal)
        print(np.shape(map_thermal))
        
        
        
        fig, ax = plt.subplots()
        plt.title("thermal electrons TA map")
        im = ax.pcolormesh(self.wl, self.tau, map_thermal.T*1000, 
                       vmin= -np.max(np.abs(map_thermal*1000)), 
                       vmax = np.max(np.abs(map_thermal*1000)),
                       rasterized=True, shading = "auto", cmap = "seismic")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel(f"Delay Time ({self.timeunit})")
        cbar = plt.colorbar(im)
        cbar.set_label(r"$\Delta$A (mOD)")
        ax.set_ylim(-0.1, 1)
        fig.show()
        
        
        map_nonthermal = self.data - map_thermal
        
        fig, ax = plt.subplots()
        plt.title("nonthermal electrons TA map")
        im = ax.pcolormesh(self.wl, self.tau, map_nonthermal.T*1000, 
                       vmin= -np.max(np.abs(map_nonthermal*1000)), 
                       vmax = np.max(np.abs(map_nonthermal*1000)),
                       rasterized=True, shading = "auto", cmap = "seismic")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel(f"Delay Time ({self.timeunit})")
        cbar = plt.colorbar(im)
        cbar.set_label(r"$\Delta$A (mOD)")
        ax.set_ylim(-0.1, 2)
        fig.show()
        
        
        
        return map_thermal
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
class TAMeasurementRow:
    
    def __init__(self, name="None", discrimination=""):
        """ 
        discrimination: how do the measurements in this row differ
        use either "sample", "fluence" or "wavelength" (wavelength is not yet integrated)
        """
        self.name = name       
        self.measurements = []
        self.discrimination = discrimination

        
        
    def add_measurements(self, measurements):
        """
        add previously imported measurement to measurement row
        import multiple measurements at a time
        """
        for measurement in measurements:
            self.measurements.append(measurement)
    
    def define_labels(self, measurement):
        """
        used to define labels for plot depending on discrimination line of measurement row
        """
        
        
        if self.discrimination == "sample":
            labeling = f"{measurement.sample}"                      
        elif self.discrimination == "fluence":
            labeling = f"{measurement.pumpfluence}" + " $\mathrm{\mu}$J/cm$^{2}$"
            #labeling = f"{measurement.pumpfluence}" + " mJ/cm$^{2}$"
        elif self.discrimination == "pump wavelength":
            labeling = f"{measurement.pumpwavelength}" + " nm"
        else:
            labeling = f"{measurement.sample}, {measurement.pumpfluence}" + " $\mathrm{\mu}$J/cm$^{2}$"
        return labeling
                                



    def get_list_of_pumpfluences(self):
        """
        returns a list of all the pumpfluences used in the measurement row
        """
        self.pumpfluence_list = [float(measurement.get_pumpfluence()) for measurement in self.measurements] 
        return self.pumpfluence_list
        
    def load_row(self):
        """
        to be done
        load a hole measurement row at once (instead of loading every measurement
        on its own and then adding it to the measurement row)
        by creating an instance of TA_data inside this method
        what is the input?
        """        
        
        
    def convert_timeunit_to_ns(self):
        for measurement in self.measurements:
            measurement.convert_timeunit_to_ns()
        
        
    def compare_dynamics(self, wl_to_be_plotted, normalize="none", save=False, invert_sign=False, 
                         xlim=[], ylim=[], show_probe_wavelength=True, figsize=[]):
        """
        compare dynamics at a certain probe wavelength for all measurements of 
        measurement row
        
        input: normalize
            no entry: not normalized
            normalize = "max": normalized to max or min of dynamics
            normalize = "last": normalized on last datapoint of dynamics
        
        input: xlim / ylim
            list with two entrys
        
        to do:
             it only evaluates if data was deleted on last measurement of row
             check if there is deleted data in any of the plotted measurements.. 
             if so, plot grey area
            
        """
        
        if type(wl_to_be_plotted) == int:
            wl_to_be_plotted = [wl_to_be_plotted]
        
        if type(wl_to_be_plotted) == float:
            wl_to_be_plotted = [wl_to_be_plotted]
        
        if invert_sign == False:
            invert_sign_value = 1
        elif invert_sign == True:    
            invert_sign_value = -1
        
        if figsize == []:
            fig, ax = plt.subplots()
        else:
            fig, ax = plt.subplots(figsize=figsize)
            
        if normalize == "max":
            for measurement in self.measurements:
                labeling = self.define_labels(measurement)
                index, dynamics = measurement.get_one_dynamic(wl_to_be_plotted)
                
                        
                if np.abs(np.nanmin(dynamics)) > np.nanmax(dynamics):
                    ax.plot(measurement.get_tau(), invert_sign_value*dynamics/np.nanmin(dynamics),
                                 label = labeling)    
                        
                else:
                    ax.plot(measurement.get_tau(), invert_sign_value*dynamics/np.nanmax(dynamics),
                             label = labeling)
            ax.set_ylabel(r"$\Delta$A (OD) (normalized)")
                    
                    
        elif normalize == "last":
             for measurement in self.measurements:
                labeling = self.define_labels(measurement)
                index, dynamics = measurement.get_one_dynamic(wl_to_be_plotted)
                ax.plot(measurement.get_tau(), invert_sign_value*dynamics/dynamics[-1],
                         label = labeling)
             ax.set_ylabel(r"$\Delta$A (OD) (normalized to long dynamic)")
            
        else:    
            for measurement in self.measurements:
                labeling = self.define_labels(measurement)
                index, dynamics = measurement.get_one_dynamic(wl_to_be_plotted)
                ax.plot(measurement.get_tau(), invert_sign_value*dynamics*1000,
                         label = labeling)
            ax.set_ylabel(r"$\Delta$A (mOD)")
        
        
        ax.grid(True)
        
        # plot grey area if data was deleted due to artifact in any one of the measurements
        isnan = []
        isnan_length = np.zeros(len(self.measurements))
        for i, measurement in enumerate(self.measurements):
            _, dynamic = measurement.get_one_dynamic(wl_to_be_plotted)
            isnan_ = np.argwhere(np.isnan(dynamic))
            isnan.append(isnan_) # create list of arrays of all isnans
            isnan_length[i] = len(isnan_) # create list containing length of arrays in isnan list of arrays
        longest_isnan_idx = np.argmax(isnan_length) # which array in isnan list is the longest (the measurement with the largest delete data-area)
        if isnan[longest_isnan_idx].size > 0:
            ax.axvspan(self.measurements[-1].tau[isnan[longest_isnan_idx][0]], 
                        self.measurements[-1].tau[isnan[longest_isnan_idx][-1]], 
                        facecolor = "0.2", alpha = 0.25)
            
        
        
        
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])
        ax.axhline(y = 0, color = "k", linewidth = 0.5)
        if show_probe_wavelength == True:
            ax.legend(title = f"Probe Wavelength: {wl_to_be_plotted[0]} nm")#, bbox_to_anchor=(1.1, 1.05))
        if show_probe_wavelength == False:
            ax.legend()
        ax.set_xlabel(f"Delay Time ({self.measurements[0].timeunit})")
        
        if self.name != "None":
            ax.set_title(f"{self.name}")
        
        if save == True:
            fig.savefig(f"{today}_compare_dynamics_{self.name}_{wl_to_be_plotted}nm.svg", bbox_inches = "tight")
        
        plt.show()
        
        return fig, ax
        
    def compare_dynamics_2(self, wl_to_be_plotted, normalize="none", save=False, invert_sign=False, 
                         xlim=[], ylim=[], av = False):
        """
        
        this needs to be merged with compare_dynamics
        
        
        
        
        has the functionality to average dynamics over a certain spectral range
        greying out deleted area (due to artifact) is currently commented out
        
        
        
        
        
        compare dynamics at a certain probe wavelength for all measurements of 
        measurement row
        
        input: normalize
            no entry: not normalized
            normalize = "max": normalized to max or min of dynamics
            normalize = "last": normalized on last datapoint of dynamics
        
        input: xlim / ylim
            list with two entrys
        
        input: av
            set to True if you want to averaged data over a certain spectral range
            give spectral range as list instead of wl_to_be_plotted
        
        
        
        to do:
             it only evaluates if data was deleted on last measurement of row
             check if there is deleted data in any of the plotted measurements.. 
             if so, plot grey area
            
        """
        
        if invert_sign == False:
            invert_sign_value = 1
        elif invert_sign == True:    
            invert_sign_value = -1
        
        fig, ax = plt.subplots()
        
        for measurement in self.measurements:
            labeling = self.define_labels(measurement)
            if av == False:
                index, dynamics = measurement.get_one_dynamic(wl_to_be_plotted)
            elif av == True:
                dynamics = measurement.get_dynamics_averaged(wl_to_be_plotted)
                
            if normalize == "max":
                
                    if np.abs(np.nanmin(dynamics)) > np.nanmax(dynamics):
                        ax.plot(measurement.get_tau(), invert_sign_value*dynamics/np.nanmin(dynamics),
                                 label = labeling)
                    else:
                        ax.plot(measurement.get_tau(), invert_sign_value*dynamics/np.nanmax(dynamics),
                                 label = labeling)
                    ax.set_ylabel(r"$\Delta$A (OD) (normalized)")
                        
                        
            elif normalize == "last":
                 
                    ax.plot(measurement.get_tau(), invert_sign_value*dynamics/dynamics[-1],
                             label = labeling)
                    ax.set_ylabel(r"$\Delta$A (OD) (normalized to long dynamic)")
                
            else:    
                
                    ax.plot(measurement.get_tau(), invert_sign_value*dynamics*1000,
                             label = labeling)
                    ax.set_ylabel(r"$\Delta$A (mOD)")
        
        
        ax.grid(True)
        
        # # plot grey area if data was deleted due to artifact in any one of the measurements
        # isnan = []
        # isnan_length = np.zeros(len(self.measurements))
        # for i, measurement in enumerate(self.measurements):
        #     _, dynamic = measurement.get_one_dynamic(wl_to_be_plotted)
        #     isnan_ = np.argwhere(np.isnan(dynamic))
        #     isnan.append(isnan_) # create list of arrays of all isnans
        #     isnan_length[i] = len(isnan_) # create list containing length of arrays in isnan list of arrays
        # longest_isnan_idx = np.argmax(isnan_length) # which array in isnan list is the longest (the measurement with the largest delete data-area)
        # if isnan[longest_isnan_idx].size > 0:
        #     plt.axvspan(self.measurements[-1].tau[isnan[longest_isnan_idx][0]], 
        #                 self.measurements[-1].tau[isnan[longest_isnan_idx][-1]], 
        #                 facecolor = "0.2", alpha = 0.25)
            
        
        
        
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])
        ax.axhline(y = 0, color = "k", linewidth = 0.5)
        if av == False:
            ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title = f"Probe Wavelength: {wl_to_be_plotted} nm")
        elif av == True:
            ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title = f"Probe Range: {wl_to_be_plotted[0]}-{wl_to_be_plotted[1]} nm")
        ax.set_xlabel(f"Delay Time ({self.measurements[0].timeunit})")
        

        plt.show()
        
        
        return fig, ax

    
    def compare_dynamics_at_different_probe_wavelengths(self, wl_list, normalize="none", save=False, invert_sign=False, 
                         xlim=[], ylim=[], probed_feature=""):
        """
        
        compare dynamics of different measurements at one specific probe wavelength per measurement
        in case you want to compare two samples where a feature is slightly shifted
        
        greying out deleted area (due to artifact) is currently commented out
        
        
        
        input: normalize
            no entry: not normalized
            normalize = "max": normalized to max or min of dynamics
            normalize = "last": normalized on last datapoint of dynamics
        
        input: xlim / ylim
            list with two entrys
        
        input: av
            set to True if you want to averaged data over a certain spectral range
            give spectral range as list instead of wl_to_be_plotted
        
        
        
        to do:
             it only evaluates if data was deleted on last measurement of row
             check if there is deleted data in any of the plotted measurements.. 
             if so, plot grey area
            
        """
        
        if invert_sign == False:
            invert_sign_value = 1
        elif invert_sign == True:    
            invert_sign_value = -1
        
        fig, ax = plt.subplots()
        
        for idx, measurement in enumerate(self.measurements):
            labeling = self.define_labels(measurement)
            
            index, dynamics = measurement.get_one_dynamic(wl_list[idx])
            
            if normalize == "max":
                
                    if np.abs(np.nanmin(dynamics)) > np.nanmax(dynamics):
                        ax.plot(measurement.get_tau(), invert_sign_value*dynamics/np.nanmin(dynamics),
                                 label = labeling)
                    else:
                        ax.plot(measurement.get_tau(), invert_sign_value*dynamics/np.nanmax(dynamics),
                                 label = labeling)
                    ax.set_ylabel(r"$\Delta$A (OD) (normalized)")
                        
                        
            elif normalize == "last":
                 
                    ax.plot(measurement.get_tau(), invert_sign_value*dynamics/dynamics[-1],
                             label = labeling)
                    ax.set_ylabel(r"$\Delta$A (OD) (normalized to long dynamic)")
                
            else:    
                
                    ax.plot(measurement.get_tau(), invert_sign_value*dynamics*1000,
                             label = labeling)
                    ax.set_ylabel(r"$\Delta$A (mOD)")
        
        
        ax.grid(True)
        
        # # plot grey area if data was deleted due to artifact in any one of the measurements
        # isnan = []
        # isnan_length = np.zeros(len(self.measurements))
        # for i, measurement in enumerate(self.measurements):
        #     _, dynamic = measurement.get_one_dynamic(wl_to_be_plotted)
        #     isnan_ = np.argwhere(np.isnan(dynamic))
        #     isnan.append(isnan_) # create list of arrays of all isnans
        #     isnan_length[i] = len(isnan_) # create list containing length of arrays in isnan list of arrays
        # longest_isnan_idx = np.argmax(isnan_length) # which array in isnan list is the longest (the measurement with the largest delete data-area)
        # if isnan[longest_isnan_idx].size > 0:
        #     plt.axvspan(self.measurements[-1].tau[isnan[longest_isnan_idx][0]], 
        #                 self.measurements[-1].tau[isnan[longest_isnan_idx][-1]], 
        #                 facecolor = "0.2", alpha = 0.25)
            
        
        
        
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])
        ax.axhline(y = 0, color = "k", linewidth = 0.5)

        # plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title = f"Probed Feature: {probed_feature}")
        ax.legend(title = f"Probed Feature: {probed_feature}")
        ax.set_xlabel(f"Delay Time ({self.measurements[0].timeunit})")
        
        # plt.title(f"{self.name}")
        
        
        plt.show()     

        return fig, ax

    def compare_spectra(self, tau_to_be_plotted=[1,10,100], xlim=[], ylim=[], save=False, normalize=False):
        """
        compare TA spectra at certain delay times for all measurements of 
        measurement row
        
        input: delay time(s) to be plotted in ps in a list
        e.g. measurement_row.compare_spectra([1,10])
        
        to do:
            
        """
        
        fig, ax = plt.subplots()
        # linestyles = ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^']
    
        if type(tau_to_be_plotted) == int:
            tau_to_be_plotted = [tau_to_be_plotted]
        
        if type(tau_to_be_plotted) == float:
            tau_to_be_plotted = [tau_to_be_plotted]
        ax.axhline(y = 0, color = "k", linewidth = 0.5)
        for idx, measurement in enumerate(self.measurements):
            labeling = self.define_labels(measurement)
            spectra = measurement.get_spectra(tau_to_be_plotted)
            if normalize == False:
                for i, spectrum in enumerate(spectra):
                    ax.plot(measurement.get_wl(), spectrum*1000,
                             label = labeling)
                ax.set_ylabel(r"$\Delta$A (mOD)")    
            if normalize == True:
                for i, spectrum in enumerate(spectra):
                    ax.plot(measurement.get_wl(), -spectrum/np.min(spectrum),
                             label = labeling)
                ax.set_ylabel(r"$\Delta$A (OD) (normalized)")    
                    
                    
        # ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title = f"Delay Time: {tau_to_be_plotted[0]} {measurement.timeunit}")
        ax.legend(title = f"Delay Time: {tau_to_be_plotted[0]} {measurement.timeunit}", bbox_to_anchor=(1.1, 1.05))
        
        ax.set_xlabel("Wavelength (nm)")
        ax.grid(True)
        
        if self.name != "None":
            ax.set_title(f"{self.name}")
        
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
            
        if ylim != []:
            ax.set_ylim(ylim[0], ylim[1])    
        
        
        # plot grey area if data was deleted due to scattered pump light in any one of the measurements
        isnan = []
        isnan_length = np.zeros(len(self.measurements))
        for i, measurement in enumerate(self.measurements):
            spectrum = measurement.get_spectra(tau_to_be_plotted)
            isnan_ = np.argwhere(np.isnan(spectrum))
            isnan.append(isnan_) # create list of arrays of all isnans
            isnan_length[i] = len(isnan_) # create list containing length of arrays in isnan list of arrays
        longest_isnan_idx = np.argmax(isnan_length) # which array in isnan list is the longest (the measurement with the largest delete data-area)
        if isnan[longest_isnan_idx].size > 0:
            ax.axvspan(self.measurements[-1].wl[isnan[longest_isnan_idx][0,1]], 
                        self.measurements[-1].wl[isnan[longest_isnan_idx][-1,1]], 
                        facecolor = "0.2", alpha = 0.25)
        
        
        
        if save == True:
            taus = ""
            for i in tau_to_be_plotted: # build string of all taus for filename
                taus = taus + "_" + str(i)
            fig.savefig(f"{today}_compare_spectra_{self.name}{taus}{self.measurements[0].timeunit}.svg", bbox_inches = "tight")
        plt.show()     

        return fig, ax

    def compare_spectra_2(self, tau_to_be_plotted=[1,10,100], xlim=[], ylim=[], 
                          save=False, normalize=False, av = False):
        """
        
        
        this needs to be merged with compare_spectra
          
        
        
        
        
        has the functionality to average spectra over a certain delay time range
        
        if av = True:
            give timerange as tau_to_be_plotted (list with two entries)
        
        
        compare TA spectra at certain delay times for all measurements of 
        measurement row
        
        input: delay time(s) to be plotted in ps in a list
        e.g. measurement_row.compare_spectra([1,10])
        
        to do:
            
        """
        
        plt.figure()
        linestyles = ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^']
    
        if type(tau_to_be_plotted) == int:
            tau_to_be_plotted = [tau_to_be_plotted]
        
        if type(tau_to_be_plotted) == float:
            tau_to_be_plotted = [tau_to_be_plotted]
      
        for idx, measurement in enumerate(self.measurements):
            labeling = self.define_labels(measurement)
            if av == False:
                spectra = measurement.get_spectra(tau_to_be_plotted)
            elif av == True:
                spectra = [0]
                spectra[0] = measurement.get_spectrum_averaged(tau_to_be_plotted)
            
            
            if normalize == False:
                for i, spectrum in enumerate(spectra):
                    plt.plot(measurement.get_wl(), spectrum*1000, linestyles[idx],
                             label = labeling)
                plt.ylabel(r"$\Delta$A (mOD)")    
            if normalize == True:
                for i, spectrum in enumerate(spectra):
                    plt.plot(measurement.get_wl(), -spectrum/np.min(spectrum), linestyles[idx],
                             label = labeling)
                plt.ylabel(r"$\Delta$A (OD) (normalized)")    
                    
                    
        plt.legend(title = f"Delay Time: {tau_to_be_plotted[0]}-{tau_to_be_plotted[1]} {measurement.timeunit}")
        plt.axhline(y = 0, color = "k", linewidth = 0.5)
        plt.xlabel("Wavelength (nm)")
        plt.grid(True)
        plt.title(f"{self.name}")
        
        if xlim != []:
            plt.xlim(xlim[0], xlim[1])
            
        if ylim != []:
            plt.ylim(ylim[0], ylim[1])    
        
        
        # plot grey area if data was deleted due to scattered pump light in any one of the measurements
        isnan = []
        isnan_length = np.zeros(len(self.measurements))
        for i, measurement in enumerate(self.measurements):
            spectrum = measurement.get_spectra(tau_to_be_plotted)
            isnan_ = np.argwhere(np.isnan(spectrum))
            isnan.append(isnan_) # create list of arrays of all isnans
            isnan_length[i] = len(isnan_) # create list containing length of arrays in isnan list of arrays
        longest_isnan_idx = np.argmax(isnan_length) # which array in isnan list is the longest (the measurement with the largest delete data-area)
        if isnan[longest_isnan_idx].size > 0:
            plt.axvspan(self.measurements[-1].wl[isnan[longest_isnan_idx][0,1]], 
                        self.measurements[-1].wl[isnan[longest_isnan_idx][-1,1]], 
                        facecolor = "0.2", alpha = 0.25)
        
        
        
        if save == True:
            taus = ""
            for i in tau_to_be_plotted: # build string of all taus for filename
                taus = taus + "_" + str(i)
            plt.savefig(f"{today}_compare_spectra_{self.name}{taus}{self.measurements[0].timeunit}.svg", bbox_inches = "tight")
        plt.show()     






    def compare_bleach_shift(self, xlim = [], fit = False):
        """
        compare bleach shift from self.wl_bleach
        fit = False for bleach analysis without fit (just reading out max values)
        fit = True for bleach analysis fitted with polynomial function
        """
        
        
        fig, ax = plt.subplots()

        for idx, measurement in enumerate(self.measurements):
            if fit == False:
                ax.plot(measurement.tau, measurement.wl_bleach, 'o',
                     label = f"{measurement.sample}, {measurement.pumpfluence}" + " $\mathrm{\mu}$J/cm$^{2}$")
            elif fit == True:
                ax.plot(measurement.tau, measurement.wl_bleach_fit, 'o', 
                     label = f"{measurement.sample}, {measurement.pumpfluence}" + " $\mathrm{\mu}$J/cm$^{2}$")
        ax.legend()
        ax.grid(True)
        ax.set_xlabel("Delay Time (ps)")
        ax.set_ylabel("Wavelength of Bleach Maximum (nm)")
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        fig.show()

        return fig, ax
     
        
        
    def compare_bleach_max(self, xlim = [], normalize = ""):
        """
        compare bleach max from self.deltaA_bleach
        as oppossed to just looking at the dynamics at one wavelength
        fit = False for bleach analysis without fit (just reading out max values)
        fit = True for bleach analysis fitted with polynomial function
        """
        # colors=['#90baefff', '#1b63c0ff', '#1b24c1ff', '#11177bff', '#61b14fff', '#fba656ff', '#ea6666ff', '#b53030ff']
        
        colors=['#b53030ff', '#ea6666ff', '#fba656ff', '#61b14fff', '#11177bff', '#1b24c1ff', '#1b63c0ff', '#90baefff']
        
        
        
        fig, ax = plt.subplots()

        for idx, measurement in enumerate(self.measurements[::-1]):
            labeling = self.define_labels(measurement)
            if normalize == "max":
                ax.plot(measurement.tau, 
                        - (np.array(measurement.deltaA_bleach)/np.max(np.abs(np.array(measurement.deltaA_bleach)))), 
                        label = labeling, color=colors[idx])
                ax.set_ylabel("$\Delta$A$_{max}$ (OD) (normalized)")
            else:
                ax.plot(measurement.tau, -1000 * np.array(measurement.deltaA_bleach), 
                        label = labeling, color=colors[idx])
                ax.set_ylabel("$\Delta$A$_{max}$ (mOD)")
                
        ax.legend()
        ax.grid(True)
        ax.set_xlabel(f"Delay Time ({self.measurements[0].timeunit})")
        
        if xlim != []:
            ax.set_xlim(xlim[0], xlim[1])
        fig.show()

        return fig, ax




    def mean_bleach_at_long_delays(self, t_long = 4000, guess_A=2):
        """
        run analyse_max_bleach on each measurement before using this method (to take care of bleach shifts) 
        
        
        calculate mean bleach at long delay time for poisson statistics
        see J. Phys. Chem. Lett. 2020, 11, 5361−5366 and supporting information thereof
        
        
        t_long: delay time at which auger processes are certainly over
        guess_A: approximate value at which bleach on long delays saturates in mOD




        adding y-intercept to poisson equation makes fit much better (but does not have any physical meaning)
        what values to use for y-error? --> weighted fitting
        calculate std for absorption cross section 

        
        """
        
        
        
        
        mean_deltaAs = []
        std_deltaAs = []        
        pumpfluences = np.array(self.get_list_of_pumpfluences())
        
        for measurement in self.measurements:
    
            
            idx = next(i for i,v in enumerate(measurement.tau) if v > t_long)
            mean_deltaA = np.mean(measurement.deltaA_bleach[idx:])
            mean_deltaAs.append(mean_deltaA)

            std_deltaA = np.std(measurement.deltaA_bleach[idx:])*1000
            std_deltaAs.append(std_deltaA)
        yerr = std_deltaAs
        
        # yerr = 0.1*np.array(mean_deltaAs)*1000
        
        # how to determine standart deviation of value that is not stable
        
        def poisson(j, c, A):
            return A * (1 - np.exp(-c*j))

        guess = [guess_A*10**-3,0.01]
        
        

        popt, pcov = curve_fit(poisson, 
                               pumpfluences, 
                               np.array(mean_deltaAs), 
                               p0 = guess,
                               maxfev = 5000)
                               # sigma = yerr,
                               # absolute_sigma = True)
        
        x_axis = np.arange(0, np.max(pumpfluences))
        
        fit = poisson(x_axis, *popt)

        self.c = popt[0]
        self.A = popt[1]
        
       


        # plot figure with second x-axis (<N>) and second y-axis (P(N>=1)))

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twiny()
        ax3 = ax1.twinx()
        
        new_tick_labelsy = np.array([0, 0.25, 0.5, 0.75, 1])
        
        def tick_function_y(y, A):
            V = A * y
            return [z*1000 for z in V]
        
        
        ax3.set_yticklabels(new_tick_labelsy)
        ax3.set_yticks(tick_function_y(new_tick_labelsy, self.A))
        ax3.set_ylabel("P(N$\geq$1)")
        
        
        xerr = np.zeros(len(pumpfluences))
        xerr.fill(5)
        
        ax1.plot(pumpfluences, np.array(mean_deltaAs) * 1000,   'o', color = "tab:orange", alpha = 0.5)
        ax1.errorbar(pumpfluences, np.array(mean_deltaAs)*1000,xerr=xerr, yerr=yerr, fmt='none', color="tab:orange")
        ax1.plot(x_axis, fit*1000, color = "tab:orange")
        
        ax1.set_xlabel(r"Pumpfluence ($\mathrm{\mu}$J/cm$^{2}$)")
        ax1.set_ylabel(r"Mean $\Delta$A on long Delays (mOD)")
        
        
        single_exc_regime = [0.3] # <N> that is still regarded as single-exciton regime
        single_exc_regime_2 = [0.4]
        new_tick_labels = np.array([0,1,2,3,4,5,6,7,8,9,10])
        
        
        def tick_function(x):
            V = x/self.c
            return [z for z in V]
        
        x_tick_locations = tick_function(new_tick_labels)
        vertical_line_single_exc_regime_x_location = tick_function(single_exc_regime)
        vertical_line_single_exc_regime_x_location_2 = tick_function(single_exc_regime_2)
        # print(x_tick_locations)
        
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(x_tick_locations)
        ax2.set_xticklabels(new_tick_labels)
        ax2.set_xlabel("<N>")
        ax2.set_xlim(ax1.get_xlim())
        ax2.axvline(vertical_line_single_exc_regime_x_location[0], ymin=0.05, ymax=0.95, ls='--')
        ax2.axvline(vertical_line_single_exc_regime_x_location_2[0], ymin=0.05, ymax=0.95, ls='--')
        # ax3.set_ylim(0,1.2)
        ax3.set_ylim(ax2.get_ylim())
        ax2.grid(True)
        ax3.grid(True)
        
        plt.show()
        
        

        # print(f"Single-exciton regime (<N> = {single_exc_regime[0]}) ends at {np.round(vertical_line_single_exc_regime_x_location[0])}", r" mu J/cm^2")

        return fig, [ax1, ax2, ax3]

    def get_crosssection(self, pump_wavelength):
        """ 
        calculate absorption cross section of particles from C from mean_bleach_at_long_delays analysis
        """
        

        # calculate photon absorption cross section
        # sigma = E(photon) * C
        # E = h*c/lambda
        
        E = sc.h * sc.speed_of_light / (pump_wavelength * 10**-9)
        
        
        
        self.absorption_crosssection = E * self.c*10**6 # c is in cm2/muJ  
        
        print(f"Absorption cross section is {self.absorption_crosssection} cm^2")

    def merge_dynamics(self, wl_to_be_plotted, xlim=[], log = False):
        """
        merge dynamics of fs-TA and mus-TA
        """
        if len(self.measurements) == 2 and self.measurements[0].ta_type == "HELIOS" and self.measurements[1].ta_type == "EOS" and self.measurements[1].timeunit != "ns":
            
            _, fs = self.measurements[0].get_one_dynamic(wl_to_be_plotted)
            fs_timescale = self.measurements[0].tau
            _, mus = self.measurements[1].get_one_dynamic(wl_to_be_plotted)
            mus_timescale = self.measurements[1].tau
            
            # calculate scaling factor
            t = fs_timescale[-1]
            idx_mus = next(i for i,v in enumerate(mus_timescale*10**6) if v > t)
            od_mus = mus[idx_mus]
            od_fs = fs[-1]
            
            scaling_factor = od_fs/od_mus
            
            self.merged_dyn = np.append(fs, scaling_factor*mus[idx_mus:])
            fs_timescale_in_ns = fs_timescale*10**-3
            mus_timescale_in_ns = mus_timescale*10**3
            self.merged_dyn_tau = np.append(fs_timescale_in_ns, mus_timescale_in_ns[idx_mus:])
            
            
            
            
            plt.figure()
            plt.plot(fs_timescale*10**-3, fs*10**3, label = "fs-TA")
            plt.plot(mus_timescale*10**3, scaling_factor*mus*10**3, label = "$\mathrm{\mu}$s-TA" + f", scaled {scaling_factor}")
            plt.ylabel("$\Delta$A (mOD)")
            plt.xlabel("Delay Time (ns)")
            plt.grid(True)
            plt.legend()
            
            if xlim != []:
                plt.xlim(xlim[0], xlim[1])
            if log == True:
                plt.xscale("log")
            
            plt.show()     
            
            
            plt.figure()
            plt.plot(self.merged_dyn_tau, self.merged_dyn*10**3, label = f"Probe Wavelength: {wl_to_be_plotted} nm")
            plt.ylabel("$\Delta$A (mOD)")
            plt.xlabel("Delay Time (ns)")
            plt.grid(True)
            plt.legend()
            
            if xlim != []:
                plt.xlim(xlim[0], xlim[1])
            if log == True:
                plt.xscale("log")
            
            plt.show()     
            
            
            
        else:
            print("Measurement Row has to contain two measurements.")
            print("The first measurement has to be the HELIOS measurement.")
            print("The second measurement has to be the EOS measurement.")
            print("Do not convert timeunit to ns before using this method.")




    def compare_absolute_bleach(self, delaytime_range, wavelength_range):
        """
        generate figure to compare absolute bleach at specified delay time and wavelength for different 
        measurement of one measurement row
        
        """
        pumpfluences = np.array(self.get_list_of_pumpfluences())
        
        list_deltaA = []
        list_deltaA_std = []
        for meas in self.measurements:
            list_deltaA.append(meas.get_deltaA_value_av(delaytime_range, wavelength_range)[0])
            list_deltaA_std.append(meas.get_deltaA_value_av(delaytime_range, wavelength_range)[1])
            
            
            
        fig, ax = plt.subplots()
        ax.errorbar(pumpfluences, 1000*np.array(list_deltaA), yerr = 1000*np.array(list_deltaA_std), fmt = 'o')
        ax.set_ylabel("$\Delta$A (mOD)")
        ax.set_xlabel("Pump Fluence ($\mathrm{\mu}$J/cm$^{2}$)")
        ax.grid(True)
        plt.show()

        return fig, ax



























