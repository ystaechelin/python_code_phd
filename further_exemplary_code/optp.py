# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 10:05:50 2022

@author: staechelin
"""


##############################################################################
#                                                                            #
#                                                                            #
#               object-oriented programming approach to                      #
#                         OPTP data                                          #
#                         08.02.2022                                         #
#                                                                            #
#                                                                            #
##############################################################################


### import section ###########################################################
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.constants as sc
from scipy.optimize import minimize, curve_fit
from datetime import date
from optp_models import drude, drude_smith, lorentz_drude, comb_model
today = date.today()
today = today.strftime("%y%m%d")
##############################################################################


### to do ####################################################################
# make every method work if only one THz-delay or OP-Delay was measured




# fit Drude etc.
# calculate refractive index from conductivity
# if save == True: plt.savefig einbauen für die einzelnen Abbildungen... mit Datum speichern
# calculate frequency range in which "enough" THz is present.. 5 % of maximum intensity?
# quality of fit
# revise figure of fit results 
##############################################################################





class OPTP():
    
    def __init__(self, path, measurement, name  = "", zero_point_in_mm = 0, x_name = "", y_name = ""):
        
        """
        initialize object by loading the already averaged measurement data
        """
        
        self.path = path
        self.measurement = measurement
        self.name = name 
        self.data, self.x_mm, self.y = self.load_data(self.path, self.measurement, x_name, y_name)
        # self.data = np.array(self.data)
        # self.data = np.reshape(self.data, (1, self.data.size))
        self.axis_to_ps(zero_point_in_mm)
        self.stepwidth = self.x_mm[1] - self.x_mm[0]
        
        
    def load_data(self, path, measurement, x_name = "", y_name = ""):
        
        """ 
        loading the data (already averaged) of a single sample
        """
        
        load_data_from = path + "\\" + measurement
        if measurement[-4:] == ".dat":
            self.data = np.genfromtxt(load_data_from, delimiter = "\t")
        else:
            self.data = np.genfromtxt(load_data_from + "_av.txt", delimiter = "\t")
        self.data = -self.data
        # self.data = np.reshape(self.data, (self.data.size, ))
        # self.data = self.data.T
        if x_name == "":
            self.x_mm = np.genfromtxt(load_data_from[:-4] + "x.dat", delimiter = "\t")
        else: 
            self.x_mm = np.genfromtxt(path + "\\" + x_name, delimiter = "\t")
        if y_name == "":    
            self.y = np.genfromtxt(load_data_from[:-4] + "y.dat", delimiter = "\t")
        else:
            self.y = np.genfromtxt(path + "\\" + y_name, delimiter = "\t")
        return self.data, self.x_mm, self.y
        
    
    
    def axis_to_ps(self, zero_point_in_mm):
        """ 
        convert x and y axis from mm delay stage movement to ps
        """
        
        self.x = (2 * self.x_mm * 10**-3)/sc.speed_of_light 
        # self.x = np.array(self.x)
        self.x = self.x + np.abs(self.x[0])  # sets beginning of measurement to time = 0
        self.y = self.y - zero_point_in_mm
        self.y = (2 * self.y * 10**-3)/sc.speed_of_light
        self.y = np.array(self.y) # change to array (in case its only one OP delay time, it would otherwise be just a float)
        self.y = np.reshape(self.y, (self.y.size, ))
        
        
    def correct_zeropoint(self, zero_point):
        """
        correct zero point of optical-pump-delay, input zero point in ps
        """
        self.y = self.y - zero_point*10**-12
        
        
           
    def cut_reflected_pulse(self, delay_time_at_which_reflection_starts_in_ps):
        
        """
        cut off reflection
        input: delay time at which reflection starts
        
        just cut data or zero-filling?
        """        
        
        idx = next(i for i,v in enumerate(self.x) if v > delay_time_at_which_reflection_starts_in_ps * 10**-12)
        if len(self.y) == 1:
            self.data[idx:] = 0
        else:
            self.data[:, idx:] = 0
        return self.data
    

    
    def plot_map(self, xlim = [], ylim = [], save = False):
        """ 
        plot map (only makes sense if whole map was measured)
        add option to save figure to .svg
        add option to crop data to certain ranges
        first plot is always way too smooth
        """        
        
        # fig, ax = plt.subplots()
        # cf = ax.contourf(self.x*10**12, self.y*10**12, self.data, 50, cmap="seismic")
        # cbar = fig.colorbar(cf, ax=ax)
        # cbar.set_label(r"$\Delta$E (t, $\tau$) (a.u.)")
        # ax.set_ylabel(r'Delay Time ($\tau$) (ps)')
        # ax.set_xlabel('Time (t) (ps)')
        # plt.show()
    
        fig, ax = plt.subplots()
        plt.title(self.name)
        plt.set_cmap("seismic")
        plt.pcolormesh(self.x*10**12, self.y*10**12, self.data*10**3, 
                        vmin= -np.max(np.abs(self.data*10**3)), 
                        vmax = np.max(np.abs(self.data*10**3)),
                        rasterized=True, shading = "auto")
        plt.xlabel('Time (t) (ps)')
        plt.ylabel(r'Delay Time ($\tau$) (ps)')
        if xlim != []:
            plt.xlim(xlim[0], xlim[1])
        if ylim != []:
            plt.ylim(ylim[0], ylim[1])
        cbar = plt.colorbar()
        # plt.gca().set_aspect('equal')
        cbar.set_label(r"$\Delta$E (t, $\tau$)")
        if save == True:
            plt.savefig(f"{today}_{self.name}_map.svg", bbox_inches = "tight")
        plt.show()


    def plot_map_pulses_dyn(self, thz_delay_to_be_plotted, op_delay_to_be_plotted):
        """ 
        plot map (only makes sense if whole map was measured)
        add option to save figure to .svg
        add option to crop data to certain ranges
        first plot is always way too smooth
        """        
        fig, ax = plt.subplots(2, 2, 
                       gridspec_kw={
                           'width_ratios': [3, 1.5],
                           'height_ratios': [3, 2]})
    
        # plot map
        plt.set_cmap("seismic")
        # ax[0][0].set_aspect("equal")
        ax[0][0].pcolormesh(self.x*10**12, self.y*10**12, self.data, 
                       vmin= -np.max(np.abs(self.data)), 
                       vmax = np.max(np.abs(self.data)), 
                       rasterized=True, shading = "auto")
        ax[0][0].set_xlabel('Time (t) (ps)')
        ax[0][0].set_ylabel(r'Delay Time ($\tau$) (ps)')
        
        
        
        # cbar = plt.colorbar()
        # # plt.gca().set_aspect('equal')
        # cbar.set_label(r"$\Delta$E (t, $\tau$) (a.u.)")
        
        
        
        #plot dynamics
        indices = []
        for t in thz_delay_to_be_plotted:
            index = next(i for i,v in enumerate(self.x) if v > t*10**-12)
            indices.append(index)
            
        for n in range(len(indices)):
            ax[0][1].plot(-self.data[:, indices[n]]*10**3, self.y,  
                          label = str(round(self.x[indices[n]]*10**12, 2)) + " ps")
            
        
        ax[0][1].set_xlabel(r"$\Delta$A")
        ax[0][1].legend()
        
        
        
        # plot spectra
        indices = []
        for tau in op_delay_to_be_plotted:
            index = next(i for i,v in enumerate(self.y) if v > tau*10**-12)
            indices.append(index)
            
   
        for n in range(len(indices)):
            ax[1][0].plot(self.x*10**12, self.data[indices[n], :]*10**3, 
                     label = str(round(self.y[indices[n]]*10**12, 0))[:2] + " ps")
        

        ax[1][0].set_ylabel(r"$\Delta$E")
        ax[1][0].set_xlabel(r"Time (t) (ps)")
        ax[1][0].set_xlim(0,10)
        ax[1][0].legend()
        
        
        # wl_to_be_plotted2 = [505]
        
        # indices = []
        # for wl in wl_to_be_plotted2:
        #     index = next(i for i,v in enumerate(self.wl) if v > wl)
        #     indices.append(index)
            
        # for n in range(len(indices)):
        #     ax[1][1].plot( self.tau, self.data[indices[n],:],  label = str(round(self.wl[indices[n]], 0))[:-2] + " nm")
   
        
        
        # ax[1][1].axhline(y = 0, color = "k", linewidth = 0.3)
        # ax[1][1].set_xlim(-0.25,5)
        # ax[1][1].set_ylim(-0.02,0.01)
        # ax[1][1].set_ylabel(r"$\Delta$A")
        # ax[1][1].set_xlabel(r"Delay Time (ps)")
        # ax[1][1].legend()
       
    
        plt.show()
        
        
        
        
        
        

    def get_dynamics(self, thz_delay):
        """
        return dynamics for certain THz delay
        input THz delay in ps
        """
        index = next(i for i,v in enumerate(self.x) if v > thz_delay*10**-12)
        dynamics = self.data[:, index]
        dynamics_delay = self.x[index] # 
        return dynamics, dynamics_delay
    
    def get_spectrum(self, op_delay):
        """
        return spectrum for certain optical pump delay
        input optical pump delay in ps
        """
        index = next(i for i,v in enumerate(self.y) if v > op_delay*10**-12)
        spec = self.data_fft[index][:,1]
        return spec
    
    
    def plot_dynamics(self, thz_delays, normalize = False, xlim = [], ylim = []):
        """
        plot dynamics at a certain THz-delay        
        """
    
        plt.figure()
        for thz_delay in thz_delays:
            if normalize == True:
                if np.max(self.get_dynamics(thz_delay)[0]) < np.abs(np.min(self.get_dynamics(thz_delay)[0])): 
                    plt.plot(self.y*10**12, -1*self.get_dynamics(thz_delay)[0]/np.max(np.abs(self.get_dynamics(thz_delay)[0])), 
                             label = str(round(self.get_dynamics(thz_delay)[1]*10**12, 2)) + " ps")
                else:
                    plt.plot(self.y*10**12, self.get_dynamics(thz_delay)[0]/np.max(np.abs(self.get_dynamics(thz_delay)[0])), 
                             label = str(round(self.get_dynamics(thz_delay)[1]*10**12, 2)) + " ps")
            else:
                plt.plot(self.y*10**12, self.get_dynamics(thz_delay)[0]*10**3, 
                         label = str(round(self.get_dynamics(thz_delay)[1]*10**12, 2)) + " ps")
        plt.legend(title = "Time (t):", framealpha = 1)
        plt.grid(True)
        plt.ylabel(r"$\Delta$E (a.u.)")
        plt.xlabel(r'Delay Time ($\tau$) (ps)')
        if xlim != []:
            plt.xlim(xlim[0], xlim[1])
        if ylim != []:
            plt.ylim(ylim[0], ylim[1])
        # plt.savefig(r"L:\Yannic\Besprechungen und Vorträge\220301_GrK_Talk\thz_dynamics.svg", bbox_inches = "tight")
        plt.show()

     
    def get_diff_pulse(self, optical_pump_delay):
        """
        return differential pulse data for certain optical pump delay
        input optical pump delay in ps
        """
        
        index = next(i for i,v in enumerate(self.y) if v > optical_pump_delay*10**-12)
        diff_pulse = self.data[index, :]
        diff_pulse_delay = self.y[index]
        return diff_pulse, diff_pulse_delay
    
    
    
    def plot_diff_pulse(self, optical_pump_delays, xlim = [], ylim = []):
        """
        plot differential THz pulse at a certain optical-pump-delay
        input optical pump delay in ps
        if multiple differential pulses were recorded (multiple optical pump delays): use self.get_diff_pulse
        
        if only one differential pulse (only one optical pump delay): use this one
        """
       
        plt.figure()
        if len(self.y) > 1:
            for optical_pump_delay in optical_pump_delays:
                plt.plot(self.x*10**12, self.get_diff_pulse(optical_pump_delay)[0]*10**3, 
                         label = str(round(self.get_diff_pulse(optical_pump_delay)[1]*10**12)) + " ps")
        if len(self.y) == 1:
            plt.plot(self.x*10**12, self.data*10**3, label = str(round(self.y[0]*10**12)) + " ps")
        plt.legend(title = r"Delay Time ($\tau$):", framealpha = 1)
        plt.grid(True)
        plt.ylabel(r"$\Delta$E (a.u.)")
        plt.xlabel('Time (t) (ps)')
        if xlim != []:
            plt.xlim(xlim[0], xlim[1])
        if ylim != []:
            plt.ylim(ylim[0], ylim[1])
        plt.show()
    
    def plot_tds(self):
        """
        plot TDS without OPTP data
        
        """
        plt.figure()
        plt.plot(self.tds_x_in_ps*10**12, self.tds[:,1], label = "E")
        plt.ylabel("Electric Field (a.u.)")
        plt.xlabel("Delay Time (ps)")
        plt.legend()
        plt.show()
        

    def correct_offset_optp(self, delay_from, delay_to):
        """ 
        correct offset in optp data
        input: 
            give intervall of THz-delays, in which no signal can be seen
            in ps
        """
        index_0 = next(i for i,v in enumerate(self.x) if v > delay_from * 10**-12)
        index_1 = next(i for i,v in enumerate(self.x) if v > delay_to * 10**-12)
     
        
        for i,diff_pulse in enumerate(self.data):
            offset = np.mean(diff_pulse[index_0 : index_1])
            self.data[i,:] = diff_pulse - offset            




        
    
    
    def load_corresponding_TDS(self, path_tds, measurement_tds):
        """ 
        provide TDS of the sample which belongs to the optp data
        input path and name of tds file
        """
        self.tds = np.genfromtxt(path_tds + "\\" + measurement_tds, delimiter = "\t")
     
        
        
        
        
    def plot_optp_tds(self, delaytimes_to_be_plotted):
        """
        is this one still needed?
        plot_optp_tds1 (see below) works better
        
        plot optp and tds together
        input list of optical pump delay times to be plotted
        
        
        to do: add possibility to scale up differential signal if needed
        its not yet possible to plot with x-axis converted to ps
        """
        plt.figure()
        plt.plot(self.tds_x_in_ps*10**12, self.tds[:,1], label = "E")
        for time in delaytimes_to_be_plotted:
            plt.plot(self.x*10**12, self.get_diff_pulse(time)[0]*100, 
                 label = r"$\Delta$E, " + str(round(self.get_diff_pulse(time)[1]*10**12,1)) + " ps")
        plt.legend()
        plt.ylabel("Electric Field (a.u.)")
        plt.xlabel(r"$\Delta$t (ps)")
        plt.grid(True)
        # plt.savefig(r"L:\Yannic\Besprechungen und Vorträge\220301_GrK_Talk\thz_optp_tds.svg", bbox_inches = "tight")
        plt.show()    
        
           
    
    
    def plot_optp_tds1(self, delaytime_to_be_plotted, scale_up = 1):
        """
        plot optp and tds together
        input  optical pump delay time to be plotted
        
        
        to do: add possibility to scale up differential signal if needed
        
        
        """
        
        self.bring_to_same_x_axis()
        self.tds_x_axis_to_ps()
        
        
        index = next(i for i,v in enumerate(self.y) if v > delaytime_to_be_plotted*10**-12)
        diff_pulse = self.data_to_tds_x_axis[index, :]
        diff_pulse_index = self.y[index]
        
        plt.figure()
        plt.plot(self.tds_x_in_ps*10**12, self.tds[:,1], label = "E")

        plt.plot(self.tds_x_in_ps*10**12, scale_up*diff_pulse, 
                 label = f"$\Delta$E (" + str(round(diff_pulse_index*10**12,1)) + " ps)")
        plt.legend()
        plt.ylabel("Electric Field (a.u.)")
        plt.xlabel(r"$\Delta$t (ps)")
        plt.grid(True)
        plt.show() 
    
    
    
  
        
        
    def bring_to_same_x_axis(self):
        """
        to compare the Fourier transforms of differential and TDS pulse, they need to have the same x axis
        usually you would measure a smaller intervall for the differential measurement, since it takes much longer
        than the TDS
        so we zero-fill the differential pulse to have the same length of x-axis as the TDS pulse
        
        differential pulse and TDS should be measured with the same step distance
        
        
        it works the following way
        1) find out from which to which indice the x-axis of the differential pulse goes within the x-axis of the TDS pulse
        2) create arrays with zeros which have the length that is needed to fill the optp data
        3) concatenate the zeros-arrays with the actual data
    
    
        builds self.data_to_tds_x_axis, the optp data which matches with the x axis of the TDS data 
        """
        #1
        index_0 = np.where(self.tds[:,0] == self.x_mm[0])
        index_1 = np.where(self.tds[:,0] == self.x_mm[-1])
        
        #2        
        zerof_0 = np.zeros((int(index_0[0]), self.y.size))
        zerof_1 = np.zeros((len(self.tds[:,0])-int(index_1[0])-1, self.y.size))
        
        #3
        self.data = np.reshape(self.data, (self.y.size, -1))
        self.data_to_tds_x_axis = np.c_[zerof_0.T, self.data]
        self.data_to_tds_x_axis = np.c_[self.data_to_tds_x_axis, zerof_1.T]
        
    def zero_fill_optp(self, number_of_datapoints = 500):
        """
        zero fill optp independently of TDS
        """
        zeros = np.zeros((number_of_datapoints, len(self.y)))
        stepwidth_x = self.x[1]-self.x[0]
        new_steps_x = np.arange(0, number_of_datapoints)*stepwidth_x + self.x[-1]
        self.x = np.hstack((self.x, new_steps_x))
        self.data = np.c_[self.data, zeros.T]
        
        
    def tds_x_axis_to_ps(self):
        """
        transform x axis of TDS data to ps
        """
        self.tds_x_in_ps = (2 * self.tds[:,0] * 10**-3)/sc.speed_of_light 
        self.tds_x_in_ps = self.tds_x_in_ps + np.abs(self.tds_x_in_ps[0])  # sets beginning of measurement to time = 0
        
    
 
    
    def additional_zero_filling(self, number_of_datapoints=500):
        """
        add additional zero filling to tds  data; might help the phase-unwrapping
        adds the number of datapoints end of data
        
        if you fourier-transform afterwards, the optp data will have the same extension of zero-filling
        """
        
        # TDS-data
        zerof = np.zeros((number_of_datapoints,2))
        zerof[:,0] = np.arange(1, number_of_datapoints + 1) * self.stepwidth + self.tds[-1, 0]
        
      
        self.tds = np.concatenate((self.tds, zerof))
        
        
        

    
    def fourier(self):
        """ 
        prepare fourier transformation by bringing TDS and differential pulse to same x-axis
        and converting x-axis to ps
        and fourier transform using the fft method
        fourier transforms 1. the tds data and 2. the optp data (each optical pump delay gets fourier transformed)
        builds self.data_fft, a list containing arrays containing the fft-results
        """
        self.bring_to_same_x_axis()
        self.tds_x_axis_to_ps()
        self.tds_fft = self.fft(self.tds[:,1], self.tds_x_in_ps)
    
        self.data_fft = [None] * len(self.y)
        for i, diff_pulse_at_certain_op_delay in enumerate(self.data_to_tds_x_axis):
            fft_v = self.fft(diff_pulse_at_certain_op_delay, self.tds_x_in_ps)
            self.data_fft[i] = fft_v
    

    
    def fourier_optp(self):
        """
        fourier transform only differential pulses (and not tds)
        """
        
        # if len(self.y) == 1:
        #     self.data_fft = self.fft(self.data, self.x)
        #     self.data_fft = self.data_fft.reshape(len(self.data_fft),1)
            
        # else:
            # self.data_fft = [None] * len(self.y)
            # for i, diff_pulse_at_certain_op_delay in enumerate(self.data):
            #     fft_v = self.fft(diff_pulse_at_certain_op_delay, self.x)
            #     self.data_fft[i] = fft_v
    
    
    
        self.data_fft = [None] * len(self.y)
        if len(self.y) > 1:
            for i, diff_pulse_at_certain_op_delay in enumerate(self.data):
                fft_v = self.fft(diff_pulse_at_certain_op_delay, self.x)
                self.data_fft[i] = fft_v
        
        elif len(self.y) == 1:
            self.data_fft[0] = self.fft(self.data, self.x)
    
    def fft(self, y, x):
        """
        does the actual job of fourier transforming the data
        input y: the data you want to fourier transform
        """        
        timestep = x[2]-x[1]
    
        sp_synthetic = np.fft.ifft(y)#*self.tds_x_in_ps.size
        #sp_synthetic = np.fft.fft(y)
    
        freq = np.fft.fftfreq(x.size, d = timestep)
    
    
        np.data = np.zeros((freq.size, 5))
        np.data[:,0] = freq
        np.data[:,1] = abs(sp_synthetic)
        np.data[:,2] = np.real(sp_synthetic)
        np.data[:,3] = np.imag(sp_synthetic)
        np.data[:,4] = np.angle(sp_synthetic)
    
        #cut out negativ frequency range
        nonzero = np.where(np.data[:,0]>0)[0]
        np.data_reduced = np.data[nonzero,:]
    
    
        return np.data_reduced
           
    
    def plot_diff_spec(self, optical_delay_times, xlim=[0,2], log=True, normalize=False, save=False):
        """
        plot amplitude of differential spectrum at given delay times
        """
        
   
        plt.figure()
        
        
        x_indices = []
        for x_value in xlim:
            idx = next(i for i,v in enumerate(self.data_fft[0][:,0]) if v > x_value * 10**12)
            x_indices.append(idx)
            
            
        for optical_delay_time in optical_delay_times:
            index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
            
            if normalize == True:
                plt.plot(self.data_fft[index][x_indices[0]:x_indices[1],0]*10**-12, 
                         self.data_fft[index][x_indices[0]:x_indices[1],1]/np.max(self.data_fft[index][x_indices[0]:x_indices[1],1]),  
                         label = str(round(self.y[index]*10**12,1)) + " ps")
                
                
            else:      
                plt.plot(self.data_fft[index][x_indices[0]:x_indices[1],0]*10**-12, 
                         self.data_fft[index][x_indices[0]:x_indices[1],1], 
                         label = str(round(self.y[index]*10**12,1)) + " ps")
            
            
        plt.ylabel("Ampitude of differential Spectrum")
        plt.xlabel("Frequency (THz)")
        plt.grid(True)
        plt.legend()
        if log == True:
            plt.yscale("log")
        
        if save == True:
            plt.savefig(f"{today}_{self.name}_absolute_spectra.svg", bbox_inches = "tight")
        
        
        plt.show()
        
    def get_dynamics_freq(self, frequency):
        """
        get dynamics at certain THz frequency
        """
        
        idx = next(i for i,v in enumerate(self.data_fft[0][:,0]) if v > frequency * 10**12)
        
        dyn = []
        for n, op_delay_time in enumerate(self.y):
            dyn.append(self.data_fft[n][idx,1])
        
        return np.array(dyn)
        
    def plot_dynamics_amplitude(self, frequencies, xlim=[], ylim=[], normalize=False):
        """
        plot dynamics of amplitude of differential spectrum at given thz frequencies
        """
        plt.figure()
        
        if type(frequencies) == int:
            frequencies = [frequencies]
        
        if type(frequencies) == float:
            frequencies = [frequencies]
        
        for frequency in frequencies:
        
            idx = next(i for i,v in enumerate(self.data_fft[0][:,0]) if v > frequency * 10**12)
        
            dyn = []
            for n, op_delay_time in enumerate(self.y):
                dyn.append(self.data_fft[n][idx,1])
            
            if normalize == True:
                #dyn = dyn/np.max(dyn)
                dyn = dyn/dyn[-1]
            plt.plot(self.y*10**12, dyn, 'o', label = f"{frequency} THz")
                
        plt.ylabel("Amplitude of differential Spectrum")
        if normalize == True:
            plt.ylabel("Amplitude of differential Spectrum (normalized)")
        plt.xlabel("Delay Time (ps)")
        plt.grid(True)
        plt.legend()
    
        if xlim!=[]:
            plt.xlim(xlim[0], xlim[1])
    
        
        if ylim!=[]:
            plt.ylim(ylim[0], ylim[1])
    
        plt.show()
            
    
    def plot_absolute_phase(self, optical_delay_times, xlim=[0,2], log=True, save=False):
        """
        plot the absolute and the phase of the fourier transformed spectra for both tds and optp at a certain delay time
        input: list of optical delay times in units of ps
        provide xlim as list with lower and upper lim in units of THz
        log: if set to True, y scale is logarithmic, default : True
        """

        # plot absolute
        plt.figure()
        if xlim != []:
            x_indices = []
            for x_value in xlim:
                idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
                x_indices.append(idx)
            plt.plot(self.tds_fft[x_indices[0]:x_indices[1], 0]*10**-12, 
                     self.tds_fft[x_indices[0]:x_indices[1], 1], label = r"E")
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                plt.plot(self.data_fft[index][x_indices[0]:x_indices[1],0]*10**-12, 
                         self.data_fft[index][x_indices[0]:x_indices[1],1], 
                         label = str(round(self.y[index]*10**12,1)) + " ps")
            
        else:
            plt.plot(self.tds_fft[:,0]*10**-12, self.tds_fft[:,1], label = r"E")
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                plt.plot(self.data_fft[index][:,0]*10**-12, self.data_fft[index][:,1], 
                         label = str(round(self.y[index]*10**12,1))  + " ps")
            
            
        plt.ylabel("Absolute of Spectrum")
        plt.xlabel("Frequency (THz)")
        plt.grid(True)
        plt.legend()
        if log == True:
            plt.yscale("log")
        
        if save == True:
            plt.savefig(f"{today}_{self.name}_absolute_spectra.svg", bbox_inches = "tight")
        
        
        plt.show()
        
        
        
        
        # plot phase
        
        plt.figure()
        if xlim != []:
            x_indices = []
            for x_value in xlim:
                idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
                x_indices.append(idx)
            plt.plot(self.tds_fft[x_indices[0]:x_indices[1],0]*10**-12, 
                    np.unwrap(self.tds_fft[x_indices[0]:x_indices[1],4]), label = r"E")
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                plt.plot(self.data_fft[index][x_indices[0]:x_indices[1],0]*10**-12, 
                         np.unwrap(self.data_fft[index][x_indices[0]:x_indices[1],4]), 
                         label = str(round(self.y[index]*10**12,1)) + " ps")
            
        else:
            plt.plot(self.tds_fft[:,0]*10**-12, np.unwrap(self.tds_fft[:,4]), label = r"E")
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                plt.plot(self.data_fft[index][:,0]*10**-12, np.unwrap(self.data_fft[index][:,4]), 
                          label = str(round(self.y[index]*10**12,1))  + " ps")
            
            
        plt.ylabel("Phase of Spectrum")
        plt.xlabel("Frequency (THz)")
        plt.grid(True)
        plt.legend()
        
        if save == True:
            plt.savefig(f"{today}{self.name}_phase_spectra.svg", bbox_inches = "tight")
        
        plt.show()
        
  
        
    def plot_absolute_phase_in_one(self, optical_delay_times, xlim=[0,2], log=True, save=False):
        """
        plot the absolute and the phase of the fourier transformed spectra for both tds and optp at a certain delay time
        input: list of optical delay times in units of ps
        provide xlim as list with lower and upper lim in units of THz
        log: if set to True, y scale is logarithmic, default : True
        """

        # plot absolute


        fig, ax = plt.subplots(2, 1)

    
        if xlim != []:
            x_indices = []
            for x_value in xlim:
                idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
                x_indices.append(idx)
            ax[0].plot(self.tds_fft[x_indices[0]:x_indices[1], 0]*10**-12, 
                     self.tds_fft[x_indices[0]:x_indices[1], 1], label = r"E")
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                ax[0].plot(self.data_fft[index][x_indices[0]:x_indices[1],0]*10**-12, 
                         self.data_fft[index][x_indices[0]:x_indices[1],1], 
                         label = str(round(self.y[index]*10**12,1)) + " ps")
            
        else:
            ax[0].plot(self.tds_fft[:,0]*10**-12, self.tds_fft[:,1], label = r"E")
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                ax[0].plot(self.data_fft[index][:,0]*10**-12, self.data_fft[index][:,1], 
                         label = str(round(self.y[index]*10**12,1))  + " ps")
            
            
        ax[0].set_ylabel("Absolute of Spectrum")
        ax[0].set_xlabel("Frequency (THz)")
        # plt.grid(True)
        ax[0].legend()
        if log == True:
            ax[0].set_yscale("log")
        
        # if save == True:
        #     plt.savefig(f"{today}_{self.name}_absolute_spectra.svg", bbox_inches = "tight")
        
        
        
        if xlim != []:
            x_indices = []
            for x_value in xlim:
                idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
                x_indices.append(idx)
            ax[1].plot(self.tds_fft[x_indices[0]:x_indices[1],0]*10**-12, 
                    np.unwrap(self.tds_fft[x_indices[0]:x_indices[1],4]), label = r"E")
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                ax[1].plot(self.data_fft[index][x_indices[0]:x_indices[1],0]*10**-12, 
                         np.unwrap(self.data_fft[index][x_indices[0]:x_indices[1],4]), 
                         label = str(round(self.y[index]*10**12,1)) + " ps")
            
        else:
            ax[1].plot(self.tds_fft[:,0]*10**-12, np.unwrap(self.tds_fft[:,4]), label = r"E")
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                ax[1].plot(self.data_fft[index][:,0]*10**-12, np.unwrap(self.data_fft[index][:,4]), 
                          label = str(round(self.y[index]*10**12,1))  + " ps")
            
            
        ax[1].set_ylabel("Phase of Spectrum")
        ax[1].set_xlabel("Frequency (THz)")
        
        # if save == True:
        #     plt.savefig(f"{today}{self.name}_phase_spectra.svg", bbox_inches = "tight")
        
        plt.show()
            
    
    
    
    def build_transferfunction(self):
        """
        build transfer function for every delay time
        
        
                             FFT of differential pulse
        transfer_function = ---------------------------
                                  FFT of TDS
                                  
                                  
        see H.J. Joyce et al., Semicond. Sci. Technol. 31(10), 103003 (2016)                          
                                  
        """
        
        self.data_transfer = [None] * len(self.y)
        spectrum_tds = self.tds_fft[:,2] + 1j * self.tds_fft[:,3]
        for i in range(len(self.y)):
            spectrum = self.data_fft[i][:,2] + 1j * self.data_fft[i][:,3] # real part + imaginary part
            #spectrum = self.data_fft[i][:, 2] * np.exp(1j * phase) # wrapped und unwrapped reingeben
            # paper von peter uhd jepsen phsae retriveal in thz tds measurements how to tutorial
            # transfer = spectrum / (spectrum_tds+spectrum)
            transfer = spectrum / spectrum_tds
            # transfer = (spectrum + spectrum_tds) / spectrum_tds # from Spies 2020 Perspective
            self.data_transfer[i] = transfer
            
        self.data_transfer = np.array(self.data_transfer)
        
    
      
    
    
    # def build_transferfunction_using_polar(self):
    #     """
    #     build transfer function for every delay time
        
        
    #                          FFT of differential pulse
    #     transfer_function = ---------------------------
    #                               FFT of TDS
                                  
                                  
    #     see H.J. Joyce et al., Semicond. Sci. Technol. 31(10), 103003 (2016)                          
                                  
    #     """
        
    #     self.data_transfer = [None] * len(self.y)
    #     # spectrum_tds = self.tds_fft[:,2] + 1j * self.tds_fft[:,3]
    #     spectrum_tds = self.tds_fft[:,1] * np.exp(1j * self.tds_fft[:,4])
    #     for i in range(len(self.y)):
    #         # spectrum = self.data_fft[i][:,2] + 1j * self.data_fft[i][:,3] # real part + imaginary part
    #         spectrum = self.data_fft[i][:, 1] * np.exp(1j * self.data_fft[i][:,4]) # wrapped und unwrapped reingeben
    #         # paper von peter uhd jepsen phsae retriveal in thz tds measurements how to tutorial
    #         transfer = spectrum / spectrum_tds
    #         
    #         self.data_transfer[i] = transfer
    #     self.data_transfer = np.array(self.data_transfer)
        
    
            
        
    
    def calc_conductivity(self, substrate="teflon"):
        """
        calculate conductivity from transfer function for every delay time
        sheet conductivity
        
                
        sigma = epsilon_0 * c * (1 + n_substrate) * transfer_function
        
        
        see H.J. Joyce et al., Semicond. Sci. Technol. 31(10), 103003 (2016)
        
        
        n_substrate only for teflon for now
        
        to do (eventually):
            calculate conductivity if sheet thickness is known
        """
        if substrate == "teflon":
            n = 1.44 + 1j * 0
        
        elif substrate == "silicon":
            n = 3.418 + 1j * 0
        
        self.build_transferfunction()
        self.data_conductivity = -sc.epsilon_0 * sc.speed_of_light * (1 + n) * self.data_transfer # approximation if DeltaT/T << 1
        # self.data_conductivity = -sc.epsilon_0 * sc.speed_of_light * (1 + n_teflon) * self.data_transfer * (1/(1+self.data_transfer))
        # self.data_conductivity = sc.epsilon_0 * sc.speed_of_light * (1 + n_teflon) * ((1/self.data_transfer)-1) # from spies 2020 perspective
   
    def plot_conductivity(self, optical_delay_times, xlim = [], ylim_r = [], ylim_i = [], save = False):
        """
        plot the calculated real and imaginary conductivity for specific optical delay time

        input optical delay times list in units of ps
        """
        
        
        # # real part
        # plt.figure()
        # if xlim != []:
        #     x_indices = []
        #     for x_value in xlim:
        #         idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
        #         x_indices.append(idx)
        #     for optical_delay_time in optical_delay_times:
        #         index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
        #         plt.plot(self.tds_fft[x_indices[0]:x_indices[1],0]*10**-12, 
        #                  np.real(self.data_conductivity[index, x_indices[0]:x_indices[1]]), 'o',
        #                  label = str(round(self.y[index]*10**12,1)) + " ps")
        
        
        
        # else:
        #     for optical_delay_time in optical_delay_times:
        #             index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
        #             plt.plot(self.tds_fft[:,0]*10**-12, np.real(self.data_conductivity[index, :]), 'o',
        #                      label = str(round(self.y[index]*10**12,1)) + " ps")
        
    
        
    
        # if ylim_r != []:
        #     plt.ylim(ylim_r[0], ylim_r[1])
            
        # plt.ylabel("Real part of conducitivty")
        # plt.xlabel("Frequency (THz)")
        # plt.grid(True)
        # plt.legend()
        # plt.title(self.name)
        # if save == True:
        #     plt.savefig(F"{today}_{self.name}_thz_con_real.svg", bbox_inches = "tight")
        # plt.show()
        
        
        # # imaginary part
        # plt.figure()
        # if xlim != []:
        #     x_indices = []
        #     for x_value in xlim:
        #         idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
        #         x_indices.append(idx)
        #     for optical_delay_time in optical_delay_times:
        #         index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
        #         plt.plot(self.tds_fft[x_indices[0]:x_indices[1],0]*10**-12, 
        #                  np.imag(self.data_conductivity[index, x_indices[0]:x_indices[1]]),  'o',
        #                  label = str(round(self.y[index]*10**12,1)) + " ps")
        
        
        
        # else:
        #     for optical_delay_time in optical_delay_times:
        #             index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
        #             plt.plot(self.tds_fft[:,0]*10**-12, np.imag(self.data_conductivity[index, :]), 'o',
        #                      label = str(round(self.y[index]*10**12,1)) + " ps")
        
    
        # if ylim_i != []:
        #     plt.ylim(ylim_i[0], ylim_i[1])
    
            
        # plt.ylabel("Imaginary part of conducitivty")
        # plt.xlabel("Frequency (THz)")
        # plt.grid(True)
        # plt.legend()
        # plt.title(self.name)
        # if save == True:
        #     plt.savefig(F"{today}_{self.name}_thz_con_imag.svg", bbox_inches = "tight")
        # plt.show()
       
    
    
    
    
        # plot both in one figure
        
        plt.figure()
        if xlim != []:
            x_indices = []
            for x_value in xlim:
                idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
                x_indices.append(idx)
            for optical_delay_time in optical_delay_times:
                index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                plt.plot(self.tds_fft[x_indices[0]:x_indices[1],0]*10**-12, 
                         np.real(self.data_conductivity[index, x_indices[0]:x_indices[1]]), 'o',
                         label = str(round(self.y[index]*10**12,1)) + " ps, real")
                plt.plot(self.tds_fft[x_indices[0]:x_indices[1],0]*10**-12, 
                         np.imag(self.data_conductivity[index, x_indices[0]:x_indices[1]]),  'o',
                         label = str(round(self.y[index]*10**12,1)) + " ps, imag")
    
        else:
            for optical_delay_time in optical_delay_times:
                    index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                    plt.plot(self.tds_fft[:,0]*10**-12, np.real(self.data_conductivity[index, :]), 'o',
                             label = str(round(self.y[index]*10**12,1)) + " ps, real")        
            for optical_delay_time in optical_delay_times:
                    index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
                    plt.plot(self.tds_fft[:,0]*10**-12, np.imag(self.data_conductivity[index, :]), 'o',
                             label = str(round(self.y[index]*10**12,1)) + " ps, imag")
                    
        if ylim_r != []:
            plt.ylim(ylim_r[0], ylim_r[1])
            
        plt.ylabel("Conducitivty")
        plt.xlabel("Frequency (THz)")
        plt.grid(True)
        plt.legend()
        plt.title(self.name)
        if save == True:
            plt.savefig(F"{today}_{self.name}_thz_con_real.svg", bbox_inches = "tight")
        plt.show()            
                    
    def fit(self, optical_delay_time, xlim, guess, model, plot = True):
        
        """
        
        not a hundred percent sure how the fitting works here.. taken from michaels script
        save fit results to self.variable and print them to screen
        introduce boundaries for fit (c from -1 to 0)
        
        in a perfect world, one should get identical results for the fit for the real and imag part?
        make this hole thing applicable to other models (drude, lorentz, generalized drude)
        
        
        the fit heavily depends on your guesses
        
        would one not have to force the two fits to have the same results? 
        --> "simultaneously fitting re and im part of conductivity" ?
        global fitting has to implemented, if I am right?
        
        
        clean up here.. there are currently three fit routines built in. which one is the correct one? neither?
        """        
        
        
        x_indices = []
        for x_value in xlim:
            idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
            x_indices.append(idx)
        
        index = next(i for i,v in enumerate(self.y) if v > optical_delay_time*10**-12)
        
        data = self.data_conductivity[index,x_indices[0]:x_indices[1]]
        omega = self.tds_fft[x_indices[0]:x_indices[1],0]
        
        
        
        if model == drude_smith: # define bounds for c_1 in case of drude_smith fit
            boundz = ((-np.inf, np.inf), (-np.inf, np.inf), (-1, 0))
        elif model == comb_model:
            boundz = ((-np.inf, np.inf), (-np.inf, np.inf), (-1, 0), (-np.inf, np.inf), (-np.inf, np.inf), (-np.inf, np.inf))
        else:
            boundz = None
     
       
        
        
        # least-square method
        # residual is difference of real part of model and real part of data
        def residuals_global(a, omega, data):
            bla1 = ((np.real(model(omega, a))-np.real(data))**2).cumsum()[-1]
            bla2 = ((np.imag(model(omega, a))-np.imag(data))**2).cumsum()[-1]
            res  = bla1 + bla2 
            return res
        
       
        
        
        fit_global = minimize(residuals_global, 
                              x0 = guess, 
                              bounds = boundz, 
                              args = (omega, data), 
                              method = "Nelder-Mead", 
                              tol = 1E-10,
                              options = {"disp": True, "maxiter": 1E5})
        
        fit_values = fit_global.x
        
        # options={'xtol': 0.00001, 'ftol': 0.00000001}
        # options to change the tolerance for convergence
        # .x since minimize gives an object, self.x of this object is the fit parameters
        # what else self.variable are in there, that might be useful? save them
        
        print(fit_values)
        
      
        
        if plot == True:
            plt.figure()
            plt.title(f"{self.name}, {model.__name__}, optical pump delay: " + str(round(self.y[index]*10**12,1)) + " ps" )
            plt.plot(omega*10**-12, np.real(data), 'o', label = r'$\sigma_{real}$', color = "tab:blue")
            plt.plot(omega*10**-12, np.real(model(omega, fit_values)), color = "tab:blue")
            plt.plot(omega*10**-12, np.imag(data), 'o', label = r'$\sigma_{imag}$', color = "tab:orange")
            plt.plot(omega*10**-12, np.imag(model(omega, fit_values)), color = "tab:orange")
            plt.ylabel("Conductivity")
            plt.xlabel("Frequency (THz)")
            plt.legend()
            plt.grid()
            plt.show()
        
        
        # calculate residual sum of squares for fit
        
        bla1 = ((np.real(model(omega, fit_values))-np.real(data))**2).cumsum()[-1]
        bla2 = ((np.imag(model(omega, fit_values))-np.imag(data))**2).cumsum()[-1]
        rss  = bla1 + bla2 
        
        
        
        # find variance of fit values
        # explanation see https://stackoverflow.com/questions/42388139/how-to-compute-standard-deviation-errors-with-scipy-optimize-least-squares
        
        # no jac for nelder-mead method
        
        
        # fit_j = fit_global.jac
        # cov = np.linalg.inv(fit_j.T.dot(fit_j))
        # var = np.sqrt(np.diagonal(cov))
        
        
        
        return fit_values,  rss, fit_global
        
    
    def cond_dynamics_at_fixed_thz(self, thz_delay):
        return
    
    
    
    def calc_cond_dynamics(self, xlim=[0.25, 1.25]):
        """
        calculate the integrated (real and imaginary) conductivity for every pump-probe delay 
        in the THz region defined by xlim
        
        does not need to be explicitely called, since it is part of self.plot_cond_dynamics
        """
        x_indices = []
        for x_value in xlim:
            idx = next(i for i,v in enumerate(self.tds_fft[:,0]) if v > x_value * 10**12)
            x_indices.append(idx)
        
        
        
        data = self.data_conductivity[:,x_indices[0]:x_indices[1]]
        data_re = np.real(data)
        data_im = np.imag(data)
        
        
        summed_cond_re = np.zeros(len(self.y))
        for i, v in enumerate(data_re):
            summed_cond_re[i] = np.sum(v)
            
        summed_cond_im = np.zeros(len(self.y))
        for i, v in enumerate(data_im):
            summed_cond_im[i] = np.sum(v)
            
        return summed_cond_re, summed_cond_im

    
        
    def plot_cond_dynamics(self, xlim=[0.25, 1.25], cond="real", normalize=False):
        """
        plot dynamics from the integrated conductivity within the frequency range given by xlim
        choose to plot real, imag or both 
        """

        
        real = self.calc_cond_dynamics(xlim)[0]
        imag = self.calc_cond_dynamics(xlim)[1]
        
        if cond == "real":
            
            plt.figure()
            plt.title(f"integrated real conductivity, {self.name}")
            plt.plot(self.y*10**12, real, 'o')
            plt.xlabel("Delay Time (ps)")
            plt.ylabel(r"$\int\sigma_{1} \enspace (\Omega^{-1})$")
            plt.grid(True)
            plt.show()
            
        
        elif cond == "imag":
            
            plt.figure()
            plt.title(f"integrated imaginary conductivity, {self.name}")
            plt.plot(self.y*10**12, imag, 'o')
            plt.xlabel("Delay Time (ps)")
            plt.ylabel(r"$\int\sigma_{2} \enspace (\Omega^{-1})$")
            plt.grid(True)
            plt.show()
            
        
        elif cond == "both" and normalize == False:
            plt.figure()
            plt.title(f"integrated conductivities, {self.name}")
            plt.plot(self.y*10**12, real, 'o', label = r"\sigma_{1}")
            plt.plot(self.y*10**12, imag, 'o', label = r"\sigma_{2}")
            plt.xlabel("Delay Time (ps)")
            plt.ylabel(r"$\int\sigma \enspace (\Omega^{-1})$")
            plt.grid(True)
            plt.show()
        elif cond == "both" and normalize == True:
            plt.figure()
            plt.title(f"integrated conductivities, {self.name}")
            plt.plot(self.y*10**12, np.abs(real)/np.max(np.abs(real)), 'o', label = r"\sigma_{1}")
            plt.plot(self.y*10**12, np.abs(imag)/np.max(np.abs(imag)), 'o', label = r"\sigma_{2}")
            plt.xlabel("Delay Time (ps)")
            plt.ylabel(r"$\int\sigma \enspace (\Omega^{-1})$")
            plt.grid(True)
            plt.show()
    
        
    def fit_cond_dynamics(self, xlim):
        """
        Fit exponential decay to integrated conductivity
        """
        
        
        real = self.calc_cond_dynamics(xlim)[0]
        idx_max = np.argmax(real)
        real_cut = real[idx_max:]
        y = self.y[idx_max:]
        
        def exp_decay(t,A, d):
            return A * np.exp(-t/d)
        
        guess = [0.1, 20*10**-12]
        popt, pcov = curve_fit(exp_decay, y, real_cut, p0 = guess)
        
        
        plt.figure()
        plt.plot(self.y*10**12, real, 'o', color = "tab:blue",  alpha = 0.5, label = str(round(popt[1]*10**12, 0)) + " ps")
        plt.plot(y*10**12, exp_decay(y, *popt), color = "tab:blue")
        plt.xlabel("Delay Time (ps)")
        plt.ylabel(r"$\int\sigma_{1} \enspace (\Omega^{-1})$")
        plt.legend()
        plt.grid(True)
        plt.show()
        
        
        return popt
        
    

        
        
        
        
        
        
        
        
        
        
        
        