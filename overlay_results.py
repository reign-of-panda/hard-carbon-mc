"""
This code plots the csv results from the hard carbon monte carlo simulations.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.spines as sp
import numpy as np
import os
from scipy import integrate

"""
File directories
"""

HEC1_dir_sites = "/HEC1_27th_June/vary_sites" # Add "400.csv" etc to access CSV file 
file_sites = "/Na_monte_carlo_results_sites"
arr_sites = ['400.csv',  '1200.csv',  '5000.csv', '10000.csv'] # '800.csv', '1600.csv',

HEC1_dir_sps = "/HEC1_27th_June/vary_sps"
file_sps = "/Na_monte_carlo_results_sps"
arr_sps = ["400.csv",  "800.csv",  "2000.csv"] #"600.csv", "1000.csv",

HEC1_dir_sps2 = "/HEC1_27th_June/vary_sps2"
file_sps2 = "/Na_monte_carlo_results_sps2_"
arr_sps2 = ["400.csv", "800.csv", "2000.csv"] #"600.csv", "800.csv", "1000.csv",

pickle_testing = "/pickle_testing/results"
file_pickle = "/Na_monte_carlo_results_sps2_"
arr_pickle = ['400.csv', '2000.csv', '5000.csv', '2000_run4.csv']

# The following three variables are used to input the correct directory for the csv
# These variables are used in the first 4 lines of code in display()
directory = HEC1_dir_sps
file = file_sps
array_element = arr_sps

# Use this to switch off certain plots (True to plot them)
pick_plot = dict([
    ('voltage', True),
    ('dS/dx', True),
    ('dQ/dV', True),
    ('dH/dx', True),
    ('S', True),
    ('d/dx(dH/dx)', True)
    ])

experimental_plot = False # Use True to plot experimental data
simulation_plot = True # Use True to plot simulation data

# Define number of plots in the image (used under create plot and formats)
num_col = 2
num_row = 3

class DisplayResults:

    def __init__(self):
        pass

    def display(self):
        """
        This is the main function to plot the results.
        
        Some changes:
        
        The experimental_data.csv file is now in current working directory (cwd)
        I have removed the results folder, and opted for a different folder organisation
        I have entered multiple files into dataframes variable
        
        
        dataframes: Data from the simulation 
        dfE: Experimental data
        """
        
        plt.rcParams['font.size'] = 14
        plt.rcParams['axes.linewidth'] = 1.2  # 1.2 for single plot, 0.5 for all 6
        plt.rcParams['lines.linewidth'] = 20.0 # Aah, this doesn't work because line width is changed later on

        cwd = os.getcwd()  # Gets current working directory.
        cwd = cwd.replace('\\', '/')
        path = cwd + directory  # This is the folder all the results are stored in.
        
        if type(array_element) == str:
            dataframes = [file + array_element] # This is to pass a single csv file
        else:
            dataframes = [file + i for i in array_element] # This is a list so you can pass multiple csv files to be overlayed on the same plot.

        colours = ['black', 'darkred', 'darkmagenta', 'darkturquoise', 'saddlebrown']  # Array of colours for the lines.

        dfE = pd.read_csv(cwd + "/experimental_data.csv")  # Reads in the experimental data as a pandas dataframe.

        # Rescale the x-axis of the experimental data.
        ratio_of_capacities = 272.4 / 338.313338  # experimental maximum capacity / theoretical maximum capacity
        dfE["x_theo"] = ratio_of_capacities * dfE["x"]
        # 'x' is the experimental x and 'x_theo' is the theoretical x.

        # Second derivative of enthalpy for experimental data. One w/ respect to the experimental x and one w/ respect to theoretical x.
        secder_enthalpy_experimental_x = np.gradient(np.array(dfE['Enthalpy dH/dx']), np.array(dfE['x']))
        secder_enthalpy_experimental_x_theo = np.gradient(np.array(dfE['Enthalpy dH/dx']), np.array(dfE['x_theo']))
        dfE['secder enthalpy x'] = secder_enthalpy_experimental_x
        dfE['secder enthalpy x theo'] = secder_enthalpy_experimental_x_theo

        # vertical shift on p.m. entropy for vibrational effect
        vibrational_shift = 0.0108  # eV K  this includes being multiplied by the ratio of capacities.
        dfE["Entropy dS/dx"] = (dfE["Entropy dS/dx"]) - vibrational_shift

        # Integrates the p.m. entropy
        entropy_list_experimental = integrate.cumtrapz(dfE['Entropy dS/dx'], dfE['x'],
                                                       initial=0)  # Contains the entropy values
        dfE['Entropy'] = entropy_list_experimental

        dfE['x_new'] = ((dfE['x_theo'] - dfE['x_theo'].iloc[0]) * dfE['x_theo'][73]) / (dfE['x_theo'][73] - dfE['x_theo'].iloc[0])  # Rescales the line so that the experimental data starts at 0.
        dfE['x'] = ((dfE['x'] - dfE['x'].iloc[0]) * dfE['x'][73]) / (dfE['x'][73] - dfE['x'].iloc[0])  # Same as above but for experimental x axis.

        # Calculates the analytical solution
        points = 1000
        x_pos = np.linspace(0, 1, points)  # x for p.m. entropy
        y_pos = np.linspace(0, 1, points)  # y for p.m. etropy
        s_x = np.linspace(0, 1, points)  # x for entropy
        s_y = np.linspace(0, 1, points)  # y for entropy
        l = 0.329217689  # This must be the same as what was used in the main script
        R = -0.0000862  # eV/K.Site
        T = 288  # K
        for index, x in enumerate(x_pos):
            if x < l:
                s_y[index] = (R * (x * np.log(x / l) - (x - l) * np.log((l - x) / l))) * T
                y_pos[index] = T * R * (np.log(x / l) - np.log((l - x) / l))
            else:
                s_y[index] = (R * l * (
                            (x / l - 1) * np.log(x / l - 1) + (1 - x) / l * np.log((1 - x) / l) - (1 - l) / l * np.log(
                        (1 - l) / l))) * T
                y_pos[index] = T * R * (np.log(x / l - 1) - np.log(1 / l - x / l))

        #  Calculates the single solid state entropy
        x_ent = np.linspace(0, 1, points)
        y_ent = np.linspace(0, 1, points)
        for index, x in enumerate(x_ent):
            y_ent[index] = T * R * (x * np.log(x) + (1-x) * np.log(1-x))
        
        """
        #
        #
        # Create plot and formats
        #
        #
        """
        
        fig, axes = plt.subplots(nrows=num_row, ncols=num_col, constrained_layout=True, squeeze=False)
        # squeeze=False is needed to prevent errors when plotting a single subplot
        plt.rc('legend', fontsize=13, handlelength=1)
        plt.rc('tick')
        lw = 1.5  # Line width
        
        plt.tick_params(bottom=True, top=True, left=True, right=True)
        plt.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
        plt.tick_params(direction='in', width=1.2, length=4.5, pad=3) # For single plot
        # plt.tick_params(direction='in', width=1, length=4.5, pad=3) # For multiple plots

        marker_list = ['v', '^', 'p', 'o']
        mark_size = 3 #0.7 for 6 plots
        
        colours = ['#176ba0', '#af4bce', 'orangered', '#48a11b', '#3caea3'] #'#af4bce'
        common_legend = ['400 Averaging Steps', '800 Averaging Steps', '2000 Averaging Steps']
        
        if num_col==2 and num_row==3: # This will work when using the original axes dimensions (3 rows, 2 columns)
            placement = dict([
                ('voltage', axes[0, 0]),
                ('dS/dx', axes[0, 1]),
                ('dQ/dV', axes[1, 0]),
                ('dH/dx', axes[1, 1]),
                ('S', axes[2, 0]),
                ('d/dx(dH/dx)', axes[2, 1])
                ])
        else: # If axes dimensions are different, I'm probably trying to plot one graph
            """
            If plotting more than one graph, the position on the plot in the subplot can be adjusted
            by appropriately altering the axes[] parameter. For the graphs that are not being plotted, 
            leave their position as axes[0, 0].
            """
            placement = dict([
                ('voltage', axes[0, 0]),
                ('dS/dx', axes[0, 0]),
                ('dQ/dV', axes[0, 0]),
                ('dH/dx', axes[0, 0]),
                ('S', axes[0, 0]),
                ('d/dx(dH/dx)', axes[0, 0])
                ])
        
        # Plots all of the experimental data
        if experimental_plot == True:
            if pick_plot['voltage'] == True:
                dfE.plot(linestyle='-', color='darkgreen', lw=lw, ax=placement['voltage'], x='x_new', y='OCV')
                dfE.plot(linestyle='-', color='darkblue', lw=lw, ax=placement['voltage'], x='x', y='OCV')
    
            if pick_plot['dS/dx'] == True:
                ax2 = dfE.plot(linestyle='-', color='darkgreen', lw=lw, ax=placement['dS/dx'], x='x_new', y='Entropy dS/dx')
                dfE.plot(linestyle='-', color='darkblue', lw=lw, ax=placement['dS/dx'], x='x', y='Entropy dS/dx')
            
            if pick_plot['dQ/dV'] == True:
                dfE.plot(linestyle='-', color='darkgreen', lw=lw, ax=placement['dQ/dV'], x='OCV', y='dQdV') 
    
            if pick_plot['dH/dx'] == True:
                dfE.plot(linestyle='-', color='darkgreen', lw=lw, ax=placement['dH/dx'], x='x_new', y='Enthalpy dH/dx')
                dfE.plot(linestyle='-', color='darkblue', lw=lw, ax=placement['dH/dx'], x='x', y='Enthalpy dH/dx')
            
            if pick_plot['S'] == True:
                ax5 = dfE.plot(linestyle='-', color='darkgreen', lw=lw, ax=placement['S'], x='x_new', y='Entropy')
    
            if pick_plot['d/dx(dH/dx)'] == True:
                dfE.plot(linestyle='-', color='darkgreen', lw=lw, ax=placement['d/dx(dH/dx)'], x='x_new', y='secder enthalpy x theo')
                dfE.plot(linestyle='-', color='darkblue', lw=lw, ax=placement['d/dx(dH/dx)'], x='x', y='secder enthalpy x')

        # Iterate through all the data to be plotted
        if simulation_plot == True:
            for count, df in enumerate(dataframes):
                df1 = pd.read_csv(path + df)  # reads file into a dataframe.
    
                df1 = df1.replace(0, np.nan).dropna(axis=0, how='all')  # For the rows with all '0' entries they are replaced with 'nan' and then these rows are dropped.
                df1 = df1.replace(np.nan, 0)  # As some legitimate 0 entries such as 0 volts we flip back the remaining from 'nan' to 0.
    
                # Integrates the p.m. entropy
                entropy_list = integrate.cumtrapz(df1['Partial molar entropy'], df1['Total mole fraction'],
                                                  initial=0)  # Contains the entropy values
                df1['Entropy'] = entropy_list
    
                # Rescale voltage profile and p.m. enthalpy by the chain rule.
                df1["adjusted voltage"] = df1["Chemical potential"] * ratio_of_capacities
                df1["adjusted enthalpy"] = df1["Partial molar enthalpy"] * ratio_of_capacities
                df1["adjusted entropy"] = df1["Partial molar entropy"] * ratio_of_capacities
                df1["adjusted dq/de"] = df1["dq/de"] * (1/ratio_of_capacities)**2
    
                # Differentiate the p.m. enthalpy to get the second derivative.
                pm_enthalpy = np.array(df1['adjusted enthalpy'])
                mole_fraction = np.array(df1['Total mole fraction'])
                secder_enthalpy = np.gradient(pm_enthalpy, mole_fraction)
                df1['secder enthalpy'] = secder_enthalpy
    
                if pick_plot['voltage'] == True:
                    ax1 = df1.plot(linestyle='-', color=colours[count], lw=lw, marker=marker_list[count], markeredgecolor=colours[count],
                                    markersize=mark_size, ax=placement['voltage'], x='Total mole fraction', y='adjusted voltage')
                    ax1.set_xlim([0, 1])
                    ax1.set_xlabel('Na content $[x]$')
                    ax1.set_ylabel('Voltage $[V]$')
                    ax1.legend(common_legend)                    
                    # ax1.legend(['Experimental data (Adjusted x)', 'Raw experimental data', 'Monte Carlo data'])
                    
                if pick_plot['dS/dx'] == True:
                    ax2 = df1.plot(linestyle='-', color=colours[count], lw=lw, marker=marker_list[count], markeredgecolor=colours[count],
                                    markersize=mark_size, ax=placement['dS/dx'], x='Total mole fraction', y='adjusted entropy')
                    # ax2.plot(x_pos, y_pos, linewidth=lw, color='red')  # Plots the ideal p.m. entropy
                    ax2.set_xlim([0, 1])
                    ax2.set_xlabel('Na content $[x]$')
                    ax2.set_ylabel('$\\frac{dS}{dx}$ $[eV K/site]$')
                    ax2.legend(common_legend)                    
                    # ax2.legend(['Experimental data (Adjusted x)', 'Raw experimental data', 'Monte Carlo data',  'Analytical solution'])
                    
                if pick_plot['dQ/dV'] == True:
                    ax3 = df1.plot(linestyle='-', color=colours[count], lw=lw, marker=marker_list[count], markeredgecolor=colours[count],
                                    markersize=mark_size, ax=placement['dQ/dV'], x='Chemical potential', y='adjusted dq/de') 
                    ax3.set_xlim([-0.1, 1])
                    ax3.set_xlabel('Voltage $[V]$')
                    ax3.set_ylabel('$\\frac{dQ}{dV}$ [$\mathregular{eV^{-1}}$]')
                    ax3.legend(common_legend)
                    # ax3.legend(['Experimental data', 'Monte Carlo Data'])
                    
                if pick_plot['dH/dx'] == True:
                    ax4 = df1.plot(linestyle='-', color=colours[count], lw=lw, marker=marker_list[count], markeredgecolor=colours[count],
                                    markersize=mark_size, ax=placement['dH/dx'], x='Total mole fraction', y='adjusted enthalpy')
                    ax4.set_xlim([0, 1])
                    ax4.set_xlabel('Na content $[x]$')
                    ax4.set_ylabel('$\\frac{dH}{dx}$ $[eV/site]$')
                    ax4.legend(common_legend)                    
                    # ax4.legend(['Experimental data (Adjusted x)', 'Raw experimental data', 'Monte Carlo data'])
                    
                if pick_plot['d/dx(dH/dx)'] == True:
                    ax5 = df1.plot(linestyle='-', color=colours[count], lw=lw, marker=marker_list[count], markeredgecolor=colours[count],
                                    markersize=mark_size, ax=placement['d/dx(dH/dx)'], x='Total mole fraction', y='secder enthalpy')
                    ax5.set_xlim([0, 1])
                    ax5.set_ylim([0, 6])
                    ax5.set_xlabel('Na content $[x]$')
                    ax5.set_ylabel('$\\frac{d^2H}{dx^2}$ $[eV/site]$')
                    ax5.legend(common_legend)
                    
                    # ax5.legend(['Experimental data (Adjusted x)', 'Raw experimental data', 'Monte Carlo data'])
                    
                if pick_plot['S'] == True:
                    ax6 = df1.plot(linestyle='-', color=colours[count], lw=lw, marker=marker_list[count], markeredgecolor=colours[count],
                                    markersize=mark_size, ax=placement['S'], x='Total mole fraction', y='Entropy')
        
                    # ax6.plot(s_x, s_y, linewidth=lw, color='red')  # Plots the entropy for l=0.32...
                    # ax6.plot(x_ent, y_ent, linewidth=lw, color='grey')  # Plots the entropy for solid state solution.
                    ax6.set_xlim([0, 1])
                    ax6.set_xlabel('Na content $[x]$')
                    ax6.set_ylabel('S $[eV K/site]$')
                    ax6.legend(common_legend)
                    # ax6.legend(['Experimental data', 'Monte Carlo data', 'Analytical solution', 'Solid state solution'], loc='upper right', bbox_to_anchor=(0.75, 0.5))
            
                

        # parameter_file = open(path + "/Input_arguments_" + uid + ".txt", "w")
        # parameter_file.write(str(self.args))
        # parameter_file.close()

        # manager = plt.get_current_fig_manager()
        # # manager.resize(*manager.window.maxsize())
        # # fig_path = cwd + "/Na_plot_results.png"
        # # plt.savefig(path + "/Na_monte_carlo_plot_" + uid + ".png")
        # plt.show()
        
        plt.savefig("Varying sps Overlaid Plots - dQ_dV", dpi = 300)

        plt.show()

    def display_averaging(self):
        """
        This can be used to display the averages by passing through the two average files for U and N. To get these files run the code with monitor_averaging = True.
        """

        cwd = os.getcwd()
        path = cwd + "/results"
        df1 = pd.read_csv(path + "/average_U.csv")  # black line
        df2 = pd.read_csv(path + "/average_N.csv")  # green line
        chem = 25  # from 0 to 35

        s1 = df1.iloc[chem]
        s1.plot()

        plt.show()

    def test_uniform(self):
        """
        Shows the pdf for a uniform distribution
        """

        s = np.random.uniform(-1.35, 0.5, 5000)
        plt.hist(s, 30, density=False)
        plt.xlabel('Interlayer point energy [eV]')
        plt.ylabel('Frequency')
        plt.show()

    def test_normal(self):
        """
        Shows the pdf for a normal distribution
        """
        s = np.random.normal(-0.42, 0.55, 5000)
        plt.hist(s, 30, density=False)
        plt.xlabel('Interlayer point energy [eV]')
        plt.ylabel('Frequency')
        plt.show()

    def test_triangular(self):
        """
        Shows the pdf for a triangular distribution
        """
        s = np.random.triangular(-1.65, 0.08, 0.08, 5000)
        plt.hist(s, bins=30, density=False)
        plt.xlabel('Interlayer point energy [eV]')
        plt.ylabel('Frequency')
        plt.show()

    def test_power(self):
        """
        Shows the pdf for a power distribution
        """
        a = 6  # shape
        samples = 5000
        max = -0.06
        min = -3.3
        s = np.random.power(a, samples) * -1 * (min - max) + min
        plt.hist(s, bins=30, density=False)
        plt.xlabel('Interlayer point energy [eV]')
        plt.ylabel('Frequency')
        plt.show()


if __name__ == '__main__':
    dr = DisplayResults()
    dr.display()
    











