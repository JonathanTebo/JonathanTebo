import csv
import numpy as np
from numpy.polynomial import polynomial as P
from scipy.integrate import simps, trapz
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import readline, glob
import os
import traceback
import sys

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]


def RtoT(res):
    expo=[]
    n_terms=0
    room_res = 77350
    quo = res/room_res
    ln = np.log(quo)
    counter = 3
    avars = [.2700714,.07323249,1.0219809,-1.9558035,3.4889894,-4.4175955,3.676856,-1.7871035,.3821843]
    first_terms = (avars[0]) + (avars[1]*(ln**0.25)) + (avars[2]*(ln))
    #a3 to a8
    for i in range(2,8):
        expo.append(ln ** i)
    #good up to here

    for i in range(len(expo)):
        n_terms += (expo[i] * avars[counter])
        if counter <=7:
            counter += 1
        else:
            break
    right_side = first_terms + n_terms
    thesis = 1/(right_side ** 4)
    return thesis


def FindEdge(A,B):
    D = []
    E = []
    F = []
    start = 0
    edge = 0
    for n in range(1,len(A)):
        if (A[n] - A[n-1]) > 0.01:
            edge = n
            break
    T_max=np.amax(B[:edge-1])
    T_min=np.amin(B[:edge-1])
    u_T=abs(T_max-T_min)/2   
    T_i = T_min + u_T
    if edge < 10:
        start = 1
    else:
        start = edge - 10
    for n in range(start,len(A)):
        D.append(A[n])
        E.append((B[n]-T_i))
        F.append((n-start)/10.000)
    return D, E, F, T_i, u_T


def FindPulse(A,B):
    i = 0
    f = 0
    inst_p = []
    r = 160
    edge = 0
    for n in range(1,len(A)):
        if (A[n] - A[n-1]) > 0.001:
            edge = n
            i = B[n]
            inst_p.append(A[n-1]*A[n-1]/r)
            break
    for n in range(edge,len(A)):
        inst_p.append(A[n]*A[n]/r)
        if A[n] < 0.001:
            f = B[n]
            break
    nrg = simps(inst_p,dx=0.1)
    nrg2 = trapz(inst_p,dx=0.1)
    u_nrg = abs(nrg-nrg2)
    return i, f, nrg, u_nrg


def func(x,a,b,c):
    return a * np.exp(-b * x) + c


def func2(x,a,b):
    return a*x**b


def CalcDeltaT(A,B,T):
    #Initialize Needed Objects
    e_extended = []
    emax = []
    emin = []
    lower_integral = []
    upper_integral= []
    umax_int = []
    umin_int = []
    upper_integral_data = []
    umax_int_data = []
    umin_int_data = []
    integral_difference = []
    int_max_diff = []
    int_min_diff = []
    start = A.index(np.amax(A))
    pulse_start_index = B.index(T)

    #Fill Arrays and Curve Fit
    trimmed_data = np.array(A[start+1:])
    trimmed_time = np.array(B[start+1:])
    popt,pcov = curve_fit(func,trimmed_time,trimmed_data,[0.5,0.05,0.05],maxfev=1500)
    perr = np.sqrt(np.diag(pcov))
    for i in range(len(B)):
        e_extended.append(func(B[i],popt[0],popt[1],popt[2]))
        emax.append(func(B[i],popt[0]-perr[0],popt[1]-perr[1],popt[2]-perr[2]))
        emin.append(func(B[i],popt[0]+perr[0],popt[1]+perr[1],popt[2]+perr[2]))        

    #Prep Arrays for Integration
    integral_time = B[pulse_start_index:]
    integral_time2 = B[pulse_start_index:]
    lower_integral_data = A[pulse_start_index:]
    upper_limit = e_extended[pulse_start_index:]
    umax = emax[pulse_start_index:]
    umin = emin[pulse_start_index:]

    #Do integration? - coded by Alex Barajas
    for i in range(len(upper_limit)):
        upper_integral_data.append(upper_limit[i]-lower_integral_data[i])
    integral_time2.reverse()
    for i in range(1,len(integral_time)):
        lower_integral.append(simps(lower_integral_data[:i+1],integral_time[:i+1]))
        upper_integral.append(abs(simps(upper_integral_data[:i+1], integral_time2[:i+1])))
    for i in range(len(upper_integral)):
        integral_difference.append(abs(upper_integral[i]-lower_integral[i]))
    eq_area_index = integral_difference.index(np.amin(integral_difference))

    #Find uncertainty in calculation
    #Take popt and add/subtract perr to get absolute differences
    #Evaluate as above for each case
    #Take the difference / 2 to get uncertainty in delta T
    for i in range(len(umax)):
        umax_int_data.append(umax[i]-lower_integral_data[i])
    integral_time2.reverse()
    for i in range(1,len(integral_time)):
        umax_int.append(abs(simps(umax_int_data[:i+1],integral_time2[:i+1])))
    for i in range(len(umax_int)):
        int_max_diff.append(abs(umax_int[i]-lower_integral[i]))
    max_eq_area_index = int_max_diff.index(np.amin(int_max_diff))
    max_Delta_T = umax[max_eq_area_index]

    for i in range(len(umin)):
        umin_int_data.append(umin[i]-lower_integral_data[i])
    integral_time2.reverse()
    for i in range(1,len(integral_time)):
        umin_int.append(abs(simps(umin_int_data[:i+1],integral_time2[:i+1])))
    for i in range(len(umin_int)):
        int_min_diff.append(abs(umin_int[i]-lower_integral[i]))
    min_eq_area_index = int_min_diff.index(np.amin(int_min_diff))
    min_Delta_T = umin[min_eq_area_index]

    u_Delta_T = abs(min_Delta_T - max_Delta_T) / 2

    return upper_limit[eq_area_index], u_Delta_T


def CalcUncC(n,u_n,t,u_t):
    unc_C = np.sqrt((u_n/t)**2 + (n*u_t/(t*t))**2)
    return unc_C


def format_exponent(ax, axis='y'):

    # Change the ticklabel format to scientific format
    ax.ticklabel_format(axis=axis, style='sci', scilimits=(-2, 2))

    # Get the appropriate axis
    if axis == 'y':
        ax_axis = ax.yaxis
        x_pos = 0.0
        y_pos = 1.0
        horizontalalignment='left'
        verticalalignment='bottom'
    else:
        ax_axis = ax.xaxis
        x_pos = 1.0
        y_pos = -0.05
        horizontalalignment='right'
        verticalalignment='top'

    plt.tight_layout()

    # Get the offset value
    offset = ax_axis.get_offset_text().get_text()

    if len(offset) > 0:
        # Get that exponent value and change it into latex format
        minus_sign = u'\u2212'
        expo = np.float(offset.replace(minus_sign, '-').split('e')[-1])
        offset_text = r'x$\mathregular{10^{%d}}$' %expo

        # Turn off the offset text that's calculated automatically
        ax_axis.offsetText.set_visible(False)

        # Add in a text box at the top of the y axis
        ax.text(x_pos, y_pos, offset_text, transform=ax.transAxes,
               horizontalalignment=horizontalalignment,
               verticalalignment=verticalalignment)
    return ax


def main(my_file):
    VFull = []      #Voltage array
    RFull = []      #Thermometer array
    RTFull = []     #Calculated temperature array from thesis
    VPulse = []     #Voltage array of the pulse onwards
    RPulse = []     #Resistance Temperature array from pulse onwards
    Time = []       #Time array

    with open(my_file) as tsv:
        csvReader = csv.reader(tsv, dialect="excel-tab")
        for line in csvReader:
            lineNum = csvReader.line_num
            if lineNum != 1:
                V = abs(float(line[0]))
                R = float(line[1])
                time = (lineNum - 1) / 10.000
                VFull.append(V)
                RFull.append(R)

    #CREATE FILENAME FOR OUTPUT FILES IN THE FORM "SAMPLE_DATE_TEMP"#
    m_f = my_file.rsplit('.', 1)
    m_f2 = m_f[0].rsplit('\\',1)
    sample_name = m_f2[1].split('_')
    date = sample_name[1] + '_' + sample_name[2] + '_' + sample_name[3]

    # CALCULATE RuO2 TEMPERATURE FROM RESISTANCE #
    for i in RFull:
        temp1 = RtoT(i)
        RTFull.append(temp1)

    VPulse, RPulse, Time, T_i, unc_T = FindEdge(VFull, RTFull)

    # CALCULATE PULSE TIME AND POWER #
    pulse_i, pulse_f, input_energy, unc_energy = FindPulse(VPulse, Time)

    #CALCULATE DELTA TEMP#
    delta_temp, unc_dT = CalcDeltaT(RPulse,Time,pulse_i)

    # CALCULATE C(T) AND k(T) #
    C = input_energy/delta_temp
    unc_C = CalcUncC(input_energy,unc_energy,delta_temp,unc_dT)
    return date + '_' + sample_name[4] , round(T_i,3), unc_T, C, unc_C


try:
    if __name__== "__main__":
        readline.set_completer_delims(' \t\n;')
        readline.parse_and_bind("tab: complete")
        readline.set_completer(complete)
        dirname = input('Folder name: ')
        q = int(input('Type 1 for Heat Capacity, Type 0 for Molar Specific Heat: '))
        if q==0:
            m = float(input('Mass of Sample: '))
            Mm = float(input('Molar Weight of Sample (Ca2Y2Cu5O10=735.688g/mol): '))
            moles = m/Mm
        else:
            moles=1

        Name = []
        Temp = []
        uncTemp = []
        Cap = []
        uncCap = []

        for filename in os.listdir(dirname):
            print(filename)
            name, T, uT, C, uC = main(dirname + '\\' + filename)
            C2=(C-func2(T,4.716E-8,2.7687799))*1000./moles
            uC2 = uC*1000./moles
            C3=C2/T
            uC3 = C3 *np.sqrt((uC2/C2)**2 + (uT/T)**2)
            Name.append(name)
            Temp.append(T)
            uncTemp.append(uT)
            Cap.append(C3)
            uncCap.append(uC3)
            print('T=', T, '+/-', uT*100/T, 'C=', C3, '+/-', uC3*100/C3)

        fig, ax = plt.subplots()
        plt.xlabel("Temperature (K)")
        plt.xlim(5.0,40.0)
        plt.minorticks_on()
        plt.grid()
        ax.errorbar(Temp, Cap, xerr=uncTemp, yerr=uncCap,ls='none', marker='.', color='red', ecolor='black', capsize=2, label = 'Data')

        if q==0:
            fig.suptitle('Sample Molar Specific Heat')
            plt.ylabel('C/T (mJ/mol K**2)')
            plt.ylim(0.0, 1E3)
            plt.savefig('OutputFigures/cvT.png')
        else:
            fig.suptitle('Heat Capacity')
            plt.ylabel('C/T (mJ/K**2)')
            #plt.ylim(0.0, 1E3)
            ax = format_exponent(ax, axis='y')
            BX = np.array(range(40))
            ax.plot(BX,func2(BX,4.716E-8,2.7687799), color = 'green',label='Addenda')
            plt.legend(loc='upper left')
            plt.savefig('OutputFigures/CvT.png')

        plt.close()

        with open('CvT_data.txt','w') as f:
            f.write('Filename' + '\t' + 'Temperature' + '\t' + 'Uncertainty T' + '\t' 'Heat Capacity' + '\t' + 'Uncertainty Heat Cap' + '\n')
            for n in range(len(Temp)):
                f.write(Name[n] + '\t' + str(Temp[n]) + '\t' + str(uncTemp[n]) + '\t' + str(Cap[n]) + '\t' + str(uncCap[n]) + '\n')
        input("Press Enter to continue...")


except Exception as err:
    traceback.print_exc(file=sys.stdout)
    input("Press Enter to continue...")