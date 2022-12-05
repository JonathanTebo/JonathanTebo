# This program calculates the heat capacity of a sample given the data
#    output by {Main_V2.vi}. Those data are heater voltage, RuO2 ther-
#    mometer resistance, base temperature, and clock time (clock time
#    available only after June 30, 2021, and only shows the second
#    hand's reading only). This program calculates all available heat
#    capacity data available.


import csv
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import curve_fit as cf
from scipy import integrate as intg


def main():
    """Sets up program, organizes functions, and displays plot."""

    # Setup.
    import_directory, config, molar_or_total, m, M, units, system_name \
        = program_setup()

    # Calculates heat capacity for all samples and creates plot.
    calculated_data = data_analysis(import_directory, molar_or_total, units,
                                    system_name, config)

    # Calculates specific heat for pure samples and creates plots.
    if config == "Pure Sample":
        pure_sample_Cv(calculated_data, M, m, config, import_directory,
                       system_name, molar_or_total)


def data_analysis(import_directory, molar_or_total, units, system_name,
                  config):
    """Calculates heat capacity for all samples and creates plots."""

    calculated_data = []  # [[[together temp_set], [together C_set]],
    #  [[addenda  temp_set], [addenda  C_set]]]

    for directory in range(len(import_directory)):
        # Read addenda data from file if already have.
        if import_directory[directory] == "Blank Addenda":
            read_or_not = input("\nEnter 'A' to calculate addenda data "
                                "anew. Enter anything else to read data "
                                "from a file: ")
            read_or_not = read_or_not.upper()
            if read_or_not != "A":
                temp_set, C_set = read_addenda()
                calculated_data.append([temp_set, C_set])
                continue

        # Enforce rule promulgated at beginning of program_setup().
        if len(os.listdir(import_directory[directory])) == 0:
            print(f"\n\nThere are no data files in the folder "
                  f"{import_directory[directory]}.")
            temp_set = np.array([])
            C_set = np.array([])
            calculated_data.append([temp_set, C_set])
            continue

        # Calculation and its progress.
        temp_set, C_set, blank_files = calculation_and_progress(
            import_directory[directory], molar_or_total, units)

        # Identifies blank data files.
        blank_data_files(blank_files)

        # Plotting.
        plot_data(import_directory[directory], system_name, molar_or_total,
                  temp_set, C_set)

        # Save data in memory for pure sample heat capacity calculation.
        calculated_data.append([temp_set, C_set])

        # Save all data for troubleshooting.
        dirname = os.path.dirname(__file__)
        save_directory = "Storage for Testing"
        save_file = save_directory + "/" + import_directory[directory] + \
                    ".txt"
        export_file = os.path.join(dirname, save_file)
        heading = ["Temperature (K)", "Heat Capacity (J/K)"]

        with open(export_file, "w", newline="") as \
                output:
            writer = csv.writer(output, delimiter="\t")
            writer.writerow(heading)  # Creates column headers.
            for i in range(len(C_set)):  # Populates file with data.
                writer.writerow([temp_set[i], C_set[i]])

    return calculated_data


def pure_sample_Cv(calculated_data, M, m, config, import_directory,
                   system_name, molar_or_total):
    """Calculates pure sample's heat capacity (if applicable)."""

    together_temp_set = calculated_data[0][0]
    together_C_set = calculated_data[0][1]
    addenda_C_and_T_sets = calculated_data[1]

    # Subtract addenda Cv from sample + addenda Cv and plots result.
    usable_addenda_C_set = addenda_fit(addenda_C_and_T_sets,
                                       together_temp_set)
    sample_C_set = together_C_set - usable_addenda_C_set

    if molar_or_total == "M":
        # Calculates specific heat from heat capacity.
        usable_addenda_C_set *= (M / m)

    # Plots final result.
    plot_data(import_directory, system_name, molar_or_total,
              together_temp_set, sample_C_set)


def program_setup():
    """Allows user to set program settings to calculate different
    experimental configurations."""

    print("\nPlease place all addenda w/ sample temperature response data "
          "in the 'Addenda + Sample' folder\nand all blank addenda "
          "temperature response data in the 'Blank Addenda' folder.")
    input("\nPress 'Enter' to continue...")

    system_name = input("\nPlease enter the sample or system's name (use "
                        "TeX for subscripts): ")

    experiment_config = input("\nEnter 'A' to calculate the addenda's heat "
                              "capacity. Enter 'B' to calculate the "
                              "addenda-sample\nsystem's heat capacity. "
                              "Enter anything else if you have data for "
                              "both and would like to cal-\nculate the "
                              "pure sample's heat capacity: ")
    experiment_config = experiment_config.upper()
    if experiment_config == "A":
        import_directory = ["Blank Addenda"]
        config = "Addenda"
        molar_or_total, m, M, units = None, None, None, None
    elif experiment_config == "B":
        import_directory = ["Addenda + Sample"]
        config = "Sample w/ Addenda"
        molar_or_total, m, M, units = molar()
    else:
        import_directory = ["Addenda + Sample", "Blank Addenda"]
        config = "Pure Sample"
        molar_or_total, m, M, units = molar()

    return import_directory, config, molar_or_total, m, M, units, \
           system_name


def molar():
    molar_or_total = input("\nEnter 'M' to calculate molar heat capacity. "
                           "Enter anything else to calculate total heat "
                           "capacity:  ")
    molar_or_total = molar_or_total.upper()
    if molar_or_total == "M":
        m = float(input("\nPlease enter the sample's mass in g: "))
        M = float(input("\nPlease enter the system's molar mass in g "
                        "mol⁻¹: "))
        units = "J mol⁻¹ K⁻¹"
    else:
        m, M = 1, 1
        units = "J K⁻¹"

    return molar_or_total, m, M, units


def calculation_and_progress(import_dir, molar_or_total, units):
    # Calculation progress tracker.
    directory_size = 0
    current_index = 0
    for file in os.listdir(import_dir):
        directory_size += 1

    # Organization.
    temp_set = np.array([])
    C_set = np.array([])
    blank_files = []
    print(f"\n\nNow calculating heat capacity for "
          f"{import_dir.lower()}:\n\n\t\t\t\t   Progress\t\t\t\t\t\t"
          f"\t\t(\t\t T\t\t\t,\t\t\tCv          )\t\t \t\t\tSource")
    for file in os.listdir(import_dir):
        measurement_C, measurement_temp = heat_capacity(file, units,
                                                        import_dir)
        if not measurement_C:  # Takes care of blank data files.
            blank_files.append(file)
            directory_size -= 1
            continue
        else:
            temp_set = np.append(temp_set, measurement_temp)
            C_set = np.append(C_set, measurement_C)
            current_index += 1
            print(f"Calculation for Cv datum "
                  f"{current_index} out of {directory_size} completed. \t\t"
                  f"({temp_set[current_index - 1]}\t,  "
                  f"{C_set[current_index - 1]}) \t\t from {file}")

    return temp_set, C_set, blank_files


def blank_data_files(blank_files):
    """Identifies blank data files."""

    if blank_files:
        print("\nThe following files contain no data:")
        for i in range(len(blank_files)):
            print(f"  - {blank_files[i]}")


def heat_capacity(file, units, directory):
    """Calls power, pulse duration, & ΔT, asks user for sample's mass &
    system's molar mass, and calculates specific heat capacity."""

    P, t_d = power_and_duration(file, directory)  # Power during & duration
    #    of pulse.
    dT, temp = DeltaT(file, directory)  # Temperature change after
    #    pulse.

    if not t_d:
        return None, None
    else:
        # c = C/n = (P*t_d / dT) / (m/M) = (P*t_d*M) / (dT*m)
        Cv = (P * t_d) / dT

    return Cv, temp


def get_data(file, directory):
    """Collects all relevant data from subdirectories of the director in
    which {Tebo_Cv_Full.py} lives."""

    with open(directory + "\\" + file, newline="") as all_data:
        data = (csv.reader(all_data, delimiter="\t"))

        V = np.array([])  # Heater voltage.
        R = np.array([])  # Bridge resistance.
        t = np.array([])  # Time of measurement.
        time = np.array([])  # Temporary storage of time value.
        minute = 0  # Indexes minutes in time adjustment.
        index = 0  # Indexes reading nr.

        for column in data:
            if (data.line_num > 1):

                # Sample voltage. Voltage during the pulse is assumed to
                #    be constant.
                # Some data files have two lines of header by some mis-
                #    take, so we check for that and correct if needed.
                try:
                    V = np.append(V, float(column[0]))
                except ValueError:
                    continue

                # Bridge resistance.
                R = np.append(R, float(column[1]) / 1_000)

                # Time of datum measurement. The time attached to each
                #    resistance datum in the LabView output .txt doc-
                #    uments is the computer clock time's second hand.
                #    This function looks for the top-of-the-minute tran-
                #    sitions and adds 60 seconds to all of the time data
                #    after each transition.
                try:  # For data taken after June 30, 2021.
                    time = np.append(time, float(column[3]))
                    if len(t) < 1:
                        t = np.append(t, time[len(t)])
                    # Finds 59.### --> 00.### and adds 60 after the transition
                    else:
                        if time[len(t) - 1] > time[len(t)]:
                            minute += 1
                        t = np.append(t, (60 * minute) + time[len(t)])
                except IndexError:  # For data taken before July 1, 2021.
                    index += 1
                    t = np.append(t, .118343465045593 * index)
                    # The assumption that each measurement is .1 s apart is inaccurate, and gives a Cv error of
                    #    -15.5%. Gregorio Ponti's CvT.py makes this assumption.
                    # The Cv vs. T sets that I calculated in Sum '21 also use that inaccurate assumption.

    return V, R, t


def power_and_duration(file, directory):
    """Imports voltage data, calls time data, and calculates power
    delivered during pulse."""

    # Sample voltage. Voltage during the pulse is assumed to be constant.
    V = get_data(file, directory)[0]

    # Pulse duration. Voltage plot is assumed to be a step function,
    #    i.e. the voltage instantly goes from zero to its maximum value,
    #    and then goes instantly from its maximum value back to zero.
    t = get_data(file, directory)[2]
    in_pulse = np.array([])
    just_did_this_step = False
    for i in range(len(V)):
        if ((max(V) - V[i]) / max(V)) <= .1:
            in_pulse = np.append(in_pulse, t[i])
        # Voltage measurements are taken to ..ॱ⁻ॱ⁻ॱ⁻.. i.e. we assume the
        #    voltage stays at the previous measured value until we
        #    actually measure it again.
        if in_pulse.size > 0 and not just_did_this_step and t[i - 1] == \
                in_pulse[-1]:
            in_pulse = np.append(in_pulse, t[i])
            just_did_this_step = True
    try:
        pulse_duration = in_pulse[-1] - in_pulse[0]  # duration = end - beg.
    except IndexError:
        return None, None

    P = (max(V)) ** 2 / 160  # P = V²/R. The heater's resistance is
    #    160. Ω for the entire temperature
    #    range under which it operates.

    return P, pulse_duration


def DeltaT(file, directory):
    """Calls normalized temperature data & normalized fit curve,
    performs ΔT calculation, returns result."""

    t = get_data(file, directory)[2]
    norm_T, norm_fit, temp = decay_fit(file, directory)
    area_difference = []

    # Finds index for which the integral under the data to the left of
    #    that index is equal to the integral between the fit curve and
    #    the data to the right of that index.
    for i in range(1, len(t)):
        area_under = intg.simps(norm_T[:i], t[:i])
        area_between = intg.simps(norm_fit[i:] - norm_T[i:], t[i:])
        area_difference.append(np.abs(area_under - area_between))

    # Since our data and fit curve are given discretely and the inte-
    #    grals are numerical, it is impossible to find an index for
    #    which the integrals exactly equal each other. Instead we search
    #    for the index at which the integrals are the closest.
    try:
        dividing_line = area_difference.index(min(area_difference))
    except ValueError:
        return None, None

    return norm_fit[dividing_line], temp


def decay_fit(file, directory):
    """Fits data with curve_fit and subtracts-off baseline temperature."""

    t = get_data(file, directory)[2]  # Time coordinate of each T datum.
    T = T_calculation(file, directory)  # Calculation of T data.
    try:
        max_at = np.argmax(T)  # Finds index of peak T value.
    except ValueError:
        return None, None, None
    fit_T = T[max_at:]  # Excludes data before peak T.
    fit_t = t[max_at:]  # Excludes data before peak T.

    opt, covar = cf(temp_fitting_fct, fit_t, fit_T, bounds=(0, 100_000))
    # opt == [tau A B]
    tau = opt[0]  # Time constant.
    A = opt[1]  # Scale factor.
    B = opt[2]  # Baseline temperature.

    # Normalizes fitted curve.
    norm_T = T - B
    norm_fit = temp_fitting_fct(t, tau, A, B) - B

    return norm_T, norm_fit, B


def addenda_fit(addenda_data, combined_sys_T):
    """Fits a curve to the addenda data so that it can be subtracted
    from the addenda + sample data to get the pure sample's heat
    capacity."""

    addenda_T = addenda_data[0]  # Temperature coordinate of each C datum.
    addenda_C = addenda_data[1]  # Addenda heat capacities for many temps.

    opt, covar = cf(addenda_fitting_fct, addenda_T, addenda_C)
    # opt == [A B C D]
    A = opt[0]
    B = opt[1]
    C = opt[2]
    D = opt[3]

    # Extrapolates addenda heat capacity to T values for which there is
    #    addenda + sample heat capacity data.
    usable_addenda_C = addenda_fitting_fct(combined_sys_T, A, B, C, D)

    return usable_addenda_C


def temp_fitting_fct(x, tau, A, B):
    """The fitting function we want to use with curve_fit. Comes from
    the theory presented in the Chabot thesis."""

    y = (A * np.exp(- x / tau)) + B    # Should I make this Ae^x/t + bx + c?

    return y


def addenda_fitting_fct(x, A, B, C, D):
    """The fitting function we want to use with addenda_fit. Arbitrary –
    we don't really care what the form of the curve is as long as it
    fits the data well."""

    y = (A * (x ** 3)) + (B * (x ** 2)) + (C * x) + D

    return y


def T_calculation(file, directory):
    """Imports raw data, performs T calculation, returns result."""

    R = get_data(file, directory)[1]  # Input data (bridge resistance).
    T_data = np.array([])  # Output data (thermometer temperature).
    for i in range(len(R)):
        # Calibration constants given in RMC manual.
        a1 = .2700714
        a2 = .07323249
        a3 = 1.0219809
        a4 = -1.9558035
        a5 = 3.4889894
        a6 = -4.4175955
        a7 = 3.6976856
        a8 = -1.7871035
        a9 = .3821843
        R0 = 77.40

        # Formula for T given in RMC manual.
        T_datum = float(1. / (a1
                              + a2 * (np.log(R[i] / R0)) ** .25
                              + a3 * (np.log(R[i] / R0)) ** 1.
                              + a4 * (np.log(R[i] / R0)) ** 2.
                              + a5 * (np.log(R[i] / R0)) ** 3.
                              + a6 * (np.log(R[i] / R0)) ** 4.
                              + a7 * (np.log(R[i] / R0)) ** 5.
                              + a8 * (np.log(R[i] / R0)) ** 6.
                              + a9 * (np.log(R[i] / R0)) ** 7.) ** 4.)

        T_data = np.append(T_data, T_datum)

    return T_data


def plot_data(import_directory, system_name, molar_or_total,
              temp_set, C_set):
    """Formats and displays graph."""

    if import_directory == "Blank Addenda":
        plt.plot(temp_set, C_set, "o", label="Addenda")
        plt.title(f"$C_v$ vs. $T$ for Addenda",fontsize=24)
        wanna_save = input("\nEnter 'A' to save the addenda's Cv data. "
                          "Enter anything else to forgo saving: ")
        wanna_save = wanna_save.upper()
        if wanna_save == "A":
            write_addenda(temp_set, C_set)
    elif import_directory == "Addenda + Sample":
        plt.plot(temp_set, C_set, "o", label=f"{system_name} w/ Addenda")
        plt.title(f"$C_v$ vs $T$ for {system_name} + Addenda",fontsize=24)
        wanna_save = input("\nEnter 'A' to save the sample + addenda Cv "
                           "data. Enter anything else to forgo saving: ")
        wanna_save = wanna_save.upper()
        if wanna_save == "A":
            write_addenda(temp_set, C_set)
    else:
        plt.plot(temp_set, C_set, "o", label=f"{system_name}")
        plt.title(f"$C_v$ vs $T$ for Pure {system_name}",fontsize=24)
        wanna_save = input("\nEnter 'A' to save the pure sample's Cv data."
                           " Enter anything else to forgo saving: ")
        wanna_save = wanna_save.upper()
        if wanna_save == "A":
            write_addenda(temp_set, C_set)
    if molar_or_total == "M" and type(import_directory) == list:
        plt.ylabel("Molar Heat Capacity (J mol$^{-1}$ K$^{-1}$)",fontsize=16)
    else:
        plt.ylabel("Heat Capacity (J K$^{-1}$)",fontsize=16)
    plt.xlabel("Temperature (K)",fontsize=16)
    plt.xticks(np.arange(0, max(temp_set) + 2, 1))
    # plt.yscale("log")
    plt.grid(which="major")
    plt.grid(which="minor", axis="y", linestyle=":")
    plt.legend(fontsize=12)
    plt.show()


def write_addenda(temp_set, C_set):
    dirname = os.path.dirname(__file__)
    save_directory = "Permanent Storage/" + \
                     input("\nEnter the name of the folder in Permanent "
                           "Storage to which you would like to save the "
                           "addenda's Cv data: ")
    save_file = save_directory + "/" + \
                input("\nEnter the name that you would like to call the "
                      "save file: ")
    export_file = os.path.join(dirname, save_file)
    heading = ["Temperature (K)", "Heat Capacity (J/K)"]

    with open(export_file, "w", newline="") as \
            output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(heading)  # Creates column headers.
        for i in range(len(C_set)):  # Populates file with data.
            writer.writerow([temp_set[i], C_set[i]])


def read_addenda():
    temp_set = np.array([])
    C_set = np.array([])

    dirname = os.path.dirname(__file__)
    read_directory = "Permanent Storage/" + input("\nEnter the name of the "
                                                  "folder in Permanent "
                                                  "Storage in which the "
                                                  "addenda's Cv data file "
                                                  "is found: ")
    read_file = read_directory + "/" + input("\nEnter the name of the file "
                                             "containing the addenda's Cv "
                                             "data: ")
    addenda_data = os.path.join(dirname, read_file)

    with open(addenda_data, newline="") as all_data:
        data = (csv.reader(all_data, delimiter="\t"))
        print(data)
        for row in data:
            if (data.line_num == 1):
                continue
            else:
                temp_set = np.append(temp_set, float(row[0]))
                C_set = np.append(C_set, float(row[1]))

    return temp_set, C_set


main()
