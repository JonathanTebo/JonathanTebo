import csv
import numpy as np
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit as cf



def main():
    """Import from second dataset using DictReader method."""

    # dirname = os.path.dirname(__file__)
    read_directory = "C:\\Users\\jonat\\Desktop\\Past Schoolwork\\Undergrad\\PHY 371C (S5)\\Heat Capacity " \
                     "Measurement\\Storage for Testing"
    for file in os.listdir(read_directory):
        xset2 = np.array([])
        yset2 = np.array([])
        print("-----------------------------------------------------------")
        print("In for loop.")
        print(read_directory)
        print(file)
        filename = read_directory + "\\" + file
        with open(filename, newline="") as all_data:
            data = csv.DictReader(all_data, delimiter="\t")
            for column in data:
                x2 = float(column["Temperature (K)"])
                y2 = float(column["Heat Capacity (J/K)"])
                xset2 = np.append(xset2, x2)
                yset2 = np.append(yset2, y2)
            print("input values (2):", xset2, "\n output values (2):", yset2)
                # Prints in the same format as for dataset 1, but with more
                # concise code.
        plt.plot(xset2, yset2, ".", label=file)
        plt.xlabel("Temperature (K)",fontsize=16)
        plt.ylabel("Heat Capacity (J/K)",fontsize=16)
        plt.legend()
    plt.title("Cv vs. T for Blank Addendum, Sample + Addendum, and Pure "
              "Sample",fontsize=24)
    plt.grid(which="major")
    plt.grid(which="minor", axis="y", linestyle=":")
    plt.show()


main()
