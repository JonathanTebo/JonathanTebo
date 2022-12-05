# This program calculates the temperature of the ruthenium(IV) oxide
#    thermometer attached to the addenda given the resistance measured
#    in it. The resistance can be read off the SRS SIM921 AC Resistance
#    Bridge display, and must be entered by the user manually into this
#    program.

import numpy as np


def main():

    # Startup message.
    print("This program calculates the temperature of the ruthenium(IV) "
          "oxide thermometer given the resistance measured through it.")
    print("The program will continue to ask for resistance inputs until the"
          " user manually terminates it.")
    print("To terminate the program at any time, enter the word 'EXIT' when"
          " asked to enter a resistance value.")
    print("-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   "
          "-   -   -   -   -   -   -   -   -   -   -   -   -   -   -")

    # Ruthenium(IV) oxide thermometer resistance.
    R = input("Please enter the resistance measured by the SRS SIM921 AC "
              "Resistance Bridge in kΩ: ")

    # Potential exit.
    if(R == "EXIT"):
        exit()
    else:
        R = float(R)

    # If too warm out of calibration range:
    while(R < 86):
        print("The thermometer only works for resistances between 86 kΩ and"
              " 246 kΩ.")
        print("This corresponds to a temperature range of 40 K down to 2 "
              "K.")
        print("Please wait for the thermometer to cool down before trying "
              "again.")
        print("-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   "
              "-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   "
              "-")
        R = input("Please enter the resistance measured by the SRS SIM921 "
                  "AC Resistance Bridge in kΩ: ")
        if(R == "EXIT"):       # Potential exit.
            exit()
        else:
            R = float(R)

    # It'll never get too cold inside of our apparatus, so don't worry
    #    about the lower bound of the thermometer's calibration range.

    # Temperature calculation for temperatures inside thermometer's
    #    calibration range. Allows user to submit valid resistances for
    #    calculation repeatedly.
    while(type(R) == float):
        T = 1. / (0.2700714
              + .07323249 * (np.log(.01291989664082687 * R)) ** .25
              + 1.0219809 * (np.log(.01291989664082687 * R)) ** 1.
              - 1.9558035 * (np.log(.01291989664082687 * R)) ** 2.
              + 3.4889894 * (np.log(.01291989664082687 * R)) ** 3.
              - 4.4175955 * (np.log(.01291989664082687 * R)) ** 4.
              + 3.6976856 * (np.log(.01291989664082687 * R)) ** 5.
              - 1.7871035 * (np.log(.01291989664082687 * R)) ** 6.
              + 0.3821843 * (np.log(.01291989664082687 * R)) ** 7.) ** 4.
        if(T > 10):
            T = round(T, 3)     # Rounds output temperature to the correct
        else:                   #    number of significant digits.
            T = round (T, 4)
        print("The temperature of the ruthenium(IV) oxide thermometer is",T,
              "K.")
        print("-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   "
              "-   -   -   -   -   -   -   -   -   -   -   -   -   -   -   "
              "-")
        R = input("Please enter the resistance measured by the SRS SIM921 "
                  "AC Resistance Bridge in kΩ: ")
        if(R == "EXIT"):       # Potential exit.
            exit()
        else:
            R = float(R)
        while(R < 86):     # Just in case the user tries an invalid R.
            print("The thermometer only works for resistances between 86 kΩ"
                  " and 246 kΩ.")
            print("This corresponds to a temperature range of 40 K down to "
                  "2 K.")
            print("Please wait for the thermometer to cool down before "
                  "trying again.")
            print("-   -   -   -   -   -   -   -   -   -   -   -   -   -   "
                  "-   -   -   -   -   -   -   -   -   -   -   -   -   -   "
                  "-   -   -")
            R = input("Please enter the resistance measured by the SRS "
                      "SIM921 AC Resistance Bridge in kΩ: ")
            if(R == "EXIT"):       # Potential exit.
                exit()
            else:
                R = float(R)


main()
