import os 
import pandas as pd
import sys
import csv
import shutil

curDir = os.getcwd()

for file in os.listdir(curDir):
	if file.startswith('test_'):
		shutil.rmtree(file)
	elif file.startswith('BMNS_params_'):
		os.remove(file)

#for file in os.listdir(curDir):
	#if file.endswith('-mc.csv'):
atom = "U38N3_wtTAR_35C" 

f = open('test/2-state_Formatted_LocalFits_{}.csv'.format(atom))
csv_reader = csv.reader(f)
next(csv_reader)
col = next(csv_reader)

df = pd.DataFrame(columns=col)
i = 0
for row in csv_reader:
	if row[0] == 'Red. Chi-sq':
		row.append(0)
	df.loc[i] = row
	i += 1

best_fit_pb = float(df[df['%s Pars' % atom]=='pB (%)']['Fit Value'].values[0])
best_fit_dw = float(df[df['%s Pars' % atom]=='dwB (ppm)']['Fit Value'].values[0])
best_fit_kex = float(df[df['%s Pars' % atom]=='kexAB (s^-1)']['Fit Value'].values[0])
#rchi2 = df[df['param']=='rchi2']['fitval'].values[0]

# larmor frequency
if str(sys.argv[1]) == 'nysbc700C':
    lf = 176.082556
elif str(sys.argv[1]) == '600C':
    lf = 150.784627
elif str(sys.argv[1]) == '700N':
    lf = 70.960783
elif str(sys.argv[1]) == '600N':
    lf = 60.810645 
# pb
BM_file1 = '''
##################################################################################
# Run the BMNS fitting program:
# > python BMNS.py -fit [BM Parameter Input File] [R1rho Data Directory] (Optional Output directory)
##################################################################################
# Define fitting setup.
# FitType: can be 'global' or 'local' or 'brute'
#          This is for global or local optimizations, not shared parameter fits.
#          'Brute' designates brute-force fixed calculations of the range of parameter
#                   space designated by lower/upper bounds on parameters.
#          - 'brutep' will generate plots at each increment point.
#             WARNING: This can take a LONG time.
#          'Bruteint' brute-forces parameter space by fitting intensities instead of
#                     R1p values
#
#          'Local' uses Levenberg-Marquardt semi-gradient descent/Gauss-Newton algorithm
#          - 'localint' fits intensities directly rather than R1p
#          'Global' uses the "Adaptive Memory Programming for Global Optimizations"
#                   algorithm, with the local 'L-BFGS-B' function, and polishes the
#                   global optimum with L-M.
# FitEqn: fit equation, "BM" for Bloch-McConnell or "Lag" for Laguerre 2-/3-state
# NumFits: is number of fit minima to find (ie. loop over fitting algorithm)
# RandomFitStart : can be 'Yes' or 'No'
#                  if 'Yes', randomly selects initial guess from parameter bounds
##################################################################################
+
FitType local
FitEqn BM
NumFits 1
RandomFitStart No

##################################################################################
# Define fit parameter data, data names, base freqs,
#  initial parameter guesses, and paramter lower and upper bounds.
#
# Add '+' to read in an additional set of parameters with given 'Name XYZ'
#   The 'Name' must match a .csv data file in given directory of the same name.
#
# Rows for parameters are as follows:
#  [Par name] [initial value] [lower bounds] [upper bounds] ([optional brute force number])
#
# If both lower and upper bounds are not given, they will be set to large values.
# '!' designates a fixed parameter that will not change throughout the fit.
# '*' designates a shared parameter that will be fitted for all data sets
#     also containing the 'x' flag, in a shared manner.
# '@' designates linear brute-force over parameter range of low to upper bounds
# '$' designates log brute-force over parameter range of low to upper bounds
#
# If R1b/c or R2b/c are fixed to 0, they will be shared with R1 / R2
#  e.g. "R1b! = 0.0" will be interpreted as "R1b = R1"
#
# lf = Larmor frequency (MHz) of the nucleus of interest
#      15N:   60.76302 (600) or  70.960783 (700)
#      13C: 150.784627 (600) or 176.090575 (700)
#
# (optional) rnddel = Fraction of data to be randomly deleted before fit
#                     e.g 'rnddel 0.1' would randomly delete 10pct of data
#
# Temp [Celsius or Kelvin] : Define temperature to calculate free energies
#
# AlignMag [Auto/Avg/GS]
#          Auto : calculates kex/dw and aligns mag depending on slow (gs) vs. fast (avg)
#          Avg : Aligns magnetization/projects along average effective field of GS/ESs
#          GS : Aligns magnetization along ground-state
#
# x-axis Lower Upper (Hz): Sets lower and upper x-axis limits for both plots
#   if not given, will automatically set them
#
# y-axis Lower Upper : Sets lower and upper y-axis limits for both plots
#   if not given, will automatically set them
#
# Trelax increment Tmax (seconds) : sets the increment delay and maximum relaxation
#  delay to simulate R1rho at.
#  Use caution with this flag, recommended that is remains commented out.
#  Array of delays is given as a linear spacing from 0 - Tmax in Tmax/Tinc number of points
#  If not defined, the program will calculate the best Tmax from the experimental
#   R1rho data.
##################################################################################

+
Name {}
lf {}
Temp 35.0
AlignMag auto
#Trelax 0.0005 0.5
x-axis -2000 2000
#y-axis 0 50
pB! {} 1e-4 0.5
pC! 0.0 1e-6 0.5
dwB {} {} {}
dwC! 0.0 -80 80
kexAB 1000.0 1.0 500000.0
kexAC! 0.0 1.0 500000.0
kexBC! 0.0 1.0 500000.0
R1 2.5 1e-6 20.
R2 16.0 1e-6 200.
R1b! 0.0
R2b! 0.0
R1c! 0.0
R2c! 0.0

'''
if best_fit_dw < 0:
	pb_dwb_min = best_fit_dw - 5
	pb_dwb_max = 0
else:
	pb_dwb_min = 0
	pb_dwb_max = best_fit_dw + 5

# pB range
for i in [0.1, 0.33, 0.4,0.5,0.6,0.7,0.8,0.9, 1,1.1,1.2,1.3,1.4,1.5, 3, 10]:
	with open ("BMNS_params_pb_{}.txt".format(i), "w") as f:
		f.write(BM_file1.format(atom, lf, i*best_fit_pb/100, best_fit_dw, pb_dwb_min, pb_dwb_max))
	os.system('python /home/sg4109/scripts/BMNS-master/BMNS.py -fit BMNS_params_pb_{0}.txt . test_pb_{0}'.format(i))

# dw
BM_file2 = '''
##################################################################################
# Run the BMNS fitting program:
# > python BMNS.py -fit [BM Parameter Input File] [R1rho Data Directory] (Optional Output directory)
##################################################################################
# Define fitting setup.
# FitType: can be 'global' or 'local' or 'brute'
#          This is for global or local optimizations, not shared parameter fits.
#          'Brute' designates brute-force fixed calculations of the range of parameter
#                   space designated by lower/upper bounds on parameters.
#          - 'brutep' will generate plots at each increment point.
#             WARNING: This can take a LONG time.
#          'Bruteint' brute-forces parameter space by fitting intensities instead of
#                     R1p values
#
#          'Local' uses Levenberg-Marquardt semi-gradient descent/Gauss-Newton algorithm
#          - 'localint' fits intensities directly rather than R1p
#          'Global' uses the "Adaptive Memory Programming for Global Optimizations"
#                   algorithm, with the local 'L-BFGS-B' function, and polishes the
#                   global optimum with L-M.
# FitEqn: fit equation, "BM" for Bloch-McConnell or "Lag" for Laguerre 2-/3-state
# NumFits: is number of fit minima to find (ie. loop over fitting algorithm)
# RandomFitStart : can be 'Yes' or 'No'
#                  if 'Yes', randomly selects initial guess from parameter bounds
##################################################################################
+
FitType local
FitEqn BM
NumFits 1
RandomFitStart No

##################################################################################
# Define fit parameter data, data names, base freqs,
#  initial parameter guesses, and paramter lower and upper bounds.
#
# Add '+' to read in an additional set of parameters with given 'Name XYZ'
#   The 'Name' must match a .csv data file in given directory of the same name.
#
# Rows for parameters are as follows:
#  [Par name] [initial value] [lower bounds] [upper bounds] ([optional brute force number])
#
# If both lower and upper bounds are not given, they will be set to large values.
# '!' designates a fixed parameter that will not change throughout the fit.
# '*' designates a shared parameter that will be fitted for all data sets
#     also containing the 'x' flag, in a shared manner.
# '@' designates linear brute-force over parameter range of low to upper bounds
# '$' designates log brute-force over parameter range of low to upper bounds
#
# If R1b/c or R2b/c are fixed to 0, they will be shared with R1 / R2
#  e.g. "R1b! = 0.0" will be interpreted as "R1b = R1"
#
# lf = Larmor frequency (MHz) of the nucleus of interest
#      15N:   60.76302 (600) or  70.960783 (700)
#      13C: 150.784627 (600) or 176.090575 (700)
#
# (optional) rnddel = Fraction of data to be randomly deleted before fit
#                     e.g 'rnddel 0.1' would randomly delete 10pct of data
#
# Temp [Celsius or Kelvin] : Define temperature to calculate free energies
#
# AlignMag [Auto/Avg/GS]
#          Auto : calculates kex/dw and aligns mag depending on slow (gs) vs. fast (avg)
#          Avg : Aligns magnetization/projects along average effective field of GS/ESs
#          GS : Aligns magnetization along ground-state
#
# x-axis Lower Upper (Hz): Sets lower and upper x-axis limits for both plots
#   if not given, will automatically set them
#
# y-axis Lower Upper : Sets lower and upper y-axis limits for both plots
#   if not given, will automatically set them
#
# Trelax increment Tmax (seconds) : sets the increment delay and maximum relaxation
#  delay to simulate R1rho at.
#  Use caution with this flag, recommended that is remains commented out.
#  Array of delays is given as a linear spacing from 0 - Tmax in Tmax/Tinc number of points
#  If not defined, the program will calculate the best Tmax from the experimental
#   R1rho data.
##################################################################################

+
Name {}
lf {}
Temp 35.0
AlignMag auto
#Trelax 0.0005 0.5
#x-axis -2000 2000
#y-axis 0 50
pB 0.01 1e-4 0.5
pC! 0.0 1e-6 0.5
dwB! {} {} {}
dwC! 0.0 -80 80
kexAB 1000.0 1.0 500000.0
kexAC! 0.0 1.0 500000.0
kexBC! 0.0 1.0 500000.0
R1 2.5 1e-6 20.
R2 16.0 1e-6 200.
R1b! 0.0
R2b! 0.0
R1c! 0.0
R2c! 0.0

'''

# dw range
dw = [best_fit_dw - 2, best_fit_dw - 1, best_fit_dw - 0.5, best_fit_dw + 0.5, best_fit_dw + 1, best_fit_dw + 2]
min_dw = min(dw) - 5
max_dw = max(dw) + 5
for i in dw:
	with open ("BMNS_params_dw_{}.txt".format(i), "w") as f:
		f.write(BM_file2.format(atom, lf, i, min_dw, max_dw))
	os.system('python /home/sg4109/scripts/BMNS-master/BMNS.py -fit BMNS_params_dw_{0}.txt . test_dw_{0}'.format(i))


# kex
BM_file3 = '''
##################################################################################
# Run the BMNS fitting program:
# > python BMNS.py -fit [BM Parameter Input File] [R1rho Data Directory] (Optional Output directory)
##################################################################################
# Define fitting setup.
# FitType: can be 'global' or 'local' or 'brute'
#          This is for global or local optimizations, not shared parameter fits.
#          'Brute' designates brute-force fixed calculations of the range of parameter
#                   space designated by lower/upper bounds on parameters.
#          - 'brutep' will generate plots at each increment point.
#             WARNING: This can take a LONG time.
#          'Bruteint' brute-forces parameter space by fitting intensities instead of
#                     R1p values
#
#          'Local' uses Levenberg-Marquardt semi-gradient descent/Gauss-Newton algorithm
#          - 'localint' fits intensities directly rather than R1p
#          'Global' uses the "Adaptive Memory Programming for Global Optimizations"
#                   algorithm, with the local 'L-BFGS-B' function, and polishes the
#                   global optimum with L-M.
# FitEqn: fit equation, "BM" for Bloch-McConnell or "Lag" for Laguerre 2-/3-state
# NumFits: is number of fit minima to find (ie. loop over fitting algorithm)
# RandomFitStart : can be 'Yes' or 'No'
#                  if 'Yes', randomly selects initial guess from parameter bounds
##################################################################################
+
FitType local
FitEqn BM
NumFits 1
RandomFitStart No

##################################################################################
# Define fit parameter data, data names, base freqs,
#  initial parameter guesses, and paramter lower and upper bounds.
#
# Add '+' to read in an additional set of parameters with given 'Name XYZ'
#   The 'Name' must match a .csv data file in given directory of the same name.
#
# Rows for parameters are as follows:
#  [Par name] [initial value] [lower bounds] [upper bounds] ([optional brute force number])
#
# If both lower and upper bounds are not given, they will be set to large values.
# '!' designates a fixed parameter that will not change throughout the fit.
# '*' designates a shared parameter that will be fitted for all data sets
#     also containing the 'x' flag, in a shared manner.
# '@' designates linear brute-force over parameter range of low to upper bounds
# '$' designates log brute-force over parameter range of low to upper bounds
#
# If R1b/c or R2b/c are fixed to 0, they will be shared with R1 / R2
#  e.g. "R1b! = 0.0" will be interpreted as "R1b = R1"
#
# lf = Larmor frequency (MHz) of the nucleus of interest
#      15N:   60.76302 (600) or  70.960783 (700)
#      13C: 150.784627 (600) or 176.090575 (700)
#
# (optional) rnddel = Fraction of data to be randomly deleted before fit
#                     e.g 'rnddel 0.1' would randomly delete 10pct of data
#
# Temp [Celsius or Kelvin] : Define temperature to calculate free energies
#
# AlignMag [Auto/Avg/GS]
#          Auto : calculates kex/dw and aligns mag depending on slow (gs) vs. fast (avg)
#          Avg : Aligns magnetization/projects along average effective field of GS/ESs
#          GS : Aligns magnetization along ground-state
#
# x-axis Lower Upper (Hz): Sets lower and upper x-axis limits for both plots
#   if not given, will automatically set them
#
# y-axis Lower Upper : Sets lower and upper y-axis limits for both plots
#   if not given, will automatically set them
#
# Trelax increment Tmax (seconds) : sets the increment delay and maximum relaxation
#  delay to simulate R1rho at.
#  Use caution with this flag, recommended that is remains commented out.
#  Array of delays is given as a linear spacing from 0 - Tmax in Tmax/Tinc number of points
#  If not defined, the program will calculate the best Tmax from the experimental
#   R1rho data.
##################################################################################

+
Name {}
lf {}
Temp 35.0
AlignMag auto
#Trelax 0.0005 0.5
#x-axis -2000 2000
#y-axis 0 50
pB 0.01 1e-4 0.5
pC! 0.0 1e-6 0.5
dwB {} {} {}
dwC! 0.0 -80 80
kexAB {} {} {}
kexAC! 0.0 1.0 500000.0
kexBC! 0.0 1.0 500000.0
R1 2.5 1e-6 20.
R2 16.0 1e-6 200.
R1b! 0.0
R2b 16.0 1e-6 200.
R1c! 0.0
R2c! 0.0

'''
if best_fit_dw < 0:
	pb_dwb_min = best_fit_dw - 5
	pb_dwb_max = 0
else:
	pb_dwb_min = 0
	pb_dwb_max = best_fit_dw + 5

# kex range
for i in [0.1, 0.33, 1, 3, 10]:
	with open ("BMNS_params_kex_{}.txt".format(i), "w") as f:
		f.write(BM_file3.format(atom, lf, best_fit_dw, pb_dwb_min, pb_dwb_max, i*best_fit_kex, i*best_fit_kex-1, i*best_fit_kex+1))
	os.system('python /home/sg4109/scripts/BMNS-master/BMNS.py -fit BMNS_params_kex_{0}.txt . test_kex_{0}'.format(i))
