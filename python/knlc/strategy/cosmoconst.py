from . import constants as const

"""
The units of H are s-1km/Mpcs-1km/MPc. Kilometers and Megaparsecs are both units of length, so their ratio is a dimensionless number. If N is the number of kilometers in a megaparsec then the ratio is just 1/N. To convert the Hubble constant to units of per second-1 just divide it by const.mpc.
"""



h = 0.7
H0 = h * 100.0 #km/sec/mpc
H0_sec = H0/const.mpc # hubble constant in sec
