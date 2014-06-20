""" This is a package for calculating electronic quantities in LiFeAs using
the tight binding model from this paper:
Y. Wang, et. al. "Superconducting gap in LiFeAs from three-dimensional 
spin-fluctuation pairing calculations", Phys. Rev. B 88, 174516 (2014)

The methods included in this package are listed below:
defineHamiltonian = build Hamiltonian using a k-space mesh
getGreens = calculate the Green's functions from the Hamiltonian for a given energy
getSpectralFunction = calculate the spectral function from the Green's function
fermiSurfaceExample = An example process using the above three codes to plot
 					the Fermi surface in the 1-Fe unit cell Brillouin Zone
					on a 201 x 201 k-space mesh. """

import numpy as np
import scipy as sci
from matplotlib.pyplot import imshow as im

cos = np.cos
sin = np.sin




def defineHamiltonian(n,npi):
	# n is the number of pixels per 2*pi/a (Number of pixels in 1st Bz)
	# npi is the integer number of pi/a's you would like to calculate out to in k-space.

	# A note on sub- and super-script notation from the paper:

	#First Fe in unit cell
	# 1 = dxy
	# 2 = dx^2-y^2
	# 3 = idxz
	# 4 = idyz
	# 5 = dz^2

	#Second Fe in unit cell
	# 6 = dxy
	# 7 = dx^2-y^2
	# 8 = -idxz
	# 9 = -idyz
	# 10 = dz^2


	## Define on-site energies
	e1 = 0.02
	e2 = -0.2605
	e3 = -0.0075
	# e4 = e3
	e5 = -0.3045

	## Define hoppings and hybridizations

	# 2D Part (Note: numbers before _ are superscript. Numbers after are subscript)
	t11_11 = 0.030
	t10_16 = -0.0185
	t20_11 = -0.010
	t21_16 = 0.0035
	t11_13 = -0.0635*1j
	t10_18 = 0.155*1j
	t11_15 = -0.090

	t10_27 = -0.2225
	t11_22 = 0.070
	t10_29 = -0.1925*1j
	t11_23 = -0.010*1j
	t10_210 = 0.1615

	t11_33 = 0.152
	t10_38 = 0.050
	t20_33 = -0.004
	t21_38 = 0.040
	t02_33 = -0.051
	t10_49 = 0.210
	t22_33 = -0.005

	t21_49 = -0.053
	t11_34 = 0.090
	t10_410 = 0.0995*1j

	t11_35 = 0.1005*1j

	# 3D Part
	t101_16 = -0.004
	t001_11 = 0.0105
	t111_11 = 0.
	t201_11 = 0.
	t201_14 = 0.
	t001_33 = -0.003
	t201_33 = 0.
	t021_33 = 0.0105
	t121_16 = 0.
	t101_18 = 0.
	t101_19 = 0.
	t211_19 = 0.
	t101_38 = 0.0115
	t121_38 = 0.
	t101_39 = 0.
	t101_49 = 0.
	t121_49 = 0.

	## Build k-space mesh.
	if np.mod(n,2)==1: # Make sure n is even
		n = n+1

	nx = n*npi+1 # Total number of points in k-space
	ny = n*npi+1

	kxv = np.linspace(-npi*np.pi,npi*np.pi,nx)
	kyv = np.linspace(-npi*np.pi,npi*np.pi,ny)
	kx,ky = np.meshgrid(kxv,kyv) # k-space mesh

	kz = np.pi/2
	k1 = kx+ky
	k2 = -kx+ky

	# Empty Hamiltonian. Dimensions are such that the first two correspond to the
	# 10 x 10 Hamiltonian, and the last two dimensions are for each point in k-space.
	H0 = np.empty((10,10,nx,ny),float)*(0+0j) # Allow for complex numbers

	## Define matrix elements in the Hamiltonian

	# 2D part
	Hpp_11 = e1+2*t11_11*(cos(k1)+cos(k2))+2*t20_11* \
		(cos(2*kx)+cos(2*ky))
	Hpp_12 = 0*kx
	Hpp_13 = 2*1j*t11_13*(sin(k1)-sin(k2))
	Hpp_14 = 2*1j*t11_13*(sin(k1)+sin(k2))
	Hpp_15 = 2*t11_15*(cos(k1)-cos(k2))

	Hpp_22 = e2+2*t11_22*(cos(k1)+cos(k2))
	Hpp_23 = 2*1j*t11_23*(sin(k1)+sin(k2))
	Hpp_24 = 2*1j*t11_23*(-sin(k1)+sin(k2))
	Hpp_25 = 0*kx

	Hpp_33 = e3+2*t11_33*(cos(k1)+cos(k2))+2*t20_33*cos(2*kx)+ \
		2*t02_33*cos(2*ky)+4*t22_33*cos(2*kx)*cos(2*ky)
	Hpp_34 = 2*t11_34*(cos(k1)-cos(k2))
	Hpp_35 = 2*1j*t11_35*(sin(k1)+sin(k2))

	Hpp_44 = e3+2*t11_33*(cos(k1)+cos(k2))+2*t02_33*cos(2*kx)+ \
		2*t20_33*cos(2*ky)+4*t22_33*cos(2*kx)*cos(2*ky)
	Hpp_45 = 2*1j*t11_35*(sin(k1)-sin(k2))

	Hpp_55 = e5*np.ones((nx,ny))     
	        
	Hpm_16 = 2*t10_16*(cos(kx)+cos(ky))+ \
		2*t21_16*((cos(k1)+cos(k2))*(cos(kx)+cos(ky))- \
		sin(k1)*(sin(kx)+sin(ky))+sin(k2)*(sin(kx)-sin(ky)))
	Hpm_17 = 0*kx
	Hpm_18 = 2*1j*t10_18*sin(kx)
	Hpm_19 = 2*1j*t10_18*sin(ky)
	Hpm_110 = 0*kx

	Hpm_27 = 2*t10_27*(cos(kx)+cos(ky))
	Hpm_28 = -2*1j*t10_29*sin(ky)
	Hpm_29 = 2*1j*t10_29*sin(kx)
	Hpm_210 = 2*t10_210*(cos(kx)-cos(ky))

	Hpm_38 = 2*t10_38*cos(kx)+2*t10_49*cos(ky)+ \
		2*t21_38*((cos(k1)+cos(k2))*cos(kx)-(sin(k1)-sin(k2))*sin(kx))+ \
		2*t21_49*((cos(k1)+cos(k2))*cos(ky)-(sin(k1)+sin(k2))*sin(ky))
	Hpm_39 = 0*kx
	Hpm_310 = 2*1j*t10_410*sin(ky)

	Hpm_49 = 2*t10_49*cos(kx)+2*t10_38*cos(ky)+ \
		2*t21_49*((cos(k1)+cos(k2))*cos(kx)-(sin(k1)-sin(k2))*sin(kx))+ \
		2*t21_38*((cos(k1)+cos(k2))*cos(ky)-(sin(k1)+sin(k2))*sin(ky))
	Hpm_410 = 2*1j*t10_410*sin(kx)

	Hpm_510 = 0.*kx

	# 3D Part
	Hpp_11 = Hpp_11+(2*t001_11+4*t111_11*(cos(k1)+cos(k2))+ \
		4*t201_11*(cos(2*kx)+cos(2*ky)))*cos(kz)
	Hpp_13 = Hpp_13-4*t201_14*sin(2*ky)*sin(kz)
	Hpp_14 = Hpp_14-4*t201_14*sin(2*kx)*sin(kz)

	Hpp_33 = Hpp_33+(2*t001_33+4*t201_33*cos(2*kx)+4*t021_33*cos(2*ky))*cos(kz)

	Hpp_44 = Hpp_44+(2*t001_33+4*t021_33*cos(2*kx)+4*t201_33*cos(2*ky))*cos(kz)

	Hpm_16 = Hpm_16+4*t101_16*(cos(kx)+cos(ky))*cos(kz)+ \
		2*t121_16*((cos(k1+ky)+cos(k1+kx))*np.exp(1j*kz)+ \
		(cos(k2+ky)+cos(k2-kx))*np.exp(-1j*kz));
	Hpm_18 = Hpm_18+4*1j*t101_18*sin(kx)*cos(kz)-4*t101_19*sin(ky)*sin(kz)+ \
		2*1j*t211_19*(sin(k1+ky)*np.exp(1j*kz)-sin(k2+ky)*np.exp(-1j*kz))
	Hpm_19 = Hpm_19+4*1j*t101_18*sin(ky)*cos(kz)-4*t101_19*sin(kx)*sin(kz)+ \
		2*1j*t211_19*(sin(k1+kx)*np.exp(1j*kz)+sin(k2-kx)*np.exp(-1j*kz))

	Hpm_38 = Hpm_38+4*(t101_38*cos(kx)+t101_49*cos(ky))*cos(kz)+ \
		2*t121_38*(cos(k1+kx)*np.exp(1j*kz)+cos(k2-kx)*np.exp(-1j*kz))+ \
		2*t121_49*(cos(k1+ky)*np.exp(1j*kz)+cos(k2+ky)*np.exp(-1j*kz))
	Hpm_39 = Hpm_39+4*1j*t101_39*(cos(kx)+cos(ky))*sin(kz)

	Hpm_49 = Hpm_49+4*(t101_49*cos(kx)+t101_38*cos(ky))*cos(kz)+ \
		2*t121_49*(cos(k1+kx)*np.exp(1j*kz)+cos(k2-kx)*np.exp(-1j*kz))+ \
		2*t121_38*(cos(k1+ky)*np.exp(1j*kz)+cos(k2+ky)*np.exp(-1j*kz))

	Hpp = np.array([[Hpp_11, Hpp_12, Hpp_13, Hpp_14, Hpp_15],
		[np.conj(Hpp_12), Hpp_22, Hpp_23, Hpp_24, Hpp_25],
		[np.conj(Hpp_13), np.conj(Hpp_23), Hpp_33, Hpp_34, Hpp_35],
		[np.conj(Hpp_14), np.conj(Hpp_24), np.conj(Hpp_34), Hpp_44, Hpp_45],
		[np.conj(Hpp_15), np.conj(Hpp_25), np.conj(Hpp_35), np.conj(Hpp_45), Hpp_55]])
	Hpm = np.array([[Hpm_16, Hpm_17, Hpm_18, Hpm_19, Hpm_110],
		[Hpm_17, Hpm_27, Hpm_28, Hpm_29, Hpm_210],
		[Hpm_18, Hpm_28, Hpm_38, Hpm_39, Hpm_310],
		[Hpm_19, Hpm_29, Hpm_39, Hpm_49, Hpm_410],
		[Hpm_110, Hpm_210, Hpm_310, Hpm_410, Hpm_510]])

	## Build full Hamiltonian 10 x 10 Matrix
	H0[0:5,0:5,:,:]=Hpp
	H0[5:10,0:5,:,:]=np.conj(Hpm)
	H0[0:5,5:10,:,:]=Hpm
	H0[5:10,5:10,:,:]=np.conj(Hpp)


	return H0

def getGreens(H0,e0,gamma):
	# Calculate bare Green's function from Hamiltonian
	# e0 corresponds to the energy at which to calculate the Green's function (in eV)
	# gamma corresponds to broadening of the bare Green's function.

	nx, ny = H0.shape[-2:] 
	greens = np.empty((nx,ny,10),float)*(0+0j)

	# e0 and gamma have to be turned into matrices for subtracting the Hamiltonian.
	e0Mat = e0*np.eye(10)
	gammaMat = 1j*gamma*np.eye(10)

	for i in range(nx): # At every point in k-space, calculate the Green's function.
		for j in range(ny):	# This takes a while due to 10x10 matrix inversion at every point.
			greens[i,j,:] = np.diag(np.linalg.inv(e0Mat + gammaMat - \
				np.squeeze(H0[:,:,i,j])))

	return greens

def getSpectralFunction(greens):
	# Calculate the spectral function (what ARPES sees) from the Green's function
	nx, ny = greens.shape[-2:]

	specFun=np.zeros((nx,ny))
	
	specFun = np.sum(greens,axis=2)

	specFun = -np.imag(specFun)
	return specFun

def fermiSurfaceExample():
	# This example calculates and plots the Fermi surface for the 1st Brillouin Zone.
	n = 200
	npi = 1
	e0 = 0
	gamma = .003

	ham = defineHamiltonian(n,npi)
	greens = getGreens(ham,e0,gamma)
	specFun = getSpectralFunction(greens)

	im(specFun)