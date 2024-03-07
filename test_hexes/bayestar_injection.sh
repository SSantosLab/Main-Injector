#!/bin/bash
lalapps_inspinj \
	`# writes output to named xml file` \
	--output $4/inj.xml \
	`# lower frequency bound (Hz)` \
	--f-lower 25 \
	--source-file $3 \
	`# type of waveform model used in sim` \
	--waveform TaylorF2threePointFivePN \
	\
	`# Time Distribution: One injection is made every time step until gps time ends.` \
	--t-distr uniform --time-step 1 --gps-start-time $1 --gps-end-time $2 \
	\
	`# Mass Distribution: Choices from totalMass (uniform in total mass), componentMass (uniform in m1/m2), gaussian, log (log dist in component mass), totalMassRatio (uniform in total and mass ratio), logTotalMassUniformMassRatio (log total, uniform mass ratio), totalMassFraction (uniform in total and m1/(m1+m2)), fixMasses (fixed values)` \
	`# --m-distr fixMasses --fixed-mass1 1.4 --fixed-mass2 1.4` \
	`# --m-distr componentMass --min-mass1 1.2 --max-mass1 1.4 --min-mass2 1.2 --max-mass2 1.4` \
	`# --m-distr totalMassRatio --min-mtotal 2.2 --max-mtotal 2.8 --min-mratio 0.85 --max-mratio 1.0` \
	--m-distr gaussian --mean-mass1 1.4 --stdev-mass1 0.1 --mean-mass2 1.4 --stdev-mass2 0.1 --min-mass1 1.2 --max-mass1 1.6 --min-mass2 1.2 --max-mass2 1.6 \
	\
	`# Distance Distribution: Choices from uniform (d), distancesquared (d^2), volume (d^3), log10 (log10(d)), or source file. Source file format unclear. Distances in kpc.` \
	`# --d-distr volume --min-distance 1 --max-distance 600e3` \
	--d-distr source \
	\
	`# Localization Distribution: Choices from random (isotropic), fixed (set using --longitude (ra) and --latitude (dec) in degrees), or from source file. Source file format unclear.` \
	--l-distr source \
	`# --l-distr fixed --longitude 8.294086e-01 --latitude -8.197192e-01` \
	\
	`# Inclination Distribution: Choices from uniform (arccos(i)), gaussian in i, or fixed i. Angles in degrees.` \
	--i-distr uniform \
	`# max angle in degrees, for use with uniform distr` \
	`# --max-inc` \
	`# --i-distr fixed --fixed-inc 20` \
	`# --i-distr gaussain --incl-std 2` \
	`# --polarization 30` \
	\
	`# Spin Distribution: Different from other distr args, choice between guassian or uniform.` \
	`# --disable-spin` \
	`# aligned forces the spins to align to the orbital angular momentum` \
	--enable-spin --min-spin1 0.0 --max-spin1 0.02 --min-spin2 0.0 --max-spin2 0.02 --aligned \
	`# --enable-spin --spin-gaussian --mean-spin1 0.5 --stdev-spin1 0.1 --mean-spin2 0.5 --stdev-spin2 0.1` \
	--disable-milkyway

bayestar-sample-model-psd \
	`# Write output to psd.xml` \
	-o $4/psd.xml \
	`# Specify noise models for desired detectors. The ones used here are for O4 design sensitivity, based on LIGO tech report T1800545.` \
	--H1=aLIGO175MpcT1800545 --L1=aLIGO175MpcT1800545 --V1=aLIGOAdVO4T1800545

bayestar-realize-coincs \
	`# Write output to coinc.xml` \
	-o $4/coinc.xml \
	`# Use the injections and noise PSDs that we generated` \
	$4/inj.xml --reference-psd $4/psd.xml \
	`# Specify which detectors are in science mode` \
	--detector H1 L1 V1 \
	`# Optionally, add Gaussian noise (rather than zero noise)` \
	--measurement-error gaussian-noise \
	`# Optionally, adjust the detection threshold: single-detector SNR, network SNR, and minimum number of detectors above threshold to form a coincidence.` \
	--snr-threshold 2.0 \
	--net-snr-threshold 8.0 \
	--min-triggers 1 \
	`# Optionally, save triggers that were below the single-detector threshold` \
	`# --keep-subthreshold` \
    -l CRITICAL

bayestar-localize-coincs -o $4 $4/coinc.xml -l CRITICAL