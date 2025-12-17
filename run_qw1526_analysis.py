#!/usr/bin/env python3
"""
QW-1526: RUNNER SCRIPT
======================
Executes the reanalysis of GW150914 with free amplitude scaling parameter.

Usage:
    python3 run_qw1526_analysis.py --label fin_test_run --mock

Arguments:
    --mock : Use simulated data (Gaussian noise + Injection) instead of real LIGO data.
             Use this for verification if you don't have GWOSC access.
"""

import argparse
import numpy as np
import bilby
from fin_waveform_model import fin_amplitude_scaled_source

# Parse Args
parser = argparse.ArgumentParser()
parser.add_argument("--label", type=str, default="qw1526_test", help="Output label")
parser.add_argument("--mock", action="store_true", help="Run with mock data injection")
parser.add_argument("--npoints", type=int, default=500, help="Number of live points")
args = parser.parse_args()

# Setup Logger
bilby.core.utils.setup_logger(outdir=args.label, label=args.label)

# ==============================================================================
# 1. SETUP PRIORS
# ==============================================================================

priors = bilby.gw.prior.BBHPriorDict(filename='GW150914.prior')

# Override for QW-1526 specific configuration
# Standard 4s priors
priors['geocent_time'] = bilby.core.prior.Uniform(
    minimum=1126259462.3, maximum=1126259462.5, name='geocent_time', latex_label='$t_c$', unit='$s$')

# ADD NEW PARAMETER: n_scaling
priors['n_scaling'] = bilby.core.prior.Uniform(
    minimum=0.5, maximum=1.5, name='n_scaling', latex_label='$n$')

# ==============================================================================
# 2. WAVEFORM GENERATOR
# ==============================================================================

waveform_arguments = dict(
    waveform_approximant='IMRPhenomD',
    reference_frequency=50.,
    minimum_frequency=20.
)

waveform_generator = bilby.gw.WaveformGenerator(
    duration=4.,
    sampling_frequency=4096.,
    frequency_domain_source_model=fin_amplitude_scaled_source,
    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
    waveform_arguments=waveform_arguments
)

# ==============================================================================
# 3. DATA SETUP (INTERFEROMETERS)
# ==============================================================================

ifos = bilby.gw.detector.InterferometerList(['H1', 'L1'])

if args.mock:
    # MOCK DATA INJECTION
    bilby.core.utils.logger.info("Running in MOCK mode with Injection")
    
    # Injection Parameters (FIN Theory Hypothesis)
    injection_parameters = dict(
        mass_1=36., mass_2=29., a_1=0., a_2=0., tilt_1=0., tilt_2=0.,
        phi_12=0., phi_jl=0., theta_jn=0.4, psi=2.659, phase=1.3,
        geocent_time=1126259462.42, ra=1.375, dec=-1.2108,
        luminosity_distance=500., # Mpc
        n_scaling=0.66            # FIN SCALING
    )
    
    ifos.set_strain_data_from_power_spectral_densities(
        sampling_frequency=4096., duration=4.,
        start_time=1126259462.42 - 2)
        
    ifos.inject_signal(waveform_generator=waveform_generator,
                       parameters=injection_parameters)
else:
    # REAL DATA (Requires GWOSC access or local frame files)
    bilby.core.utils.logger.info("Running with REAL GW150914 Data")
    
    # Trigger time for GW150914
    trigger_time = 1126259462.4
    
    # Fetch open data
    for ifo in ifos:
        ifo.set_strain_data_from_gwosc(
            overlap=0,
            duration=4.,
            start_time=trigger_time - 2,
            end_time=trigger_time + 2,
            channel_name=f'{ifo.name}:GWOSC-4KHZ_R1_STRAIN' # Adjust channel if needed
        )

# ==============================================================================
# 4. LIKELIHOOD & SAMPLER
# ==============================================================================

likelihood = bilby.gw.GravitationalWaveTransient(
    interferometers=ifos,
    waveform_generator=waveform_generator
)

# Run Sampler
result = bilby.run_sampler(
    likelihood=likelihood,
    priors=priors,
    sampler='dynesty',
    npoints=args.npoints,
    injection_parameters=injection_parameters if args.mock else None,
    outdir=args.label,
    label=args.label,
    conversion_function=bilby.gw.conversion.generate_all_bbh_parameters
)

# Plot
result.plot_corner()
