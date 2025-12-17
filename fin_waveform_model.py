"""
QW-1526: FIN Theory Waveform Model Wrapper for Bilby
====================================================

This module defines a source model function compatible with Bilby.
It wraps standard LALSimulation waveforms (e.g., IMRPhenomD)
and applies the FIN amplitude scaling correction:
    h_fin = h_gr * (1 / DL^(n-1))

Since h_gr already contains 1/DL, calculating h_gr(DL=1) * 1/DL^n
is equivalent to taking standard output and multiplying by:
    Correction = (1/DL^n) / (1/DL) = DL^(1-n) = 1 / DL^(n-1)

Usage in Bilby:
    waveform_generator = bilby.gw.WaveformGenerator(
        frequency_domain_source_model=fin_amplitude_scaled_source,
        ...
    )
"""

import numpy as np
import bilby
import lalsimulation as lalsim
import lal

def fin_amplitude_scaled_source(frequency_array, mass_1, mass_2, luminosity_distance,
                                a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl,
                                theta_jn, phase, n_scaling=1.0, **kwargs):
    """
    Frequency domain source model with modified amplitude scaling.
    
    Parameters
    ----------
    n_scaling : float
        The FIN scaling exponent. 
        n=1.0 corresponds to GR.
        n=0.66 corresponds to Active Vacuum.
    
    All other parameters are standard Bilby/LAL source parameters.
    """
    
    # 1. Generate Standard GR Waveform
    # We use IMRPhenomD as the baseline (robust, available)
    # Note: We pass all parameters except 'n_scaling' to the standard model
    
    # Use standard Bilby wrapper for IMRPhenomD
    # We call the lal_binary_black_hole source model from bilby
    
    waveform_kwargs = kwargs.copy()
    
    # Generate GR polarizations (dict: 'plus', 'cross')
    # This assumes intrinsic distance scaling of 1/luminosity_distance
    gr_polarizations = bilby.gw.source.lal_binary_black_hole(
        frequency_array=frequency_array,
        mass_1=mass_1,
        mass_2=mass_2,
        luminosity_distance=luminosity_distance,
        a_1=a_1, tilt_1=tilt_1, phi_12=phi_12,
        a_2=a_2, tilt_2=tilt_2, phi_jl=phi_jl,
        theta_jn=theta_jn, phase=phase,
        **waveform_kwargs
    )
    
    if gr_polarizations is None:
        return None

    # 2. Apply FIN Scaling Correction
    # GR Amplitude is proportional to 1/DL
    # FIN Amplitude is proportional to 1/DL^n
    # Correction Factor = (1/DL^n) / (1/DL) = DL^(1-n)
    
    # Avoid division by zero if DL is weird, though priors handle bounds
    if luminosity_distance <= 0:
        return gr_polarizations # Should not happen
        
    correction_factor = luminosity_distance**(1.0 - n_scaling)
    
    fin_polarizations = {}
    for mode in gr_polarizations:
        fin_polarizations[mode] = gr_polarizations[mode] * correction_factor
        
    return fin_polarizations
