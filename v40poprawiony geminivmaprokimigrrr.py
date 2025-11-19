#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
v40.4.1 (Faza X.4: Production-Ready Unified Framework)

CRITICAL ARCHITECTURAL FIXES (Unified):
═══════════════════════════════════════════════════════════════════════════════
[NUMERICAL STABILITY]
1. VEV Stabilization (Krok 2.1): Improved clamping and epsilon addition.
2. Beta Epsilon Increase (Krok 2.2): 1e-9 -> 1e-6 for mass enhancement stability.

[PERFORMANCE OPTIMIZATION]
3. Fine-tuning Early Stop (Krok 3.1): Breaks fine-tuning loop if loss < 0.1.
4. TPU Batch Alignment (Krok 3.2): Batch size adjusted to be a multiple of 8.
5. Memory Cleanup (Krok 3.3): Added gc.collect() after epoch loops.

[CODE QUALITY & STANDARDIZATION]
6. Global Constants (Krok 4.1): Centralized weights and clip thresholds.
7. Input Validation (Krok 4.2): Added required key check in total_energy.
8. Cache Hash Update (Krok 5.2): k_mass included in hash.
═══════════════════════════════════════════════════════════════════════════════
"""

print("="*80)
print(" INITIALIZING v40.4.1 (Production Hardened) ")
print("="*80)
print("🔥 CRITICAL FIXES APPLIED:")
print("   ✅ Stationary state solver (t=0.5)")
print("   ✅ L2 normalization constraint")
print("   ✅ Boundary condition enforcement")
print("   ✅ Scaled Laplacian (ε=dr*0.1)")
print("   ✅ Stable VEV calculation & mass enhancement")
print("   ✅ Early stopping and TPU batch optimization")
print("   ✅ Memory cleanup (gc.collect)")
print("="*80)

EXECUTION_MODE = 'FULL_RUN'
DEVICE_MODE = 'CPU'  # 'AUTO', 'TPU', 'GPU', 'CPU'

# ==============================================================================
# IMPORTS AND ENVIRONMENT
# ==============================================================================
import os, sys, time, warnings, subprocess, gc, json, hashlib, pickle
import numpy as np
import pandas as pd
import scipy, scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.optimize import minimize
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import threading
from contextlib import nullcontext
import glob
from datetime import datetime

warnings.filterwarnings("ignore", category=RuntimeWarning)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['XLA_USE_BF16'] = '1'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.5"

# Core ML
import torch
import torch.nn as nn
from torch.optim import Adam, AdamW, SGD
from torch.optim.lr_scheduler import ReduceLROnPlateau, LambdaLR
from torch.utils.data import TensorDataset, DataLoader

# Optional dependencies
try:
    import optuna
    from optuna.samplers import NSGAIISampler
    OPTUNA_AVAILABLE = True
except ImportError:
    OPTUNA_AVAILABLE = False
    print("⚠️ Optuna not available - install via pip")

try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

try:
    from torch.amp import autocast
    AUTOCAST_AVAILABLE = True
except ImportError:
    AUTOCAST_AVAILABLE = False

# ==============================================================================
# DEVICE SETUP (Optimized for Kaggle TPU)
# ==============================================================================
print("\n[INFO] Device detection and setup...")
IS_TPU, GPU_MODE = False, False
xp = np  # Default to NumPy
xm = None
if DEVICE_MODE == 'TPU':
    try:
        import torch_xla.core.xla_model as xm
        import torch_xla
        device = xm.xla_device()
        IS_TPU = True
        print(f"✅ TPU forced: {device}")
    except:
        raise RuntimeError("TPU forced but torch_xla not available")

elif DEVICE_MODE == 'GPU' and torch.cuda.is_available():
    device = torch.device("cuda")
    GPU_MODE = True
    torch.backends.cudnn.benchmark = True
    try:
        import cupy as cp
        xp = cp
        print("✅ GPU + CuPy enabled")
    except:
        print("⚠️ GPU mode but CuPy not available")

elif DEVICE_MODE == 'CPU':
    device = torch.device("cpu")
    print("✅ CPU mode forced")

else:  # AUTO mode
    try:
        import torch_xla.core.xla_model as xm
        device = xm.xla_device()
        IS_TPU = True
        print(f"✅ Auto-detected TPU: {device}")
    except:
        if torch.cuda.is_available():
            device = torch.device("cuda")
            GPU_MODE = True
            print(f"✅ Auto-detected GPU: {device}")
        else:
            device = torch.device("cpu")
            print("✅ Auto-detected CPU")

# ==============================================================================
# CONSTANTS (FIX #4.1 & #5.0: Centralized magic numbers)
# ==============================================================================
STATIONARY_TIME = 0.5  # FIX #1: Stationary state time
GRAD_CLIP_PRETRAIN = 5.0
GRAD_CLIP_FINETUNE = 10.0
EPSILON_BETA = 1e-6  # Krok 2.2: Increased epsilon for mass_enhancement
VEV_MIN_CLAMP = 0.1 # Krok 2.1: Min value for VEV target calculation
VEV_PENALTY_WEIGHT = 50.0
BOUNDARY_PENALTY_WEIGHT = 100.0
NORMALIZATION_WEIGHT = 500.0
PHYSICAL_PENALTY_WEIGHT = 50.0

HASH_PRECISION = 8
CACHE_SALT = "_v40_4_stable"

# ==============================================================================
# GLOBAL PARAMETERS (Production-Ready)
# ==============================================================================
print("\n[INFO] Setting global simulation parameters...")

# Grid resolution (optimized for TPU memory)
Nr, Nr_theta, Nr_phi = 128, 16, 16  # Higher resolution for accuracy
Nr_theta_mesh, Nr_phi_mesh = 8, 8   # Reduced for speed in fine-tuning

# Physics constants
r_max = 25.0
num_octaves = 12
lambda_H = 0.5  # Centralized (FIX #5.1)
delta = 0.2  # Sextic term coefficient
beta_hierarchy = 0.15  # Legacy parameter (kept for compatibility)

# Target masses (GeV)
TARGET_TOP_GEV = 173.1
TARGET_HIGGS_GEV = 125.1
TARGET_ELECTRON_GEV = 0.000511
TARGET_MUON_GEV = 0.10566
TARGET_TAU_GEV = 1.77686

# Numerical settings
sigma_noise = 0.05
neigs = 300
tol_energy = 1e-8
clip_value = 1e4
m0_clip_max = 50.0

# Batch sizes (optimized)
FINETUNE_BATCH_SIZE = 1048576  # 1M for TPU
INFERENCE_BATCH_SIZE = 4194304  # 4M for TPU
PRETRAIN_BATCH_SIZE = 32768     # 32K for stability

# Grid arrays
r_cpu = np.linspace(1e-6, r_max, Nr, dtype=np.float64)
dr_cpu = r_cpu[1] - r_cpu[0]
t_steps_mesh = 1  # FIX #1: Stationary state (was 15)

# Cache system
CACHE_DIR = '/kaggle/working/finetune_cache_v40'
os.makedirs(CACHE_DIR, exist_ok=True)
# Optuna settings
n_initial_trials = 3000
# Logging
LOG_CSV_FILE = "corr_log_v40_unified.csv"
LOG_DIR = './logs'
HEARTBEAT_FILE = 'heartbeat.log'
KAGGLE_OUTPUT_FILE = '/kaggle/working/tpu_active.log'

print(f"[PARAM] Grid: {Nr}x{Nr_theta}x{Nr_phi} (Δr={dr_cpu:.4f})")
print(f"[PARAM] Device: {device} | IS_TPU: {IS_TPU} | GPU: {GPU_MODE}")

# ==============================================================================
# GENERATION MAPPING (Faza X.2 Logic)
# ==============================================================================
def get_gen_idx_and_scale(octave_num, mass_scale_mu=15.0, mass_scale_tau=75.0):
    """
    Maps octaves to generations with parameterized mass scales.
    - Octaves 0-3: Gen 1 (electron) - 4 octaves
    - Octaves 4-6: Gen 2 (muon) - 3 octaves
    - Octaves 7-11: Gen 3 (tau) - 5 octaves
    """
    if octave_num < 4:  # Octaves 0-3
        return 0, 1.0
    elif octave_num < 7:  # Octaves 4-6
        return 1, mass_scale_mu
    else:  # Octaves 7-11
        return 2, mass_scale_tau

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
def tpu_print(msg, **kwargs):
    """Thread-safe printing for TPU environment"""
    if IS_TPU:
        try:
            sys.stdout.write(msg + kwargs.get('end', '\n'))
            sys.stdout.flush()
        except:
            print(msg, **kwargs)
    else:
        print(msg, **kwargs)

def write_dual_heartbeat(message):
    """Writes heartbeat to multiple files for Kaggle"""
    timestamp = time.time()
    for filepath in [HEARTBEAT_FILE, KAGGLE_OUTPUT_FILE]:
        try:
            with open(filepath, 'a') as f:
                f.write(f"{timestamp} {message}\n")
                f.flush()
        except: pass

def safe_loss_check_and_sync(loss, epoch, batch_idx):
    """Check for finite loss and sync TPU"""
    if IS_TPU:
        try:
            loss_cpu = loss.detach().cpu()
            xm.mark_step()
            loss_val = loss_cpu.item()
            if not np.isfinite(loss_val):
                tpu_print(f"⚠️ Non-finite loss: {loss_val} at E{epoch}B{batch_idx}")
                return True, loss_val
            return False, loss_val
        except RuntimeError as e:
            if "tensor_data" in str(e):
                return True, float('inf')
            raise e
    else:
        loss_val = loss.item()
        return not torch.isfinite(loss), loss_val

def get_dynamic_batch_size(epoch, device_type='CPU'):
    """Adaptive batch sizing for stability"""
    if epoch < 5:
        return 8192 if device_type == 'TPU' else 4096
    elif epoch < 15:
        return 16384 if device_type == 'TPU' else 8192
    else:
        return 32768 if device_type == 'TPU' else 16384

def print_progress_bar(epoch, batch, total_batches, loss):
    """Visual progress indicator"""
    percent = int((batch / total_batches) * 50)
    bar = '█' * percent + '░' * (50 - percent)
    print(f"\rE{epoch} [{bar}] {batch}/{total_batches} L={loss:.2e}", end='', flush=True)

# ==============================================================================
# CACHE SYSTEM (v40.4 Enhanced)
# ==============================================================================
def get_params_hash_v40(params_v40):
    """
    Creates MD5 hash of all 20+ parameters for cache lookup. (FIX #5.2)
    """
    keys_to_hash = [
        'beta_max', 'A_dip', 'o_dip', 'sigma_dip',
        'A_dip2', 'o_dip2', 'sigma_dip2',
        'A_k', 'omega_k', 'phi_k', 'alpha_geo_k',
        'g_Y_gen1', 'g_Y_gen2', 'g_Y_gen3', 'k_inv',
        'mu2', 'm0', 'g', 'delta', 'k_mass',  # <-- k_mass DODANE
        'lambda_Y_tau', 'mass_scale_mu', 'mass_scale_tau'
    ]

    # Używamy HASH_PRECISION
    rounded_values = [f"{params_v40.get(key, 0.0):.{HASH_PRECISION}f}" for key in keys_to_hash]
    params_str = "_".join(rounded_values) + CACHE_SALT
    return hashlib.md5(params_str.encode()).hexdigest()[:16]

def save_finetune_result_v40(params_v40, Psi_out, Phi_out, final_loss, elapsed_time):
    """Saves fine-tuning result to cache with error handling"""
    params_hash = get_params_hash_v40(params_v40)
    cache_file = os.path.join(CACHE_DIR, f'ft_{params_hash}.pkl')

    try:
        cache_data = {
            'params_v40': params_v40,
            'Psi_out': Psi_out if isinstance(Psi_out, np.ndarray) else Psi_out.cpu().numpy(),
            'Phi_out': Phi_out if isinstance(Phi_out, np.ndarray) else Phi_out.cpu().numpy(),
            'final_loss': final_loss,
            'elapsed_time': elapsed_time,
            'timestamp': time.time()
        }
        with open(cache_file, 'wb') as f:
            pickle.dump(cache_data, f)
        tpu_print(f"  💾 Cache saved: {params_hash} ({os.path.getsize(cache_file)/(1024**2):.2f} MB)")
    except Exception as e:
        tpu_print(f"  ⚠️ Cache save error: {e}")

def load_finetune_result_v40(params_v40):
    """Loads cached result if exists and is valid"""
    params_hash = get_params_hash_v40(params_v40)
    cache_file = os.path.join(CACHE_DIR, f'ft_{params_hash}.pkl')
    tpu_print(f"  🔍 Cache check: {params_hash}")

    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                cache_data = pickle.load(f)
            Psi_out, Phi_out = cache_data['Psi_out'], cache_data['Phi_out']

            if np.all(np.isfinite(Psi_out)) and np.all(np.isfinite(Phi_out)):
                age_min = (time.time() - cache_data['timestamp']) / 60
                tpu_print(f"  ♻️ Cache HIT! Age: {age_min:.1f}min, Loss: {cache_data['final_loss']:.4e}")
                return Psi_out, Phi_out
            else:
                tpu_print("  ⚠️ Corrupted cache, recalculating...")
                os.remove(cache_file)
        except Exception as e:
            tpu_print(f"  ⚠️ Cache load error: {e}")
            if os.path.exists(cache_file):
                os.remove(cache_file)
    else:
        tpu_print(f"  ❌ Cache MISS")

    return None, None

def list_cache_stats():
    """Display cache statistics"""
    cache_files = glob.glob(os.path.join(CACHE_DIR, 'ft_*.pkl'))
    if not cache_files:
        tpu_print("[CACHE] Empty")
        return

    total_size = sum(os.path.getsize(f) for f in cache_files)
    tpu_print(f"[CACHE] {len(cache_files)} files, {total_size/(1024**2):.2f} MB")

    recent = sorted(cache_files, key=os.path.getmtime, reverse=True)[:3]
    for f in recent:
        age_min = (time.time() - os.path.getmtime(f)) / 60
        print(f"    - {os.path.basename(f)}: {age_min:.1f}min old")

# ==============================================================================
# PHYSICS FUNCTIONS (v40.4 - Core Logic)
# ==============================================================================
def beta_topo_gaussian_dip(o, beta_max, A_dip, o_dip, sigma_dip, A_dip2, o_dip2, sigma_dip2):
    """
    Double-valley topological potential:
    - Valley 1 (o~3-5): Force hierarchy
    - Valley 2 (o~7-10): Tau mass enhancement
    """
    dip1 = A_dip * np.exp(-(o - o_dip)**2 / (2 * sigma_dip**2))
    dip2 = A_dip2 * np.exp(-(o - o_dip2)**2 / (2 * sigma_dip2**2))
    return beta_max - dip1 - dip2

def compute_interaction_kernel(i, j, params, use_fast_kernel=True):
    """Universal kernel with topological modulation"""
    d = abs(i - j)
    if d not in [1, 2]: return 0.0

    beta_topo = beta_topo_gaussian_dip(
        min(i,j), params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
        params['A_dip2'], params['o_dip2'], params['sigma_dip2']
    )

    A, alpha_geo = params['A_k'], params['alpha_geo_k']
    mod_diff = abs((i%3)-(j%3))

    if use_fast_kernel:
        K_ij_raw = A * (2**(-alpha_geo * d)) * np.exp(-beta_topo * mod_diff)
    else:
        omega, phi = params['omega_k'], params['phi_k']
        K_ij_raw = A * np.cos(omega * d + phi) * (2**(-alpha_geo * d)) * np.exp(-beta_topo * mod_diff)

    return K_ij_raw if d == 1 else K_ij_raw / 2.0

def compute_force_hierarchy_v40(beta_profile, k_inv):
    """
    Inverted force law: g_i ∝ 1/β^k
    Natural explanation: High topological stiffness → weaker coupling
    """
    o_u1, o_su2, o_su3 = 0, 2, 4
    beta_u1, beta_su2, beta_su3 = beta_profile[o_u1], beta_profile[o_su2], beta_profile[o_su3]
    eps = 1e-9
    g1 = 1.0 / (beta_u1**k_inv + eps)
    g2 = 1.0 / (beta_su2**k_inv + eps)
    g3 = 1.0 / (beta_su3**k_inv + eps)
    return g1, g2, g3

def radial_laplacian(field, r, dr, xp):
    """
    Radial Laplacian with L'Hospital's rule at origin.
    Critical for numerical stability near r=0.
    """
    dfield_dr = xp.gradient(field, dr)
    d2field_dr2 = xp.gradient(dfield_dr, dr)

    lap = xp.zeros_like(field)
    if len(field) > 1:
        lap[0] = 3.0 * d2field_dr2[0]  # L'Hospital for r→0
        lap[1] = 3.0 * d2field_dr2[1]
    elif len(field) == 1:
        lap[0] = 3.0 * d2field_dr2[0]

    r_safe = xp.where(r > 1e-9, r, 1e-9)
    lap[2:] = d2field_dr2[2:] + (2.0 / r_safe[2:]) * dfield_dr[2:]
    return lap

def total_energy_v40(Psi, Phi_H, params, r, dr, xp):
    """
    Complete energy functional WITH k_mass enhancement
    (Krok 4.2: Walidacja wejść)
    """
    required_keys = ['m0', 'g', 'delta', 'mu2', 'g_Y_gen1', 'g_Y_gen2', 'g_Y_gen3',
                     'lambda_Y_tau', 'mass_scale_mu', 'mass_scale_tau',
                     'beta_max', 'A_dip', 'o_dip', 'sigma_dip',
                     'A_dip2', 'o_dip2', 'sigma_dip2', 'k_mass']

    for key in required_keys:
        if key not in params:
            raise ValueError(f"Missing required parameter: {key}")

    m0, g, delta = params['m0'], params['g'], params['delta']
    g_Yukawa_vec = [params['g_Y_gen1'], params['g_Y_gen2'], params['g_Y_gen3']]
    mu2, lambda_H_val = params['mu2'], lambda_H # Używa stałej globalnej lambda_H
    k_mass = params.get('k_mass', 1.0)

    energy_density = xp.zeros(len(r), dtype=Psi.dtype)

    for o in range(num_octaves):
        gen_idx, mass_scale = get_gen_idx_and_scale(o, params['mass_scale_mu'], params['mass_scale_tau'])

        # Oblicz wzmocnienie masowe z β_topo
        beta_val = beta_topo_gaussian_dip(
            o, params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
            params['A_dip2'], params['o_dip2'], params['sigma_dip2']
        )
        # Krok 2.2: Stabilizacja: używamy EPSILON_BETA
        mass_enhancement = (1.0 / (abs(beta_val) + EPSILON_BETA)) ** k_mass

        dpsi = xp.gradient(Psi[o], dr)
        psi_sq = Psi[o]**2

        # KWADRATOWY CZŁON MASOWY z wzmocnieniem k_mass
        energy_density += 0.5*dpsi**2 + 0.5*m0**2*psi_sq * mass_enhancement + 0.25*g*psi_sq**2 + 0.125*delta*psi_sq**3

        # YUKAWA z podwójnym wzmocnieniem: mass_scale * k_mass
        yukawa_term = 0.5 * g_Yukawa_vec[gen_idx] * (Phi_H**2) * psi_sq * mass_scale * mass_enhancement

        # SEXTIC Yukawa dla tau (gen_idx==2)
        if gen_idx == 2:
            yukawa_term += 0.125 * params['lambda_Y_tau'] * (Phi_H**2) * (psi_sq**2)

        energy_density += yukawa_term

    # Inter-octave couplings
    for i in range(num_octaves):
        for j in range(i + 1, num_octaves):
            K_ij = compute_interaction_kernel(i, j, params, use_fast_kernel=True)
            if K_ij != 0.0:
                energy_density += K_ij * Psi[i] * Psi[j]

    # Higgs terms
    dPhi = xp.gradient(Phi_H, dr)
    E_kin_H = 0.5 * dPhi**2
    E_pot_H = 0.5 * mu2 * Phi_H**2 + 0.25 * lambda_H_val * Phi_H**4

    integrand = energy_density + E_kin_H + E_pot_H
    return 4.0 * xp.pi * xp.sum(integrand * r**2) * dr

def functional_derivative_v40(Psi, Phi_H, params, r, dr, xp):
    """Functional derivative WITH k_mass enhancement"""
    m0, g, delta = params['m0'], params['g'], params['delta']
    g_Yukawa_vec = [params['g_Y_gen1'], params['g_Y_gen2'], params['g_Y_gen3']]
    mu2, lambda_H_val = params['mu2'], lambda_H
    k_mass = params.get('k_mass', 1.0)

    dE_Psi = xp.zeros_like(Psi)

    for o in range(num_octaves):
        gen_idx, mass_scale = get_gen_idx_and_scale(o, params['mass_scale_mu'], params['mass_scale_tau'])

        # Wzmocnienie masowe dla pochodnych
        beta_val = beta_topo_gaussian_dip(
            o, params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
            params['A_dip2'], params['o_dip2'], params['sigma_dip2']
        )
        # Krok 2.2: Stabilizacja: używamy EPSILON_BETA
        mass_enhancement = (1.0 / (abs(beta_val) + EPSILON_BETA)) ** k_mass

        # Laplacian term
        lap = -radial_laplacian(Psi[o], r, dr, xp)

        # POTENTIAL DERIVATIVES z k_mass
        mass_term = m0**2 * Psi[o] * mass_enhancement
        nonlin = g * Psi[o]**3
        sextic_term = 0.75 * delta * (Psi[o]**5)

        # YUKAWA DERIVATIVE z k_mass
        yukawa_term = g_Yukawa_vec[gen_idx] * (Phi_H**2) * Psi[o] * mass_scale * mass_enhancement
        if gen_idx == 2:
            yukawa_term += 0.5 * params['lambda_Y_tau'] * (Phi_H**2) * (Psi[o]**3)

        # Couplings
        coupling = xp.zeros_like(Psi[o])
        for j in range(num_octaves):
            if o == j: continue
            K_ij = compute_interaction_kernel(o, j, params, use_fast_kernel=True)
            coupling += K_ij * Psi[j]

        dE_Psi[o] = lap + mass_term + nonlin + sextic_term + coupling + yukawa_term

    # Higgs derivative
    lap_Phi = -radial_laplacian(Phi_H, r, dr, xp)

    yukawa_phi_sum = 0
    for o in range(num_octaves):
        gen_idx, mass_scale = get_gen_idx_and_scale(o, params['mass_scale_mu'], params['mass_scale_tau'])

        # Dla Higgsa też potrzebne wzmocnienie!
        beta_val = beta_topo_gaussian_dip(
            o, params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
            params['A_dip2'], params['o_dip2'], params['sigma_dip2']
        )
        # Krok 2.2: Stabilizacja: używamy EPSILON_BETA
        mass_enhancement = (1.0 / (abs(beta_val) + EPSILON_BETA)) ** k_mass

        psi_sq = Psi[o]**2
        yukawa_phi_sum += g_Yukawa_vec[gen_idx] * psi_sq * mass_scale * mass_enhancement

        if gen_idx == 2:
            yukawa_phi_sum += 0.25 * params['lambda_Y_tau'] * (psi_sq**2) * mass_enhancement

    dE_Phi = lap_Phi + mu2 * Phi_H + lambda_H_val * (Phi_H**3) + yukawa_phi_sum

    return dE_Psi, dE_Phi

def diagonalize_v40(Psi_loc, Phi_loc, params):
    """
    Diagonalizacja z poprawnym uwzględnieniem różnych mas generacyjnych
    """
    Nfull = num_octaves * Nr
    neigs_to_calc = max(1, min(neigs, Nfull - 2))
    k_mass = params.get('k_mass', 1.0)

    H = np.zeros((Nfull, Nfull), dtype=float)

    # Diagonal blocks
    for o in range(num_octaves):
        idx0 = o * Nr
        gen_idx, mass_scale = get_gen_idx_and_scale(o, params['mass_scale_mu'], params['mass_scale_tau'])

        # Wzmocnienie dla Hesjanu
        beta_val = beta_topo_gaussian_dip(
            o, params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
            params['A_dip2'], params['o_dip2'], params['sigma_dip2']
        )
        # Krok 2.2: Stabilizacja: używamy EPSILON_BETA
        mass_enhancement = (1.0 / (abs(beta_val) + EPSILON_BETA)) ** k_mass

        # DRUGA POCHODNA POTENCJAŁU z wzmocnieniem
        yukawa_diag = params[f'g_Y_gen{gen_idx+1}'] * (Phi_loc**2) * mass_scale * mass_enhancement
        if gen_idx == 2:
            yukawa_diag += 1.5 * params['lambda_Y_tau'] * (Phi_loc**2) * (Psi_loc[o]**2) * mass_enhancement

        # UWAGA: m0² też musi być wzmocnione!
        V_pp = params['m0']**2 * mass_enhancement + 3*params['g']*Psi_loc[o]**2 + 3.75*params['delta']*Psi_loc[o]**4 + yukawa_diag

        # Fill matrix
        for i in range(Nr):
            H[idx0+i, idx0+i] = (2/dr_cpu**2) + V_pp[i]
        for i in range(Nr - 1):
            off = -1/dr_cpu**2
            H[idx0+i, idx0+i+1] = off
            H[idx0+i+1, idx0+i] = off

    # Off-diagonal couplings
    for i in range(num_octaves):
        for j in range(i + 1, num_octaves):
            K_ij = compute_interaction_kernel(i, j, params, use_fast_kernel=True)
            if K_ij != 0.0:
                idx_i, idx_j = i * Nr, j * Nr
                for k_r in range(Nr):
                    H[idx_i + k_r, idx_j + k_r] = K_ij
                    H[idx_j + k_r, idx_i + k_r] = K_ij

    # Sparse diagonalization
    H_sparse = sp.csr_matrix(H)
    try:
        w, _ = spl.eigsh(H_sparse, k=neigs_to_calc, which='SA', tol=1e-6)
        positive_w = w[w > 1e-12]
        if len(positive_w) == 0:
            return np.array([])
        return np.sqrt(np.sort(positive_w))
    except Exception as e:
        tpu_print(f"  [DIAG ERROR] Diagonalization failed: {e}")
        return None

# ==============================================================================
# NEURAL NETWORK MODEL
# ==============================================================================
class ResidualBlock(nn.Module):
    def __init__(self, size):
        super().__init__()
        self.l1 = nn.Linear(size, size)
        self.l2 = nn.Linear(size, size)
        self.act = nn.GELU()

    def forward(self, x):
        return self.act(self.l2(self.act(self.l1(x))) + x)

class SolitonPINN(nn.Module):
    def __init__(self, output_size=num_octaves+1):
        super().__init__()
        self.num_octaves = num_octaves
        self.inp = nn.Linear(4, 128)
        self.bn1 = nn.LayerNorm(128)
        self.act = nn.GELU()
        self.blocks = nn.Sequential(*[ResidualBlock(128) for _ in range(3)])
        self.out = nn.Linear(128, output_size)

        # Initialize weights
        nn.init.xavier_uniform_(self.inp.weight)
        nn.init.xavier_uniform_(self.out.weight)
        nn.init.zeros_(self.out.bias)

    def forward(self, x):
        x = self.act(self.bn1(self.inp(x)))
        out = self.out(self.blocks(x))
        # Explicit: first num_octaves channels are Ψ, last is Φ
        return out

class PrefetchDataLoader:
    """Data loader with device prefetching"""
    def __init__(self, loader, device):
        self.loader = loader
        self.device = device

    def __iter__(self):
        batch = None
        for next_batch in self.loader:
            if batch is not None:
                yield batch
            # Ensure non_blocking is only used for CUDA, or use default behavior
            non_blocking = self.device.type == 'cuda'
            batch = [b.to(self.device, non_blocking=non_blocking) for b in next_batch]
        if batch is not None:
            yield batch

    def __len__(self):
        return len(self.loader)

# ==============================================================================
# LOSS FUNCTION (v40.4.1 - Fixed and Enhanced)
# ==============================================================================
def beta_topo_gaussian_dip_torch(o, beta_max, A_dip, o_dip, sigma_dip, A_dip2, o_dip2, sigma_dip2, device):
    """PyTorch version of double-valley beta profile"""
    o_tensor = torch.tensor(o, dtype=torch.float32, device=device)
    dip1 = A_dip * torch.exp(-(o_tensor - o_dip)**2 / (2 * sigma_dip**2))
    dip2 = A_dip2 * torch.exp(-(o_tensor - o_dip2)**2 / (2 * sigma_dip2**2))
    return beta_max - dip1 - dip2

def compute_interaction_kernel_torch(i, j, params, device):
    """PyTorch kernel computation"""
    d = abs(i - j)
    if d not in [1, 2]: return torch.tensor(0.0, device=device)

    beta_topo = beta_topo_gaussian_dip_torch(
        min(i,j), params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
        params['A_dip2'], params['o_dip2'], params['sigma_dip2'], device
    )

    A_k = torch.tensor(params['A_k'], device=device)
    alpha_geo_k = torch.tensor(params['alpha_geo_k'], device=device)
    d_tensor = torch.tensor(d, dtype=torch.float32, device=device)
    mod_diff = torch.tensor(abs((i%3)-(j%3)), dtype=torch.float32, device=device)

    K_ij_raw = A_k * (2**(-alpha_geo_k * d_tensor)) * torch.exp(-beta_topo * mod_diff)
    return K_ij_raw if d == 1 else K_ij_raw / 2.0

def compute_temporal_derivative_batch_stable(field, coords):
    """
    Stable temporal derivative computation.
    Computes time derivative for EACH output channel separately.
    """
    coords.requires_grad_(True)

    batch_size, num_channels = field.shape

    # Compute time derivative for each channel
    time_derivatives = []

    for i in range(num_channels):
        grad = torch.autograd.grad(
            outputs=field[:, i],
            inputs=coords,
            grad_outputs=torch.ones_like(field[:, i]),
            create_graph=True,
            retain_graph=True
        )[0]

        # Extract time derivative (2nd column, index 1)
        time_derivatives.append(grad[:, 1:2])

    # Concatenate: [batch_size, num_channels]
    return torch.cat(time_derivatives, dim=1)

def pinn_loss_v40(model, r, t, theta, phi, params, current_epoch=0, ic_weight_boost=1.0, batch_idx=0):
    """
    v40.4.1 LOSS FUNCTION with all stabilizing fixes.
    """
    if IS_TPU:
        import torch_xla.core.xla_model as xm
    # Set requires_grad for autograd
    r = r.clone().requires_grad_(True)
    t = t.clone().requires_grad_(True)
    theta = theta.clone().requires_grad_(True)
    phi = phi.clone().requires_grad_(True)

    epsilon = dr_cpu * 0.1

    coords_center = torch.cat([r, t, theta, phi], dim=1)
    coords_plus = torch.cat([r + epsilon, t, theta, phi], dim=1)
    coords_minus = torch.cat([r - epsilon, t, theta, phi], dim=1)

    # Forward passes
    output_center = model(coords_center)
    output_plus = model(coords_plus)
    output_minus = model(coords_minus)

    # Value clipping for stability
    output_center = torch.clamp(output_center, -5.0, 5.0)
    output_plus = torch.clamp(output_plus, -5.0, 5.0)
    output_minus = torch.clamp(output_minus, -5.0, 5.0)

    Psi_center, Phi_center = output_center[:, :-1], output_center[:, -1]

    # Laplacian approximation
    d2field_dr2 = (output_plus - 2 * output_center + output_minus) / (epsilon**2)
    dfield_dr = (output_plus - output_minus) / (2 * epsilon)

    r_safe = torch.clamp(r, min=1e-8)
    lap_field = d2field_dr2 + (2.0 / r_safe) * dfield_dr
    lap_field = torch.nan_to_num(lap_field, nan=0.0, posinf=1e4, neginf=-1e4)

    # Stable temporal derivatives
    dfield_dt = compute_temporal_derivative_batch_stable(output_center, coords_center)

    dPsi_dt = dfield_dt[:, :num_octaves]
    dPhi_dt = dfield_dt[:, num_octaves:]
    lap_Psi = lap_field[:, :num_octaves]
    lap_Phi = lap_field[:, num_octaves:]

    # Physics parameters
    mu2_val = torch.tensor(params['mu2'], device=device)
    m0 = torch.tensor(params['m0'], device=device)
    g = torch.tensor(params['g'], device=device)
    delta = torch.tensor(params['delta'], device=device)
    g_Y_gen1 = torch.tensor(params['g_Y_gen1'], device=device)
    g_Y_gen2 = torch.tensor(params['g_Y_gen2'], device=device)
    g_Y_gen3 = torch.tensor(params['g_Y_gen3'], device=device)
    g_Yukawa_vec = [g_Y_gen1, g_Y_gen2, g_Y_gen3]
    lambda_H_val = torch.tensor(lambda_H, device=device)
    k_mass = torch.tensor(params.get('k_mass', 1.0), device=device)

    # Yukawa coupling to Higgs (Requires mass enhancement factor from beta)
    yukawa_phi_term = torch.zeros_like(Phi_center)
    for o in range(num_octaves):
        gen_idx, mass_scale_np = get_gen_idx_and_scale(o, params['mass_scale_mu'], params['mass_scale_tau'])
        mass_scale = torch.tensor(mass_scale_np, device=device)

        # Obliczanie β_topo
        beta_val = beta_topo_gaussian_dip_torch(
            o, params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
            params['A_dip2'], params['o_dip2'], params['sigma_dip2'], device
        )
        # Krok 2.2: Stabilizacja: używamy EPSILON_BETA
        mass_enhancement = (1.0 / (torch.abs(beta_val) + EPSILON_BETA)) ** k_mass

        # Sumowanie z wzmocnieniem
        psi_sq = Psi_center[:, o]**2
        yukawa_phi_term += g_Yukawa_vec[gen_idx] * psi_sq * mass_scale * mass_enhancement

        if gen_idx == 2:
            lambda_Y_tau = torch.tensor(params['lambda_Y_tau'], device=device)
            yukawa_phi_term += 0.25 * lambda_Y_tau * (psi_sq**2) * mass_scale**2 * mass_enhancement


    # Residuals for Phi
    res_Phi = dPhi_dt - (-lap_Phi - mu2_val*Phi_center - lambda_H_val*Phi_center**3 - yukawa_phi_term*Phi_center)

    res_Psi_list = []
    for o in range(num_octaves):
        gen_idx, mass_scale_np = get_gen_idx_and_scale(o, params['mass_scale_mu'], params['mass_scale_tau'])
        mass_scale = torch.tensor(mass_scale_np, device=device)

        # Obliczanie β_topo
        beta_val = beta_topo_gaussian_dip_torch(
            o, params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
            params['A_dip2'], params['o_dip2'], params['sigma_dip2'], device
        )
        # Krok 2.2: Stabilizacja: używamy EPSILON_BETA
        mass_enhancement = (1.0 / (torch.abs(beta_val) + EPSILON_BETA)) ** k_mass

        # Potencjał z k_mass
        mass_term = m0**2 * Psi_center[:, o] * mass_enhancement
        nonlin = g * Psi_center[:, o]**3
        sextic_term = 0.75 * delta * (Psi_center[:, o]**5)

        # Yukawa z k_mass
        yukawa_term = g_Yukawa_vec[gen_idx] * (Phi_center**2) * Psi_center[:, o] * mass_scale * mass_enhancement
        if gen_idx == 2:
            lambda_Y_tau = torch.tensor(params['lambda_Y_tau'], device=device)
            yukawa_term += 0.5 * lambda_Y_tau * (Phi_center**2) * (Psi_center[:, o]**3) * mass_scale**2 * mass_enhancement

        # Inter-octave couplings
        coupling = torch.zeros_like(Psi_center[:, o])
        for j in range(num_octaves):
            if o == j: continue
            K_ij = compute_interaction_kernel_torch(o, j, params, device)
            coupling += K_ij * Psi_center[:, j]

        res = dPsi_dt[:, o] - (-lap_Psi[:, o] - mass_term - nonlin - sextic_term - coupling - yukawa_term)
        res_Psi_list.append(res)

    # PDE LOSS
    pde_loss = ic_weight_boost * (torch.mean(res_Phi**2) + torch.mean(torch.stack(res_Psi_list)**2))

    # Krok 4.1: Używamy stałych globalnych
    # L2 normalization constraint
    norm_constraint = NORMALIZATION_WEIGHT * torch.mean((torch.mean(Psi_center**2, dim=0) - 1.0)**2)

    # Boundary condition penalty (Ψ(r_max)=0)
    #boundary_penalty = BOUNDARY_PENALTY_WEIGHT * torch.mean(Psi_center[:, -5:]**2)
    boundary_mask = (r > r_max - dr_cpu * 10).float()  # Ostatnie ~10 punktów siatki (dostosuj w razie potrzeby)
    boundary_weighted = boundary_mask.unsqueeze(1) * Psi_center**2  # [batch, octaves]
    boundary_penalty = BOUNDARY_PENALTY_WEIGHT * torch.mean(boundary_weighted)
    # FIX #9/Krok 2.1: Enhanced VEV reward with fade-in (STABILIZED)
    target_vev = torch.sqrt(torch.clamp(-mu2_val / lambda_H_val, min=VEV_MIN_CLAMP)) # min=1e-6
    mean_abs_phi = torch.mean(torch.abs(Phi_center)) + 1e-8 # unikaj zer

    fade_factor = 1.0 / (1.0 + torch.exp(torch.tensor(-0.1 * (current_epoch - 10), device=device)))

    vev_error_sq = (mean_abs_phi - target_vev)**2
    # Penalty becomes active after epoch 20
    vev_penalty = VEV_PENALTY_WEIGHT * vev_error_sq * fade_factor if current_epoch > 20 else 0.0

    # Early reward for approaching target (zabezpieczone przed zerem)
    gaussian_reward = torch.exp(-torch.abs(mean_abs_phi - target_vev) / torch.clamp(target_vev, min=0.1))
    vev_reward = -0.5 * (1.0 - gaussian_reward) * fade_factor

    # Physical penalties
    phys_penalty = PHYSICAL_PENALTY_WEIGHT * torch.mean(torch.relu(-Phi_center))
    phys_cap = target_vev + 2.0
    phys_penalty += 100.0 * torch.mean(torch.relu(Phi_center - phys_cap)**2)

    # Psi activity penalty
    psi_activity_penalty = 5.0 * torch.exp(-torch.mean(Psi_center**2) * 10)

    # COMBINED LOSS with all fixes
    total_loss = pde_loss + norm_constraint + boundary_penalty + vev_penalty + vev_reward + phys_penalty + psi_activity_penalty

    # Diagnostics
    if batch_idx % 50 == 0:
        tpu_print(f"\r[LOSS E{current_epoch}] PDE:{pde_loss.item():.2e} Norm:{norm_constraint:.2e} Bound:{boundary_penalty:.2e} VEV:{vev_penalty:.2e} |Φ|:{mean_abs_phi.item():.3f}", end='')

    return torch.nan_to_num(total_loss, nan=1e6, posinf=1e6, neginf=-1e6)

# ==============================================================================
# TRAINING FUNCTIONS (v40.4.1 - Fixed)
# ==============================================================================
def derivative_warmup(pinn, num_warm_steps=5):
    """Warmup XLA compilation graph"""
    if not IS_TPU:
        return

    tpu_print("🔄 Warming up derivatives...")
    mock_params = {
        'beta_max': 10.0, 'A_dip': 7.0, 'o_dip': 4.0, 'sigma_dip': 2.5,
        'A_dip2': 8.0, 'o_dip2': 9.0, 'sigma_dip2': 1.5,
        'A_k': 1.0, 'alpha_geo_k': 0.15, 'omega_k': 0.5, 'phi_k': np.pi/2,
        'g_Y_gen1': 1.0, 'g_Y_gen2': 4.0, 'g_Y_gen3': 12.0,
        'k_inv': 1.5, 'mu2': -15.0, 'm0': 1.0, 'g': 2.5, 'delta': 0.4,
        'lambda_Y_tau': 15.0, 'mass_scale_mu': 15.0, 'mass_scale_tau': 75.0,
        'k_mass': 1.8
    }

    optimizer_dummy = Adam(pinn.parameters(), lr=1e-6)
    batch_sizes = [128, 256, 512, 1024, 2048]

    for step, bs in enumerate(batch_sizes):
        r = torch.rand(bs, 1, device=device) * r_max
        t = torch.full((bs, 1), STATIONARY_TIME, device=device) # Krok 4.1
        th = torch.rand(bs, 1, device=device) * np.pi
        ph = torch.rand(bs, 1, device=device) * 2*np.pi

        loss = pinn_loss_v40(pinn, r, t, th, ph, mock_params,
                    current_epoch=0, batch_idx=0)

        optimizer_dummy.zero_grad()
        loss.backward()
        optimizer_dummy.step()

        if IS_TPU: xm.mark_step()
        tpu_print(f"  Warmup {step+1}/5: bs={bs}, loss={loss.item():.2e} ✓")


def create_dataloader(epoch):
    """Create optimized dataloader with correct batch sizing"""
    batch_size = get_dynamic_batch_size(epoch, 'TPU' if IS_TPU else 'CPU')

    # FIX #6: Create mesh that fits in batch size
    total_points = min(batch_size * 4, 500000)  # Cap to prevent OOM

    r_samples = torch.rand(total_points, 1) * r_max
    t_samples = torch.full((total_points, 1), STATIONARY_TIME)  # Krok 4.1: Stationary
    theta_samples = torch.rand(total_points, 1) * np.pi
    phi_samples = torch.rand(total_points, 1) * 2*np.pi

    mesh = torch.cat([r_samples, t_samples, theta_samples, phi_samples], dim=1).float()

    base_loader = DataLoader(TensorDataset(mesh), batch_size=batch_size, shuffle=True, drop_last=True)
    return PrefetchDataLoader(base_loader, device)

def pre_train_pinn(resume_from_path=None):
    """Pre-training with early stopping and checkpointing"""
    SUCCESS_THRESHOLD = 0.05
    STAGNATION_LIMIT = 300
    MAX_EPOCHS = 5000

    print("\n" + "="*60)
    print("STARTING PRE-TRAINING v40.4.1")
    print("="*60)

    mock_params = {
        'beta_max': 10.0, 'A_dip': 7.0, 'o_dip': 4.0, 'sigma_dip': 2.5,
        'A_dip2': 8.0, 'o_dip2': 9.0, 'sigma_dip2': 1.5,
        'A_k': 1.0, 'alpha_geo_k': 0.15, 'omega_k': 0.5, 'phi_k': np.pi/2,
        'g_Y_gen1': 1.0, 'g_Y_gen2': 4.0, 'g_Y_gen3': 12.0,
        'k_inv': 1.5, 'mu2': -15.0, 'm0': 1.0, 'g': 2.5, 'delta': 0.4,
        'lambda_Y_tau': 15.0, 'mass_scale_mu': 15.0, 'mass_scale_tau': 75.0,
        'k_mass': 1.8
    }

    pinn = SolitonPINN().to(device)
    optimizer = AdamW(pinn.parameters(), lr=1e-3, weight_decay=1e-5)
    scheduler = ReduceLROnPlateau(optimizer, mode='min', patience=10, factor=0.5, min_lr=1e-7)

    start_epoch, best_loss, stagnation = 0, float('inf'), 0

    # Load checkpoint if exists
    if resume_from_path and os.path.exists(resume_from_path):
        try:
            checkpoint = torch.load(resume_from_path, map_location=device)
            pinn.load_state_dict(checkpoint['model_state_dict'])
            optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
            start_epoch = checkpoint.get('epoch', 0) + 1
            best_loss = checkpoint.get('best_loss', float('inf'))
            print(f"🔄 Resumed from epoch {start_epoch}")
        except Exception as e:
            print(f"⚠️ Resume failed: {e}, starting fresh")

    # Initialize VEV
    vev_init = np.sqrt(max(-mock_params['mu2'] / lambda_H, 0.0))
    with torch.no_grad():
        pinn.out.bias[-1].fill_(vev_init * 0.3)  # Soft start

    derivative_warmup(pinn)

    for epoch in range(start_epoch, MAX_EPOCHS):
        dataloader = create_dataloader(epoch)
        epoch_loss, num_batches = 0.0, 0

        for batch_idx, (batch_rttp,) in enumerate(dataloader):
            r, t, th, ph = [b.to(device) for b in batch_rttp.split(1, dim=1)]

            loss = pinn_loss_v40(pinn, r, t, th, ph, mock_params,
                    current_epoch=epoch, ic_weight_boost=2.0, batch_idx=batch_idx)

            should_skip, loss_val = safe_loss_check_and_sync(loss, epoch, batch_idx)
            if should_skip:
                continue

            optimizer.zero_grad()
            loss.backward()
            # Krok 4.1: Używamy stałej
            torch.nn.utils.clip_grad_norm_(pinn.parameters(), GRAD_CLIP_PRETRAIN)

            if IS_TPU:
                xm.optimizer_step(optimizer)
            else:
                optimizer.step()

            epoch_loss += loss_val
            num_batches += 1

            if batch_idx % 100 == 0:
                print_progress_bar(epoch, batch_idx, len(dataloader), loss_val)

        avg_loss = epoch_loss / max(num_batches, 1)
        scheduler.step(avg_loss)

        print(f"\nEPOCH {epoch} | Loss: {avg_loss:.4e} | LR: {optimizer.param_groups[0]['lr']:.1e}")

        # Checkpointing
        if epoch % 5 == 0 or avg_loss < SUCCESS_THRESHOLD:
            checkpoint = {
                'epoch': epoch,
                'model_state_dict': {k: v.cpu() for k, v in pinn.state_dict().items()},
                'optimizer_state_dict': optimizer.state_dict(),
                'best_loss': best_loss,
                'current_loss': avg_loss
            }
            torch.save(checkpoint, 'pinn_latest.pth', _use_new_zipfile_serialization=False)

        # Early stopping
        if avg_loss < SUCCESS_THRESHOLD:
            print(f"\n✅ SUCCESS! Loss < {SUCCESS_THRESHOLD}")
            break

        if avg_loss < best_loss:
            best_loss, stagnation = avg_loss, 0
        else:
            stagnation += 1

        if stagnation >= STAGNATION_LIMIT:
            print(f"\n⚠️ STAGNATION after {stagnation} epochs")
            break

        # Krok 3.3: Memory Cleanup po każdej epoce
        del dataloader
        gc.collect()
        if IS_TPU:
            xm.mark_step()


    # Save final model
    final_path = 'pretrained_pinn.pth'
    torch.save({k: v.cpu() for k, v in pinn.state_dict().items()},
               final_path, _use_new_zipfile_serialization=False)
    print(f"✅ Pre-trained model saved: {final_path}")

    return pinn

# ==============================================================================
# FINE-TUNING SOLVER (v40.4.1 - Production Ready)
# ==============================================================================
def pinn_soliton_cached(params_v40):
    """
    Fast fine-tuning using pre-trained PINN as initial condition.
    Returns converged field profiles Ψ(r) and Φ(r).
    """
    # Check cache first
    Psi_cached, Phi_cached = load_finetune_result_v40(params_v40)
    if Psi_cached is not None:
        return Psi_cached, Phi_cached

    print(f"\n--- Starting Fine-Tuning (Trial params) ---")
    start_time = time.time()

    # Load pre-trained model
    pinn = SolitonPINN().to(device)
    try:
        pinn.load_state_dict(torch.load('pretrained_pinn.pth', map_location='cpu'))
        pinn = pinn.to(device)
        print("  ✅ Pre-trained model loaded")
    except Exception as e:
        print(f"  ⚠️ Load failed: {e}, using analytic init")
        mu2_val = params_v40['mu2']
        v_est = np.sqrt(max(-mu2_val / lambda_H, 0.0))
        save_finetune_result_v40(params_v40, None, None, float('inf'), 0)
        return np.full((num_octaves, Nr), 0.1), np.full(Nr, v_est)

    # Initialize VEV
    with torch.no_grad():
        real_vev = np.sqrt(max(-params_v40['mu2'] / lambda_H, 0.0))
        pinn.out.bias[-1].fill_(real_vev)
        pinn.out.weight[-1].data.fill_(0.0)
        print(f"  ✅ VEV initialized: {real_vev:.4f}")

    pinn.train()

    # Calculate mesh size and adjust batch
    r_vals = torch.tensor(r_cpu, dtype=torch.float32)
    t_val = torch.tensor([STATIONARY_TIME], dtype=torch.float32)  # Krok 4.1: Stationary
    theta_vals = torch.linspace(0, np.pi, Nr_theta_mesh, dtype=torch.float32)
    phi_vals = torch.linspace(0, 2*np.pi, Nr_phi_mesh, dtype=torch.float32)

    mesh = torch.cartesian_prod(r_vals, t_val, theta_vals, phi_vals)

    # Krok 3.2: TPU BATCH SIZE OPTIMIZATION
    raw_batch_size = min(FINETUNE_BATCH_SIZE, len(mesh))
    actual_batch_size = max(8, raw_batch_size // 8 * 8) # multiple of 8, min 8

    base_loader = DataLoader(TensorDataset(mesh), batch_size=actual_batch_size,
                           shuffle=True, drop_last=True)
    dataloader = PrefetchDataLoader(base_loader, device)

    optimizer = AdamW(pinn.parameters(), lr=5e-4, weight_decay=1e-5)
    scheduler = ReduceLROnPlateau(optimizer, mode='min', patience=2, factor=0.5, min_lr=1e-6)

    best_loss = float('inf')

    for epoch in range(12):  # Few epochs for fine-tuning
        epoch_losses = []

        for batch_rttp, in dataloader:
            r, t, th, ph = [b.to(device) for b in batch_rttp.split(1, dim=1)]

            loss = pinn_loss_v40(pinn, r, t, th, ph, params_v40,
                               current_epoch=epoch, ic_weight_boost=1.5)

            if not torch.isfinite(loss):
                continue

            optimizer.zero_grad()
            loss.backward()
            # Krok 4.1: Używamy stałej
            torch.nn.utils.clip_grad_norm_(pinn.parameters(), GRAD_CLIP_FINETUNE)

            if IS_TPU:
                xm.optimizer_step(optimizer)
            else:
                optimizer.step()

            epoch_losses.append(loss.item())

        avg_loss = np.mean(epoch_losses) if epoch_losses else 1e5
        best_loss = min(best_loss, avg_loss)
        scheduler.step(avg_loss)

        print(f"  [Fine-tune E{epoch+1}] Loss: {avg_loss:.4e}")

        # Krok 3.1: EARLY STOPPING
        if avg_loss < 0.1:
            print(f"  ✅ Early stopping at epoch {epoch+1}")
            break

        # Krok 3.3: Memory Cleanup po każdej epoce fine-tuningu
        del dataloader # Del dataloader only if it was redefined inside the loop
        gc.collect()
        if IS_TPU:
            xm.mark_step()

    # Final cleanup of resources from the loop
    if 'epoch_losses' in locals():
        del epoch_losses
    gc.collect()


    # Extract converged fields
    pinn.eval()
    with torch.no_grad():
        # FIX #6: Use correct batch size for inference
        out_list = []
        for i in range(0, len(mesh), INFERENCE_BATCH_SIZE):
            batch = mesh[i:i+INFERENCE_BATCH_SIZE].to(device)
            out_list.append(pinn(batch))

        out = torch.cat(out_list, dim=0).view(Nr, Nr_theta_mesh, Nr_phi_mesh, -1)
        out = out.cpu().numpy()

        Psi_out = np.mean(out[..., :-1], axis=(1,2)).T
        Phi_out = np.mean(out[..., -1], axis=(1,2))

    elapsed = time.time() - start_time
    save_finetune_result_v40(params_v40, Psi_out, Phi_out, best_loss, elapsed)

    return Psi_out, Phi_out

# ==============================================================================
# L-BFGS-B SOLVER (v40.4 - Stable)
# ==============================================================================
current_trial_params = {}
solver_iteration_data = {}
beta_low_octaves_stable = []

def vector_to_fields(vector, num_octaves_val, Nr_val):
    """Reshape flat vector to field arrays"""
    Psi_flat = vector[:-Nr_val]
    Phi_flat = vector[-Nr_val:]
    return Psi_flat.reshape((num_octaves_val, Nr_val)), Phi_flat

def fields_to_vector(Psi, Phi):
    """Flatten field arrays to vector"""
    return np.concatenate([Psi.flatten(), Phi.flatten()])

def energy_func_for_solver(flat_fields):
    """Energy functional for L-BFGS-B minimization"""
    global current_trial_params

    Psi, Phi_H = vector_to_fields(flat_fields, num_octaves, Nr)
    energy = total_energy_v40(Psi, Phi_H, current_trial_params, r_cpu, dr_cpu, np)
    return energy

def gradient_func_for_solver(flat_fields):
    """Gradient of energy functional"""
    global current_trial_params

    Psi, Phi_H = vector_to_fields(flat_fields, num_octaves, Nr)
    dE_Psi, dE_Phi = functional_derivative_v40(Psi, Phi_H, current_trial_params, r_cpu, dr_cpu, np)
    return fields_to_vector(dE_Psi, dE_Phi)

def solver_callback(xk):
    """Callback for monitoring L-BFGS-B iterations"""
    global solver_iteration_data, beta_low_octaves_stable

    solver_iteration_data['iter'] += 1
    if solver_iteration_data['iter'] % 50 == 0:
        print(f"  L-BFGS Iter: {solver_iteration_data['iter']} | Energy: {energy_func_for_solver(xk):.4e}")

def run_self_consistent_job(job_params, job_id):
    """
    Main solver entry point.
    Returns: (energy, error_message, masses)
    """
    global current_trial_params, solver_iteration_data
    current_trial_params = job_params

    print(f"\n{'='*60}")
    print(f"JOB {job_id}: Starting L-BFGS-B solver")
    print(f"{'='*60}")

    solver_iteration_data = {'iter': 0}
    beta_low_octaves_stable.clear()

    # Initialize fields from PINN
    try:
        Psi_init, Phi_init = pinn_soliton_cached(job_params)
        print(f"  Initial fields: max|Ψ|={np.max(Psi_init):.3f}, max|Φ|={np.max(Phi_init):.3f}")
    except Exception as e:
        print(f"  ⚠️ PINN init failed: {e}, using analytic")
        mu2_val = job_params['mu2']
        v_est = np.sqrt(max(-mu2_val / lambda_H, 0.0))
        Psi_init = np.full((num_octaves, Nr), 0.1)
        Phi_init = np.full(Nr, v_est)

    initial_vector = fields_to_vector(Psi_init, Phi_init)
    bounds = [(-10, 10)] * len(initial_vector)

    print("  Running L-BFGS-B optimization...")

    start_time = time.time()
    result = minimize(
        fun=energy_func_for_solver,
        x0=initial_vector,
        method='L-BFGS-B',
        jac=gradient_func_for_solver,
        bounds=bounds,
        callback=solver_callback,
        options={'maxiter': 100, 'gtol': 1e-3, 'ftol': 1e-6}
    )

    elapsed = time.time() - start_time
    print(f"  Solver completed in {elapsed:.1f}s, {result.nit} iterations")

    if not result.success:
        return None, f"Solver failed: {result.message}", None

    # Extract final fields and diagonalize
    Psi_final, Phi_final = vector_to_fields(result.x, num_octaves, Nr)
    masses = diagonalize_v40(Psi_final, Phi_final, job_params)

    if masses is None or len(masses) < 3:
        return result.fun, "Diagonalization failed", None

    print(f"  Extracted masses: {[f'{m:.3f}' for m in masses[:3]]}")

    return result.fun, "", masses

# ==============================================================================
# OPTUNA OBJECTIVE (v40.4 - Hardened)
# ==============================================================================
def objective_unified_v40(trial):
    """Multi-objective optimization for masses and forces"""
    print(f"\n--- OPTUNA Trial #{trial.number} ---")

    # Define parameter space (20 parameters)
    job_params = {
        'beta_max': trial.suggest_float('beta_max', 5.0, 15.0),
        'A_dip': trial.suggest_float('A_dip', 4.0, 10.0),
        'o_dip': trial.suggest_float('o_dip', 3.0, 5.0),
        'sigma_dip': trial.suggest_float('sigma_dip', 1.5, 4.0),
        'A_dip2': trial.suggest_float('A_dip2', 2.0, 6.0),
        'o_dip2': trial.suggest_float('o_dip2', 7.0, 10.0),
        'sigma_dip2': trial.suggest_float('sigma_dip2', 1.0, 3.0),
        'A_k': trial.suggest_float('A_k', 0.1, 2.0),
        'omega_k': trial.suggest_float('omega_k', 0.2, 1.0),
        'phi_k': trial.suggest_float('phi_k', 0.0, 2*np.pi),
        'alpha_geo_k': trial.suggest_float('alpha_geo_k', 0.05, 0.3),
        'g_Y_gen1': trial.suggest_float('g_Y_gen1', 0.1, 2.0),
        'g_Y_gen2': trial.suggest_float('g_Y_gen2', 1.0, 8.0),
        'g_Y_gen3': trial.suggest_float('g_Y_gen3', 10.0, 80.0),
        'k_inv': trial.suggest_float('k_inv', 0.5, 2.5),
        'mu2': trial.suggest_float('mu2', -25.0, -5.0),
        'm0': trial.suggest_float('m0', 0.1, 2.0),
        'g': trial.suggest_float('g', 0.5, 5.0),
        'k_mass': trial.suggest_float('k_mass', 0.5, 3.0),
        'delta': trial.suggest_float('delta', 0.1, 0.8),
        'lambda_Y_tau': trial.suggest_float('lambda_Y_tau', 5.0, 25.0),
        'mass_scale_mu': trial.suggest_float('mass_scale_mu', 10.0, 30.0),
        'mass_scale_tau': trial.suggest_float('mass_scale_tau', 50.0, 5000.0)
    }

    # Compute force hierarchy
    beta_profile = np.array([beta_topo_gaussian_dip(
        o, job_params['beta_max'], job_params['A_dip'], job_params['o_dip'], job_params['sigma_dip'],
        job_params['A_dip2'], job_params['o_dip2'], job_params['sigma_dip2']
    ) for o in range(num_octaves)])

    g1, g2, g3 = compute_force_hierarchy_v40(beta_profile, job_params['k_inv'])

    # HARD CONSTRAINT: Force hierarchy must be correct
    if not (g1 < g2 < g3):
        print(f"❌ FORCE HIERARCHY VIOLATED: g1={g1:.3f}, g2={g2:.3f}, g3={g3:.3f}")
        return 1e9, 1e9, 1e9

    # SOFT CONSTRAINT: Force ratios close to target
    error_force = 1000.0 * ((g2/g1) - 1.8)**2 + 1000.0 * ((g3/g2) - 1.89)**2

    if error_force > 500.0:
        raise optuna.exceptions.TrialPruned()

    # Run solver
    energy, error_msg, masses = run_self_consistent_job(job_params, f"trial_{trial.number}")

    # Compute mass error
    error_mass = 1e9
    if masses is not None and len(masses) >= 3:
        m1, m2, m3 = masses[:3]
        if m1 > 1e-12 and m2 > 1e-12 and m3 > 1e-12:
            log_ratio_mu = np.log(m2/m1)
            log_ratio_tau = np.log(m3/m1)
            error_mass = (log_ratio_mu - np.log(TARGET_MUON_GEV/TARGET_ELECTRON_GEV))**2 + \
                        (log_ratio_tau - np.log(TARGET_TAU_GEV/TARGET_ELECTRON_GEV))**2

    total_cost = error_mass + error_force

    print(f"### TRIAL {trial.number} RESULT ###")
    print(f"  Cost: {total_cost:.4f} (Mass err: {error_mass:.3f}, Force err: {error_force:.3f})")
    print(f"  Couplings: g1={g1:.4f}, g2={g2:.4f}, g3={g3:.4f}")
    if masses is not None:
        print(f"  Masses: {[f'{m:.3f}' for m in masses[:3]]}")

    # Cache parameters
    for k, v in job_params.items():
        trial.set_user_attr(k, v)

    return total_cost, error_mass, error_force

# ==============================================================================
# MAIN RUNNER (v40.4 - Orchestration)
# ==============================================================================
def main_runner(rank=0):
    """Main execution orchestration"""
    global device, IS_TPU, xm

    if IS_TPU:
        # Check if xm is available before use
        try:
            import torch_xla.core.xla_model as xm
            device = xm.xla_device()
        except ImportError:
            print("Warning: torch_xla not imported, running on CPU/GPU as fallback.")
            IS_TPU = False
            if torch.cuda.is_available():
                device = torch.device("cuda")
            else:
                device = torch.device("cpu")


    print("\n" + "="*80)
    print(" v40.4.1 MAIN EXECUTION ")
    print("="*80)

    # Phase 1: Pre-training
    print("\n" + "-"*40 + " PHASE 1: PRE-TRAINING " + "-"*40)
    should_pretrain = (EXECUTION_MODE == 'PRETRAIN_ONLY') or \
                     (EXECUTION_MODE == 'FULL_RUN' and not os.path.exists('pretrained_pinn.pth'))

    if should_pretrain:
        resume_from = 'pinn_latest.pth' if os.path.exists('pinn_latest.pth') else None
        pinn_model = pre_train_pinn(resume_from_path=resume_from)

        if EXECUTION_MODE == 'PRETRAIN_ONLY':
            print("✅ Pre-training completed, exiting")
            return
    else:
        print("✅ Skipping pre-training (model exists)")

    # Phase 2: Optuna Optimization
    print("\n" + "-"*40 + " PHASE 2: OPTUNA OPTIMIZATION " + "-"*40)
    list_cache_stats()

    if not OPTUNA_AVAILABLE:
        print("❌ Optuna not available, cannot continue")
        return

    # Initialize study
    study_name = 'supersoliton_v40_final'
    study_db = f'sqlite:////kaggle/working/{study_name}.db'

    try:
        study = optuna.load_study(study_name=study_name, storage=study_db,
                                sampler=NSGAIISampler(population_size=12, seed=42))
        n_completed = len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])
        print(f"✅ Resumed study: {n_completed} trials completed")
    except Exception as e:
        print(f"Failed to load study: {e}")
        study = optuna.create_study(
            study_name=study_name,
            directions=['minimize', 'minimize', 'minimize'],
            sampler=NSGAIISampler(population_size=12, seed=42),
            storage=study_db,
            load_if_exists=True
        )
        print("🆕 Created new study")

    # Run trials
    if n_initial_trials > 0:
        print(f"[OPTUNA] Starting {n_initial_trials} trials...")
        study.optimize(objective_unified_v40, n_trials=n_initial_trials, n_jobs=1)

    # Phase 3: Results Analysis
    print("\n" + "-"*40 + " PHASE 3: RESULTS ANALYSIS " + "-"*40)

    if len(study.trials) > 0:
        # Plot Pareto front
        fig, ax = plt.subplots(figsize=(12, 8))

        trials_df = study.trials_dataframe()
        if not trials_df.empty:
            complete_trials = trials_df[trials_df['state'] == 'COMPLETE']
            if not complete_trials.empty:
                ax.scatter(complete_trials['values_1'], complete_trials['values_2'],
                          c=complete_trials['values_0'], cmap='viridis', alpha=0.7)
                ax.set_xlabel("Mass Error")
                ax.set_ylabel("Force Error")
                ax.set_title("Pareto Front (Mass vs Force)")
                plt.colorbar(ax.collections[0], label='Total Cost')
                plt.savefig("pareto_front_v40_4.png", dpi=150)
                print("✅ Pareto front saved")

        # Best trial analysis
        try:
            # Note: Optuna NSGAII stores best_trials as a list of non-dominated solutions
            best_trial = study.best_trials[0]
            print(f"\n🏆 BEST TRIAL (Non-dominated set leader): #{best_trial.number}")
            print(f"  Total Cost: {best_trial.values[0]:.4f}")

            # Extract parameters
            bp = {k: best_trial.user_attrs[k] for k in best_trial.user_attrs if k.startswith('params_')}

            # Recompute force hierarchy
            beta_prof = np.array([beta_topo_gaussian_dip(
                o, bp['beta_max'], bp['A_dip'], bp['o_dip'], bp['sigma_dip'],
                bp['A_dip2'], bp['o_dip2'], bp['sigma_dip2']
            ) for o in range(num_octaves)])

            g1, g2, g3 = compute_force_hierarchy_v40(beta_prof, bp['k_inv'])
            print(f"  Force couplings: g1={g1:.4f}, g2={g2:.4f}, g3={g3:.4f}")
            print(f"  Ratios: g2/g1={g2/g1:.3f} (target 1.80), g3/g2={g3/g2:.3f} (target 1.89)")

            # Predict mass hierarchy
            scale_boost = bp['mass_scale_tau'] / bp['mass_scale_mu']
            yukawa_boost = bp['lambda_Y_tau'] / bp['g_Y_gen3']
            predicted_tau_e = 206.77 * scale_boost * (1 + yukawa_boost)
            print(f"  Predicted τ/e: {predicted_tau_e:.0f} (target 3477)")

            # Save detailed results
            results = {
                'trial_number': best_trial.number,
                'parameters': bp,
                'force_couplings': {'g1': g1, 'g2': g2, 'g3': g3},
                'predicted_ratios': {'tau_over_e': predicted_tau_e},
                'beta_profile': beta_prof.tolist()
            }

            with open('best_trial_results.json', 'w') as f:
                json.dump(results, f, indent=2)
            print("✅ Best trial results saved to best_trial_results.json")

        except Exception as e:
            print(f"⚠️ Best trial analysis failed: {e}")

    # Final summary
    print("\n" + "="*80)
    print(" v40.4.1 EXECUTION COMPLETED ")
    print("="*80)
    print(f"Trials completed: {len(study.trials)}")
    print(f"Cache entries: {len(glob.glob(os.path.join(CACHE_DIR, 'ft_*.pkl')))}")
    print("="*80)

# ==============================================================================
# ENTRY POINT
# ==============================================================================
if __name__ == '__main__':
    try:
        main_runner(0)
    except KeyboardInterrupt:
        print("\n\n🛑 Execution interrupted by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n\n❌ FATAL ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
