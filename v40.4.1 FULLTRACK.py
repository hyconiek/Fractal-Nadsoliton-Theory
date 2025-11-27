#!/usr/bin/env python3
# Author: Krzysztof Żuchowski
# Version: v40.4.1-FAST (Kaggle Optimized)

"""
v40.4.1 FAST TRACK (Faza X.4: Rapid Verification)

CRITICAL ARCHITECTURAL FIXES (Unified & Optimized):
═══════════════════════════════════════════════════════════════════════════════
[FAST TRACK SETTINGS]
1. Grid Size: Reduced to 64 (spatial) x 8 (angular) for speed.
2. Epochs: Pre-train capped at 300, Fine-tune at 5.
3. Trials: Optuna set to 50 trials for rapid landscape scanning.

[NUMERICAL STABILITY]
1. VEV Stabilization: Improved clamping and epsilon addition.
2. Beta Epsilon Increase: 1e-6 for mass enhancement stability.
═══════════════════════════════════════════════════════════════════════════════
"""

print("="*80)
print(" INITIALIZING v40.4.1 (FAST TRACK EDITION) ")
print("="*80)

EXECUTION_MODE = 'FULL_RUN'
DEVICE_MODE = 'AUTO'  # 'AUTO', 'TPU', 'GPU', 'CPU'

# ==============================================================================
# IMPORTS AND ENVIRONMENT
# ==============================================================================
import os, sys, time, warnings, subprocess, gc, json, hashlib, pickle
import numpy as np
import pandas as pd
import scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import glob

warnings.filterwarnings("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
os.environ['XLA_USE_BF16'] = '1'

# Install Optuna if missing
try:
    import optuna
    from optuna.samplers import NSGAIISampler
    OPTUNA_AVAILABLE = True
except ImportError:
    print("⚠️ Installing Optuna...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "optuna"])
    import optuna
    from optuna.samplers import NSGAIISampler
    OPTUNA_AVAILABLE = True

import torch
import torch.nn as nn
from torch.optim import AdamW
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data import TensorDataset, DataLoader

# ==============================================================================
# DEVICE SETUP (Optimized for Kaggle)
# ==============================================================================
print("\n[INFO] Device detection...")
IS_TPU = False
try:
    import torch_xla
    import torch_xla.core.xla_model as xm
    device = xm.xla_device()
    IS_TPU = True
    print(f"✅ TPU DETECTED: {device}")
except:
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print(f"✅ GPU DETECTED: {device}")
    else:
        device = torch.device("cpu")
        print(f"⚠️ CPU FALLBACK (Slow)")

# ==============================================================================
# GLOBAL PARAMETERS (FAST TRACK)
# ==============================================================================
# Grid resolution (Optimized for Speed)
Nr = 256          # Was 128
Nr_theta = 8     # Was 16
Nr_phi = 8       # Was 16

# Physics constants
r_max = 25.0
num_octaves = 12
lambda_H = 0.5
delta = 0.2
beta_hierarchy = 0.15

# Target masses (GeV) - For Reference/Loss
TARGET_ELECTRON_GEV = 0.000511
TARGET_MUON_GEV = 0.10566
TARGET_TAU_GEV = 1.77686

# Constants
STATIONARY_TIME = 0.5
GRAD_CLIP = 5.0
EPSILON_BETA = 1e-6
VEV_MIN_CLAMP = 0.1

# Optimization Config
PRETRAIN_EPOCHS = 5000    # Was 5000
FINETUNE_EPOCHS = 5      # Was 12
OPTUNA_TRIALS = 50       # Was 3000

# Batch sizes
BATCH_SIZE = 4096 if not IS_TPU else 16384

# Grid arrays
r_cpu = np.linspace(1e-6, r_max, Nr, dtype=np.float64)
dr_cpu = r_cpu[1] - r_cpu[0]

# Paths
CACHE_DIR = 'finetune_cache_fast'
os.makedirs(CACHE_DIR, exist_ok=True)

print(f"[CONFIG] Grid: {Nr}x{Nr_theta}x{Nr_phi}")
print(f"[CONFIG] Pretrain Epochs: {PRETRAIN_EPOCHS} | Trials: {OPTUNA_TRIALS}")

# ==============================================================================
# PHYSICS FUNCTIONS
# ==============================================================================
def get_gen_idx_and_scale(octave_num, mass_scale_mu=15.0, mass_scale_tau=75.0):
    if octave_num < 4: return 0, 1.0
    elif octave_num < 7: return 1, mass_scale_mu
    else: return 2, mass_scale_tau

def beta_topo_gaussian_dip(o, beta_max, A_dip, o_dip, sigma_dip, A_dip2, o_dip2, sigma_dip2):
    dip1 = A_dip * np.exp(-(o - o_dip)**2 / (2 * sigma_dip**2))
    dip2 = A_dip2 * np.exp(-(o - o_dip2)**2 / (2 * sigma_dip2**2))
    return beta_max - dip1 - dip2

def beta_topo_gaussian_dip_torch(o, beta_max, A_dip, o_dip, sigma_dip, A_dip2, o_dip2, sigma_dip2, dev):
    o_tensor = torch.tensor(o, dtype=torch.float32, device=dev)
    dip1 = A_dip * torch.exp(-(o_tensor - o_dip)**2 / (2 * sigma_dip**2))
    dip2 = A_dip2 * torch.exp(-(o_tensor - o_dip2)**2 / (2 * sigma_dip2**2))
    return beta_max - dip1 - dip2

def compute_interaction_kernel(i, j, params, use_fast_kernel=True):
    d = abs(i - j)
    if d not in [1, 2]: return 0.0

    beta_topo = beta_topo_gaussian_dip(
        min(i,j), params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
        params['A_dip2'], params['o_dip2'], params['sigma_dip2']
    )

    A, alpha_geo = params['A_k'], params['alpha_geo_k']
    mod_diff = abs((i%3)-(j%3))

    K_ij_raw = A * (2**(-alpha_geo * d)) * np.exp(-beta_topo * mod_diff)
    return K_ij_raw if d == 1 else K_ij_raw / 2.0

def compute_interaction_kernel_torch(i, j, params, dev):
    d = abs(i - j)
    if d not in [1, 2]: return torch.tensor(0.0, device=dev)

    beta_topo = beta_topo_gaussian_dip_torch(
        min(i,j), params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
        params['A_dip2'], params['o_dip2'], params['sigma_dip2'], dev
    )

    A = torch.tensor(params['A_k'], device=dev)
    alpha_geo = torch.tensor(params['alpha_geo_k'], device=dev)

    K_ij_raw = A * (2**(-alpha_geo * d)) * torch.exp(-beta_topo * abs((i%3)-(j%3)))
    return K_ij_raw if d == 1 else K_ij_raw / 2.0

def compute_force_hierarchy(beta_profile, k_inv):
    # Force hierarchy logic from theory
    g1 = 1.0 / (beta_profile[0]**k_inv + 1e-9)
    g2 = 1.0 / (beta_profile[2]**k_inv + 1e-9)
    g3 = 1.0 / (beta_profile[4]**k_inv + 1e-9)
    return g1, g2, g3

def diagonalize_v40(Psi_loc, Phi_loc, params):
    """Diagonalization with mass enhancement logic"""
    Nfull = num_octaves * Nr
    neigs_to_calc = 10 # Only need lowest for mass ratios
    k_mass = params.get('k_mass', 1.0)

    H = np.zeros((Nfull, Nfull), dtype=float)

    for o in range(num_octaves):
        idx0 = o * Nr
        gen_idx, mass_scale = get_gen_idx_and_scale(o, params['mass_scale_mu'], params['mass_scale_tau'])

        beta_val = beta_topo_gaussian_dip(
            o, params['beta_max'], params['A_dip'], params['o_dip'], params['sigma_dip'],
            params['A_dip2'], params['o_dip2'], params['sigma_dip2']
        )
        mass_enhancement = (1.0 / (abs(beta_val) + EPSILON_BETA)) ** k_mass

        # Potential terms
        yukawa = params[f'g_Y_gen{gen_idx+1}'] * (Phi_loc**2) * mass_scale * mass_enhancement
        if gen_idx == 2:
            yukawa += 1.5 * params['lambda_Y_tau'] * (Phi_loc**2) * (Psi_loc[o]**2) * mass_enhancement

        V_pp = params['m0']**2 * mass_enhancement + 3*params['g']*Psi_loc[o]**2 + yukawa

        # Laplacian discretization
        for i in range(Nr):
            H[idx0+i, idx0+i] = (2/dr_cpu**2) + V_pp[i]
        for i in range(Nr - 1):
            off = -1/dr_cpu**2
            H[idx0+i, idx0+i+1] = off
            H[idx0+i+1, idx0+i] = off

    # Off-diagonal couplings
    for i in range(num_octaves):
        for j in range(i + 1, num_octaves):
            K_ij = compute_interaction_kernel(i, j, params)
            if K_ij != 0:
                idx_i, idx_j = i * Nr, j * Nr
                for k in range(Nr):
                    H[idx_i + k, idx_j + k] = K_ij
                    H[idx_j + k, idx_i + k] = K_ij

    # Solve
    try:
        vals, _ = spl.eigsh(sp.csr_matrix(H), k=neigs_to_calc, which='SA', tol=1e-4)
        return np.sqrt(np.sort(vals[vals > 1e-12]))
    except:
        return np.array([])

# ==============================================================================
# NEURAL NETWORK (PINN)
# ==============================================================================
class SolitonPINN(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(4, 128), nn.LayerNorm(128), nn.GELU(),
            nn.Linear(128, 128), nn.GELU(),
            nn.Linear(128, 128), nn.GELU(),
            nn.Linear(128, num_octaves + 1)
        )
        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)

    def forward(self, x):
        return self.net(x)

def pinn_loss(model, coords, params):
    # Enable gradient for inputs
    coords.requires_grad_(True)
    output = model(coords)

    Psi = output[:, :num_octaves]
    Phi = output[:, -1:]

    # Calculate derivatives (simplified for speed)
    # Laplacian via finite difference approximation on random batch is hard
    # Using Autograd for time/space derivatives

    grads = torch.autograd.grad(output, coords, torch.ones_like(output), create_graph=True)[0]
    # coords: [r, t, theta, phi]
    # grads: [batch, 4] -> THIS IS WRONG for vector output.
    # Need per-channel derivative. Optimized approach:

    # Approximation: Minimize Energy Functional directly instead of PDE residual for speed
    # E = Kinetic + Potential

    # Only radial gradient for now (dominant)
    dr_Psi = torch.autograd.grad(Psi.sum(), coords, create_graph=True)[0][:, 0:1]
    dr_Phi = torch.autograd.grad(Phi.sum(), coords, create_graph=True)[0][:, 0:1]

    # Potential Energy construction (Batched)
    # This approximates the PDE solution by minimizing energy

    loss_energy = torch.mean(0.5 * dr_Psi**2) + torch.mean(0.5 * dr_Phi**2)

    # Self-interaction & Yukawa (Simplified for PINN pre-training)
    rho = Psi**2
    loss_pot = torch.mean(0.25 * params['g'] * rho**2 + 0.5 * params['mu2'] * Phi**2 + 0.25 * lambda_H * Phi**4)

    # Coupling
    loss_coup = 0
    # Simple approximation for coupling

    # Boundary conditions (r->max => Psi->0)
    r = coords[:, 0]
    mask_bound = (r > r_max * 0.9).float()
    loss_bc = torch.mean(mask_bound * (Psi**2 + (Phi - 246.0)**2)) # VEV anchor

    return loss_energy + loss_pot + 10.0 * loss_bc

# ==============================================================================
# TRAINING & OPTIMIZATION
# ==============================================================================
def pre_train():
    print("\n[PRE-TRAIN] Starting PINN warmup...")
    model = SolitonPINN().to(device)
    optim = AdamW(model.parameters(), lr=1e-3)

    # Mock params for pretraining
    mock_p = {
        'm0': 1.0, 'g': 1.0, 'mu2': -100.0, 'delta': 0.1,
        'g_Y_gen1': 0.1, 'g_Y_gen2': 1.0, 'g_Y_gen3': 10.0,
        'lambda_Y_tau': 10.0, 'k_mass': 1.0,
        'beta_max': 10.0, 'A_dip': 5.0, 'o_dip': 4.0, 'sigma_dip': 2.0,
        'A_dip2': 5.0, 'o_dip2': 8.0, 'sigma_dip2': 2.0,
        'A_k': 1.0, 'alpha_geo_k': 0.2, 'omega_k': 0.5, 'phi_k': 0.0,
        'mass_scale_mu': 15.0, 'mass_scale_tau': 75.0
    }

    for epoch in range(PRETRAIN_EPOCHS):
        # Generate random batch
        r = torch.rand(BATCH_SIZE, 1, device=device) * r_max
        t = torch.zeros(BATCH_SIZE, 1, device=device)
        th = torch.rand(BATCH_SIZE, 1, device=device) * np.pi
        ph = torch.rand(BATCH_SIZE, 1, device=device) * 2*np.pi
        coords = torch.cat([r, t, th, ph], dim=1)

        loss = pinn_loss(model, coords, mock_p)

        optim.zero_grad()
        loss.backward()
        if IS_TPU: xm.optimizer_step(optim)
        else: optim.step()

        if epoch % 50 == 0:
            l_val = loss.item()
            print(f"  Epoch {epoch}: Loss = {l_val:.4e}")
            if l_val < 0.1: break

    torch.save(model.state_dict(), 'pinn_base.pth')
    print("[PRE-TRAIN] Model saved.")
    return model

def objective(trial):
    # Suggest params
    p = {
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

    # 1. Check Force Hierarchy (Fast check)
    beta_prof = np.array([beta_topo_gaussian_dip(o, **{k:v for k,v in p.items() if k in ['beta_max','A_dip','o_dip','sigma_dip','A_dip2','o_dip2','sigma_dip2']}) for o in range(num_octaves)])
    g1, g2, g3 = compute_force_hierarchy(beta_prof, p['k_inv'])

    if not (g1 < g2 < g3):
        raise optuna.exceptions.TrialPruned("Force hierarchy violated")

    # 2. Solve System (Simplified: Start from flat, optimize energy)
    # We skip full PINN fine-tuning in FAST mode and use direct minimization if possible,
    # but here we'll do a very short PINN run or direct L-BFGS on CPU

    # Initialize fields analytically
    v_est = np.sqrt(max(-p['mu2']/lambda_H, 0.01))
    Psi_init = np.full((num_octaves, Nr), 0.1)
    Phi_init = np.full(Nr, v_est)

    # Quick Diagonalization check
    masses = diagonalize_v40(Psi_init, Phi_init, p)

    if len(masses) < 3:
        return 1e9, 1e9

    m1, m2, m3 = masses[:3]

    # Loss calculation
    err_mass = (np.log(m2/m1) - np.log(TARGET_MUON_GEV/TARGET_ELECTRON_GEV))**2 + \
               (np.log(m3/m1) - np.log(TARGET_TAU_GEV/TARGET_ELECTRON_GEV))**2

    err_force = ((g2/g1) - 1.8)**2 + ((g3/g2) - 1.89)**2

    total_cost = err_mass + err_force * 10.0

    return total_cost, err_mass

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================
if __name__ == "__main__":
    print("\n[EXEC] Starting Fast Track Optimization...")

    # 1. Pretrain base model
    # p_model = pre_train() # Skip in fast check, rely on analytic init

    # 2. Run Optimization
    study = optuna.create_study(
        directions=['minimize', 'minimize'],
        sampler=NSGAIISampler(population_size=10)
    )

    print(f"[EXEC] Running {OPTUNA_TRIALS} trials...")
    try:
        study.optimize(objective, n_trials=OPTUNA_TRIALS)
    except KeyboardInterrupt:
        print("Stopped by user.")

    # 3. Report
    print("\n" + "="*60)
    print(" OPTIMIZATION RESULTS ")
    print("="*60)

    if len(study.best_trials) > 0:
        best = study.best_trials[0]
        print(f"🏆 Best Trial #{best.number}")
        print(f"   Cost: {best.values[0]:.4f}")
        print(f"   Mass Error: {best.values[1]:.4f}")
        print("\nParameters:")
        print(json.dumps(best.params, indent=2))

        # Save
        with open('best_params_fast.json', 'w') as f:
            json.dump(best.params, f)

    print("="*60)
    print("DONE.")
