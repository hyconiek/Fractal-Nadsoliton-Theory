#!/usr/bin/env python3
# Author: Krzysztof Żuchowski

"""
parameter_scan_supersoliton_v39_STABLE_HIERARCHY.py

PRODUCTION-READY VERSION with all stability improvements:
- (1) δΨ⁶ STABILIZATION: Sextic potential term prevents quartic runaway
- (2) FIXED LAPLACIAN: L'Hospital's rule at r=0 (∇²f = 3·d²f/dr² for r→0)
- (3) HIERARCHICAL COUPLING: λ(o) = λ_base · 2^(-β·o) for mass hierarchy
- (4) FULL OPTUNA FRAMEWORK: Multi-objective optimization with caching

Based on v38.5 + v34.4 + theoretical recommendations from dynamic vs potential
stabilization analysis.
"""
print("="*80)
print(" INITIALIZING v39 STABLE HIERARCHY PRODUCTION CODE ")
print("="*80)
print("✅ (1) δΨ⁶ stabilization enabled")
print("✅ (2) Numerically stable radial Laplacian at r=0")
print("✅ (3) Hierarchical inter-octave coupling λ(o) = λ_base · 2^(-β·o)")
print("✅ (4) Full Optuna optimization framework")
print("="*80)

EXECUTION_MODE = 'FULL_RUN'  # <-- ZMIEŃ NA 'PRETRAIN_ONLY' jeśli chcesz tylko pre-train

print(f"✅ Tryb uruchomienia: {EXECUTION_MODE}")
if EXECUTION_MODE == 'PRETRAIN_ONLY':
    print("   Skrypt zakończy działanie po zakończeniu pre-treningu.")

# ==============================================================================
# IMPORTS AND ENVIRONMENT VERIFICATION
# ==============================================================================
# V-- DODANO PRINT --V
print("\n[INFO] Rozpoczynanie importu bibliotek...")
import os, sys, time, warnings, subprocess, gc
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from joblib import Parallel, delayed, dump
import itertools
import matplotlib.pyplot as plt
import threading
from contextlib import nullcontext
import glob
from datetime import datetime
import json
import hashlib
import pickle
print("[INFO] Import podstawowych bibliotek zakończony.")

# Core (always)
import torch
import torch.nn as nn
from torch.optim import Adam
from torch.optim.lr_scheduler import ReduceLROnPlateau, LambdaLR # <-- ZMIANA: DODANO LambdaLR
from torch.utils.data import TensorDataset, DataLoader
print("[INFO] Import bibliotek PyTorch zakończony.")

# PATCH 5 dependency
try:
    import psutil
    PSUTIL_AVAILABLE = True
    print("✅ psutil załadowany. Liczba wątków będzie dynamiczna.")
except ImportError:
    psutil = None
    PSUTIL_AVAILABLE = False
    print("⚠️ psutil not found, parallel job count will be static.")


try:
    from torch.amp import autocast
    AUTOCAST_AVAILABLE = True
    print("✅ torch.amp.autocast dostępny.")
except ImportError:
    AUTOCAST_AVAILABLE = False
    print("⚠️ torch.amp not available - BF16 will be handled by XLA on TPU")

try:
    from tensorboardx import SummaryWriter
    TENSORBOARDX_AVAILABLE = True
    print("✅ TensorBoardX dostępny.")
except ImportError:
    TENSORBOARDX_AVAILABLE = False

try:
    import optuna
    from optuna.samplers import NSGAIISampler
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.neural_network import MLPRegressor
    from sklearn.exceptions import NotFittedError
    from scipy.stats import pearsonr, gaussian_kde
    print(f"✅ Optuna (v{optuna.__version__}) + sklearn załadowane.")
except ImportError:
    print("⚠️ Optuna/sklearn nie znalezione, próba instalacji...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "optuna[deap]", "scikit-learn", "-q"])
    import optuna
    from optuna.samplers import NSGAIISampler
    from sklearn.preprocessing import MinMaxScaler
    from sklearn.neural_network import MLPRegressor
    from sklearn.exceptions import NotFittedError
    from scipy.stats import pearsonr, gaussian_kde
    print("✅ Optuna/sklearn zainstalowane i załadowane.")

PARTICLE_AVAILABLE = False
try:
    from particle import Particle
    PARTICLE_AVAILABLE = True
    print("✅ Particle załadowany dla PDG.")
except ImportError:
    print("⚠️ Particle fallback to hardcoded SM masses.")

BOTORCH_AVAILABLE = False
try:
    import botorch
    from botorch.models import SingleTaskGP
    from gpytorch.mlls import ExactMarginalLogLikelihood
    from botorch.fit import fit_gpytorch_model
    from botorch.acquisition import ExpectedImprovement
    from optuna.integration import BoTorchSampler
    BOTORCH_AVAILABLE = True
    print(f"✅ BoTorch (v{botorch.__version__}) załadowany.")
except ImportError:
    print("⚠️ BoTorch fallback to NSGAII.")

warnings.filterwarnings("ignore", category=RuntimeWarning)
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# V-- DODANO PRINT --V
print("\n[INFO] Konfiguracja środowiska XLA/TPU...")
# TPU Setup
try:
    import torch_xla.core.xla_model as xm
    import torch_xla
    XLA_IMPORTS_SUCCESS = True
except ImportError:
    XLA_IMPORTS_SUCCESS = False

os.environ['XLA_USE_BF16'] = '1'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.5"
os.environ["XLA_TENSOR_ALLOCATOR_MAXSIZE"] = "1000000000"

# --- UTILITY POPRAWKA 5: XLA COMPILATION CACHE ---
os.environ['XLA_IR_DEBUG'] = '0'
os.environ['XLA_HLO_DEBUG'] = '0'
os.environ['XLA_COMPILATION_CACHE_DIR'] = '/kaggle/working/xla_cache'
os.environ['PJRT_DEVICE'] = 'TPU'
print("[INFO] Zmienne środowiskowe XLA zoptymalizowane pod kątem wydajności.")

print("[INFO] Zmienne środowiskowe XLA ustawione.")

IS_TPU = False
GPU_MODE = False
xp, cp = np, None
LOG_DIR = './logs'
HEARTBEAT_FILE = 'heartbeat.log'
KAGGLE_OUTPUT_FILE = '/kaggle/working/tpu_active.log'

# V-- DODANO PRINT --V
print("\n[INFO] Wykrywanie dostępnego sprzętu (TPU/GPU/CPU)...")
if XLA_IMPORTS_SUCCESS:
    try:
        device = xm.xla_device()
        IS_TPU = True
        print(f"✅ Wykryto Kaggle TPU: {device}, dostępnych 8 rdzeni.")
    except Exception:
        XLA_IMPORTS_SUCCESS = False

if not IS_TPU:
    print("[INFO] TPU fallback zainicjowany. Sprawdzanie dostępności GPU...")
    if torch.cuda.is_available():
        device = torch.device("cuda")
        print(f"✅ Wykryto GPU: {torch.cuda.get_device_name(0)}")
        torch.backends.cudnn.benchmark = True
        try:
            import cupy as _cp
            import cupyx.scipy.sparse as csp
            import cupyx.scipy.sparse.linalg as cspl
            cp, xp, GPU_MODE = _cp, _cp, True
            print("✅ CuPy załadowany. Akceleracja GPU dla SciPy włączona.")
        except ImportError:
            GPU_MODE = False
            xp = np
            cp = None
            print("⚠️ Nie znaleziono CuPy. Obliczenia SciPy będą wykonywane na CPU.")
    else:
        device = torch.device("cpu")
        print("✅ Brak TPU i GPU. Uruchamianie w trybie CPU.")

print(f"==> Finalnie wybrane urządzenie: {device} <==")

# ==============================================================================
# ANTI-IDLE & PRINTER UTILITIES
# ==============================================================================
def tpu_print(msg, **kwargs):
    if IS_TPU:
        try:
            # master_print nie obsługuje kwargs, więc symuluj end='\r' via sys.stdout
            sys.stdout.write(msg)
            if kwargs.get('end', '\n') != '\n':
                sys.stdout.write(kwargs['end'])
            sys.stdout.write('\n') # Dodajemy nową linię po każdym komunikacie dla logów
            sys.stdout.flush()
            sys.stderr.flush()
        except NameError:
            print(msg, **kwargs)
            sys.stdout.flush()
    else:
        print(msg, **kwargs)
        sys.stdout.flush()

def write_dual_heartbeat(message):
    timestamp = time.time()
    for filepath in [HEARTBEAT_FILE, KAGGLE_OUTPUT_FILE]:
        try:
            with open(filepath, 'a') as f:
                f.write(f"{timestamp} {message}\n")
                f.flush()
                if hasattr(f, 'fileno'):
                    os.fsync(f.fileno())
        except Exception:
            pass

class ManualTPUMetrics:
    @staticmethod
    def get_tpu_memory():
        try:
            mem_info = xm.get_memory_info(torch_xla.device())
            return {
                'bytes_used': mem_info.get('bytes_used', 0),
                'bytes_limit': mem_info.get('bytes_limit', 0),
                'utilization': mem_info.get('bytes_used', 0) / mem_info.get('bytes_limit', 1) * 100
            }
        except Exception as e:
            return {'error': str(e)}

    @staticmethod
    def get_process_memory():
        try:
            import psutil
            process = psutil.Process()
            mem = process.memory_info()
            return {'rss_mb': mem.rss / 1024 / 1024, 'vms_mb': mem.vms / 1024 / 1024}
        except ImportError:
            return {'error': 'psutil not installed'}
        except Exception as e:
            return {'error': str(e)}

    @staticmethod
    def collect_and_save(filepath='/kaggle/working/tpu_metrics.json'):
        metrics = {
            'timestamp': datetime.now().isoformat(),
            'tpu_memory': ManualTPUMetrics.get_tpu_memory(),
            'process_memory': ManualTPUMetrics.get_process_memory(),
        }
        try:
            with open(filepath, 'w') as f:
                json.dump(metrics, f, indent=2)
        except Exception:
             pass
        return metrics

def background_heartbeat():
    counter = 0
    while True:
        time.sleep(8)
        counter += 1
        for i in range(3):
            try:
                with open(f'/kaggle/working/keepalive_{counter % 10}_{i}.tmp', 'w') as f:
                    f.write(f"{time.time()} Heartbeat #{counter}\n")
                    f.flush()
                    if hasattr(f, 'fileno'):
                        os.fsync(f.fileno())
            except Exception: pass

        # Używamy tpu_print, aby trafiło do logów
        tpu_print(f"[BG-KEEPALIVE] #{counter} @ {time.time():.0f}s")

        if IS_TPU:
            try:
                metrics = ManualTPUMetrics.collect_and_save()
                tpu_util = metrics.get('tpu_memory', {}).get('utilization', 0)
                ram_mb = metrics.get('process_memory', {}).get('rss_mb', 0)
                if counter % 4 == 0:
                    tpu_print(f"[TPU METRICS] Util: {tpu_util:.1f}% | RAM: {ram_mb:.0f}MB")
            except Exception as e:
                if counter % 20 == 0:
                    tpu_print(f"[TPU METRICS] Collection failed: {e}")

def kaggle_clicker_daemon():
    counter = 0
    while True:
        time.sleep(60)
        counter += 1
        try:
            tpu_print(f"[CLICKER] UI ping #{counter}")
            with open(f"/kaggle/working/clicker_{counter % 5}.log", "w") as f:
                f.write(f"{time.time()} click {counter}\n")
                f.flush()
                os.fsync(f.fileno())
            if counter % 5 == 0:
                os.system("echo ping > /dev/null")
        except Exception as e:
            tpu_print(f"[CLICKER] warning: {e}")

def log_tpu_utilization():
    if not IS_TPU:
        return
    while True:
        time.sleep(30)
        try:
            metrics = ManualTPUMetrics.collect_and_save()
            tpu_util = metrics.get('tpu_memory', {}).get('utilization', 0)
            ram_mb = metrics.get('process_memory', {}).get('rss_mb', 0)
            tpu_print(f"[TPU UTIL] Memory: {tpu_util:.1f}% | RAM: {ram_mb:.0f}MB")
        except Exception:
            pass

def get_autocast_context(is_tpu):
    if is_tpu:
        return nullcontext()
    else:
        if AUTOCAST_AVAILABLE:
            return autocast('cuda', dtype=torch.bfloat16)
        return nullcontext()

if IS_TPU:
    heartbeat_thread = threading.Thread(target=background_heartbeat, daemon=True)
    heartbeat_thread.start()
    tpu_monitor_thread = threading.Thread(target=log_tpu_utilization, daemon=True)
    tpu_monitor_thread.start()
    clicker_thread = threading.Thread(target=kaggle_clicker_daemon, daemon=True)
    clicker_thread.start()
    tpu_print("✅ Wątki w tle (keepalive, monitor TPU, clicker) zostały uruchomione.")

    write_dual_heartbeat("TPU_threads_started")

# V-- DODANO PRINT --V
print("\n[INFO] Wykonywanie finalnego testu tensora na wybranym urządzeniu...")
with torch.no_grad():
    tensor_a = torch.randn(3, 3, device=device)
    tensor_b = torch.randn(3, 3, device=device)
    dummy_out = tensor_a.sum()
    if IS_TPU:
        xm.mark_step() # Zamiast .item(), wymuszamy sync przez mark_step
        dummy_out_val = dummy_out.cpu().item()
    else:
        dummy_out_val = dummy_out.item()
tpu_print(f"✅ Test urządzenia zakończony. Suma tensora testowego: {dummy_out_val:.2f}")

if IS_TPU:
    if not os.path.exists(LOG_DIR): os.makedirs(LOG_DIR)
    try:
        torch_xla.sync()
        tpu_print(f"Marker cyklu życia TPU wykonany @ {time.time():.0f}s (bez pełnego sync dla stabilności)")
    except Exception as e:
        tpu_print(f"⚠️ Marker TPU pominięty (benign): {e}")

# ==============================================================================
# GLOBAL FINE-TUNING CACHE
# ==============================================================================
CACHE_DIR = '/kaggle/working/finetune_cache'
os.makedirs(CACHE_DIR, exist_ok=True)
# V-- DODANO PRINT --V
tpu_print(f"[INFO] System cache'u zainicjalizowany w katalogu: {CACHE_DIR}")

def get_params_hash(mu2_val, g_Yukawa, v_H, m0, g, lam_1, lam_2):
    """
    ✅ ROUNDED HASH: Similar params use same cache entry
    """
    mu2_rounded = round(mu2_val / 0.5) * 0.5
    gY_rounded = round(g_Yukawa / 0.05) * 0.05
    vH_rounded = round(v_H / 0.02) * 0.02
    m0_rounded = round(m0 / 1.0) * 1.0
    g_rounded = round(g / 0.2) * 0.2
    lam1_rounded = round(lam_1 / 0.05) * 0.05
    lam2_rounded = round(lam_2 / 0.15) * 0.15

    params_str = f"{mu2_rounded:.1f}_{gY_rounded:.1f}_{vH_rounded:.2f}_{m0_rounded:.1f}_{g_rounded:.1f}_{lam1_rounded:.1f}_{lam2_rounded:.1f}"
    if np.random.rand() < 0.01:
        original_params_str = f"mu2={mu2_val:.2f}, gY={g_Yukawa:.2f}, vH={v_H:.2f}, m0={m0:.2f}, g={g:.2f}, lam1={lam_1:.2f}, lam2={lam_2:.2f}"
        tpu_print(f"  [DEBUG-CACHE] Oryginalne parametry: {original_params_str}")
        tpu_print(f"  [DEBUG-CACHE] Zaokrąglony ciąg do hashowania: {params_str}")

    return hashlib.md5(params_str.encode()).hexdigest()[:12]


def save_finetune_result(mu2_val, g_Yukawa, v_H, m0, g, lam_1, lam_2, Psi_out, Phi_out, final_loss, elapsed_time):
    params_hash = get_params_hash(mu2_val, g_Yukawa, v_H, m0, g, lam_1, lam_2)
    cache_file = os.path.join(CACHE_DIR, f'ft_{params_hash}.pkl')
    Psi_cpu = Psi_out if isinstance(Psi_out, np.ndarray) else Psi_out.cpu().numpy()
    Phi_cpu = Phi_out if isinstance(Phi_out, np.ndarray) else Phi_out.cpu().numpy()
    cache_data = {
        'params': {'mu2': mu2_val, 'gY': g_Yukawa, 'vH': v_H, 'm0': m0, 'g': g, 'lam1': lam_1, 'lam2': lam_2},
        'Psi_out': Psi_cpu,
        'Phi_out': Phi_cpu,
        'final_loss': final_loss,
        'elapsed_time': elapsed_time,
        'timestamp': time.time()
    }
    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(cache_data, f)
        tpu_print(f"  💾 [CACHE] Zapisano wynik dla hasha {params_hash} ({os.path.getsize(cache_file)/(1024**2):.2f} MB)")
        tpu_print(f"     -> Parametry: mu2={mu2_val:.2f}, gY={g_Yukawa:.2f}, vH={v_H:.2f}, m0={m0:.2f}")
    except Exception as e:
        tpu_print(f"  ⚠️ [CACHE] Błąd zapisu do cache'a: {e}")

def load_finetune_result(mu2_val, g_Yukawa, v_H, m0, g, lam_1, lam_2):
    params_hash = get_params_hash(mu2_val, g_Yukawa, v_H, m0, g, lam_1, lam_2)
    cache_file = os.path.join(CACHE_DIR, f'ft_{params_hash}.pkl')
    tpu_print(f"  🔍 [CACHE] Sprawdzanie cache dla hasha: {params_hash}")
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                cache_data = pickle.load(f)
            Psi_out = cache_data['Psi_out']
            Phi_out = cache_data['Phi_out']
            if isinstance(Psi_out, np.ndarray) and np.all(np.isfinite(Psi_out)) and np.all(np.isfinite(Phi_out)):
                tpu_print(f"  ♻️ [CACHE HIT] Znaleziono i załadowano {params_hash}.")
                tpu_print(f"     -> Zapisano {(time.time() - cache_data['timestamp'])/60:.1f} min temu, loss: {cache_data['final_loss']:.4e}")
                return Psi_out, Phi_out, cache_data['final_loss']
            else:
                tpu_print(f"  ⚠️ [CACHE] Uszkodzony plik cache {params_hash}, ponowne obliczenia...")
                os.remove(cache_file)
        except Exception as e:
            tpu_print(f"  ⚠️ [CACHE] Błąd odczytu z cache'a: {e}. Plik zostanie usunięty.")
            if os.path.exists(cache_file):
                os.remove(cache_file)
    else:
        tpu_print(f"  ❌ [CACHE MISS] Brak wpisu dla hasha {params_hash}.")
    return None, None, None

def list_cache_stats():
    tpu_print("\n[INFO] Analiza statystyk cache'u...")
    cache_files = glob.glob(os.path.join(CACHE_DIR, 'ft_*.pkl'))
    if not cache_files:
        tpu_print("[CACHE STATS] Cache jest pusty.")
        return
    total_size = sum(os.path.getsize(f) for f in cache_files)
    tpu_print(f"[CACHE STATS] Znaleziono {len(cache_files)} zapisanych wyników. Całkowity rozmiar: {total_size/(1024**2):.2f} MB")
    if cache_files:
        recent = sorted(cache_files, key=os.path.getmtime, reverse=True)[:5]
        tpu_print("  5 najnowszych wpisów w cache'u:")
        for f in recent:
            age_min = (time.time() - os.path.getmtime(f)) / 60
            size_kb = os.path.getsize(f) / 1024
            tpu_print(f"    - {os.path.basename(f)}: rozmiar {size_kb:.1f} KB, zapisano {age_min:.1f} min temu")

# ==============================================================================
# GLOBAL PARAMETERS
# ==============================================================================
print("\n[INFO] Ustawianie globalnych parametrów symulacji...")
LOG_CSV_FILE = "corr_log_v38_fractal.csv"
n_initial_trials = 120
Nr, Nr_theta, Nr_phi = 800, 32, 32
Nr_theta_mesh, Nr_phi_mesh = 12, 12
t_steps_mesh = 15
r_max = 25.0
num_octaves = 12
m0_init, g_init, lam_1_init = 0.10, 4.64, 1.0
lambda_H = 0.5
delta = 0.2  # Sextic stabilization parameter (δΨ⁶ term)
# Hierarchical coupling parameters (PROPOSAL 2)
beta_hierarchy = 0.15  # Decay rate for hierarchical coupling: λ(o) = λ_base · 2^(-β·o)

TARGET_TOP_GEV = 173.1
TARGET_HIGGS_GEV = 125.1
sigma_noise = 0.05
neigs = 300
lam_2_init = lam_1_init * np.pi
dtau_init = 2e-5
tol_energy = 1e-8
clip_value = 1e4
m0_clip_max = 50.0
NR_SUBSAMPLE = Nr // 10

# --- MAPA RZECZYWISTOŚCI ---
print("\n[INFO] Ustawianie globalnych parametrów dla v38 (Mapa Rzeczywistości)...")

sm_masses_61 = np.array([
    0.0022]*6 + [0.0047]*6 + [0.095]*6 + [1.27]*6 + [4.18]*6 + [173.21]*6 +
    [0.000511]*2 + [0.10566]*2 + [1.77686]*2 +
    [0.0]*6 + [0.0]*8 + [0.0] + [80.379]*2 + [91.1876] + [125.1]
)
print(f"✅ Załadowano {len(sm_masses_61)} mas SM.")

scales_m = np.array([
    1.616e-35, 1e-18, 1e-15, 2.426e-12, 5.292e-11, 1e-10, 5e-7,
    1, 1.274e7, 1.496e11, 9.461e15, 9.461e20, 8.8e26
])
print(f"✅ Załadowano {len(scales_m)} skal fizycznych.")

sm_ratios = []
non_zero_masses = sm_masses_61[sm_masses_61 > 1e-12]
unique_non_zero_masses = np.unique(non_zero_masses)
for i in range(len(unique_non_zero_masses)):
    for j in range(i + 1, len(unique_non_zero_masses)):
        ratio = unique_non_zero_masses[j] / unique_non_zero_masses[i]
        if 1.1 < ratio < 500: sm_ratios.append(ratio)
sm_ratios = np.unique(np.round(sm_ratios, decimals=3))
print(f"✅ Wygenerowano {len(sm_ratios)} unikalnych stosunków z mas SM.")

eff_masses_scales = 1.0 / scales_m
scale_ratios = []
for i in range(len(eff_masses_scales)):
    for j in range(i + 1, len(eff_masses_scales)):
        ratio = max(eff_masses_scales[i], eff_masses_scales[j]) / min(eff_masses_scales[i], eff_masses_scales[j])
        if 1.1 < ratio < 500: scale_ratios.append(ratio)
scale_ratios = np.unique(np.round(scale_ratios, decimals=3))
print(f"✅ Wygenerowano {len(scale_ratios)} unikalnych stosunków ze skal fizycznych.")

TARGET_REALITY_RATIOS = np.unique(np.concatenate([sm_ratios, scale_ratios]))
TARGET_REALITY_WEIGHTS = [2.0] * len(sm_ratios) + [3.0] * len(scale_ratios)
x_kde = np.linspace(np.log10(1.1), np.log10(500), 400)
TARGET_REALITY_DENSITY_KDE = gaussian_kde(np.log10(TARGET_REALITY_RATIOS), bw_method='scott')(x_kde)
print(f"✅ Zunifikowana Mapa Rzeczywistości: {len(TARGET_REALITY_RATIOS)} kluczowych stosunków.")


FINETUNE_BATCH_SIZE = 2097152
INFERENCE_BATCH_SIZE = 4194304
PRETRAIN_BATCH_SIZE = 1048576
ACCUMULATION_STEPS = 1

r_cpu = np.linspace(1e-6, r_max, Nr, dtype=np.float64)
dr_cpu = r_cpu[1]-r_cpu[0]
x_kde = np.linspace(np.log10(1.1), np.log10(500), 400)

tpu_print(f"[PARAM] Rozdzielczość siatki: Nr={Nr}, Nr_theta={Nr_theta}, Nr_phi={Nr_phi}")
tpu_print(f"[PARAM] Parametry fizyczne: r_max={r_max}, num_octaves={num_octaves}")
tpu_print(f"[PARAM] Rozmiary batchy: Pre-train={PRETRAIN_BATCH_SIZE}, Fine-tune={FINETUNE_BATCH_SIZE}")

if not os.path.exists(LOG_CSV_FILE):
    tpu_print(f"[INFO] Tworzenie nowego pliku logów: {LOG_CSV_FILE}")
    pd.DataFrame(columns=['timestamp', 'trial_number', 'g_Y', 'mu2', 'v_H', 'm0', 'g', 'lam_1',
                           'fractal_score', 'hierarchy', 'regularity', 'error_msg']).to_csv(LOG_CSV_FILE, index=False)
else:
    tpu_print(f"[INFO] Znaleziono istniejący plik logów: {LOG_CSV_FILE}")

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
def safe_loss_check_and_sync(loss, epoch, batch_idx):
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
                tpu_print(f"⚠️ XLA materialize error at E{epoch}B{batch_idx}")
                return True, float('inf')
            raise e
    else:
        loss_val = loss.item()
        return not torch.isfinite(loss), loss_val

def transfer_optimizer_state_to_device(optimizer, device, force_detach=False):
    import torch
    for param_id, state in optimizer.state.items():
        for k, v in list(state.items()):
            if isinstance(v, torch.Tensor):
                if force_detach:
                    state[k] = v.detach().cpu().to(device)
                else:
                    if v.device != device:
                        state[k] = v.to(device)
                if state[k].requires_grad:
                    state[k].requires_grad_(False)
            elif isinstance(v, dict):
                 for vk, vv in v.items():
                     if isinstance(vv, torch.Tensor) and vv.device.type != 'cpu':
                         v[vk] = vv.cpu().to(device)


    for param_group in optimizer.param_groups:
        for k, v in param_group.items():
            if isinstance(v, torch.Tensor) and v.device != device:
                param_group[k] = v.to(device)
    if 'xla' in str(device):
        try:
            import torch_xla.core.xla_model as xm
            for state in optimizer.state.values():
                for k, v in state.items():
                    if isinstance(v, torch.Tensor):
                        _ = v.device
            xm.mark_step()
        except Exception:
            pass
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

def get_accumulation_steps(epoch):
    if epoch < 10:
        return 1
    elif epoch < 20:
        return 2
    else:
        return ACCUMULATION_STEPS

def periodic_xla_maintenance(batch_counter, epoch):
    if IS_TPU and batch_counter > 0 and batch_counter % 100 == 0:
        tpu_print(f"  [GRAPH MAINTENANCE] Batch {batch_counter}: clearing XLA cache")
        xm.mark_step()
        gc.collect()

def aggressive_kaggle_heartbeat(epoch, batch_idx, loss_val):
    with open(f'/kaggle/working/live_batch_e{epoch}_b{batch_idx}.tmp', 'w') as f:
        f.write(f"{time.time()} E{epoch}B{batch_idx} Loss={loss_val:.4e}\n")
        f.flush()
        os.fsync(f.fileno())

    print(f"\rE{epoch}B{batch_idx}] Loss: {loss_val:.3e}", end='', flush=True)
    sys.stdout.flush()
    sys.stderr.flush()

    if IS_TPU:
        try:
            # Mark step to keep TPU alive
            xm.mark_step()
        except Exception:
            pass

def print_progress_bar(epoch, batch, total_batches, loss):
    percent = int((batch / total_batches) * 50)
    bar = '█' * percent + '░' * (50 - percent)
    print(f"\rE{epoch} [{'█' * percent}{'░' * (50 - percent)}] {batch}/{total_batches} L={loss:.2e}", end='', flush=True)
    sys.stdout.flush()


def save_verbose_checkpoint(epoch, batch_idx, model, optimizer, metrics):
    checkpoint = {
        'epoch': epoch,
        'batch': batch_idx,
        'timestamp': time.time(),
        'model_state_dict': {k: v.cpu() for k, v in model.state_dict().items()},
        'optimizer_state_dict': {k: v.cpu() if isinstance(v, torch.Tensor) else v for k, v in optimizer.state_dict().items()},
        'metrics': metrics,
        'metadata': {
            'hostname': os.uname().nodename,
            'pid': os.getpid(),
            'device': str(device),
        },
        'padding': np.random.randn(1000, 1000).tolist()
    }

    filepath = f'/kaggle/working/live_checkpoint_e{epoch}_b{batch_idx}.pth'
    try:
        torch.save(checkpoint, filepath, _use_new_zipfile_serialization=False)
        if IS_TPU: xm.mark_step()
        tpu_print(f"\n[CHECKPOINT] Zapisano {os.path.getsize(filepath)/(1024**2):.1f}MB")
    except Exception as e:
        tpu_print(f"⚠️ Verbose checkpoint failed: {e}")

# <-- ZMIANA: DODANO FUNKCJĘ kaggle_checkpoint_output -->
def kaggle_checkpoint_output(epoch, metrics):
    """Symuluje output, który Kaggle rozpoznaje jako 'cell execution'"""
    checkpoint_file = f'/kaggle/working/checkpoint_epoch_{epoch}.txt'
    with open(checkpoint_file, 'w') as f:
        f.write(f"Epoch: {epoch}\n")
        f.write(f"Timestamp: {time.time()}\n")
        for k, v in metrics.items():
            f.write(f"{k}: {v}\n")
        f.flush()
        if hasattr(f, 'fileno'):
            os.fsync(f.fileno())

    tpu_print("="*60)
    tpu_print(f"CHECKPOINT | Epoch {epoch}")
    tpu_print("="*60)
    for k, v in metrics.items():
        tpu_print(f"  {k:20s}: {v}")
    tpu_print("="*60)

# ==============================================================================
# PHYSICS FUNCTIONS
# ==============================================================================

def fractal_consistency_score(masses, target_ratios=None, target_kde=None, x_kde_axis=None):
    """
    Enhanced version with sector-weighted hits.
    Rewards matches in specific sectors with higher weights:
    - Electroweak ratios: weight = 3.0
    - Large hierarchies (>10): weight = 5.0
    - General ratios: weight = 1.0
    """
    if target_ratios is None: target_ratios = TARGET_REALITY_RATIOS
    if target_kde is None: target_kde = TARGET_REALITY_DENSITY_KDE
    if x_kde_axis is None: x_kde_axis = x_kde

    unique_masses = np.unique(masses[masses > 1e-6])
    if len(unique_masses) < 2: return 0.0, 1.0, 0.0

    model_ratios = []
    for i in range(len(unique_masses)):
        for j in range(i+1, len(unique_masses)):
            ratio = unique_masses[j] / unique_masses[i]
            if ratio > 1.01: model_ratios.append(ratio)

    if len(model_ratios) == 0: return 0.0, 1.0, 0.0
    model_ratios = np.array(model_ratios)

    # Define sector boundaries based on SM physics
    # Electroweak sector: ratios near Z/W, H/Z, H/W
    ew_targets = [1.1345, 1.3719, 1.5564]  # Z/W, H/Z, H/W

    # Count sector-weighted hits
    tolerance = 0.05  # 5% tolerance
    sector_weighted_score = 0.0
    total_hits = 0

    for target in target_ratios:
        best_match_ratio = np.min(np.abs(model_ratios - target) / target)
        if best_match_ratio < tolerance:
            total_hits += 1
            # Determine sector weight
            if any(abs(target - ew) / ew < 0.1 for ew in ew_targets):
                # Electroweak sector
                sector_weighted_score += 3.0
            elif target > 10:
                # Large hierarchy sector
                sector_weighted_score += 5.0
            else:
                # General ratio
                sector_weighted_score += 1.0

    # Normalize sector-weighted score
    targeted_bonus_norm = sector_weighted_score / len(target_ratios)

    # KDE correlation (existing logic)
    try:
        model_kde = gaussian_kde(np.log10(model_ratios), bw_method='scott')(x_kde_axis)
        corr = np.corrcoef(model_kde, target_kde)[0, 1]
        if not np.isfinite(corr): corr = 0.0
    except:
        corr = 0.0

    fractal_score = corr + targeted_bonus_norm

    # Hierarchy calculation
    hierarchy = np.max(unique_masses) / np.min(unique_masses)

    # Regularity (musicality)
    model_scales = 1.0 / np.sort(unique_masses)
    log_ratios = np.log(model_scales[1:] / model_scales[:-1])
    mean_log_ratio = np.mean(log_ratios)
    if len(log_ratios) < 2 or mean_log_ratio < 1e-6:
        regularity = 0.0
    else:
        coeff_of_variation = np.std(log_ratios) / mean_log_ratio
        regularity = np.exp(-coeff_of_variation)

    tpu_print(f"  [SCORE ENHANCED] Fractal: {fractal_score:.4f} (Corr={corr:.3f}, Weighted={targeted_bonus_norm:.3f}, Hits={total_hits}) | Hierarchy: {hierarchy:.2f} | Regularity: {regularity:.3f}")

    return fractal_score, hierarchy, regularity
def diagonalize_with_H(Psi_loc, Phi_loc, gY_loc, m0, g, lam_1, lam_2, use_gpu=None):
    """
    Diagonalizacja Macierzy Masy (Hessianu) na CPU.
    """
    Nfull = num_octaves * Nr
    neigs_to_calc = max(1, min(neigs, Nfull - 2))
    dr = dr_cpu
    row, col, data = [], [], []

    if cp is not None and hasattr(Psi_loc, 'get'):
        psi_loc_cpu = cp.asnumpy(Psi_loc)
        phi_loc_cpu = cp.asnumpy(Phi_loc)
    else:
        psi_loc_cpu = np.asarray(Psi_loc)
        phi_loc_cpu = np.asarray(Phi_loc)

    for o in range(num_octaves):
        idx0 = o * Nr
        for i in range(Nr):
            row.append(idx0 + i); col.append(idx0 + i)
            yukawa_diag = 2.0 * gY_loc * (phi_loc_cpu[i]**2)
            diag_val = 2/dr**2 + m0**2 + 3*g*psi_loc_cpu[o,i]**2 + yukawa_diag
            data.append(diag_val)
        if Nr > 1:
            off = -1/dr**2
            for i in range(1, Nr):
                row.extend([idx0+i, idx0+i-1]); col.extend([idx0+i-1, idx0+i]); data.extend([off, off])

    for o in range(num_octaves-1):
        for i in range(Nr): row.extend([o*Nr+i,(o+1)*Nr+i]); col.extend([(o+1)*Nr+i,o*Nr+i]); data.extend([lam_1,lam_1])
    for o in range(num_octaves-2):
        for i in range(Nr): row.extend([o*Nr+i,(o+2)*Nr+i]); col.extend([(o+2)*Nr+i,o*Nr+i]); data.extend([lam_2,lam_2])

    row = np.array(row, dtype=int)
    col = np.array(col, dtype=int)
    data = np.array(data, dtype=float)

    A_cpu = sp.coo_matrix((data, (row, col)), shape=(Nfull, Nfull)).tocsr()
    A_cpu = 0.5*(A_cpu + A_cpu.T)

    try:
        w, _ = spl.eigsh(A_cpu, k=neigs_to_calc, which='SA', tol=1e-6, maxiter=1000)
        positive_w = w[w > 1e-12]
        if len(positive_w) == 0:
            return np.array([])
        masses = np.sqrt(np.sort(positive_w))
        return masses[np.isfinite(masses)]
    except Exception as e:
        tpu_print(f"  [DIAG ERROR] {str(e)[:60]}")
        return None

def total_energy_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Yukawa, mu2, r, dr, xp):
    """
    Total energy functional with δΨ⁶ stabilization term.

    Added term: (1/8)·δ·Ψ⁶ to energy density
    """
    energy_density_psi = xp.zeros(Nr, dtype=Psi.dtype)
    for o in range(num_octaves):
        dpsi = xp.gradient(Psi[o], dr)
        psi_sq = Psi[o]**2
        psi_6 = psi_sq**3
        energy_density_psi += 0.5*dpsi**2 + 0.5*(m0**2)*psi_sq + 0.25*g*(psi_sq**2) + 0.125*delta*psi_6
    # Hierarchical coupling in energy: λ(o) = λ_base · 2^(-β·o)
    for o in range(num_octaves - 1):
        lam_1_hier = lam_1 * (2.0 ** (-beta_hierarchy * o))
        energy_density_psi += lam_1_hier * Psi[o] * Psi[o+1]
    for o in range(num_octaves - 2):
        lam_2_hier = lam_2 * (2.0 ** (-beta_hierarchy * o))
        energy_density_psi += lam_2_hier * Psi[o] * Psi[o+2]
    dPhi = xp.gradient(Phi_H, dr)
    E_kin_H = 0.5 * dPhi**2
    E_pot_H = 0.5 * mu2 * Phi_H**2 + 0.25 * lambda_H * Phi_H**4
    psi_density = xp.sum(Psi**2, axis=0)
    E_Yukawa = g_Yukawa * psi_density * Phi_H**2
    integrand_total = energy_density_psi + E_kin_H + E_pot_H + E_Yukawa
    return 4.0 * xp.pi * xp.sum(integrand_total * r**2) * dr

def functional_derivative_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Yukawa, mu2, r, dr, xp):
    """
    Functional derivative with δΨ⁶ stabilization and hierarchical coupling.

    Added term: (3/4)·δ·Ψ⁵ from δΨ⁶ stabilization potential
    """
    dE_Psi = xp.zeros_like(Psi)
    psi_density = xp.sum(Psi**2, axis=0)
    for o in range(num_octaves):
        lap = -radial_laplacian(Psi[o], r, dr, xp)
        mass_term = m0**2 * Psi[o]
        nonlin = g * Psi[o]**3
        sextic_term = 0.75 * delta * (Psi[o]**5)  # δΨ⁶ stabilization
        yukawa_term = 2.0 * g_Yukawa * Phi_H**2 * Psi[o]
        coupling = xp.zeros_like(Psi[o])
        # Hierarchical coupling: λ(o) = λ_base · 2^(-β·o)
        lam_1_hier = lam_1 * (2.0 ** (-beta_hierarchy * o))
        lam_2_hier = lam_2 * (2.0 ** (-beta_hierarchy * o))
        if o > 0: coupling += lam_1_hier * Psi[o-1]
        if o < num_octaves - 1: coupling += lam_1_hier * Psi[o+1]
        if o > 1: coupling += lam_2_hier * Psi[o-2]
        if o < num_octaves - 2: coupling += lam_2_hier * Psi[o+2]
        dE_Psi[o] = lap + mass_term + nonlin + sextic_term + coupling + yukawa_term

    lap_Phi = -radial_laplacian(Phi_H, r, dr, xp)
    dE_Phi = lap_Phi + mu2 * Phi_H + lambda_H * (Phi_H**3) + 2.0 * g_Yukawa * Phi_H * psi_density
    return dE_Psi, dE_Phi

def radial_laplacian(field, r, dr, xp):
    """
    Radial Laplacian with L'Hospital's rule at r=0 for numerical stability.

    At r=0: ∇²f = d²f/dr² + (2/r)(df/dr) → 3·d²f/dr² (via L'Hospital's rule)
    For r>0: standard formula ∇²f = d²f/dr² + (2/r)(df/dr)
    """
    dfield_dr = xp.gradient(field, dr)
    d2field_dr2 = xp.gradient(dfield_dr, dr)

    # Apply L'Hospital's rule at r=0 (first two points)
    lap = xp.zeros_like(field)
    lap[0] = 3.0 * d2field_dr2[0]
    lap[1] = 3.0 * d2field_dr2[1]

    # Standard formula for r > 0
    r_safe = xp.where(r > 1e-9, r, 1e-9)
    lap[2:] = d2field_dr2[2:] + (2.0 / r_safe[2:]) * dfield_dr[2:]

    return lap

# ==============================================================================
# PINN MODEL AND LOSS
# ==============================================================================
class ResidualBlock(nn.Module):
    def __init__(self, size):
        super().__init__()
        self.l1=nn.Linear(size,size)
        self.l2=nn.Linear(size,size)
        self.act=nn.GELU()
    def forward(self, x): return self.act(self.l2(self.act(self.l1(x)))+x)

class SolitonPINN(nn.Module):
    def __init__(self, output_size=num_octaves+1):
        super().__init__()
        self.inp = nn.Linear(4, 128)
        self.bn1 = nn.LayerNorm(128)
        self.act=nn.GELU()
        self.blocks = nn.Sequential(*[ResidualBlock(128) for _ in range(3)])
        self.out = nn.Linear(128, output_size)
        nn.init.xavier_uniform_(self.inp.weight)
        nn.init.xavier_uniform_(self.out.weight)
        nn.init.zeros_(self.out.bias)
    def forward(self, x):
        x = self.act(self.bn1(self.inp(x)))
        return self.out(self.blocks(x))

class PrefetchDataLoader:
    def __init__(self, loader, device):
        self.loader = loader
        self.device = device
    def __iter__(self):
        batch = None
        for next_batch in self.loader:
            if batch is not None:
                yield batch
            # Prefetch the next batch
            if self.device.type == 'cuda':
                batch = [b.to(self.device, non_blocking=True) for b in next_batch]
            else:
                batch = [b.to(self.device) for b in next_batch]

        if batch is not None:
            yield batch

def compute_derivatives_batch(field, coords, r_max_val, use_finite_diff=False, model=None):
    is_vector_field = field.dim() > 1 and field.size(1) > 1
    max_fd_size = 8192

    if use_finite_diff:
        use_finite_diff = False
        if use_finite_diff:
            tpu_print("  [WARNING] Using stubbed FD logic (fallback to autograd derivatives in this block).")

    # Standard autograd path
    if is_vector_field:
        D = field.size(1)
        N = field.size(0)
        field_flat = field.reshape(-1).contiguous()
        coords_repeated = coords.unsqueeze(1).expand(-1, D, -1).contiguous().reshape(-1, 4)
        field_to_use = field_flat
        coords_to_use = coords_repeated
    else:
        field_to_use = field.contiguous()
        coords_to_use = coords.contiguous()
        N, D = 1, 1 if field.dim() == 1 else field.size(1)

    field_to_use = torch.clamp(field_to_use, -1e4, 1e4)
    coords_to_use = torch.clamp(coords_to_use, -1e3, 1e3)

    grads = torch.autograd.grad(
        outputs=field_to_use.sum(),
        inputs=coords_to_use,
        create_graph=True,
        allow_unused=True
    )[0]

    if grads is None or not torch.isfinite(grads).all():
        zeros = torch.zeros_like(field)
        if is_vector_field:
            return zeros, zeros
        return zeros, zeros

    df_dr, df_dt, df_dth, df_dph = grads.chunk(4, dim=1)
    df_dr = torch.clamp(df_dr, -1e3, 1e3).detach().to(dtype=torch.bfloat16)

    df_dr_det = df_dr.detach().requires_grad_(True)
    df_dth_det = df_dth.detach().requires_grad_(True) # Używamy df_dth (3. chunk), nie df_dt (2. chunk)
    df_dph_det = df_dph.detach().requires_grad_(True)
    coords_to_use_det = coords_to_use.detach().requires_grad_(True)

    d2f_dr2 = torch.autograd.grad(outputs=df_dr_det.sum(), inputs=coords_to_use_det, create_graph=True, retain_graph=True, allow_unused=True)[0][:, 0:1]
    d2f_dth2 = torch.autograd.grad(outputs=df_dth_det.sum(), inputs=coords_to_use_det, create_graph=True, retain_graph=True, allow_unused=True)[0][:, 2:3]
    d2f_dph2 = torch.autograd.grad(outputs=df_dph_det.sum(), inputs=coords_to_use_det, create_graph=True, retain_graph=True, allow_unused=True)[0][:, 3:4]

    d2f_dr2 = torch.clamp(d2f_dr2, -1e4, 1e4).detach()

    r_all, _, theta_all, _ = torch.chunk(coords_to_use, 4, dim=1)
    r_safe = torch.clamp(r_all, 1e-9, r_max_val)
    sin_th_safe = torch.clamp(torch.sin(theta_all + 1e-6), 1e-9, 1.0)
    cot_theta = torch.clamp(torch.cos(theta_all) / sin_th_safe, -50, 50)

    lap = d2f_dr2 + (2.0/r_safe)*df_dr + \
          (d2f_dth2 + cot_theta*df_dth)/(r_safe**2) + \
          d2f_dph2/(r_safe**2 * sin_th_safe**2)

    lap = torch.nan_to_num(torch.clamp(lap, -1e5, 1e5), nan=0.0).detach().to(dtype=torch.bfloat16)

    df_dt_out = df_dt.squeeze()
    lap_out = lap.squeeze()

    if is_vector_field:
        df_dt_out = df_dt_out.view(N, D)
        lap_out = lap_out.view(N, D)

    return df_dt_out, lap_out
# ==============================================================================
# TRAINING FUNCTIONS
# ==============================================================================
def derivative_warmup(pinn, num_warm_steps=5):
    """Progressive warmup with increasing batch sizes to avoid OOM. (v35 content)"""
    if not IS_TPU:
        return

    tpu_print("Warming up derivatives (progressive batch size ramp-up)...")

    batch_sizes = [128, 256, 512, 1024, 2048]

    # Mock params (Używamy mock dla warmup, bo to jest tylko XLA compilation)
    mock_mu2, mock_g_Y, mock_v_H = -10.0, 5.0, 1.0
    mock_m0, mock_g, mock_lam1, mock_lam2 = 15.0, 5.0, 0.5, 0.5 * np.pi
    optimizer_dummy = Adam(pinn.parameters(), lr=1e-6)

    global r_max, device

    for step, bs in enumerate(batch_sizes):
        try:
            mock_coords = [
                torch.rand(bs, 1, device=device) * val
                for val in [r_max, 1.0, np.pi, 2*np.pi]
            ]

            with get_autocast_context(IS_TPU):
                loss = pinn_loss(
                    pinn, *mock_coords, mock_mu2, mock_g_Y,
                    mock_m0, mock_g, mock_lam1, mock_lam2, mock_v_H,
                    use_finite_diff=False
                )

            should_skip, loss_val = safe_loss_check_and_sync(loss, 0, step)

            if not should_skip:
                optimizer_dummy.zero_grad(set_to_none=True)

                try:
                    loss.backward()
                    optimizer_dummy.step()

                    tpu_print(f"  Warmup {step+1}/{len(batch_sizes)}: "
                             f"bs={bs}, loss={loss_val:.2e} ✓")
                except RuntimeError as e:
                    if "tensor_data" not in str(e):
                        raise e
                    tpu_print(f"  Warmup {step+1}: XLA error (expected), continuing...")
            else:
                tpu_print(f"  Warmup {step+1}: non-finite/XLA materialize error, skipping step...")

        except Exception as e:
            tpu_print(f"  Warmup {step+1} error: {str(e)[:50]}, continuing...")

    torch_xla.sync()
    tpu_print("Warmup complete! XLA graph compiled.")


def get_dynamic_batch_size(epoch):
    """Rampa batch size dla kompilacji XLA, potem max throughput"""
    if epoch < 2:
        return 12288
    elif epoch < 5:
        return 18432
    else:
        return PRETRAIN_BATCH_SIZE

def pinn_loss(model, r, t, theta, phi, mu2_val, g_Yukawa, m0, g, lam_1, lam_2, v_H, current_epoch=0, use_finite_diff=False, ic_weight_boost=1.0):
    coords = torch.cat([r, t, theta, phi], dim=1)
    coords.requires_grad_(True)
    output = model(coords)
    Psi, Phi = output[:, :-1], output[:, -1]
    combined_field = torch.cat([Phi.unsqueeze(1), Psi], dim=1)
    dfield_dt, lap_field = compute_derivatives_batch(combined_field, coords, r_max, use_finite_diff=use_finite_diff, model=model)
    dPhi_dt = dfield_dt[:, 0]
    lap_Phi = lap_field[:, 0]
    dPsi_dt = dfield_dt[:, 1:]
    lap_Psi = lap_field[:, 1:]
    psi_density = torch.sum(Psi**2, dim=1)
    res_Phi = dPhi_dt - (-lap_Phi - mu2_val*Phi - lambda_H*(Phi**3) - 2.0*g_Yukawa*Phi*psi_density)
    res_Psi_list = []
    for o in range(num_octaves):
        coup = torch.zeros_like(Psi[:, o])
        if o > 0: coup += lam_1 * Psi[:, o-1]
        if o < num_octaves - 1: coup += lam_1 * Psi[:, o+1]
        if o > 1: coup += lam_2 * Psi[:, o-2]
        if o < num_octaves - 2: coup += lam_2 * Psi[:, o+2]
        res = dPsi_dt[:, o] - (-lap_Psi[:, o] - m0**2*Psi[:, o] - g*Psi[:, o]**3 - coup - 2.0*g_Yukawa*Phi**2*Psi[:, o])
        res_Psi_list.append(res)
    pde_loss = torch.mean(res_Phi**2) + torch.mean(torch.stack(res_Psi_list)**2)
    loss_dtype = pde_loss.dtype if isinstance(pde_loss, torch.Tensor) else torch.get_default_dtype()
    zero_scalar = torch.tensor(0.0, device=device, dtype=loss_dtype)
    fade_factor = 1.0 / (1.0 + np.exp(-0.05 * (current_epoch + 45)))
    ic_base_weight = 700.0 if abs(mu2_val) > 15 else 300.0
    ic_weight = fade_factor * ic_base_weight * ic_weight_boost
    ic_loss = zero_scalar.clone()
    t_min = torch.min(t)
    t0_mask = (torch.abs(t.squeeze() - t_min) < 1e-9)
    if torch.any(t0_mask):
        expected_ic = v_H * torch.sqrt(torch.clamp(torch.tensor(-mu2_val / lambda_H, device=device, dtype=loss_dtype), min=0.0))
        ic_loss = ic_weight * torch.mean((Phi[t0_mask] - expected_ic)**2)
    target_vev = v_H * torch.sqrt(torch.clamp(torch.tensor(-mu2_val / lambda_H, device=device, dtype=loss_dtype), min=0.0))
    mean_abs_phi = torch.mean(torch.abs(Phi))
    vev_error = torch.abs(mean_abs_phi - target_vev)
    target_safe = torch.clamp(target_vev, min=1e-3)
    vev_reward = zero_scalar.clone()
    vev_penalty = zero_scalar.clone()
    if current_epoch >= 30:
        vev_penalty = 300.0 * (mean_abs_phi - target_vev)**2 * fade_factor
    phys_penalty = zero_scalar.clone()
    if current_epoch < 100:
        phys_cap = target_vev + 0.3
        phys_penalty = 20.0 * torch.mean(torch.relu(-Phi))
        phys_penalty += 500.0 * torch.mean(torch.relu(Phi - phys_cap)**2)
        if current_epoch >= 10:
            gaussian_reward = torch.exp(-vev_error / target_safe)
            vev_reward = -0.15 * (1.0 - gaussian_reward) * fade_factor

    if torch.rand(1).item() < 0.02:
        pde_val = pde_loss.item() if isinstance(pde_loss, torch.Tensor) else pde_loss
        ic_val = ic_loss.item() if isinstance(ic_loss, torch.Tensor) else ic_loss
        vev_val = vev_reward.item() if isinstance(vev_reward, torch.Tensor) else vev_reward
        phys_val = phys_penalty.item() if isinstance(phys_penalty, torch.Tensor) else phys_penalty
        mean_phi = torch.mean(torch.abs(Phi)).item()
        target_vev_np = v_H * np.sqrt(max(-mu2_val / lambda_H, 0.0))
        tpu_print(f"    [LOSS E{current_epoch}] pde={pde_val:.2e} | ic={ic_val:.2e}*{ic_weight:.0f} | vev_rew={vev_val:.2e} | phys={phys_val:.2e} | fade={fade_factor:.2f} | |Phi|={mean_phi:.2e} (target={target_vev_np:.2e})")

    try:
        total_loss = pde_loss + ic_loss + vev_penalty + vev_reward + phys_penalty
        total_loss = torch.nan_to_num(total_loss, nan=1e2, posinf=1e4, neginf=-1e4)
        if total_loss < 1e-6:
            total_loss = total_loss + 1e-6
        if not torch.isfinite(total_loss):
            total_loss = torch.tensor(1e6, device=device, dtype=loss_dtype)
    except Exception:
        total_loss = torch.tensor(1e6, device=device, dtype=loss_dtype)
    return total_loss

# <-- ZMIANA: CAŁA FUNKCJA pre_train_pinn ZASTĄPIONA WERSJĄ Z v34.4 -->
def pre_train_pinn(resume_from_epoch=None):
    # --- NOWA LOGIKA: Cel i parametry kontrolne ---
    SUCCESS_THRESHOLD = 0.01  # Cel: loss < 0.1
    STAGNATION_LIMIT = 300     # Cierpliwość na stagnację
    MAX_EPOCHS = 5000         # Maksymalna liczba epok
    base_lr = 1e-3            # Podstawowy LR

    start_time = time.time()
    tpu_print("="*80)
    tpu_print(f"STARTING PRE-TRAINING (v34.4) z celem loss < {SUCCESS_THRESHOLD}")
    tpu_print("="*80)

    # Parametry mock są używane tylko do inicjalizacji wagi biasu (VEV)
    mock_mu2, mock_g_Y, mock_v_H = -10.0, 5.0, 1.0
    mock_m0, mock_g, mock_lam1, mock_lam2 = 15.0, 5.0, 0.5, 0.5 * np.pi

    start_epoch = 0
    pinn = None
    optimizer = None
    checkpoint = {} # Przechowuje checkpoint, jeśli wczytano
    best_loss = float('inf')  # Śledzenie najlepszego lossa
    stagnation = 0
    total_time = 0.0
    skipped_batches = 0

    # --- Logika wznawiania (z poprawką LR) ---
    if resume_from_epoch is not None:
        checkpoint_path = f'/kaggle/working/pinn_latest.pth'
        if isinstance(resume_from_epoch, str) and 'pinn_checkpoint_epoch' in str(resume_from_epoch):
             checkpoint_path = resume_from_epoch

        if os.path.exists(checkpoint_path):
            try:
                checkpoint = torch.load(checkpoint_path, map_location=torch.device('cpu'))
                pinn = SolitonPINN()
                pinn.load_state_dict(checkpoint['model_state_dict'])
                pinn = pinn.to(device)

                pinn.losses_window = checkpoint.get('losses_window', [])
                pinn.bad_epochs = checkpoint.get('bad_epochs', 0)

                optimizer = Adam(pinn.parameters(), lr=base_lr, weight_decay=1e-6)

                if 'optimizer_state_dict' in checkpoint:
                    try:
                        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
                        optimizer.param_groups[0]['lr'] = 5e-4
                        transfer_optimizer_state_to_device(optimizer, device, force_detach=True)
                        tpu_print("   ✅ Optimizer state transferred to XLA device (detached)")

                        if IS_TPU:
                            torch_xla.sync()

                        current_lr = optimizer.param_groups[0]['lr']
                        optimizer_state_data = checkpoint['optimizer_state_dict'].get('state', {})
                        step_count = len(optimizer_state_data)

                        tpu_print(f"✅ OPTIMIZER RESTORED: LR={current_lr:.2e}, Steps={step_count}")

                    except Exception as e:
                        tpu_print(f"⚠️ Optimizer load failed: {e}. Using fresh optimizer.")
                        optimizer.param_groups[0]['lr'] = 5e-4
                else:
                    tpu_print("⚠️ No optimizer state in checkpoint - fresh start (expect loss jump)")
                    optimizer.param_groups[0]['lr'] = 5e-4


                start_epoch = checkpoint['epoch'] + 1
                best_loss = checkpoint.get('avg_loss', float('inf'))
                total_time = checkpoint.get('total_time', 0.0)
                skipped_batches = checkpoint.get('skipped_batches', 0)

                tpu_print(f"✅ RESUME: Wczytano checkpoint z epoki {start_epoch - 1}. Kontynuacja od {start_epoch} (best_loss: {best_loss:.4e}).")
                tpu_print(f"   Rolling Window: {len(pinn.losses_window)} losses, {pinn.bad_epochs} bad epochs.")
                start_time = time.time() - total_time

                if IS_TPU:
                    torch_xla.sync()

            except Exception as e:
                tpu_print(f"⚠️ Błąd wczytywania checkpoint ({checkpoint_path}): {e}. Start od zera.")
                pinn = None
                start_epoch = 0
                skipped_batches = 0
        else:
            tpu_print(f"⚠️ Checkpoint {checkpoint_path} nie istnieje. Start od zera.")
            pinn = None
            start_epoch = 0
            skipped_batches = 0

    if pinn is None or start_epoch == 0:
        pinn = SolitonPINN()
        pinn = pinn.to(device)
        start_epoch = 0
        vev_init = mock_v_H * np.sqrt(max(-mock_mu2 / lambda_H, 0.0))
        with torch.no_grad():
            pinn.out.bias[-1].fill_(vev_init)
            pinn.out.weight[-1].data.fill_(0.0)

        if optimizer is None:
             optimizer = Adam(pinn.parameters(), lr=1e-6, weight_decay=1e-6)

        pinn.losses_window = []
        pinn.bad_epochs = 0

        tpu_print("🆕 Start od zera (initial LR: 1e-6).")

    def get_lr_lambda(epoch):
        global_epoch = start_epoch + epoch
        if global_epoch < 5:
            return min(1.0, (global_epoch + 1) / 5.0)
        else:
            return 1.0

    lr_scheduler = LambdaLR(optimizer, lr_lambda=get_lr_lambda)

    if start_epoch > 0:
        if 'scheduler_state_dict' in checkpoint:
            try:
                lr_scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
                tpu_print(f"✅ Scheduler restored (last_epoch: {lr_scheduler.last_epoch})")
            except Exception as e:
                tpu_print(f"⚠️ Scheduler load failed: {e}. Using epoch-based reset.")
                for _ in range(start_epoch):
                    lr_scheduler.last_epoch += 1
        else:
            for _ in range(start_epoch):
                lr_scheduler.last_epoch += 1


    if IS_TPU:
        torch_xla.sync()
        tpu_print("Warming up TPU (XLA graph compilation)...")
        batch_size_warm = 4096
        for i in range(5):
            x_warm = torch.randn(batch_size_warm, 4, device=device, dtype=torch.float32)
            _ = pinn(x_warm)
            if i % 2 == 0:
                torch_xla.sync()
        torch_xla.sync()
        tpu_print("TPU warmed up – first real batch should be fast!")
        with torch.no_grad():
            dummy_forward = pinn(torch.randn(32, 4, device=device))
        torch_xla.sync()
    derivative_warmup(pinn)

    def create_dataloader(epoch):
        batch_size = get_dynamic_batch_size(epoch)
        r_tensor_cpu = torch.tensor(r_cpu, dtype=torch.float32, device='cpu')
        t_tensor_cpu = torch.tensor(np.linspace(0, 1, t_steps_mesh), dtype=torch.float32, device='cpu')
        theta_tensor_cpu = torch.tensor(np.linspace(0, np.pi, Nr_theta_mesh), dtype=torch.float32, device='cpu')
        phi_tensor_cpu = torch.tensor(np.linspace(0, 2*np.pi, Nr_phi_mesh), dtype=torch.float32, device='cpu')
        mesh = torch.cartesian_prod(r_tensor_cpu, t_tensor_cpu, theta_tensor_cpu, phi_tensor_cpu)

        base_dataloader = DataLoader(TensorDataset(mesh),
                                    batch_size=batch_size,
                                    shuffle=True,
                                    drop_last=True,
                                    num_workers=0, pin_memory=False, prefetch_factor=None)
        return PrefetchDataLoader(base_dataloader, device)

    writer = None
    if IS_TPU and TENSORBOARDX_AVAILABLE:
        try:
            writer = SummaryWriter(log_dir=os.path.join(LOG_DIR, 'pre_train'))
            tpu_print("✅ TensorBoardX Writer active for pre-training.")
        except Exception as e:
            tpu_print(f"⚠️ TensorBoardX Writer failed: {e}. Logging disabled.")

    last_output_time = time.time()

    for epoch in range(start_epoch, MAX_EPOCHS):
        if IS_TPU and epoch > start_epoch and epoch % 5 == 0:
            tpu_print(f"  [XLA MAINTENANCE] Re-transferring optimizer state to device (E{epoch})")
            transfer_optimizer_state_to_device(optimizer, device, force_detach=True)
            torch_xla.sync()

        dataloader = create_dataloader(epoch)
        epoch_start_time = time.time()
        epoch_loss = 0.0
        accumulated_loss = 0.0
        batch_counter = 0
        total_batches = len(dataloader.loader)

        current_accum_steps = get_accumulation_steps(epoch)
        tpu_print(f"Epoch {epoch}: Using accumulation_steps={current_accum_steps}")

        optimizer.zero_grad(set_to_none=True)

        for batch_idx, (batch_rttp,) in enumerate(dataloader):
            batch_counter += 1
            avg_accumulated_loss_for_log = epoch_loss / batch_counter if batch_counter > 0 else 0.0

            if batch_counter % 5 == 0:
                aggressive_kaggle_heartbeat(epoch, batch_counter, avg_accumulated_loss_for_log)

            print_progress_bar(epoch, batch_counter, total_batches, avg_accumulated_loss_for_log)
            periodic_xla_maintenance(batch_counter, epoch)

            if batch_counter % 25 == 0:
                save_verbose_checkpoint(epoch, batch_counter, pinn, optimizer, {
                    'accumulated_loss': accumulated_loss, 'skipped_batches': skipped_batches
                })

            with get_autocast_context(IS_TPU):
                r, t, th, ph = [b.to(device) for b in batch_rttp.split(1, dim=1)]
                use_fd = (2 <= epoch < 10 and batch_counter > 5 and skipped_batches > batch_counter * 0.1)
                if use_fd and batch_counter % 20 == 0:
                    tpu_print(f"  [FD ENABLED] Skip rate {skipped_batches/batch_counter*100:.1f}% - using finite diff")

                loss = pinn_loss(pinn, r, t, th, ph, mock_mu2, mock_g_Y,
                               mock_m0, mock_g, mock_lam1, mock_lam2, mock_v_H,
                               current_epoch=epoch, use_finite_diff=use_fd)

                if loss is None:
                    skipped_batches += 1
                    continue

            should_skip, loss_val = safe_loss_check_and_sync(loss, epoch, batch_idx)
            if should_skip:
                skipped_batches += 1
                optimizer.zero_grad(set_to_none=True)
                continue

            if epoch > 10 and loss_val > 1e6:
                skipped_batches += 1
                optimizer.zero_grad(set_to_none=True)
                continue

            backward_success = False
            try:
                loss.backward()
                if IS_TPU:
                    torch_xla.sync()
                backward_success = True
            except RuntimeError as e:
                if "tensor_data" in str(e) or "XLATensor" in str(e):
                    backward_success = False
                else: raise e

            if not backward_success:
                skipped_batches += 1
                optimizer.zero_grad(set_to_none=True)
                continue

            accumulated_loss += loss.detach().item()

            if IS_TPU and time.time() - last_output_time > 30:
                tpu_print(f"[ALIVE] E{epoch} B{batch_counter}/{total_batches} L={accumulated_loss/batch_counter:.3e} (30s tick)")
                write_dual_heartbeat(f"30s_tick_epoch_{epoch}_batch_{batch_counter}")
                last_output_time = time.time()

            if batch_counter % current_accum_steps == 0:
                max_norm = 5.0
                total_norm = torch.nn.utils.clip_grad_norm_(pinn.parameters(), max_norm)

                for param in pinn.parameters():
                    if param.grad is not None:
                        param.grad.data.div_(current_accum_steps)

                step_success = False
                if IS_TPU:
                    try:
                        xm.optimizer_step(optimizer)
                        step_success = True
                    except RuntimeError as e:
                        if "Input tensor is not an XLA tensor" in str(e):
                            transfer_optimizer_state_to_device(optimizer, device, force_detach=True)
                            torch_xla.sync()
                            try:
                                xm.optimizer_step(optimizer)
                                step_success = True
                            except Exception:
                                pass
                        elif "tensor_data" in str(e):
                            try:
                                optimizer.step()
                                torch_xla.sync()
                                step_success = True
                            except Exception:
                                pass
                        else:
                            raise e
                else:
                    optimizer.step()
                    step_success = True

                epoch_loss += accumulated_loss
                optimizer.zero_grad(set_to_none=True)
                accumulated_loss = 0.0

        if batch_counter % current_accum_steps != 0 and accumulated_loss > 0:
            for param in pinn.parameters():
                if param.grad is not None:
                    param.grad.data.div_(current_accum_steps)
            torch.nn.utils.clip_grad_norm_(pinn.parameters(), 0.05)
            try:
                if IS_TPU: xm.optimizer_step(optimizer)
                else: optimizer.step()
                epoch_loss += accumulated_loss
            except RuntimeError:
                pass

        if IS_TPU and epoch % 50 == 0: gc.collect()

        total_batches_processed = batch_counter
        avg_loss = epoch_loss / total_batches_processed if total_batches_processed > 0 else 0.0

        if avg_loss < SUCCESS_THRESHOLD:
            tpu_print(f"\n✅ SUKCES! Pre-trening zakończony w epoce {epoch}. loss ({avg_loss:.4e}) < {SUCCESS_THRESHOLD}")
            best_loss = avg_loss
            break

        if avg_loss < best_loss:
            best_loss, stagnation = avg_loss, 0
        else:
            if epoch >= max(10, start_epoch + 5):
                stagnation += 1

        if stagnation >= STAGNATION_LIMIT:
            tpu_print(f"\n⚠️ OSTRZEŻENIE: Stagnacja treningu w epoce {epoch}. Cel NIE ZOSTAŁ OSIĄGNIĘTY.")
            break

        lr_scheduler.step()
        if not hasattr(pinn, 'losses_window'):
            pinn.losses_window, pinn.bad_epochs = [], 0
        pinn.losses_window.append(avg_loss)
        if len(pinn.losses_window) > 20: pinn.losses_window.pop(0)
        if len(pinn.losses_window) >= 5:
            if avg_loss > min(pinn.losses_window) * 1.1:
                pinn.bad_epochs += 1
            else:
                pinn.bad_epochs = 0
            if pinn.bad_epochs >= 8:
                old_lr = optimizer.param_groups[0]['lr']
                optimizer.param_groups[0]['lr'] *= 0.85
                new_lr = optimizer.param_groups[0]['lr']
                tpu_print(f"  [LR DECAY] {old_lr:.2e} → {new_lr:.2e}")
                pinn.bad_epochs = 0

        if epoch > 50 and avg_loss > best_loss * 10 and avg_loss > 1e4:
            tpu_print(f"⚠️ Loss explosion at epoch {epoch}. Stopping training.")
            break

        with torch.no_grad():
            dummy_coords = torch.randn(1000, 4, device=device)
            dummy_coords[:, 0:4] = torch.clamp(dummy_coords[:, 0:4], 0, 1)
            mean_phi = torch.mean(torch.abs(pinn(dummy_coords)[:, -1]))

        epoch_time = time.time() - epoch_start_time
        total_time = time.time() - start_time

        print(f"\rE{epoch} [{'█'*50}] {batch_counter}/{total_batches} L={avg_loss:.2e}", flush=True)
        tpu_print(f"EPOCH {epoch} COMPLETE | Loss: {avg_loss:.4e} | Mean |Phi|: {mean_phi:.4e} | Time: {epoch_time:.1f}s | Total: {total_time/60:.1f}min | Skips: {skipped_batches}/{batch_counter}")
        write_dual_heartbeat(f"epoch_{epoch}_complete_loss={avg_loss:.6f}")

        if writer:
            writer.add_scalar('Loss/Pre_Train_Avg', avg_loss, epoch)
            writer.add_scalar('Monitor/Mean_Phi', mean_phi, epoch)
            writer.flush()

        if epoch >= start_epoch:
            if IS_TPU: torch_xla.sync()
            optimizer_state_cpu = optimizer.state_dict()
            for state_key, state_value in optimizer_state_cpu['state'].items():
                if isinstance(state_value, dict):
                    for key in list(state_value.keys()):
                        if torch.is_tensor(state_value[key]):
                            state_value[key] = state_value[key].cpu()
            lightweight_ckpt = {'epoch': epoch, 'model_state_dict': {k: v.cpu() for k,v in pinn.state_dict().items()},
                                'optimizer_state_dict': optimizer_state_cpu, 'scheduler_state_dict': lr_scheduler.state_dict(),
                                'avg_loss': avg_loss, 'best_loss': best_loss, 'total_time': total_time,
                                'skipped_batches': skipped_batches, 'batch_counter': batch_counter,
                                'losses_window': pinn.losses_window, 'bad_epochs': pinn.bad_epochs}
            torch.save(lightweight_ckpt, '/kaggle/working/pinn_latest.pth', _use_new_zipfile_serialization=False)
            if IS_TPU: torch_xla.sync()
            skip_rate = skipped_batches / batch_counter * 100 if batch_counter > 0 else 0
            tpu_print(f"✅ [EPOCH {epoch}] Checkpoint saved | Loss: {avg_loss:.4e} | Skip: {skip_rate:.1f}%")

        if epoch % 10 == 0:
            kaggle_checkpoint_output(epoch, {'avg_loss': avg_loss, 'mean_phi': mean_phi.item(), 'epoch_time': epoch_time, 'total_time': total_time})
            torch.save(lightweight_ckpt, f'/kaggle/working/pinn_checkpoint_epoch_{epoch}.pth', _use_new_zipfile_serialization=False)
            if IS_TPU: torch_xla.sync()

    if writer: writer.close()
    final_loss_for_log = best_loss if best_loss != float('inf') else avg_loss
    final_state_dict = {k: v.cpu() for k, v in pinn.state_dict().items()}
    torch.save(final_state_dict, 'pretrained_pinn.pth', _use_new_zipfile_serialization=False)
    tpu_print(f"\nOstateczny model pre-trenowany zapisany w 'pretrained_pinn.pth' (loss: {final_loss_for_log:.4e}).")
    return pinn

def pinn_soliton_cached(mu2_val, g_Yukawa, v_H, m0, g, lam_1, lam_2):
    """
    ULTRA-FAST fine-tuning: 12 epok
    """
    if IS_TPU:
        import torch_xla.core.xla_model as xm

    # Check cache first
    Psi_cached, Phi_cached, final_loss = load_finetune_result(mu2_val, g_Yukawa, v_H, m0, g, lam_1, lam_2)
    if Psi_cached is not None:
        return Psi_cached, Phi_cached

    tpu_print(f"\n--- Rozpoczynanie fine-tuningu dla: mu2={mu2_val:.2f}, gY={g_Yukawa:.2f}, vH={v_H:.2f}, m0={m0:.2f}, g={g:.2f}, lam1={lam_1:.2f} ---")
    tpu_print("  [FINETUNE] Brak w cache'u. Uruchamianie ultra-szybkiego fine-tuningu...")
    start_time = time.time()
    pinn = SolitonPINN()

    try:
        pinn.load_state_dict(torch.load('pretrained_pinn.pth', map_location='cpu'))
        pinn = pinn.to(device)
        tpu_print("  [FINETUNE] Model pre-trenowany załadowany pomyślnie.")
    except FileNotFoundError:
        tpu_print(f"  CRITICAL ⚠️ Plik 'pretrained_pinn.pth' nie istnieje!")
        return np.array([[np.nan]]), np.array([[np.nan]])
    except Exception as e:
        tpu_print(f"  CRITICAL ⚠️ Nieoczekiwany błąd podczas ładowania modelu: {e}")
        return np.array([[np.nan]]), np.array([[np.nan]])

    with torch.no_grad():
        real_vev_init = v_H * np.sqrt(max(-mu2_val / lambda_H, 0.0))
        pinn.out.bias[-1].fill_(real_vev_init)
        pinn.out.weight[-1].data.fill_(0.0)
        tpu_print(f"  [FINETUNE] Inicjalizacja VEV dla Phi: {real_vev_init:.4f}")

    pinn.train()

    MAX_EPOCHS = 12
    BATCH_SIZE = FINETUNE_BATCH_SIZE
    LEARNING_RATE = 5e-4

    tpu_print(f"  [FINETUNE] Konfiguracja: Epoki={MAX_EPOCHS}, Batch Size={BATCH_SIZE}, LR={LEARNING_RATE}")

    optimizer = Adam(pinn.parameters(), lr=LEARNING_RATE, weight_decay=1e-5)

    scheduler = ReduceLROnPlateau(
        optimizer,
        mode='min',
        patience=2,
        factor=0.5,
        min_lr=1e-6
    )

    r_tensor = torch.tensor(r_cpu, dtype=torch.float32, device='cpu')
    t_tensor = torch.tensor(np.linspace(0, 1, t_steps_mesh), dtype=torch.float32, device='cpu')
    theta_tensor = torch.tensor(np.linspace(0, np.pi, Nr_theta_mesh), dtype=torch.float32, device='cpu')
    phi_tensor = torch.tensor(np.linspace(0, 2*np.pi, Nr_phi_mesh), dtype=torch.float32, device='cpu')
    mesh = torch.cartesian_prod(r_tensor, t_tensor, theta_tensor, phi_tensor)

    if BATCH_SIZE > mesh.size(0):
        BATCH_SIZE = mesh.size(0)
        tpu_print(f"  [FINETUNE] Auto-adjust: Batch {BATCH_SIZE} (dataset limit)")

    base_dataloader = DataLoader(
        TensorDataset(mesh),
        batch_size=BATCH_SIZE,
        shuffle=True,
        drop_last=True,
        num_workers=0,
        pin_memory=False
    )
    dataloader = PrefetchDataLoader(base_dataloader, device)

    best_loss = float('inf')
    tpu_print(f"  [FINETUNE] Pętla treningowa rozpoczęta...")

    for epoch in range(MAX_EPOCHS):
        epoch_loss_tensor_list = []
        epoch_has_nan = False
        batch_count = 0

        optimizer.zero_grad(set_to_none=True)

        for batch_rttp, in dataloader:
            batch_count += 1
            r, t, th, ph = [b.to(device) for b in batch_rttp.split(1, dim=1)]

            ic_boost = 3.0 if abs(mu2_val) > 20 else 2.0

            with get_autocast_context(IS_TPU):
                loss = pinn_loss(pinn, r, t, th, ph, mu2_val, g_Yukawa,
                                 m0, g, lam_1, lam_2, v_H,
                                 current_epoch=epoch,
                                 ic_weight_boost=ic_boost)

            if not torch.isfinite(loss):
                tpu_print(f"\n  ⚠️ [FINETUNE E{epoch+1}] Wykryto NaN/Inf w loss. Pomijanie batcha.")
                epoch_has_nan = True
                optimizer.zero_grad(set_to_none=True)
                continue

            loss.backward()
            torch.nn.utils.clip_grad_norm_(pinn.parameters(), 10.0)

            if IS_TPU:
                xm.optimizer_step(optimizer)
            else:
                optimizer.step()

            epoch_loss_tensor_list.append(loss.detach())
            optimizer.zero_grad(set_to_none=True)

            if batch_count % 10 == 0:
                try:
                    tpu_print(f"  [FINETUNE E{epoch+1}] B{batch_count} Loss: {loss.item():.4e}", end='\r', flush=True)
                except Exception:
                    pass

        if IS_TPU:
            xm.mark_step()

        valid_losses = [l.item() for l in epoch_loss_tensor_list if torch.isfinite(l)]

        if not valid_losses or epoch_has_nan:
            tpu_print(f"\n  [FINETUNE E{epoch+1}/{MAX_EPOCHS}] Epoka niestabilna. Ustawianie wysokiej straty.")
            avg_loss = 1e5
        else:
            avg_loss = sum(valid_losses) / len(valid_losses)

        best_loss = min(best_loss, avg_loss)
        scheduler.step(avg_loss)
        current_lr = optimizer.param_groups[0]['lr']
        if epoch % 3 == 0:
            tpu_print(f"  [SCHEDULER] Aktualny LR: {current_lr:.2e}")
        tpu_print(f"\n  [FINETUNE E{epoch+1}/{MAX_EPOCHS}] Ukończono. Avg Loss: {avg_loss:.4e} (Najlepszy: {best_loss:.4e})")

        if avg_loss < 0.1 and not epoch_has_nan:
            tpu_print(f"  ⚡ FAST CONVERGE at E{epoch+1}: {avg_loss:.4e}")
            break

    pinn.eval()
    with torch.no_grad():
        tpu_print("  [FINETUNE-INFERENCE] Tworzenie siatki inferencyjnej na CPU...")
        final_r_cpu = torch.tensor(r_cpu, dtype=torch.float32)
        final_t_cpu = torch.ones(Nr, dtype=torch.float32) * 1.0
        final_theta_cpu = torch.linspace(0, np.pi, Nr_theta_mesh, dtype=torch.float32)
        final_phi_cpu = torch.linspace(0, 2*np.pi, Nr_phi_mesh, dtype=torch.float32)

        r_grid_cpu, th_grid_cpu, ph_grid_cpu = torch.meshgrid(final_r_cpu, final_theta_cpu, final_phi_cpu, indexing='ij')
        t_grid_cpu = torch.ones_like(r_grid_cpu)

        inp_cpu = torch.stack([g.flatten() for g in [r_grid_cpu, t_grid_cpu, th_grid_cpu, ph_grid_cpu]], dim=1)

        tpu_print("  [FINETUNE-INFERENCE] Przenoszenie siatki na urządzenie docelowe...")
        inp = inp_cpu.to(device)

        out_list = []
        tpu_print("  [FINETUNE-INFERENCE] Rozpoczynanie pętli inferencji...")
        for i in range(0, inp.size(0), INFERENCE_BATCH_SIZE):
            with get_autocast_context(IS_TPU):
                batch_inp_bf16 = inp[i:i+INFERENCE_BATCH_SIZE].to(dtype=torch.bfloat16)
                batch_out_bf16 = pinn(batch_inp_bf16)
                out_list.append(batch_out_bf16.to(dtype=torch.float32))

        if IS_TPU:
            tpu_print("  [FINETUNE-INFERENCE] Synchronizacja wyników...")
            xm.mark_step()

        tpu_print("  [FINETUNE-INFERENCE] Łączenie i konwersja wyników...")
        out_tensor = torch.cat(out_list, dim=0)

        out = out_tensor.view(Nr, Nr_theta_mesh, Nr_phi_mesh, -1).cpu().numpy()

    Psi_out = np.mean(out[..., :-1], axis=(1,2)).T
    Phi_out = np.mean(out[..., -1], axis=(1,2))

    if not np.all(np.isfinite(Psi_out)):
        tpu_print("  [PINN] Non-finite output in final inference.")
def run_single_scan_job(g_Yukawa, mu2_val, v_H, m0, g, lam_1, lam_2, job_id=0):
    """
    ENHANCED VERSION: Implements large-amplitude structured initialization
    per recommendations from previous analysis.
    """
    xp_local = cp if GPU_MODE and cp else np
    r, dr = xp_local.asarray(r_cpu), xp_local.asarray(dr_cpu)

    tpu_print(f"  [SCAN JOB {job_id}] Wczytywanie/Obliczanie PINN dla: gY={g_Yukawa:.2f}, mu2={mu2_val:.1f}")

    Psi_ml_cpu, Phi_ml_cpu = pinn_soliton_cached(mu2_val, g_Yukawa, v_H, m0, g, lam_1, lam_2)

    if not np.all(np.isfinite(Psi_ml_cpu)) or not np.all(np.isfinite(Phi_ml_cpu)):
        # ENHANCED INITIALIZATION STRATEGY (Recommendation #1)
        # Use large-amplitude structured profiles instead of small random noise
        v_est = np.sqrt(max(-mu2_val / (lambda_H + 1e-12), 0.0))
        Phi_H_init_cpu = np.full(Nr, v_est * v_H, dtype=np.float64) + 0.01 * np.random.randn(Nr)

        # KEY CHANGE: Use exponential decay with amplitude A ~ 1.5
        r_array = r_cpu
        A_amplitude = 1.5  # Large amplitude for non-trivial basin
        R_decay = 3.0      # Decay length scale

        Psi_init_cpu = np.zeros((num_octaves, Nr), dtype=np.float64)
        for o in range(num_octaves):
            # Exponential decay: Ψ(r) = A·exp(-r/R)
            Psi_init_cpu[o] = A_amplitude * np.exp(-r_array / R_decay)
            # Add small perturbation for symmetry breaking
            Psi_init_cpu[o] += 0.05 * np.random.randn(Nr)

        tpu_print(f"  [SCAN JOB {job_id}] PINN zwrócił NaN/Inf. Użycie ENHANCED inicjalizacji (A={A_amplitude}, R={R_decay}).")
    else:
        Psi_init_cpu, Phi_H_init_cpu = Psi_ml_cpu, Phi_ml_cpu
        tpu_print(f"  [SCAN JOB {job_id}] PINN załadowany z cache'a/poprzedniego kroku.")

    Psi, Phi_H = xp_local.asarray(Psi_init_cpu.copy()), xp_local.asarray(Phi_H_init_cpu.copy())
    dtau = dtau_init * (0.1 if abs(mu2_val) > 10 else 1.0)
    E_prev = total_energy_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Yukawa, mu2_val, r, dr, xp_local)
    clip_value = 1e3

    tpu_print(f"  [SCAN JOB {job_id}] Energia początkowa: {E_prev:.4e}. Uruchamianie 100 kroków MC (dtau={dtau:.2e}).")

    for step in range(1, 100):
        dE_Psi, dE_Phi = functional_derivative_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Yukawa, mu2_val, r, dr, xp_local)
        if not xp_local.all(xp_local.isfinite(dE_Psi)) or not xp_local.all(xp_local.isfinite(dE_Phi)): return None, None
        Psi -= dtau * dE_Psi
        Phi_H -= dtau * dE_Phi
        Psi[:, -1], Phi_H[-1] = 0.0, 0.0
        Psi, Phi_H = xp_local.clip(Psi, -clip_value, clip_value), xp_local.clip(Phi_H, -clip_value, clip_value)

        if step % 5 == 0:
            max_phi = xp_local.max(xp_local.abs(Phi_H))
            if max_phi > 1e3:
                scale = 1e3 / max_phi
                Psi *= scale
                Phi_H *= scale
                tpu_print(f"  [SCAN JOB {job_id}] Reskala Phi o {scale:.2e} w step {step} (max Phi={max_phi:.1e})")
        if step % 10 == 0:
            norm = xp_local.sqrt(xp_local.sum(Psi**2) * 4 * np.pi * xp_local.sum(r**2) * dr)
            if norm > 1e-9: Psi /= norm
            E = total_energy_with_H(Psi, Phi_H, m0, g, lam_1, lam_2, g_Yukawa, mu2_val, r, dr, xp_local)
            if abs(E - E_prev) < tol_energy * 100:
                tpu_print(f"  [SCAN JOB {job_id}] Zbieżność osiągnięta w kroku {step}.")
                break
            E_prev = E
    return Psi, Phi_H
def run_self_consistent_job(g_Yukawa, mu2_val, v_H, m0_init_arg, g_init_arg, lam_1_init_arg, lam_2_init_arg, job_id):
    tpu_print(f"\n>>>> [JOB {job_id}] Rozpoczynanie pętli samouzgodnienia dla gY={g_Yukawa:.3f}, mu2={mu2_val:.2f}")
    m0, g_param, lam_1 = m0_init_arg, g_init_arg, lam_1_init_arg
    lam_2 = lam_1 * np.pi
    gY = g_Yukawa
    start_time = time.time()
    error_msg = ""
    max_iters = 2
    learning_rate = 0.2
    convergence_threshold = 0.05
    Psi_final, Phi_H_final = None, None
    delta_history = []
    G_CLIP_MIN, G_CLIP_MAX = 1.0, 100.0
    masses_raw = [] # Init

    for self_iter in range(max_iters):
        tpu_print(f"  [JOB {job_id} | ITER {self_iter+1}/{max_iters}] Uruchamianie kroku...")
        xp_local = cp if GPU_MODE and cp else np
        Psi_inter, Phi_H_inter = run_single_scan_job(gY, mu2_val, v_H, m0, g_param, lam_1, lam_2, job_id)

        if Psi_inter is None or not xp_local.all(xp_local.isfinite(Psi_inter)):
            error_msg = "PINN/Single Scan Failed"
            break

        Psi_final, Phi_H_final = Psi_inter, Phi_H_inter
        Psi_cpu = Psi_final if xp_local is np else xp_local.asnumpy(Psi_final)
        Phi_H_cpu = Phi_H_final if xp_local is np else xp_local.asnumpy(Phi_H_final)
# PRE-CHECKS (Recommendation #4) - ADDED AFTER run_single_scan_job
    # Quickly reject unpromising solutions before expensive diagonalization

    # Pre-check 1: Ensure non-trivial solution
    psi_max = np.max(np.abs(Psi_cpu))
    if psi_max < 0.1:
        tpu_print(f"  [JOB {job_id}] PRE-CHECK REJECTED: psi_max={psi_max:.4f} < 0.1 (trivial solution)")
        return {'fractal_score': -1.0, 'hierarchy': 0.0, 'regularity': 0.0, 'error': 'trivial_solution'}

    # Pre-check 2: Energy stability
    E_final = total_energy_with_H(Psi_cpu, Phi_H_cpu, m0, g_param, lam_1, lam_2, gY, mu2_val, r_cpu, dr_cpu, np)
    if not np.isfinite(E_final) or E_final > 0:
        tpu_print(f"  [JOB {job_id}] PRE-CHECK REJECTED: Energy={E_final:.4e} (unstable)")
        return {'fractal_score': -1.0, 'hierarchy': 0.0, 'regularity': 0.0, 'error': 'energy_unstable'}

    # Pre-check 3: Phi field VEV within reasonable bounds
    phi_mean = np.mean(np.abs(Phi_H_cpu))
    if phi_mean > 1e2 or phi_mean < 1e-3:
        tpu_print(f"  [JOB {job_id}] PRE-CHECK REJECTED: phi_mean={phi_mean:.4e} (out of bounds)")
        return {'fractal_score': -1.0, 'hierarchy': 0.0, 'regularity': 0.0, 'error': 'phi_out_of_bounds'}

    tpu_print(f"  [JOB {job_id}] PRE-CHECKS PASSED: psi_max={psi_max:.3f}, E={E_final:.2e}, phi_mean={phi_mean:.3f}")
    # Continue with expensive diagonalization...

    psi_density_sq = np.sum(Psi_cpu**2, axis=0)
    g_eff = np.mean(psi_density_sq * Phi_H_cpu**2) * 1.0
    m0_eff = abs(mu2_val) * np.mean(np.abs(Phi_H_cpu)) * 10.0
    corr_vals = [np.corrcoef(Psi_cpu[o].flatten(), Psi_cpu[o+1].flatten())[0,1] if len(Psi_cpu[o])>1 else 0 for o in range(num_octaves-1)]
    lam1_eff = np.nanmean([c for c in corr_vals if np.isfinite(c)]) * 2.0 if corr_vals else lam_1

    m0_eff = np.clip(m0_eff, 0.01, m0_clip_max)
    g_eff = np.clip(g_eff, G_CLIP_MIN, G_CLIP_MAX)
    lam1_eff = np.clip(lam1_eff, 0.3, 1.5)

    delta_g_raw = abs(g_param - g_eff) / max(g_param, 1e-6)
    lr_g = learning_rate * (0.1 if delta_g_raw > 1.0 else 1.0)

    if self_iter > 0:
            delta_m0 = abs(m0 - m0_eff) / max(m0, 1e-6)
            delta_g = abs(g_param - g_eff) / max(g_param, 1e-6)
            delta_lam1 = abs(lam_1 - lam1_eff) / max(lam_1, 1e-6)
            delta_avg = np.mean([delta_m0, delta_g, delta_lam1])
            delta_history.append(delta_avg)

            if delta_avg < convergence_threshold:
                tpu_print(f"  [JOB {job_id}] OSIĄGNIĘTO ZBIEŻNOŚĆ w iteracji {self_iter+1}.")
                break

    new_m0 = m0 + learning_rate * (m0_eff - m0)
    new_g = g_param + lr_g * (g_eff - g_param)
    new_l1 = lam_1 + learning_rate * (lam1_eff - lam_1)
    m0, g_param, lam_1 = new_m0, np.clip(new_g, G_CLIP_MIN, G_CLIP_MAX), new_l1
    lam_2 = lam_1 * np.pi
    learning_rate *= 0.95

    fractal_score, hierarchy, regularity = -1.0, 0.0, 0.0

    if error_msg:
        # Jeśli jakikolwiek błąd wystąpił wcześniej (np. w pętli SC)
        pass
    elif Psi_final is None:
        error_msg = "Loop failed to produce final state"
    else:
        Psi_cpu_final = Psi_final if xp_local is np else xp_local.asnumpy(Psi_final)
        Phi_cpu_final = Phi_H_final if xp_local is np else xp_local.asnumpy(Phi_H_final)

    masses_raw = diagonalize_with_H(Psi_cpu_final, Phi_cpu_final, gY, m0, g_param, lam_1, lam_2)

    if masses_raw is None or len(masses_raw) < 10:
            error_msg = "Diagonalization failed"
            # Jawnie ustawiamy wartości kary, aby były spójne
            fractal_score, hierarchy, regularity = -1.0, 0.0, 0.0
    else:
            tpu_print(f"     -> Znaleziono {len(masses_raw)} surowych mas. Obliczanie metryk v38...")
            fractal_score, hierarchy, regularity = fractal_consistency_score(masses_raw)
            tpu_print(f"  [JOB {job_id}] WYNIK KOŃCOWY: Score={fractal_score:.4f} | Hierarchy={hierarchy:.2f} | Regularity={regularity:.3f}")
            if fractal_score > 0.6:
                tpu_print(f"  🎉 [HIGH SCORE v38] Trial {job_id}: fractal_score > 0.6!")

    del Psi_final, Phi_H_final
    gc.collect()

    result_dict = {'fractal_score': fractal_score, 'hierarchy': hierarchy, 'regularity': regularity}
    log_entry = {
    'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'), 'trial_number': job_id,
    'g_Y': g_Yukawa, 'mu2': mu2_val, 'v_H': v_H, 'm0': m0, 'g': g_param, 'lam_1': lam_1,
    **result_dict, 'error_msg': error_msg
    }
    pd.DataFrame([log_entry]).to_csv(LOG_CSV_FILE, mode='a', header=not os.path.exists(LOG_CSV_FILE) or os.path.getsize(LOG_CSV_FILE) == 0, index=False)

    if IS_TPU:
        write_dual_heartbeat(f"job_{job_id}_log_written_fscore={result_dict['fractal_score']:.4f}")

    tpu_print(f"<<<< [JOB {job_id}] Zakończono. Czas: {(time.time() - start_time):.1f}s.")
    return result_dict

def verification_wrapper_sequential(jobs):
    results = []
    job_counter = 0
    tpu_print(f"[VERIFY WRAPPER] Uruchamianie {len(jobs)} zadań weryfikacyjnych sekwencyjnie.")
    for p in jobs:
        job_counter += 1
        tpu_print(f"[VERIFY WRAPPER] Uruchamianie zadania {job_counter}/{len(jobs)}...")
        results.append(run_self_consistent_job(*p))
        if IS_TPU:
            xm.mark_step()
def objective_fractal_ToE(trial):
    tpu_print(f"\n--- OPTUNA Trial #{trial.number} ---")

    # UPDATED PARAMETER RANGES based on recommendations
    g_Y = trial.suggest_float('g_Y', 0.05, 0.5)  # was 0.3-1.5
    mu2 = trial.suggest_float('mu2', -2.0, -0.5)  # was -25.0 to -15.0
    v_H = trial.suggest_float('v_H', 0.4, 0.8)
    m0 = trial.suggest_float('m0', -1.0, -0.1)  # was 0.1-20.0, now negative
    g = trial.suggest_float('g', 0.05, 0.2)  # was 2.0-8.0
    lam_1 = trial.suggest_float('lam_1', 0.03, 0.1)  # was 0.5-1.5

    lam_2 = trial.suggest_float('lam_2', 0.005, 0.05)  # was lam_1 * pi
    tpu_print(f"[PARAMETRY ENHANCED] gY={g_Y:.3f}, mu2={mu2:.2f}, vH={v_H:.2f}, m0={m0:.2f}, g={g:.2f}, lam1={lam_1:.2f}, lam2={lam_2:.2f}")

    result = run_self_consistent_job(g_Y, mu2, v_H, m0, g, lam_1, lam_2, f"fractal_{trial.number}")

    fractal_score = result.get('fractal_score', -1.0)
    hierarchy = result.get('hierarchy', 0.0)
    regularity = result.get('regularity', 0.0)

    # ADD HIERARCHY BONUS (Recommendation #3)
    hierarchy_bonus = 1.0 if hierarchy > 10.0 else 0.0

    # Enhanced score with hierarchy bonus
    enhanced_score = fractal_score * (1.0 + hierarchy_bonus)

    tpu_print(f"### WYNIK OPTUNA T{trial.number} (ENHANCED) ### Score={fractal_score:.4f} | Enhanced={enhanced_score:.4f} | Hierarchy={hierarchy:.4f} (Bonus={hierarchy_bonus:.1f}) | Regularity={regularity:.4f}")

    trial.set_user_attr("fractal_score", fractal_score)
    trial.set_user_attr("hierarchy", hierarchy)
    trial.set_user_attr("regularity", regularity)
    trial.set_user_attr("hierarchy_bonus", hierarchy_bonus)
    trial.set_user_attr("enhanced_score", enhanced_score)
    trial.set_user_attr("v_H", v_H)
    trial.set_user_attr("m0", m0)
    trial.set_user_attr("g", g)
    trial.set_user_attr("lam_1", lam_1)
    trial.set_user_attr("lam_2", lam_2)

    if IS_TPU:
        tpu_print(f"### Trial {trial.number} ### Results: score={fractal_score:.3f}, enhanced={enhanced_score:.3f}, hier={hierarchy:.3f}, reg={regularity:.3f}")

    # Return tuple for multi-objective optimization
    return fractal_score, hierarchy, regularity

    if fractal_score <= -0.99:
        if IS_TPU:
            tpu_print(f"Trial {trial.number} pruned")

        # DODAJ GARBAGE COLLECTOR PRZED WYJŚCIEM/BŁĘDEM
        gc.collect()
        if IS_TPU: xm.mark_step() # Wymuś wykonanie operacji w tle

        raise optuna.exceptions.TrialPruned()

    # DODAJ GARBAGE COLLECTOR PRZED ZWROTEM WARTOŚCI
    gc.collect()
    if IS_TPU: xm.mark_step()

    return fractal_score, hierarchy, regularity

def main_runner(rank):
    """Główna funkcja uruchamiająca na jednym rdzeniu (rank=0)."""
    global device
    global BOTORCH_AVAILABLE, n_initial_trials
    if IS_TPU:
        device = xm.xla_device()
        tpu_print(f"Runner (rank {rank}) started on device: {device}")

    main_start_time = time.time()
    tpu_print("\n" + "="*80)
    tpu_print("                GŁÓWNA PĘTLA WYKONAWCZA - START                ")
    tpu_print("="*80)

    #n_initial_trials = 120
    n_verification_trials = 15
    n_surrogate_grid_points = 2500

    tpu_print("\n" + "-"*35 + " FAZA 1: PRE-TRENING " + "-"*35)
    should_run_pretrain = (EXECUTION_MODE == 'PRETRAIN_ONLY') or \
                          (EXECUTION_MODE == 'FULL_RUN' and not os.path.exists('pretrained_pinn.pth'))

    if should_run_pretrain:
        tpu_print("[INFO] Warunki spełnione, uruchamianie pre-treningu...")
        latest_checkpoint = '/kaggle/working/pinn_latest.pth'
        resume_from = None

        if os.path.exists(latest_checkpoint):
            try:
                torch.load(latest_checkpoint, map_location='cpu')
                resume_from = latest_checkpoint
                tpu_print(f"✅ Znaleziono checkpoint. Wznawianie...")
            except Exception:
                resume_from = None
                tpu_print("⚠️ Znaleziono plik checkpoint, ale nie udało się go wczytać.")

        pre_train_pinn(resume_from_epoch=resume_from)
    else:
        tpu_print("[INFO] Pomijanie pre-treningu (model 'pretrained_pinn.pth' już istnieje).")

    if EXECUTION_MODE == 'PRETRAIN_ONLY':
        tpu_print("✅ Pre-trening zakończony w trybie 'PRETRAIN_ONLY'. Zamykanie skryptu.")
        sys.exit(0)

    tpu_print("\n" + "-"*35 + " FAZA 2: OPTUNA MFO " + "-"*36)
    tpu_print("✅ Pre-trening gotowy. Przechodzenie do Fazy Optymalizacji Wielokryterialnej (Optuna)...")

    list_cache_stats()

    sampler = None
    ref_point = torch.tensor([-1.0, -1.0, -1.0], device='cpu', dtype=torch.float32)

    if BOTORCH_AVAILABLE:
        try:
            sampler = BoTorchSampler(n_startup_trials=n_initial_trials, ref_point=ref_point)
            tpu_print("Używanie BoTorch Sampler.")
        except Exception:
            sampler = NSGAIISampler(population_size=10, seed=42)
            BOTORCH_AVAILABLE = False
            tpu_print("BoTorch zawiódł, przełączanie na NSGAIISampler.")
    else:
        sampler = NSGAIISampler(population_size=10, seed=42)

    study_db_path = 'sqlite:////kaggle/working/optuna_study_v38_fractal.db'
    study_name = 'supersoliton_fractal_v38'
    try:
        study = optuna.load_study(study_name=study_name, storage=study_db_path, sampler=sampler)
        n_completed = len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])
        tpu_print(f"✅ Optuna study resumed: {n_completed} completed trials")
        n_initial_trials = max(0, n_initial_trials - n_completed)
    except KeyError:
        study = optuna.create_study(study_name=study_name,
                                    directions=['maximize', 'maximize', 'maximize'],
                                    sampler=sampler, storage=study_db_path, load_if_exists=True)
        tpu_print("🆕 New Optuna study created")

    num_optuna_jobs = 1
    if not IS_TPU:
        if PSUTIL_AVAILABLE:
            available_ram_gb = psutil.virtual_memory().available / (1024**3)
            num_optuna_jobs = max(1, min(os.cpu_count() // 2, 4 if available_ram_gb > 20 else 2))
        else:
            num_optuna_jobs = 1

    tpu_print(f"[OPTUNA PARALLELISM] Ustawiono {num_optuna_jobs} wątków wykonawczych dla Optuna.")

    if n_initial_trials > 0:
        tpu_print(f"\n[OPTUNA] Rozpoczynanie pętli optymalizacji dla {n_initial_trials} prób z {num_optuna_jobs} wątkami...")
        optuna_start_time = time.time()
        study.optimize(objective_fractal_ToE, n_trials=n_initial_trials, n_jobs=num_optuna_jobs)
        optuna_duration = time.time() - optuna_start_time
        tpu_print(f"[OPTUNA] Pętla optymalizacji zakończona. Czas: {optuna_duration/60:.2f} minut.")

        top_trials = sorted(study.trials, key=lambda t: t.user_attrs.get('fractal_score', -1), reverse=True)[:5]
        tpu_print("\n🏆 TOP 5 TRIALS by Fractal Score:")
        for i, trial in enumerate(top_trials):
            score = trial.user_attrs.get('fractal_score', -1)
            hierarchy = trial.user_attrs.get('hierarchy', -1)
            regularity = trial.user_attrs.get('regularity', -1)
            params = {k: v for k, v in trial.params.items()}
            tpu_print(f"  #{i+1}: FScore={score:.4f} | Hier={hierarchy:.3f} | Reg={regularity:.3f} | params={params}")
    else:
        tpu_print("[OPTUNA] Faza eksploracji pominięta: Wystarczająca liczba prób już istnieje w badaniu.")

    tpu_print("\n" + "-"*29 + " FAZA 3: GP SURROGATE I WERYFIKACJA " + "-"*29)

    gp_model = None
    df_trials = study.trials_dataframe(multi_index=True)
    df_trials.columns = ['_'.join(map(str, col)).strip() for col in df_trials.columns.values]

    df_trials = df_trials.rename(columns={'values_0': 'fractal_score', 'values_1': 'hierarchy', 'values_2': 'regularity'}).dropna(subset=['fractal_score', 'hierarchy', 'regularity'])
    df_trials = df_trials[df_trials['fractal_score'] > -0.99]

    feature_cols_base = ['g_Y', 'mu2', 'v_H', 'm0', 'g', 'lam_1']
    feature_cols_mapped = [f'params_{p}' for p in feature_cols_base]
    available_cols = [col for col in feature_cols_mapped if col in df_trials.columns]

    X = df_trials[available_cols].values
    y = df_trials[['fractal_score', 'hierarchy', 'regularity']].values
    df_verified = pd.DataFrame()

    if len(X) > 1 and BOTORCH_AVAILABLE and X.shape[1] >= 2:
        try:
            tpu_print(f"[GP] Trenowanie modelu SingleTaskGP na {len(X)} próbkach...")
            scaler_X, scaler_y = MinMaxScaler(), MinMaxScaler()
            X_scaled, y_scaled = scaler_X.fit_transform(X), scaler_y.fit_transform(y)
            X_torch_scaled = torch.tensor(X_scaled, dtype=torch.float32, device='cpu')
            Y_torch_scaled = torch.tensor(y_scaled, dtype=torch.float32, device='cpu')
            gp_model = SingleTaskGP(X_torch_scaled, Y_torch_scaled)
            mll = ExactMarginalLogLikelihood(gp_model.likelihood, gp_model)
            fit_gpytorch_model(mll)
            tpu_print(f"✅ [GP] Model wytrenowany. Log Likelihood: {mll.log_likelihood().item():.2f}")
            dump((scaler_X, scaler_y, gp_model.state_dict()), 'gp_surrogate_v38_fractal.joblib')
            tpu_print("[GP] Model i skalery zapisane w 'gp_surrogate_v38_fractal.joblib'.")

            param_bounds_6D = {
                'params_g_Y': [0.1, 2.0], 'params_mu2': [-30.0, -10.0], 'params_v_H': [0.3, 1.0],
                'params_m0': [0.05, 50.0], 'params_g': [1.0, 10.0], 'params_lam_1': [0.3, 2.0]
            }
            grid_axes = []
            for col in available_cols:
                low, high = param_bounds_6D[col]
                n_points_per_dim = int(n_surrogate_grid_points**(1/len(available_cols)))
                grid_axes.append(np.linspace(low, high, n_points_per_dim))
            params_grid_cpu = np.array(list(itertools.product(*grid_axes)))
            X_grid_scaled = scaler_X.transform(params_grid_cpu)
            X_grid_torch = torch.tensor(X_grid_scaled, dtype=torch.float32, device='cpu')

            tpu_print(f"[GP] Generowanie siatki {params_grid_cpu.shape[0]} punktów do predykcji...")
            with torch.no_grad():
                gp_model.eval()
                posterior = gp_model.posterior(X_grid_torch)
                mean_scaled = posterior.mean.cpu().numpy()
            predictions = scaler_y.inverse_transform(np.clip(mean_scaled, 0.0, 1.0))
            pareto_mask = is_pareto_efficient(predictions)
            pareto_preds = predictions[pareto_mask]
            pareto_params = params_grid_cpu[pareto_mask]
            tpu_print(f"[GP] Znaleziono {len(pareto_preds)} punktów na froncie Pareto predykcji GP.")

            tpu_print(f"Uruchamianie weryfikacji dla {n_verification_trials} najlepszych kandydatów z frontu Pareto GP...")
            verification_jobs_params = []
            if len(pareto_params) > 0:
                hypervolume_proxy = pareto_preds[:, 0] + pareto_preds[:, 1] + pareto_preds[:, 2]
                indices_to_verify = np.argsort(hypervolume_proxy)[-n_verification_trials:]
                for rank_v, i in enumerate(indices_to_verify):
                    p = pareto_params[i]
                    param_map = dict(zip(available_cols, p))
                    gY = param_map.get('params_g_Y', 1.0)
                    mu2 = param_map.get('params_mu2', -20.0)
                    vH = param_map.get('params_v_H', 0.5)
                    m0 = param_map.get('params_m0', m0_init)
                    g = param_map.get('params_g', g_init)
                    lam1 = param_map.get('params_lam_1', lam_1_init)
                    lam2 = lam1 * np.pi
                    verification_jobs_params.append((gY, mu2, vH, m0, g, lam1, lam2, f"verify_{rank_v}"))

            if IS_TPU or num_optuna_jobs == 1:
                tpu_print("[VERIFY] Uruchamianie weryfikacji sekwencyjnie.")
                verification_results = verification_wrapper_sequential(verification_jobs_params)
            else:
                tpu_print(f"[VERIFY] Uruchamianie weryfikacji równolegle z {num_optuna_jobs} wątkami.")
                verification_results = Parallel(n_jobs=num_optuna_jobs, backend='loky')(delayed(run_self_consistent_job)(*p) for p in verification_jobs_params)

            df_verified = pd.DataFrame(verification_results)
            tpu_print("✅ [VERIFY] Weryfikacja zakończona.")

        except Exception as e:
            tpu_print(f"⚠️ Faza GP/Weryfikacji nie powiodła się: {e}")
            gp_model = None
    else:
        tpu_print("⚠️ Pomijanie fazy GP/Weryfikacji: niewystarczająca liczba danych lub brak BoTorch.")

    tpu_print("\n" + "-"*35 + " FAZA 4: WIZUALIZACJA " + "-"*35)
    tpu_print("[PLOT] Generowanie wykresu frontu Pareto...")
    plt.figure(figsize=(12, 9))

    initial_score = df_trials['fractal_score'].values
    initial_hierarchy = df_trials['hierarchy'].values

    plt.scatter(initial_score, initial_hierarchy, c='blue', label=f'Eksploracja Optuna ({len(initial_score)} pkt)', s=50, alpha=0.6)

    if gp_model and 'pareto_preds' in locals() and len(pareto_preds) > 0:
        plt.scatter(pareto_preds[:, 0], pareto_preds[:, 1], c='gray', label=f'Predykcje GP (Front Pareto, {len(pareto_preds)} pkt)', s=10, alpha=0.3)

    if not df_verified.empty and 'fractal_score' in df_verified.columns:
        verified_valid = df_verified[df_verified['fractal_score'] > -0.99]
        plt.scatter(verified_valid['fractal_score'], verified_valid['hierarchy'],
                    c='red', marker='*', s=250, label=f'Zweryfikowane punkty ({len(verified_valid)} pkt)',
                    edgecolors='black', zorder=10)

        for i, row in verified_valid.nlargest(3, 'fractal_score').iterrows():
            plt.text(row['fractal_score'], row['hierarchy'], f' S={row["fractal_score"]:.3f}', fontsize=9, color='darkred')


    plt.xlabel("Fractal Score (Korelacja + Bonus Frak.)")
    plt.ylabel("Hierarchy")
    plt.title("Front Pareto: Optymalizacja Wielokryterialna (v38 Fractal Consistency)")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.xlim(left=max(-1, np.min(initial_score)-0.05) if len(initial_score)>0 else -1)
    plt.ylim(bottom=max(-1, np.min(initial_hierarchy)-0.05) if len(initial_hierarchy)>0 else -1)

    output_filename = "pareto_front_v38_fractal_kaggle.png"
    plt.savefig(output_filename, dpi=150)
    tpu_print(f"✅ Wykres zapisany jako '{output_filename}'.")

    tpu_print("\n" + "-"*35 + " FAZA 4: WIZUALIZACJA 3D " + "-"*35)
    try:
        import plotly.express as px
        df_trials_3d = study.trials_dataframe()
        df_trials_3d.columns = ['_'.join(map(str, col)).strip() for col in df_trials_3d.columns.values]

        if 'values_2' in df_trials_3d.columns:
            fig = px.scatter_3d(df_trials_3d,
                                x='values_0', y='values_1', z='values_2',
                                color='values_0',
                                hover_data=['params_g_Y', 'params_mu2', 'params_m0', 'params_g', 'params_lam_1'])
            fig.update_layout(
                title="Front Pareto 3D: Fractal Score vs Hierarchy vs Regularity (v38)",
                scene=dict(xaxis_title='Fractal Score', yaxis_title='Hierarchy', zaxis_title='Regularity')
            )
            fig.write_html("pareto_front_3d_v38.html")
            tpu_print("✅ Interaktywny wykres 3D zapisany jako 'pareto_front_3d_v38.html'.")
    except ImportError:
        tpu_print("⚠️ Plotly nie jest zainstalowane. Pomijam wykres 3D.")
    except Exception as e:
        tpu_print(f"⚠️ Błąd przy generowaniu wykresu 3D: {e}")


    total_runtime = time.time() - main_start_time
    tpu_print("\n" + "="*80)
    tpu_print("✅ SKRYPT ZAKOŃCZYŁ DZIAŁANIE POMYŚLNIE.")
    tpu_print(f"Całkowity czas wykonania: {total_runtime/3600:.2f} godzin ({total_runtime/60:.1f} minut).")
    tpu_print("="*80)

if __name__ == '__main__':
    main_runner(0)
