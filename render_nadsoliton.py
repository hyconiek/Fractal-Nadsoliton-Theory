import numpy as np
import imageio
import math
from tqdm import tqdm  # Do paska postępu (pip install tqdm)

# === PARAMETERS FROM PDF ===
ALPHA = 2.772589
BETA = 0.010000
OMEGA = 0.785398
PHI = 0.523599
OCTAVES = 12

# === RENDER SETTINGS ===
WIDTH = 400       # Szerokość GIF-a
HEIGHT = 400      # Wysokość GIF-a
FPS = 30          # Klatki na sekundę
DURATION = 10     # Czas trwania w sekundach
SCALES_PER_SEC = 0.000000000000001 # Szybkość zmiany skali (rzędy wielkości na sekundę)
FILENAME = "nadsoliton_descent_smooth.gif"

def get_amplitude(d):
    fd = float(d)
    num = ALPHA * np.cos(OMEGA * fd + PHI)
    den = 1.0 + BETA * fd
    return num / den

# Precompute amplitudes for speed
AMPLITUDES = [get_amplitude(i) for i in range(1, OCTAVES + 1)]

def hologram(uv, time):
    # uv is a meshgrid (Y, X)
    field_sum = np.zeros_like(uv[0])

    for i in range(OCTAVES):
        k = float(i + 1)
        amp = AMPLITUDES[i]

        # Slow wave motion (0.05 and 0.035 speed coefficients)
        modeX = np.cos(uv[1] * k + time * 0.05)
        modeY = np.cos(uv[0] * k - time * 0.035)

        field_sum += amp * (modeX * modeY)

    return field_sum

def render_frame(t, width, height):
    # 1. Setup coordinates (-0.5 to 0.5, aspect corrected)
    y = np.linspace(-0.5, 0.5, height)
    x = np.linspace(-0.5, 0.5, width)
    if width > height:
        x *= width / height
    else:
        y *= height / width

    Y, X = np.meshgrid(y, x, indexing='ij')

    # 2. Zoom Mechanism (Smooth Descent)
    # We want to descend at a rate of SCALES_PER_SEC orders of magnitude per second.
    # Start scale is arbitrary (e.g. 10^5), we go down.
    # scale = 10^(Start - Rate * t)

    start_log_scale = 15.0 # Start high enough to look "macro"
    current_log_scale = start_log_scale - (SCALES_PER_SEC * t)

    # Loop the zoom for seamless GIF? Or just one pass?
    # User asked for "fluidity" and "10 seconds length".
    # Let's make it continuous descent.

    scale = np.power(10.0, current_log_scale)

    # Apply scale
    Y *= scale
    X *= scale

    # 3. Compute Field
    field = hologram([Y, X], t)
    density = field * field

    # 4. Coloring (Vectorized)
    val = density * 0.08

    # Initialize color channels
    r = np.zeros_like(val)
    g = np.zeros_like(val)
    b = np.zeros_like(val)

    # Deep Space (Black/Blue)
    b += 0.15 * val

    # Structure (Cyan/Purple mixing logic)
    mask1 = np.clip((val - 0.2) / (0.25 - 0.2), 0.0, 1.0)
    r = r * (1 - mask1) + 0.0 * mask1
    g = g * (1 - mask1) + 0.6 * mask1
    b = b * (1 - mask1) + 1.0 * mask1

    # Energy Cores (Orange)
    mask2 = np.clip((val - 0.6) / (0.8 - 0.6), 0.0, 1.0)
    r = r * (1 - mask2) + 1.0 * mask2
    g = g * (1 - mask2) + 0.8 * mask2
    b = b * (1 - mask2) + 0.2 * mask2

    # Hotspots (White)
    mask3 = np.clip((val - 1.2) / (1.5 - 1.2), 0.0, 1.0)
    r += mask3
    g += mask3
    b += mask3

    # Clamp and convert to uint8
    img = np.stack([r, g, b], axis=-1)
    img = np.clip(img, 0.0, 1.0)
    return (img * 255).astype(np.uint8)

def create_gif():
    print(f"Generating GIF: {WIDTH}x{HEIGHT}, {FPS} fps, {DURATION}s")
    frames = []
    total_frames = FPS * DURATION

    for i in tqdm(range(total_frames)):
        t = i / FPS
        frame = render_frame(t, WIDTH, HEIGHT)
        frames.append(frame)

    imageio.mimsave(FILENAME, frames, fps=FPS, loop=0)
    print(f"Done! Saved to {FILENAME}")

if __name__ == "__main__":
    create_gif()
