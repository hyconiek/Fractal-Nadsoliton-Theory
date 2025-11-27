import numpy as np
import imageio
import math
from tqdm import tqdm

# === PARAMETERS FROM PDF ===
ALPHA = 2.772589
BETA = 0.010000
OMEGA = 0.785398
PHI = 0.523599
OCTAVES = 12

# === RENDER SETTINGS ===
WIDTH = 400
HEIGHT = 400
FPS = 20
DURATION = 5
FILENAME = "nadsoliton_scale10_fast.gif"
SCALE = 10.0
SPEED_MULTIPLIER = 5.0  # Speed up factor

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

    # Effective time scaled by speed multiplier
    t_eff = time * SPEED_MULTIPLIER

    for i in range(OCTAVES):
        k = float(i + 1)
        amp = AMPLITUDES[i]

        # Adjusted base speeds for more dynamic movement
        # Original: 0.05 and 0.08
        # We keep the base speeds low so the multiplier handles the acceleration
        modeX = np.cos(uv[1] * k + t_eff * 0.5)
        modeY = np.cos(uv[0] * k - t_eff * 0.8)

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

    # 2. Apply Fixed Scale
    Y *= SCALE
    X *= SCALE

    # 3. Compute Field
    field = hologram([Y, X], t)
    density = field * field

    # 4. Coloring (Vectorized) - Optimized for Intermediate Scale
    val = density * 0.08

    # Initialize color channels
    r = np.zeros_like(val)
    g = np.zeros_like(val)
    b = np.zeros_like(val)

    # Dark matter web (Blue/Purple)
    mask0 = np.clip((val - 0.0) / (0.3 - 0.0), 0.0, 1.0)
    r += 0.05 * mask0
    g += 0.0 * mask0
    b += 0.2 * mask0

    # Structure bridges (Magenta/Red)
    mask1 = np.clip((val - 0.2) / (0.5 - 0.2), 0.0, 1.0)
    r = r * (1 - mask1) + 0.6 * mask1
    g = g * (1 - mask1) + 0.1 * mask1
    b = b * (1 - mask1) + 0.4 * mask1

    # Energy nodes (Orange/Gold)
    mask2 = np.clip((val - 0.5) / (0.9 - 0.5), 0.0, 1.0)
    r = r * (1 - mask2) + 1.0 * mask2
    g = g * (1 - mask2) + 0.7 * mask2
    b = b * (1 - mask2) + 0.2 * mask2

    # Singularities (White)
    mask3 = np.clip((val - 0.9) / (1.5 - 0.9), 0.0, 1.0)
    r = r * (1 - mask3) + 1.0 * mask3
    g = g * (1 - mask3) + 1.0 * mask3
    b = b * (1 - mask3) + 1.0 * mask3

    # Clamp and convert to uint8
    img = np.stack([r, g, b], axis=-1)
    img = np.clip(img, 0.0, 1.0)
    return (img * 255).astype(np.uint8)

def create_gif():
    print(f"Generating GIF: {WIDTH}x{HEIGHT}, {FPS} fps, {DURATION}s, Scale {SCALE}, Speed {SPEED_MULTIPLIER}x")
    frames = []
    total_frames = FPS * DURATION

    time_points = np.linspace(0, DURATION, total_frames, endpoint=False)

    for t in tqdm(time_points):
        frame = render_frame(t, WIDTH, HEIGHT)
        frames.append(frame)

    imageio.mimsave(FILENAME, frames, fps=FPS, loop=0)
    print(f"Done! Saved to {FILENAME}")

if __name__ == "__main__":
    create_gif()
