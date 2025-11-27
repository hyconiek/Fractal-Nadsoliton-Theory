import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from scipy.fft import ifftn

# Re-implement QW-355 logic for high-res visualization
def coupling_kernel(d, alpha_geo=0.202, beta_tors=0.137, omega=2.0):
    return alpha_geo * np.cos(omega * d) / (1.0 + beta_tors * d)

def create_coupling_matrix(n_octaves):
    S = np.zeros((n_octaves, n_octaves))
    for i in range(n_octaves):
        for j in range(n_octaves):
            if i != j:
                S[i, j] = coupling_kernel(np.abs(i - j))
    return S

def generate_high_res_hologram(n_octaves=12, grid_size=256): # High resolution
    # 1. Get 1D state
    S = create_coupling_matrix(n_octaves)
    eigenvalues, eigenvectors = eigh(S)
    idx_max = np.argmax(np.abs(eigenvalues))
    state_1d = eigenvectors[:, idx_max]

    # 2. Holographic Lift (1D -> 3D) via FFT
    # Pad to grid size
    n_fft = grid_size
    # Map 1D state to 3D Fourier coefficients (radial symmetry assumption for projection)
    fourier_3d = np.zeros((grid_size, grid_size, grid_size), dtype=complex)
    
    # Simple radial mapping logic from original QW-355
    center = grid_size // 2
    for x in range(grid_size):
        for y in range(grid_size):
            for z in range(grid_size):
                # Distance from DC component (0,0,0) in frequency space
                # We use a simple octant mapping to fill k-space
                k_r = np.sqrt(x**2 + y**2 + z**2)
                if k_r < n_octaves and k_r >= 0:
                     # Map octave index to radial distance
                     idx = int(k_r)
                     if idx < len(state_1d):
                        fourier_3d[x, y, z] = state_1d[idx]

    # Enforce Hermitian symmetry for real output (simplified approach)
    # Actually, let's just take the magnitude of the complex IFFT result 
    # as we are looking for density structure |psi|^2
    
    space_3d = np.fft.ifftn(fourier_3d)
    density_3d = np.abs(space_3d)**2
    
    return density_3d

# Generate high-res data
grid_size = 256 # Higher resolution than original 16
density_3d = generate_high_res_hologram(grid_size=grid_size)

# Take a central slice for 2D visualization
slice_idx = grid_size // 2
heatmap_data = density_3d[:, :, slice_idx]

# Plotting
plt.figure(figsize=(12, 10))
plt.imshow(heatmap_data, cmap='inferno', interpolation='lanczos') # Lanczos for smooth interpolation
plt.colorbar(label='Information Density |Ψ|²')
plt.title(f"QW-355: Holographic Lift - High Resolution Slice (Grid {grid_size}x{grid_size})", fontsize=16)
plt.xlabel("X (spatial)", fontsize=12)
plt.ylabel("Y (spatial)", fontsize=12)

# Enhance visibility of "blobs" (potential nadsoliton structure)
plt.contour(heatmap_data, levels=10, colors='white', alpha=0.2, linewidths=0.5)

plt.tight_layout()
plt.savefig('histogram_output.png', dpi=150, bbox_inches='tight')
print("Image saved to: histogram_output.png")
plt.show()
