"""
Superscreen simulation pipeline for superconducting diode effect
Analyzes Si/MgO/CoFeB/Ta/Nb stack with varying CoFeB thickness
Stack orientation: Si substrate (bottom) → MgO → CoFeB → Ta → Nb (top)
Magnetic field from CoFeB propagates UPWARD through Ta into Nb
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Material parameters
# Niobium superconductor parameters
lambda_L = 39e-9  # London penetration depth (m)
xi = 38e-9  # Coherence length (m)
Tc = 9.2  # Critical temperature (K)
T_op = 4.2  # Operating temperature (K)

# Calculate temperature-dependent penetration depth
lambda_eff = lambda_L / np.sqrt(1 - (T_op/Tc)**4)

# Device geometry (must match MuMax3)
device_width = 2.048e-6  # m (2.048 μm)
device_length = 4.096e-6  # m (4.096 μm)
Nb_thickness = 20e-9  # m (20 nm)
Ta_thickness = 2e-9  # m (2 nm)

# CoFeB thicknesses to analyze (matching MuMax3)
cofeb_thicknesses = np.array([0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]) * 1e-9

# Storage for results
results = {
    'thickness_nm': [],
    'B_exchange_mT': [],
    'I_c_plus': [],
    'I_c_minus': [],
    'diode_efficiency': [],
    'screening_current_max': []
}

print("="*60)
print("Superconducting Diode Effect Analysis")
print("Si/MgO/CoFeB/Ta/Nb Stack (Nb on top)")
print("="*60)

# Physical constants
mu_0 = 4 * np.pi * 1e-7  # Permeability of free space
Ms = 1.4e6  # A/m - Saturation magnetization of CoFeB

for t_cofeb in cofeb_thicknesses:
    t_nm = t_cofeb * 1e9
    print(f"\nProcessing CoFeB thickness: {t_nm:.1f} nm")
    
    # Load magnetic field from MuMax3 output
    # Stack: CoFeB is BELOW Ta and Nb - field propagates UPWARD
    # In practice: 
    # from superscreen.io import read_ovf
    # B_field = read_ovf(f'Bfield_CoFeB_{t_nm:.1f}nm_at_interface.ovf')
    # For this example, we'll create synthetic field data based on expected physics
    
    # Exchange field propagates upward from CoFeB through Ta spacer into Nb
    # Field decays with distance from CoFeB surface
    # Simplified model: B_z(stray) ≈ μ₀ * Ms * (t_cofeb / (2 * distance_above))
    # Distance is measured UPWARD from CoFeB surface
    
    distance_to_Nb = Ta_thickness + Nb_thickness/2  # Vertical distance upward
    
    # Stray field magnitude pointing upward (+z direction)
    # Simplified dipole approximation for thin film
    B_exchange = mu_0 * Ms * (t_cofeb / (2 * distance_to_Nb))
    
    # Create synthetic spatial field profile
    # Real implementation would load from MuMax3 OVF file
    # Field points upward (+z) from CoFeB magnetization
    nx, ny = 128, 256
    x = np.linspace(-device_width/2, device_width/2, nx)
    y = np.linspace(-device_length/2, device_length/2, ny)
    X, Y = np.meshgrid(x, y)
    
    # Non-uniform field with edge effects (stronger at center, weaker at edges)
    B_z = B_exchange * (1 + 0.2 * np.exp(-((X/(device_width*0.3))**2 + 
                                           (Y/(device_length*0.3))**2)))
    
    print(f"  Exchange field magnitude: {B_exchange*1e3:.3f} mT")
    
    # Superconducting diode effect calculation
    # The diode effect arises from Rashba spin-orbit coupling + exchange field
    # Creating asymmetric potential landscape
    # Simplified model: ΔI_c/I_c ≈ α_R * B_ex / (ℏ * v_F / (2e))
    
    # Estimate critical current density (simplified Ginzburg-Landau)
    # For thin Nb films at T=4.2K
    j_c0 = 1e10  # A/m²
    
    # Asymmetry parameter (depends on Rashba coupling strength)
    # α_R typically 1e-10 eV·m for heavy metal interfaces
    alpha_R = 1e-10  # eV·m (typical Rashba parameter)
    hbar = 1.055e-34  # J·s
    v_F = 1e6  # m/s (Fermi velocity, typical for metals)
    e_charge = 1.602e-19  # C
    
    # Calculate asymmetry factor
    energy_scale = hbar * v_F / (2 * e_charge)  # Convert to energy
    asymmetry_factor = (alpha_R * e_charge * B_exchange) / energy_scale
    
    # Calculate critical currents
    I_c_base = j_c0 * device_width * Nb_thickness
    I_c_plus = I_c_base * (1 + asymmetry_factor)
    I_c_minus = I_c_base * (1 - asymmetry_factor)
    
    # Diode efficiency
    if (I_c_plus + I_c_minus) > 0:
        diode_eff = (I_c_plus - I_c_minus) / (I_c_plus + I_c_minus)
    else:
        diode_eff = 0
    
    # Maximum screening current (proportional to applied field)
    J_screen_max = B_exchange / (mu_0 * lambda_eff)
    I_screen_max = J_screen_max * device_width * Nb_thickness
    
    # Store results
    results['thickness_nm'].append(t_nm)
    results['B_exchange_mT'].append(B_exchange * 1e3)
    results['I_c_plus'].append(I_c_plus * 1e6)  # Convert to μA
    results['I_c_minus'].append(I_c_minus * 1e6)
    results['diode_efficiency'].append(diode_eff * 100)  # Convert to %
    results['screening_current_max'].append(I_screen_max * 1e6)  # μA
    
    print(f"  I_c+ = {I_c_plus*1e6:.2f} μA")
    print(f"  I_c- = {I_c_minus*1e6:.2f} μA")
    print(f"  Diode efficiency = {diode_eff*100:.3f}%")

# Convert to arrays for plotting
thickness_array = np.array(results['thickness_nm'])
B_exchange_array = np.array(results['B_exchange_mT'])
I_c_plus_array = np.array(results['I_c_plus'])
I_c_minus_array = np.array(results['I_c_minus'])
diode_eff_array = np.array(results['diode_efficiency'])
screening_array = np.array(results['screening_current_max'])

print("\n" + "="*60)
print("Results Summary:")
print("="*60)
print(f"{'Thickness (nm)':<15} {'B_ex (mT)':<12} {'I_c+ (μA)':<12} {'I_c- (μA)':<12} {'η (%)':<10}")
print("-"*60)
for i in range(len(thickness_array)):
    print(f"{thickness_array[i]:<15.1f} {B_exchange_array[i]:<12.3f} "
          f"{I_c_plus_array[i]:<12.2f} {I_c_minus_array[i]:<12.2f} "
          f"{diode_eff_array[i]:<10.3f}")
print("="*60)

# Create comprehensive plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Superconducting Diode Effect vs CoFeB Thickness\nSi/MgO/CoFeB/Ta(2nm)/Nb(20nm) - Nb on Top', 
             fontsize=14, fontweight='bold')

# Plot 1: Exchange field vs thickness
ax1 = axes[0, 0]
ax1.plot(thickness_array, B_exchange_array, 'o-', color='#2E86AB', 
         linewidth=2, markersize=8, markerfacecolor='white', markeredgewidth=2)
ax1.set_xlabel('CoFeB Thickness (nm)', fontsize=11)
ax1.set_ylabel('Exchange Field (mT)', fontsize=11)
ax1.set_title('Exchange Field Strength', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0.5, 2.1])

# Plot 2: Diode efficiency vs thickness
ax2 = axes[0, 1]
ax2.plot(thickness_array, diode_eff_array, 'o-', color='#A23B72',
         linewidth=2, markersize=8, markerfacecolor='white', markeredgewidth=2)
ax2.set_xlabel('CoFeB Thickness (nm)', fontsize=11)
ax2.set_ylabel('Diode Efficiency (%)', fontsize=11)
ax2.set_title('Superconducting Diode Efficiency', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0.5, 2.1])

# Plot 3: Critical currents vs thickness
ax3 = axes[1, 0]
ax3.plot(thickness_array, I_c_plus_array, 'o-', color='#F18F01', 
         linewidth=2, markersize=8, label='$I_c^+$', markerfacecolor='white', markeredgewidth=2)
ax3.plot(thickness_array, I_c_minus_array, 's-', color='#C73E1D',
         linewidth=2, markersize=8, label='$I_c^-$', markerfacecolor='white', markeredgewidth=2)
ax3.set_xlabel('CoFeB Thickness (nm)', fontsize=11)
ax3.set_ylabel('Critical Current (μA)', fontsize=11)
ax3.set_title('Asymmetric Critical Currents', fontsize=12, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.set_xlim([0.5, 2.1])

# Plot 4: Screening current vs thickness
ax4 = axes[1, 1]
ax4.plot(thickness_array, screening_array, 'o-', color='#6A4C93',
         linewidth=2, markersize=8, markerfacecolor='white', markeredgewidth=2)
ax4.set_xlabel('CoFeB Thickness (nm)', fontsize=11)
ax4.set_ylabel('Max Screening Current (μA)', fontsize=11)
ax4.set_title('Meissner Screening Current', fontsize=12, fontweight='bold')
ax4.grid(True, alpha=0.3)
ax4.set_xlim([0.5, 2.1])

plt.tight_layout()
plt.savefig('superconducting_diode_analysis.png', dpi=300, bbox_inches='tight')
print("\nPlot saved as 'superconducting_diode_analysis.png'")
plt.show()

# Additional analysis plot: Combined view
fig2, ax = plt.subplots(figsize=(10, 6))
ax2 = ax.twinx()

line1 = ax.plot(thickness_array, B_exchange_array, 'o-', color='#2E86AB',
                linewidth=2.5, markersize=10, label='Exchange Field', 
                markerfacecolor='white', markeredgewidth=2)
line2 = ax2.plot(thickness_array, diode_eff_array, 's-', color='#A23B72',
                 linewidth=2.5, markersize=10, label='Diode Efficiency',
                 markerfacecolor='white', markeredgewidth=2)

ax.set_xlabel('CoFeB Thickness (nm)', fontsize=12, fontweight='bold')
ax.set_ylabel('Exchange Field (mT)', fontsize=12, fontweight='bold', color='#2E86AB')
ax2.set_ylabel('Diode Efficiency (%)', fontsize=12, fontweight='bold', color='#A23B72')

ax.tick_params(axis='y', labelcolor='#2E86AB')
ax2.tick_params(axis='y', labelcolor='#A23B72')

ax.grid(True, alpha=0.3)
ax.set_xlim([0.5, 2.1])

lines = line1 + line2
labels = [l.get_label() for l in lines]
ax.legend(lines, labels, loc='upper left', fontsize=11)

plt.title('Exchange Field and Diode Efficiency vs CoFeB Thickness\nSi/MgO/CoFeB/Ta(2nm)/Nb(20nm) - Nb on Top',
          fontsize=13, fontweight='bold', pad=20)
plt.tight_layout()
plt.savefig('diode_combined_analysis.png', dpi=300, bbox_inches='tight')
print("Combined plot saved as 'diode_combined_analysis.png'")
plt.show()

# Save results to CSV
with open('diode_results.csv', 'w') as f:
    f.write('thickness_nm,B_exchange_mT,I_c_plus_uA,I_c_minus_uA,diode_efficiency_percent,screening_current_max_uA\n')
    for i in range(len(thickness_array)):
        f.write(f'{thickness_array[i]:.1f},{B_exchange_array[i]:.6f},'
                f'{I_c_plus_array[i]:.6f},{I_c_minus_array[i]:.6f},'
                f'{diode_eff_array[i]:.6f},{screening_array[i]:.6f}\n')

print("Results saved to 'diode_results.csv'")

print("\n" + "="*60)
print("Analysis Complete!")
print("="*60)