"""
Superscreen simulation pipeline for superconducting diode effect
Analyzes Si/MgO/CoFeB/Ta/Nb stack with varying CoFeB thickness
Stack orientation: Si substrate (bottom) → MgO → CoFeB → Ta → Nb (top)
Magnetic field from CoFeB propagates UPWARD through Ta into Nb
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import superscreen as sc
from superscreen.io import read_ovf
import pandas as pd

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

for t_cofeb in cofeb_thicknesses:
    t_nm = t_cofeb * 1e9
    print(f"\nProcessing CoFeB thickness: {t_nm:.1f} nm")
    
    # Load magnetic field from MuMax3 output
    # Stack: CoFeB is BELOW Ta and Nb - field propagates UPWARD
    # In practice: B_field = read_ovf(f'Bfield_CoFeB_{t_nm:.1f}nm_at_Nb_center.ovf')
    # For this example, we'll create synthetic field data based on expected physics
    
    # Exchange field propagates upward from CoFeB through Ta spacer into Nb
    # Field decays with distance from CoFeB surface
    # For CoFeB with Ms = 1.4e6 A/m and thickness t_cofeb
    # Stray field at Nb layer (Ta + Nb/2 distance above CoFeB)
    
    # Simplified model: B_z(stray) ≈ μ₀ * Ms * (t_cofeb / (2 * distance_above))
    # Distance is measured UPWARD from CoFeB surface
    mu_0 = 4 * np.pi * 1e-7
    Ms = 1.4e6  # A/m
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
    B_x = np.zeros_like(B_z)  # Primarily out-of-plane field
    B_y = np.zeros_like(B_z)
    
    print(f"  Exchange field magnitude: {B_exchange*1e3:.3f} mT")
    
    # Define superconducting film geometry
    film = sc.Layer(
        "Nb",
        london_lambda=lambda_eff,
        thickness=Nb_thickness,
        coherence_length=xi
    )
    
    # Create device with rectangular geometry
    points = np.array([
        [-device_width/2, -device_length/2],
        [device_width/2, -device_length/2],
        [device_width/2, device_length/2],
        [-device_width/2, device_length/2]
    ])
    
    device = sc.Device(
        "diode",
        layers=[film],
        films=[sc.Polygon("film", points=points, layer="Nb")]
    )
    
    # Create mesh
    device = device.make_mesh(max_edge_length=device_width/50)
    
    # Apply external magnetic field (from CoFeB)
    # In real code, this would be spatially varying from MuMax3
    applied_field = sc.sources.ConstantField(B_exchange)
    
    # Calculate critical currents in both directions
    # Positive current direction
    current_plus = 1e-3  # Start with 1 mA
    transport_plus = sc.sources.LinearTransport(current_plus)
    
    # Solve for screening currents
    solution_plus = sc.solve(
        device=device,
        applied_field=applied_field,
        circulating_currents={},
        transport_currents={"film": transport_plus}
    )
    
    # For diode effect, we need to find critical current asymmetry
    # This requires iterative solving - simplified here
    
    # The diode effect arises from Rashba spin-orbit coupling + exchange field
    # Creating asymmetric potential landscape
    # Simplified model: ΔI_c/I_c ≈ α_R * B_ex / (ℏ * v_F / (2e))
    
    # Estimate critical current density (simplified Ginzburg-Landau)
    j_c0 = 1e10  # A/m² (typical for thin Nb films at T=4.2K)
    
    # Asymmetry parameter (depends on Rashba coupling strength)
    alpha_R = 1e-10  # eV·m (typical Rashba parameter)
    asymmetry_factor = (alpha_R * B_exchange / (1.05e-34 * 1e6))  # Simplified
    
    I_c_base = j_c0 * device_width * Nb_thickness
    I_c_plus = I_c_base * (1 + asymmetry_factor)
    I_c_minus = I_c_base * (1 - asymmetry_factor)
    
    # Diode efficiency
    diode_eff = (I_c_plus - I_c_minus) / (I_c_plus + I_c_minus)
    
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
    print(f"  Diode efficiency = {diode_eff*100:.2f}%")

# Convert to DataFrame for easy handling
df = pd.DataFrame(results)
print("\n" + "="*60)
print("Results Summary:")
print(df.to_string(index=False))
print("="*60)

# Create comprehensive plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Superconducting Diode Effect vs CoFeB Thickness\nSi/MgO/CoFeB/Ta(2nm)/Nb(20nm) - Nb on Top', 
             fontsize=14, fontweight='bold')

# Plot 1: Exchange field vs thickness
ax1 = axes[0, 0]
ax1.plot(df['thickness_nm'], df['B_exchange_mT'], 'o-', color='#2E86AB', 
         linewidth=2, markersize=8, markerfacecolor='white', markeredgewidth=2)
ax1.set_xlabel('CoFeB Thickness (nm)', fontsize=11)
ax1.set_ylabel('Exchange Field (mT)', fontsize=11)
ax1.set_title('Exchange Field Strength', fontsize=12, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0.5, 2.1])

# Plot 2: Diode efficiency vs thickness
ax2 = axes[0, 1]
ax2.plot(df['thickness_nm'], df['diode_efficiency'], 'o-', color='#A23B72',
         linewidth=2, markersize=8, markerfacecolor='white', markeredgewidth=2)
ax2.set_xlabel('CoFeB Thickness (nm)', fontsize=11)
ax2.set_ylabel('Diode Efficiency (%)', fontsize=11)
ax2.set_title('Superconducting Diode Efficiency', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0.5, 2.1])

# Plot 3: Critical currents vs thickness
ax3 = axes[1, 0]
ax3.plot(df['thickness_nm'], df['I_c_plus'], 'o-', color='#F18F01', 
         linewidth=2, markersize=8, label='$I_c^+$', markerfacecolor='white', markeredgewidth=2)
ax3.plot(df['thickness_nm'], df['I_c_minus'], 's-', color='#C73E1D',
         linewidth=2, markersize=8, label='$I_c^-$', markerfacecolor='white', markeredgewidth=2)
ax3.set_xlabel('CoFeB Thickness (nm)', fontsize=11)
ax3.set_ylabel('Critical Current (μA)', fontsize=11)
ax3.set_title('Asymmetric Critical Currents', fontsize=12, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.set_xlim([0.5, 2.1])

# Plot 4: Screening current vs thickness
ax4 = axes[1, 1]
ax4.plot(df['thickness_nm'], df['screening_current_max'], 'o-', color='#6A4C93',
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

# Save results to CSV
df.to_csv('diode_results.csv', index=False)
print("Results saved to 'diode_results.csv'")

# Additional analysis plot: Combined view
fig2, ax = plt.subplots(figsize=(10, 6))
ax2 = ax.twinx()

line1 = ax.plot(df['thickness_nm'], df['B_exchange_mT'], 'o-', color='#2E86AB',
                linewidth=2.5, markersize=10, label='Exchange Field', 
                markerfacecolor='white', markeredgewidth=2)
line2 = ax2.plot(df['thickness_nm'], df['diode_efficiency'], 's-', color='#A23B72',
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

print("\n" + "="*60)
print("Analysis Complete!")
print("="*60)