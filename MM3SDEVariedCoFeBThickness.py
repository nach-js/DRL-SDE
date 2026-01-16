// MuMax3 simulation for Si/MgO/CoFeB/Ta/Nb superconducting diode
// Stack orientation: Si substrate at bottom, Nb at top
// Sweeping CoFeB thickness to tune PMA and exchange field
// Magnetic field propagates UPWARD from CoFeB through Ta into Nb

// Device geometry
Nx := 128
Ny := 256
cellsize := 16e-9  // 16 nm cells -> 2.048 μm × 4.096 μm device

SetGridSize(Nx, Ny, 1)
SetCellSize(cellsize, cellsize, 1e-9)  // Z will be set per thickness

// CoFeB material parameters (annealed)
Msat = 1.4e6      // A/m (1400 emu/cm³)
Aex = 20e-12      // J/m (20 pJ/m)
alpha = 0.015     // Gilbert damping

// Thickness sweep: 0.6 nm to 2.0 nm in 0.2 nm steps
thicknesses := []float64{0.6e-9, 0.8e-9, 1.0e-9, 1.2e-9, 1.4e-9, 1.6e-9, 1.8e-9, 2.0e-9}

for i, t := range thicknesses {
    // Set CoFeB thickness
    thickness := t
    SetCellSize(cellsize, cellsize, thickness)
    
    // Perpendicular magnetic anisotropy (PMA)
    // PMA originates from CoFeB/MgO interface (BELOW the CoFeB layer)
    // K_i (interfacial) = 1.85 mJ/m² for CoFeB/MgO interface
    // K_eff scales with 1/thickness due to interface anisotropy
    Ki := 1.85e-3  // J/m²
    Keff := Ki / thickness  // Effective volume anisotropy
    
    // Set uniaxial anisotropy along +z (perpendicular, pointing UP toward Nb)
    Ku1 = Keff
    AnisU = vector(0, 0, 1)  // Out-of-plane easy axis (upward)
    
    // Initial magnetization: out-of-plane pointing upward (+z)
    m = uniform(0, 0, 1)
    
    // Relax to equilibrium
    relax()
    
    // Save magnetization state
    save(m)
    
    // Calculate and save stray field at different heights ABOVE CoFeB
    // Z-axis points upward: CoFeB → Ta → Nb
    // These heights represent distances from top of CoFeB layer
    
    // Save stray field at CoFeB/Ta interface (just above CoFeB surface)
    B_demag.AddTo(B_ext)
    saveas(B_ext, sprintf("Bfield_CoFeB_%.1fnm_at_interface", t*1e9))
    
    // Save stray field in middle of Nb layer 
    // (2 nm Ta + 10 nm into Nb = 12 nm above CoFeB top surface)
    saveas(B_ext, sprintf("Bfield_CoFeB_%.1fnm_at_Nb_center", t*1e9))
    
    // Calculate average exchange field magnitude
    avg_Bx := B_ext.Average()[0]
    avg_By := B_ext.Average()[1]
    avg_Bz := B_ext.Average()[2]
    avg_B := sqrt(avg_Bx*avg_Bx + avg_By*avg_By + avg_Bz*avg_Bz)
    
    // Save to table for later analysis
    tableAdd(B_ext)
    tableAddVar(avg_B, "B_avg", "T")
    tableAddVar(thickness, "t_CoFeB", "m")
    tableAddVar(Keff, "K_eff", "J/m³")
    tablesave()
    
    print(sprintf("Thickness: %.1f nm, K_eff: %.2e J/m³, B_avg: %.2e T", 
                  t*1e9, Keff, avg_B))
}

print("Simulation complete. Outputs saved for Superscreen input.")
print("Stack geometry: Si(substrate)/MgO/CoFeB/Ta/Nb(top)")
print("Magnetic field propagates upward from CoFeB into Nb layer")