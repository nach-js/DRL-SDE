// MuMax3 simulation for Si/MgO/CoFeB/Ta/Nb superconducting diode
// Stack orientation: Si substrate at bottom, Nb at top
// Sweeping CoFeB thickness to tune PMA and exchange field
// Magnetic field propagates UPWARD from CoFeB through Ta into Nb

// Device geometry
Nx := 128
Ny := 256
cellsize := 16e-9  // 16 nm cells -> 2.048 μm × 4.096 μm device

// CoFeB material parameters (annealed)
Msat = 1.4e6      // A/m (1400 emu/cm³)
Aex = 20e-12      // J/m (20 pJ/m)
alpha = 0.015     // Gilbert damping

// Thickness sweep: 0.6 nm to 2.0 nm in 0.2 nm steps
t_cofeb_06 := 0.6e-9
t_cofeb_08 := 0.8e-9
t_cofeb_10 := 1.0e-9
t_cofeb_12 := 1.2e-9
t_cofeb_14 := 1.4e-9
t_cofeb_16 := 1.6e-9
t_cofeb_18 := 1.8e-9
t_cofeb_20 := 2.0e-9

// Interfacial anisotropy constant
Ki := 1.85e-3  // J/m² for CoFeB/MgO interface

// Function to run simulation for a given thickness
thickness := t_cofeb_06
SetGridSize(Nx, Ny, 1)
SetCellSize(cellsize, cellsize, thickness)

// Calculate effective anisotropy
Keff := Ki / thickness
Ku1 = Keff
AnisU = vector(0, 0, 1)  // Out-of-plane easy axis (upward)

// Initial magnetization: out-of-plane pointing upward (+z)
m = uniform(0, 0, 1)

// Relax to equilibrium
relax()

// Save magnetization state
save(m)
saveas(m, "m_CoFeB_0.6nm")

// Save demagnetization field (stray field)
saveas(B_demag, "Bfield_CoFeB_0.6nm_at_interface")

tableautosave(10e-12)
tablesave()

print("Thickness: 0.6 nm complete")

// Repeat for 0.8 nm
thickness = t_cofeb_08
SetCellSize(cellsize, cellsize, thickness)
Keff = Ki / thickness
Ku1 = Keff
m = uniform(0, 0, 1)
relax()
saveas(m, "m_CoFeB_0.8nm")
saveas(B_demag, "Bfield_CoFeB_0.8nm_at_interface")
tablesave()
print("Thickness: 0.8 nm complete")

// Repeat for 1.0 nm
thickness = t_cofeb_10
SetCellSize(cellsize, cellsize, thickness)
Keff = Ki / thickness
Ku1 = Keff
m = uniform(0, 0, 1)
relax()
saveas(m, "m_CoFeB_1.0nm")
saveas(B_demag, "Bfield_CoFeB_1.0nm_at_interface")
tablesave()
print("Thickness: 1.0 nm complete")

// Repeat for 1.2 nm
thickness = t_cofeb_12
SetCellSize(cellsize, cellsize, thickness)
Keff = Ki / thickness
Ku1 = Keff
m = uniform(0, 0, 1)
relax()
saveas(m, "m_CoFeB_1.2nm")
saveas(B_demag, "Bfield_CoFeB_1.2nm_at_interface")
tablesave()
print("Thickness: 1.2 nm complete")

// Repeat for 1.4 nm
thickness = t_cofeb_14
SetCellSize(cellsize, cellsize, thickness)
Keff = Ki / thickness
Ku1 = Keff
m = uniform(0, 0, 1)
relax()
saveas(m, "m_CoFeB_1.4nm")
saveas(B_demag, "Bfield_CoFeB_1.4nm_at_interface")
tablesave()
print("Thickness: 1.4 nm complete")

// Repeat for 1.6 nm
thickness = t_cofeb_16
SetCellSize(cellsize, cellsize, thickness)
Keff = Ki / thickness
Ku1 = Keff
m = uniform(0, 0, 1)
relax()
saveas(m, "m_CoFeB_1.6nm")
saveas(B_demag, "Bfield_CoFeB_1.6nm_at_interface")
tablesave()
print("Thickness: 1.6 nm complete")

// Repeat for 1.8 nm
thickness = t_cofeb_18
SetCellSize(cellsize, cellsize, thickness)
Keff = Ki / thickness
Ku1 = Keff
m = uniform(0, 0, 1)
relax()
saveas(m, "m_CoFeB_1.8nm")
saveas(B_demag, "Bfield_CoFeB_1.8nm_at_interface")
tablesave()
print("Thickness: 1.8 nm complete")

// Repeat for 2.0 nm
thickness = t_cofeb_20
SetCellSize(cellsize, cellsize, thickness)
Keff = Ki / thickness
Ku1 = Keff
m = uniform(0, 0, 1)
relax()
saveas(m, "m_CoFeB_2.0nm")
saveas(B_demag, "Bfield_CoFeB_2.0nm_at_interface")
tablesave()
print("Thickness: 2.0 nm complete")

print("Simulation complete. Outputs saved for Superscreen input.")
print("Stack geometry: Si(substrate)/MgO/CoFeB/Ta/Nb(top)")
print("Magnetic field propagates upward from CoFeB into Nb layer")