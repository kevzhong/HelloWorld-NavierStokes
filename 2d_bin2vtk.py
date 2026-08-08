# Read in the raw binary .dat files dumped by the code and write as rectlinear VTK files (.vtr)
# Files are dumped in single precision (float32) to be economical with memory
import numpy as np
from vtk import vtkXMLRectilinearGridWriter, vtkRectilinearGrid, vtkFloatArray
from vtk.util.numpy_support import numpy_to_vtk

tstart = 1000
tend = 24000

tskip = 1000

#fields = ['c','u','w']
fields = ['uxz']

x_coords = np.fromfile("outputdir/xm_grid.dat")
y_coords = np.fromfile("outputdir/zm_grid.dat")
z_coords = np.zeros_like(x_coords, dtype=np.float32)  # 2D grid

x_coords = np.float32(x_coords)
y_coords = np.float32(y_coords)

Nx = x_coords.size
Ny = y_coords.size

# Create the rectilinear grid
grid = vtkRectilinearGrid()
grid.SetDimensions(Nx, Ny, 1)  # 1 slice in Z-direction

# Convert the coordinate arrays to VTK format
x_vtk = numpy_to_vtk(x_coords)
y_vtk = numpy_to_vtk(y_coords)
z_vtk = numpy_to_vtk(z_coords)

# Set coordinates
grid.SetXCoordinates(x_vtk)
grid.SetYCoordinates(y_vtk)
grid.SetZCoordinates(z_vtk)

writer = vtkXMLRectilinearGridWriter()

# Loop over fields and timesteps

for fpref in fields:
    for nt in range(tstart,tend+tskip,tskip):
        print(f"Writing outputdir/{fpref}_it{nt}.vtr" )
        
        field_fname = f"outputdir/{fpref}_it{nt}.dat"
        vtr_fname = f"outputdir/{fpref}_it{nt}.vtr"

        field_data = np.fromfile(field_fname, dtype=np.float64).reshape(Nx, Ny, order='F')
        field_data = np.float32(field_data)

        field_vtk = numpy_to_vtk(field_data.ravel(order='F'))

        # Add the field data to the grid
        grid.GetPointData().SetScalars(field_vtk)

        # Create the writer and set the file name
        writer.SetFileName(vtr_fname)
        writer.SetInputData(grid)

        # Write the file in binary format
        writer.SetDataModeToBinary()
        writer.Write()

