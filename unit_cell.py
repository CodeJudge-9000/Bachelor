atomicRadius = {
    'S' : 1.80,
    'Mo' : 1.39
}

def centUnitCell(system, addL = 0, timesVDW = 1):
    """
    Resizes unit cell to encompas all atoms, then centers the cell.
    
    For optimizations be sure to add additional spacing to the unit cell.
    
    Units are in angstrom.
    
    Input:
     - ase system of atoms
     - (optional) length added to all sides of unit cell
     - (optional) amount of Van Der Waals radii to add to the unit cell in all directions. Default is 1.
    
    Output:
     - None
    """
    
    # Get farthest distances between outer particles in all directions
    # Also get the symbols of the outer atoms at the same time
    p = system[0].position
    xLow, xHigh = p[0], p[0]
    yLow, yHigh = p[1], p[1]
    zLow, zHigh = p[2], p[2]
    
    outerType = ''
    
    for a in system:
        tempPos = a.position
        tempSymb = a.symbol
        
        # First x
        if tempPos[0] < xLow:
            xLow = tempPos[0]
            outerType = tempSymb
        elif tempPos[0] > xHigh:
            xHigh = tempPos[0]
            outerType = tempSymb
        
        # Then y
        if tempPos[1] < yLow:
            yLow = tempPos[1]
            outerType = tempSymb
        elif tempPos[1] > yHigh:
            yHigh = tempPos[1]
            outerType = tempSymb
            
        # Then z
        if tempPos[2] < zLow:
            zLow = tempPos[2]
            outerType = tempSymb
        elif tempPos[2] > zHigh:
            zHigh = tempPos[2]
            outerType = tempSymb
    

    # Redefine the unit cell size
    R = atomicRadius[outerType]
    system.set_cell([xHigh - xLow + timesVDW*R + addL, yHigh - yLow + timesVDW*R + addL, zHigh - zLow + timesVDW*R + addL])
    system.center()
    
    return xLow, xHigh, yLow, yHigh, zLow, zHigh
