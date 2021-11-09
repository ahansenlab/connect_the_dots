import pandas as pd
import numpy as np


# load trajectory
df = pd.read_csv(snakemake.input.traj)

# load chromatic shift data
chromatic_shift_data = np.load(snakemake.input.chromatic_shift)

shiftX = chromatic_shift_data['chromatic_shift_X']
shiftY = chromatic_shift_data['chromatic_shift_Y']
shiftZ = chromatic_shift_data['chromatic_shift_Z']

def chromatic_aberration_correction(pos,shift_data):
    """
    Takes in a pixel position and returns the chromatic shift corrected position
    Parameters
    ----------
    pos : float or np.array of floats 
        Position in pixel values of a dot for one coordinate dimension. Thus, `pos`
        should be a Channel_1 position, as Channel_0 is taken as the reference.
        
        
    shift_data : numpy.ndarray
        Array containing two columns. The first column is the pixel position, the second
        is the shift at that pixel. By default it is assumed that the shift is calulcated 
        was Channel_1-Channel_0. 
        
    Output
    ------
    shift : float
        Chromatic shift between the channels
    """
    # retrieve file
    return np.interp(pos,shift_data[:,0],shift_data[:,1])


# overwrite the positions with the corrected, shifted positions
df['x'] = df.apply(lambda col: col.x - chromatic_aberration_correction(col.x,shiftX) 
                       if col.channel!=0 else col.x, axis=1)
df['y'] = df.apply(lambda col: col.y - chromatic_aberration_correction(col.y,shiftY) 
                       if col.channel!=0 else col.y, axis=1)
df['z'] = df.apply(lambda col: col.z - chromatic_aberration_correction(col.z,shiftZ) 
                       if col.channel!=0 else col.z, axis=1)


df.to_csv(snakemake.output.traj)