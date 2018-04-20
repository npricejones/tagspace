from datetime import datetime

def gettimestr():
   current = datetime.utcnow()
   timestr = '{0:04d}-{1:02d}-{2:02d}.{3:02d}.{4:02d}.{5:02d}.{6:05d}'.format(current.year,
                                                               current.month,
                                                            current.day,
                                                            current.hour,
                                                            current.minute,
                                                            current.second,
                                                            current.microsecond)
   return timestr

# Function to slice a h5py dataset on element number, element name, pixel number of wavelength