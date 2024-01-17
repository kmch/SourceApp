# import sys
# sys.path.append('..') # to access plotea
from plotea.mpl2d import Shade

from typing import Tuple, Dict, List #, Final, Iterator, , Optional, Union
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes._axes import Axes
import pandas as pd

# import networkx as nx
import xarray as xr
import pyproj
import geopandas as gpd
import rasterio

def convert_xy_to_ij(x, y, extent, dx, dy):
  x1, x2, y1, y2 = extent
  i = int((x - x1) / dx)
  j = int((y2 - y) / dy) # yes!! because the array is numbered from top y
  return i, j

def gauss(x, mu=0, sigma=1):
  """
  Unnormalized Gaussian distribution.
  
  Parameters
  ----------
  
  Returns
  -------
  y : type(x)
    Gaussian evaluated at x. 
  
  Notes
  -----
  Some people use alpha (1/e point)
  instead of the sigma (standard deviation)
  to define the width of the Gaussian. 
  They are related through: alpha = sigma * sqrt(2)
  
  """
  return np.exp(-((x - mu)**2) / (2 * sigma**2))
def gaussND(dims, centers, radius, **kwargs):
  """
  This should be vectorized for speed.
  
  """
  coords = np.indices(dims)
  A = np.zeros(dims)
  for center in centers: 
    distance = np.sqrt((coords[0] - center[0])**2 + (coords[1]-center[1])**2 + (coords[2]-center[2])**2) 
    a = gauss(distance, 0, radius)
    A += a
  return A

def create_coord_transform(inp_coord='latlon', out_coord='bng'):
  code = {
    'latlon': 'EPSG:4326',    # WGS84 (latitude and longitude)
    'openstreet': 'EPSG:3857', 
    'bng'   : 'EPSG:27700',   # British National Grid
  }

  # Define coordinate reference systems
  inp = pyproj.CRS(code[inp_coord])  
  out = pyproj.CRS(code[out_coord])  

  # Create a pyproj transform
  transformer = pyproj.Transformer.from_crs(inp, out, always_xy=True)
  return transformer
def log10(df):
  return df.applymap(lambda x: np.log10(x) if pd.notnull(x) else x)
def log(df): # alias
  return log10(df)
def round2accuracy(to_round, accuracy):
    return np.round(to_round / accuracy) * accuracy
def sort(datastruct, flip=True):
  
  if type(datastruct) == pd.DataFrame:
    ascending = True if flip else False
    df = datastruct
    means = df.mean()
    sorted_columns = means.sort_values(ascending=ascending).index
    df_sorted = df[sorted_columns]
    return df_sorted
  elif type(datastruct) == dict:
    return {k: datastruct[k] for k in sorted(datastruct)}
  else:
    raise TypeError('Unsupported type: %s' % type(datastruct))

class SampleNode:
  """
  For duck-typing of pyfastunmix.SampleNode
  """
  def __init__(self):
    self.name = None
    self.x = None
    self.y = None
    self.area = None
    self.total_upstream_area = None
    self.downstream_node = None
    self.my_export_rate = None
    self.my_flux = None
    self.my_total_flux = None
    self.my_total_tracer_flux = None
    self.my_tracer_flux = None
    self.my_tracer_value = None
    self.upstream_nodes = None

def read_d8_header(file_path):
  with open(file_path, 'r') as file:
    header_lines = [next(file) for _ in range(6)]
  header_dict = {}
  for line in header_lines:
      line_parts = line.strip().split()
      if len(line_parts) == 2:
          key, value = line_parts
          header_dict[key.lower()] = float(value)

  xllcorner = header_dict.get('xllcorner', 0.0)
  yllcorner = header_dict.get('yllcorner', 0.0)
  cellsize = header_dict.get('cellsize', 1.0)
  ncols = int(header_dict.get('ncols', 0))
  nrows = int(header_dict.get('nrows', 0))

  x1 = xllcorner
  x2 = xllcorner + cellsize * ncols
  y1 = yllcorner
  y2 = yllcorner + cellsize * nrows

  dx = cellsize
  dy = cellsize
  nx = ncols
  ny = nrows

  return x1, x2, y1, y2, dx, dy, nx, ny
def extract_d8_subdomain(x1, x2, y1, y2, fname_d8_all, fname_d8, pad_width=10, pad_value=0):
  # Extract, pad & save to file 
  data_array = xr.open_dataset(fname_d8_all)['z'].sel(x=slice(x1, x2), 
                                                      y=slice(y1, y2))
  data_array_vals = data_array.transpose('y', 'x').values
  data_array_vals[np.isnan(data_array_vals)] = -9999
  data_array_vals = np.flipud(data_array_vals)
  # Pad with zeroes
  data_array_vals =  np.pad(data_array_vals, pad_width, 'constant', 
                            constant_values=(pad_value))

  # Save to file (incl. padding in the header)
  dx = data_array.x[1] - data_array.x[0]
  dx = int(dx.values)
  ny = data_array.shape[1] + 2 * pad_width
  nx = data_array.shape[0] + 2 * pad_width
  newx1 = data_array.x.values[0] - pad_width * dx
  newx2 = data_array.x.values[-1] + pad_width * dx
  newy1 = data_array.y.values[0] - pad_width * dx
  newy2 = data_array.y.values[-1] + pad_width * dx
  header = f"ncols {ny}\n"
  header += f"nrows {nx}\n"
  header += f"xllcorner {newx1}\n"
  header += f"yllcorner {newy1}\n"
  header += f"cellsize {dx}\n"

  with open(fname_d8, 'w') as f:
    f.write(header)
    np.savetxt(f, data_array_vals, fmt='%d')
  data = np.loadtxt(fname_d8, skiprows=6)
  data = data[:-1, :-1] # invert yaxis?
  
  # Create coordinate arrays
  x_coords = np.linspace(newx1, newx2, data.shape[1])
  y_coords = np.linspace(newy1, newy2, data.shape[0])
  # Create xarray DataArray
  data = xr.DataArray(data, coords=[y_coords, x_coords], dims=['y', 'x'])

  return data




def exclude_thames_trunk(samples: pd.DataFrame) -> pd.DataFrame:
  return samples[samples['Water_Category_Details'] != 'Thames']
def exclude_all_trunk(samples: pd.DataFrame) -> pd.DataFrame:
  # labels copied from Excel in three steps
  string = """O2 Arena 2 3
  O2 Arena 2 2
  O2 Arena 2 1
  O2 Arena 1 3
  O2 Arena 1 2
  O2 Arena 1 1
  Bermondsey 3
  Bermondsey 2
  Bermondsey 1
  Tower Bridge 3
  Tower Bridge 2
  Tower Bridge 1
  Southwark Bridge 2
  Southwark Bridge 1
  Millenium Bridge 3
  Millenium Bridge 2
  Millenium Bridge 1
  Blackfrias Bridge 3
  Blackfrias Bridge 2
  Blackfrias Bridge 1
  Gabriels Wharf 2
  Gabriels Wharf 1
  Westminster Bridge 3
  Westminster Bridge 2
  Westminster Bridge 1
  H. Parliament3
  H. Parliament2
  H. Parliament1
  Lambeth Bridge 3
  Lambeth Bridge 2
  Lambeth Bridge 1
  Vauxhall Tower 3
  Vauxhall Tower 2
  Vauxhall Tower 1
  US Embassy 3
  US Embassy 2
  US Embassy 1
  Battersea Park 3
  Battersea Park 2
  Battersea Park 1
  Wandsworth Park 3
  Wandsworth Park 2
  Wandsworth Park 1
  Putney Bridge 3
  Putney Bridge 2
  Putney Bridge 1
  Putney Rowing 3
  Putney Rowing 2
  Putney Rowing 1
  Kew Bridge 3
  Kew Bridge 2
  Kew Bridge 1
  Kew Palace 3
  Kew Palace 2
  Kew Palace 1
  Richmond 3
  Richmond 2
  Richmond 1
  Teddington 2
  Teddington 1
  Kingston 3
  Kingston 2
  Kingston 1
  Upstream Erith 3
  Upstream Erith 2
  Upstream Erith 1
  Downstream Crossness 3
  Downstream Crossness 2
  Downstream Crossness 1
  Crossness Outfall 3
  Crossness Outfall 2
  Crossness Outfall 1
  Upstream Crossness 3
  Upstream Crossness 2
  Upstream Crossness 1
  Downstream Beckton 3
  Downstream Beckton 2
  Downstream Beckton 1 
  Beckton Outfall 3
  Beckton Outfall 2
  Beckton Outfall 1
  Erith 3
  Erith 2
  Erith 1
  Junction
  Thames
  WWTP"""
  labels = string.splitlines()
  for label in labels:
    label = label.strip()
    samples = samples[samples['Water_Category_Details'] != label]
  return samples
def find_period_with_most_data(samples: pd.DataFrame, \
  chemical: str, time_window) -> Tuple[int, pd.Timestamp, Dict[str, int], pd.DataFrame]:
  
  samples = exclude_thames_trunk(samples)

  df = select_single_chemical_subset(samples, chemical)
  
  # Exclude NaNs
  df = df.dropna(subset=[chemical])

  # print('After dropna: ', len(df))

  # Convert the 'Sampling_Date' column to pandas datetime format
  df['Sampling_Date'] = pd.to_datetime(df['Sampling_Date'])

  # Sort the DataFrame by the 'Sampling_Date' column
  df = df.sort_values('Sampling_Date')

  # Set the 'Sampling_Date' column as the index
  df = df.set_index('Sampling_Date')

  # Initialize variables to store the maximum count and the corresponding period
  max_count = 0
  max_period = None

  dates = {}
  counts = []
  best_period = None
  # Iterate over each date in the DataFrame
  for date in df.index:

      # Define the start and end dates for the period
      start_date = date
      end_date = start_date + pd.DateOffset(days=time_window-1)

      # Filter the DataFrame for the current period
      period_df = df.loc[start_date:end_date]

      # Count the number of samples in the period
      count = len(period_df)
      
      dates[date] = count
      
      # Check if the current period has a higher count than the previous maximum
      if count > max_count:
          max_count = count
          max_period = (start_date, end_date)
          best_period = period_df

  # Print the maximum count and the corresponding period
  # print("Maximum count:", max_count)
  # print("Period:", max_period)

  return max_count, max_period, dates, best_period
def select_single_chemical_subset(samples: pd.DataFrame, chemical: str):
  key = chemical
  non_chems = ['x_bng', 'y_bng', 'Sample Code', 'Location', 'Water_Category_Details', 'Year', 'Sampling_Date', 'Latitude', 'Longitude']
  chems = samples.drop(non_chems, axis=1).columns.tolist()
  return samples.drop([i for i in chems if i != key], axis=1)
def select_subset_drain(rivers, river_name):
  rivers = rivers[rivers[col1].str.contains(name)]
def select_subset(samples, rivers, name):
    col1 = 'WB_NAME'
    col2 = 'Water_Category_Details'
    if rivers is not None:
      rivers = rivers[rivers[col1].str.contains(name)]
    samples = samples[samples[col2].str.contains(name)]
    return rivers, samples
def split_dataset_into_chemicals(samples: pd.DataFrame) -> dict[str, pd.DataFrame]:
  non_chems = ['x_bng', 'y_bng', 'Sample Code', 'Location', 'Water_Category_Details', 'Year', 'Sampling_Date', 'Latitude', 'Longitude']
  chems = samples.drop(non_chems, axis=1).columns.tolist()
  di = {}
  for key in chems:
    di[key] = samples.drop([i for i in chems if i != key], axis=1)
  return di  


def set_lims(ax, river_name):
  rn = river_name
  if rn == 'lea':
    ax.set_xlim((530000,545000))
    ax.set_ylim((180000,192600))
  elif rn == 'hogsmill':
    ax.set_xlim((515000,523000))
    ax.set_ylim((162010,171000))
  elif rn == 'beverley':
    ax.set_xlim((519000,530000))
    ax.set_ylim((162010,177500))
  elif rn == 'wandle':
    ax.set_xlim((523500,532500))
    ax.set_ylim((163500,177500))
  elif rn == 'guc':
    ax.set_xlim((504300,520000))
    ax.set_ylim((175000,183000))
  else:
    raise ValueError('Unknown river name (%s)' % river_name)
  return ax
def set_lims_around_sample(ax, df, sample_id, pad):
  sdf = df[df['Sample.Code'] == sample_id]
  x, y = sdf['x_coordinate'], sdf['y_coordinate']
  x, y = float(x.iloc[0]), float(y.iloc[0])
  x1, x2 = x - pad, x + pad
  y1, y2 = y - pad, y + pad
  ax.set_xlim(x1,x2)
  ax.set_ylim(y1,y2)
  return ax



def plot_data_coverage_over_time(dates: dict, \
  max_count: int, time_window: int, chemical: str):
  plt.subplots(figsize=(15,10))

  # Convert the dictionary to separate lists for x and y values
  x = list(dates.keys())
  y = list(dates.values())

  # Plot the data
  plt.plot(x, y, '.-')

  # Customize the plot
  plt.xlabel('First day of a %s-day time-window' % time_window)
  plt.ylabel('No. of non-NaN samples in the time-window')
  plt.title('%s - max count: %s' % (chemical, max_count))
  plt.xticks(rotation=45)  # Rotate x-axis labels for better readability
  return plt.gca()
def plot_samples_in_BNG_color(samples: pd.DataFrame, c: str, ax: Axes, **kwargs) -> Axes:
  xcol = 'x_bng'
  ycol = 'y_bng'
  sc1 = ax.plot(samples[xcol].values, samples[ycol].values, 'o', color=c, **kwargs)

  return plt.gca()
def plot_samples_in_BNG(samples: pd.DataFrame, chem: str, ax: Axes, \
  col='chem', nans=False, **kwargs) -> Axes:
  xcol = 'x_bng'
  ycol = 'y_bng'
  if not nans:
    kwargs['edgecolor'] = kwargs.get('edgecolor', 'none')
  

  if col == 'chem':
    if nans:
      # kwargs['edgecolor'] = kwargs.get('edgecolor', 'r')
      # kwargs['facecolor'] = kwargs.get('facecolor', 'none')
      sc1 = ax.scatter(samples[xcol].values, samples[ycol].values, \
                       facecolor='none', edgecolor='r', **kwargs)
    else:
      kwargs['cmap'] = kwargs.get('cmap', 'viridis')
      sc1 = ax.scatter(samples[xcol].values, samples[ycol].values, 
                      c=np.log10(samples[chem].values), **kwargs)
      cbar = colorbar(sc1, ax)
      cbar.set_label('log10(concentration)')
  elif col == 'date':
    import matplotlib.dates as mdates
    samples = samples.copy()
    samples.reset_index(inplace=True)
    dates1 = pd.to_datetime(samples['Sampling_Date'])
    dates1 = mdates.date2num(dates1)
    if nans:
      # print(';kdhaf')
      # kwargs['edgecolor'] = kwargs.get('edgecolor', 'r')
      # kwargs['facecolor'] = kwargs.get('facecolor', 'none')
      sc1 = ax.scatter(samples[xcol].values, samples[ycol].values, 
                    edgecolor='r', facecolor='none', **kwargs)
    # kwargs['cmap'] = kwargs.get('cmap', 'tab20c')
    else:
      sc1 = ax.scatter(samples[xcol].values, samples[ycol].values, 
                    c=dates1, **kwargs)
    cbar = colorbar(sc1, ax, format=mdates.DateFormatter('%Y-%m-%d'))
    cbar.set_label('Sampling date')
  else:
    raise TypeError()

  return plt.gca()
def plot_samples_in_BNG_dates(samples: pd.DataFrame, chem: str, ax: Axes, **kwargs) \
  -> Axes:
  xcol = 'x_bng'
  ycol = 'y_bng'
  kwargs['edgecolor'] = kwargs.get('edgecolor', 'none')
  kwargs['cmap'] = kwargs.get('cmap', 'viridis')
  sc1 = ax.scatter(samples[xcol].values, samples[ycol].values, 
                    c=np.log10(samples[chem].values), **kwargs)
  cbar = colorbar(sc1, ax)
  cbar.set_label('log10(concentration)')
  return plt.gca()
def plot_topography(X, Y, Z, **kwargs):
  fig, ax = plt.subplots(figsize=[12,8])
  shader = Shade()
  extent=[X.min(), X.max(), Y.min(), Y.max()]
  shader.plot(Z.T, cmap='Greys_r', extent=extent, aspect='equal',
              azdeg=45, altdeg=45, label='Topography, m a.s.l.', **kwargs)
  ax.set_xlabel('X - British National Grid, m')
  ax.set_ylabel('Y - British National Grid, m')
  ax.invert_yaxis()
  return fig, ax

