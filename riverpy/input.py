"""
TODO
Define 'magic' column names such as 'id', 'x', 'y', 'lat', 'lon', 'date'
as global vars. 
"""
from autologging import logged, traced
from bs4 import BeautifulSoup
from datetime import datetime, timedelta
import json
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, NoNorm
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import numpy as np
import os
import pandas as pd
import random
import requests
from sklearn.cluster import DBSCAN
import subprocess
import urllib.request

import autocatchments as ac
from compysitional.pca import CompositionalPCA
from iogeo.api import read
from plotea.mpl2d import *
from riverpy.ioapi import *
from riverpy.plot import *
from riverpy.utils import create_coord_transform, round2accuracy,\
   extract_d8_subdomain, sort, convert_xy_to_ij
import sample_network_unmix as snu
from sewage import get_current_discharge_status

@logged
class Geospatial:
  def get_extent(self):
    """
    Works only for object having associated DataFrames.
    """
    x1 = self.df['x'].min()
    x2 = self.df['x'].max()
    y1 = self.df['y'].min()
    y2 = self.df['y'].max()
    return [x1, x2, y1, y2]
  def get_extent_zoom(self, x, y, pad=1000):
    x1, x2 = x - pad, x + pad
    y1, y2 = y - pad, y + pad
    return x1, x2, y1, y2
  def _restrict_to_catchment(self, df, catch, xcol='x', ycol='y'):
    self.__log.info('Excluding samples outside of the catchment')
    basin = catch.whole
    x_, y_ = [], []
    for x, y in zip(df[xcol], df[ycol]):
      i, j = convert_xy_to_ij(x, y, catch.extent, 50, 50)
      # NOTE the reversed index order!
      if (j >= 0 and j < basin.shape[0]) and (i >= 0 and i < basin.shape[1]) and basin[j, i]:
        x_.append(x)
        y_.append(y)
    df = df[(df[xcol].isin(x_)) & (df[ycol].isin(y_))]
    return df
  def _restrict_to_extent(self, df, extent, xcol='x', ycol='y'):
    x1, x2, y1, y2 = extent
    self.__log.debug(f'Excluding records outside the [{extent}] domain')
    df = df.loc[((df[xcol] > x1) & (df[xcol] < x2) & (df[ycol] > y1) & (df[ycol] < y2))]
    return df

@logged
class Acquisition(Geospatial): # this is superflous
  def __init__(self, data_samples):
    self.data_samples = data_samples
    self.network = data_samples.graph.graph
  def create(self, mode='from_data_samples', **kwargs):
    if mode == 'from_data_samples':
      out = self._create_from_data_samples(**kwargs)
    else:
      raise ValueError(f'Unknown mode: {mode}')
    return out
  def _create_from_data_samples(self, data_samples):   
    self.data_samples = data_samples
    self.network = data_samples.graph.graph     
@logged
class Catch(Geospatial):
  @classmethod
  def init(cls, ID, *args, **kwargs):
    if ID == 'wandle':
      subclass = CatchWandle
    elif ID == 'hogsmill':
      subclass = CatchHogsmill
    elif ID == 'thames':
      subclass = CatchThames   
    else:
      raise NotImplementedError(f'ID {ID}')
    return subclass(*args, **kwargs)
  def get_watershed(self, drainage, thresh=1e6):
    max_down = self.most_downstream_point
    self.drainage = drainage
    self.extent = drainage.extent
    self.__log.info(f'Using most-downstream point: {max_down}')
    self.__log.info(f'Snapping to drainage using {thresh} m2 threshold.')
    ds = DataSamples(pd.DataFrame({
       'id': ['1'], 
       'x': max_down[0], 
       'y': max_down[1],
       'lat': None,
       'lon': None}))
    ds.process(drainage, thresh=thresh, relocate_dict={}, use_letters=0)
    self.subcatchments = snu.get_unique_upstream_areas(ds.graph.graph)
    return self.subcatchments  
  def plot(self, ax=None, c='k', lw=1, ls='-', zorder=None):
    ax = get_ax(ax)
    extent = self.extent.copy() 
    extent[-2], extent[-1] = extent[-1], extent[-2]
    for a in self.subcatchments.values():
      ax.contour(a, extent=extent, colors=c, linewidths=lw, linestyles=ls, zorder=zorder)
    ax.set_aspect('equal') 
    return ax  
  def plot_most_downstream(self, ax=None, marker='*', s=20, **kwargs):
    ax = get_ax(ax)
    x, y = self.most_downstream_point
    ax.scatter(x, y, marker=marker, s=s, **kwargs)
    return ax
@logged
class CatchWandle(Catch):
  def __init__(self):
    self.most_downstream_point = (525550, 175100)
@logged
class CatchHogsmill(Catch):
  def __init__(self):
    self.most_downstream_point = (517750, 169170)
@logged
class CatchThames(Catch):
  def __init__(self):
    self.most_downstream_point = (552537, 169170)

@logged
class Catchment(Geospatial):
  """
  Catchment defined by its subcatchments a sample-network.
  """
  def __init__(self, network, drainage_extent):
    self.extent = drainage_extent
    self.subcatchments = snu.get_unique_upstream_areas(network)
    # whole catchment:
    self.whole = np.logical_or.reduce(list(self.subcatchments.values()))
  def plot(self, ax=None, c='k', lw=1, ls='--', zorder=None):
    ax = get_ax(ax)
    extent = self.extent.copy() 
    extent[-2], extent[-1] = extent[-1], extent[-2]
    for a in self.subcatchments.values():
      ax.contour(a, extent=extent, colors=c, linewidths=lw, linestyles=ls, zorder=zorder)
    ax.set_aspect('equal') 
    AxesFormatter()._set_xylabels()
    return ax  
@logged
class CatchmentOptimal(Catchment):
  def __init__(self, subcatchments: dict, drainage_extent):
    self.subcatchments = subcatchments
    self.extent = drainage_extent

@traced
@logged
class CSOs(Geospatial):
  def __init__(self, source='edm-api', extent=None, **kwargs):
    self.xcol = 'Easting'
    self.ycol = 'Northing'
    self.tcol = 'Total spill duration (hours)'
    if source is not None:
      self.read(source, extent=extent, **kwargs)
  def get_data(self, date_start, date_end):
    df = read('../../data/sewage-all-discharge-alerts-till-23-12-11.csv', index_col=0)
    df = df[df.AlertType == 'Start']
    df['date'] = pd.to_datetime(df['DateTime']).dt.date
    from datetime import datetime, timedelta
    start_date = datetime.strptime(date_start, '%Y-%m-%d').date()
    end_date = datetime.strptime(date_end, '%Y-%m-%d').date()

    # Iterate from start_date to end_date
    d = {}
    current_date = start_date
    while current_date <= end_date:
      d[current_date] = len(df[df['date'] == current_date]) 
      current_date += timedelta(days=1)
    return d
  def plot(self, ax=None, extent=None, marker='x', c='k', lw=2, \
    s=50, s_factor=1, cbar=False, cbar_label=None, **kwargs):
    ax = get_ax(ax)
    df = self._restrict_to_extent(self.df, extent) if extent is not None else self.df
    if c == 'duration':
      c = df[self.tcol]
      cbar_label = kwargs.get('cbar_label', self.tcol)
    if s == 'duration':
      s = df[self.tcol] * s_factor
    sc = ax.scatter(df['Easting'], df['Northing'], marker=marker, c=c, s=s, lw=lw, **kwargs)  
    if cbar:
      colorbar_flexible(sc, ax, cbar_label)
    AxesFormatter()._set_xylabels()
    return ax
  def read(self, source='edm-api', extent=None, restrict_to_catch=True,\
     **kwargs):
    if source == 'defra':
      df = self.read_from_DEFRA(**kwargs)
    elif source == 'edm-api':
      df = self.read_from_EDM()
    elif source == 'edm-starts':
      df = self.read_from_EDM_starts(**kwargs) 
    elif source == 'edm-xlsx':     
      year = kwargs['year']
      permit_column = 'EA Permit Reference\n(EA Consents Database)' \
        if year == 2021 else 'Permit No.'
      df = self.read_from_EDM_xlsx(year, permit_column)
    else:
      raise ValueError('Unkown source %s' % source)  
    self.df = self._restrict_to_extent(df, extent) if extent is not None else df
    if restrict_to_catch:
      catch = kwargs['catch']
      self.df = self._restrict_to_catchment(self.df, catch, xcol='Easting', ycol='Northing')
    return self.df
  def read_from_DEFRA(self, fname='../data/water-discharges.csv', \
    only_active=True):
    """
    Downloaded from:
    https://environment.data.gov.uk/public-register/downloads/water-discharges

    Parameters
    ----------
    fname : str, optional
        _description_, by default '../../data/water-discharges/water-discharges.csv'
    """
    df = pd.read_csv(fname)
    if only_active:
      df = df[pd.isna(df['Revocation Date'])]
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    return df
  def read_from_EDM_starts(self, fname='../data/all-EDM-discharge-starts.csv', \
    drop_duplicates=True):
    df = read(fname)
    if drop_duplicates:
      df = df.drop_duplicates(subset=['X', 'Y'])
    df = df.rename(columns={'X': self.xcol, 'Y': self.ycol})
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    return df
  def read_from_EDM(self):
    df = get_current_discharge_status()
    df = df.rename(columns={'X': self.xcol, 'Y': self.ycol})   
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x) 
    return df    
  def read_from_EDM_xlsx(self, year, permit_column: str):
    fname = f'../../data/edm-annual-return-{str(year)}.xlsx'
    skiprows = None if year >= 2021 else lambda x: x % 2 != 0
    df = pd.read_excel(fname, skiprows=skiprows)
    df.reset_index(drop=True, inplace=True)
    df = df.rename(columns=lambda x: x.rstrip())
    df.rename(columns={permit_column: 'PermitNumber'}, inplace=True)
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    if year == 2021:
      old = 'Total Duration (hrs) all spills prior to processing through 12-24h count method'
    else:
      old = 'Total Duration (hours) of all spills prior to processing through 12-24 hour counting method'
    df.rename(columns={old: self.tcol}, inplace=True)
    df[self.tcol] = pd.to_numeric(df[self.tcol], errors='coerce')
    # merge
    edmapi = self.read(source='edm-api')
    df = df.merge(edmapi, on='PermitNumber', how='left')
    return df
  def _annotate(self, ax, c, column='name', offset=None,):
    xcol, ycol = 'Easting', 'Northing'
    column = 'WWTP_NAME'
    offset = 130
    for i, row in df.iterrows():
        ax.annotate(row[column], (row[xcol]+2*offset, row[ycol]+offset), c='k', 
                    fontsize=12)      
  def _restrict_to_extent(self, df, extent):
    x1, x2, y1, y2 = extent
    xcol, ycol = self.xcol, self.ycol
    self.__log.debug(f'Excluding records outside the [{extent}] domain')
    df = df.loc[((df[xcol] > x1) & (df[xcol] < x2) & (df[ycol] > y1) & (df[ycol] < y2))]
    return df

@logged
class Domain(Geospatial):
  def __init__(self, x1, x2, y1, y2):
    self.x1 = x1
    self.x2 = x2
    self.y1 = y1
    self.y2 = y2    
    self.extent = np.array([x1, x2, y1, y2])
  def plot(self, *args, **kwargs):
    from plotea.mpl2d import plot_rect
    kwargs['c'] = kwargs.get('c', 'k')
    return plot_rect(*self.extent, **kwargs)
  def plot_patch(self, ax=None, c='k', lw=2, ls='--', fill=False, alpha=1):
    ax = get_ax(ax)
    
    x1, x2, y1, y2 = self.extent
    extent_rect = Rectangle((x1,y1), x2-x1, y2-y1, fill=fill, \
      color=c, linewidth=lw, linestyle=ls, alpha=alpha)
    ax.add_patch(extent_rect)
    return ax 
@logged
class DomainWandle(Domain):
  def __init__(self):
    super().__init__(521500, 540000, 151500, 179500)
@logged
class DomainHogsmill(Domain):
  def __init__(self):
    super().__init__(516500, 525500, 152500, 171000)
@logged
class DomainThames(Domain):
  def __init__(self):
    super().__init__(504332, 552537, 162037, 192521)

@logged
class Data(Geospatial):
  """
  Single-chemical data created from a dictionary
  and assigned 'value' chemical name in the df.
  Used to 

  Parameters
  ----------
  Geospatial : _type_
      _description_
  """
  def __init__(self, acquisition, values: dict):
    self.acquisition = acquisition
    di = acquisition.data_samples.xy
    self.xy = {k: di[k] for k in sorted(di)}
    self.dict = {k: values[k] for k in sorted(values)}
    dict1 = self.xy
    dict2 = self.dict
    df1 = pd.DataFrame.from_dict(dict1, orient='index', columns=['x', 'y'])
    df2 = pd.DataFrame.from_dict(dict2, orient='index', columns=['value'])
    df1.reset_index(inplace=True)
    df2.reset_index(inplace=True)
    self.df = pd.merge(df1, df2, on='index', suffixes=('', '_dict2'))
    self.df.rename(columns={'index': 'id'}, inplace=True)
    self.df['lat'] = None
    self.df['lon'] = None
    self.samples = DataSamples(self.df)
  # def prep(self, dataset='icl', river='wandle', max_dist=100):
  #   if dataset != 'icl':
  #     raise ValueError('dataset: %s' % dataset)
  #   dataset = DatasetICL()
  #   if river != 'wandle':
  #     raise ValueError('river: %s' % river)
  #   data_samples = dataset.select_subset(river)
  #   clusters = data_samples.cluster(max_dist=max_dist)
  #   means = clusters.export_statistic('mean')
  #   fname_d8 = drainage.fname
  #   means.process(drainage, fname_d8, {'F': [None, 170500]}, use_letters=1)
  #   stds = clusters.export_statistic('std')
  #   stds.process(drainage, fname_d8, {'F': [None, 170500]}, use_letters=1)
  #   acquisition = Acquisition(means)
  #   return acquiistion, data_samples, clusters, means, stds
  def plot(self, ax=None, **kwargs):
    ax = get_ax(ax)
    kwargs['column'] = kwargs.get('column', 'value')
    kwargs['annotate'] = kwargs.get('annotate', 'id')
    return self.samples.plot(ax=ax, **kwargs)
@logged
class DatasetICL(Geospatial):
  def __init__(self, path='../../data/ICL_London_Dataset_2019-2021.xlsx',\
    detection_threshold=3, add_to_extent=1000):
    """
    _summary_

    Parameters
    ----------
    path : str, optional
        _description_, by default '../../data/ICL_London_Dataset_2019-2021.xlsx'
    detection_threshold : float
         LLOD (ng/l)
    pad : int, optional
        _description_, by default 1000
    """
    self.path = path
    self.detection_threshold = detection_threshold
    self.df = self.read(path)
    self.clean()
    self.chems = self.get_list_of_chemicals()
    self.cecs = self.get_list_of_CECs()
    c_min = self.get_min_concentration(self.chems)
    self.__log.info('Min. measured concentration in the dataset: %s ng/l' % round(c_min, 2))
    self.fillna(fill_value=detection_threshold)
    self.extent = self._get_extent(pad=add_to_extent)
  def clean(self, key_lon='Longitude', key_lat='Latitude'):
    df = self.df
    
    self.__log.info('Removing trailing spaces from columns names...')
    df = df.rename(columns=lambda x: x.rstrip())
    
    self.__log.info('Renaming some columns...')
    df.rename(inplace=True, columns={
      'Sample Code': 'id',
      'Latitude': 'lat',
      'Longitude': 'lon',
      'Sampling_Date': 'date'
    })

    self.__log.info('Adding BNG coordinates as x,y columns...')
    transformer = create_coord_transform('latlon', 'bng')

    df['x'], df['y'] = transformer.transform(
      df['lon'].values, df['lat'].values) 

    # self.__log.info('Converting metres to km...')
    # df['x'], df['y'] = df['x'] / 1000., df['y'] / 1000.

    self.__log.info("Adding 'campaign' column...")
    df = self._add_column_campaign(df)
    

    df = df.rename(columns={
      'Mefenamic_acid': 'Mefenamic acid',
      'Salicylic_acid': 'Salicylic acid',
    })

    self.df = df

    return self.df
  def fillna(self, fill_value):
    self.__log.warning('Filling blanks with LLOD = %s ng/l' % fill_value)
    self.df[self.chems] = self.df[self.chems].fillna(fill_value)
  def get_min_concentration(self, chems):
    return np.min(self.df[chems])
  def get_list_of_chemicals(self, df=None):
    if df is None:
      df = self.df
    non_chems = ['x', 'y', 'campaign', 'id', 'Location', 'Water_Category_Details', 'Year', 'date', 'lat', 'lon']
    chems = sorted(df.drop(non_chems, axis=1).columns.tolist())
    self.chems = chems
    return self.chems
  def get_list_of_CECs(self):
    
    chems_all = self.chems
    chems = ['Salicylic acid', 'Cocaine', 'Benzoylecgonine', \
      'Imidacloprid', 'Venlafaxine', 'Acetamiprid', 
      'Diclofenac', 'Trimethoprim', 'Sulfamethoxazole',
      'Tramadol', 'Carbamazepine' # actually not recommended by Leon
    ]
    self.__log.info('Added Tramadol & Carbamazepine which are not CECs but still worth looking.')
    for ch in chems_all:
      if 'azole' in ch or 'mycin' in ch:
        chems.append(ch)

    chems = [i for i in chems if i != 'Lincomycin']
    self.__log.info('Excluding Lincomycin - no data.')
    chems = sorted(list(set(chems)))
    return chems
  def plot(self, ax=None, pad=1000, **kwargs):
    ax = get_ax(ax)
    
    x1 = min(self.df['x']) - pad
    x2 = max(self.df['x']) + pad
    y1 = min(self.df['y']) - pad
    y2 = max(self.df['y']) + pad    
    self.extent = [x1, x2, y1, y2]
    ds = DataSamples(self.df)
    return ds.plot(**kwargs)
  def read(self, file_path, **kwargs):
    self.df = read(file_path, **kwargs)
    return self.df
  def select_subset(self, river_name):
    samples = self.df
    rn = river_name
    col1 = 'Location'
    col2 = 'Water_Category_Details'
    if rn == 'lea':
      strings = ['Lea', 'Lee', 'Channelsea', 'Pymmes']
      condition = samples[col2].str.contains(strings[0])
      for string in strings[1:]:
        condition |= samples[col2].str.contains(string)
    elif rn == 'wandle':
      condition = samples[col2].str.contains('Wandle')
    elif rn == 'guc':
      strings = ['Brent', 'Grand_Union_Canal', 'Paddington_Arm', \
        'Frays_River', 'Slough_Arm']
      condition = samples[col2].str.contains(strings[0])
      for string in strings[1:]:
        condition |= samples[col2].str.contains(string)
      # extra condition
      condition |= samples[col1].str.contains('GUC')
    elif rn == 'hogsmill':
      strings = ['Hogsmill','Teddington_Lock_Canal', 'Teddington_Weir',\
        'Teddington_Kingston',  'Kingston', 'Teddington_Lock_Canal']
      condition = samples[col2].str.contains(strings[0])
      for string in strings[1:]:
        condition |= samples[col2].str.contains(string)
        # extra condition
      condition |= samples[col1].str.contains('Hogsmill')
    elif rn == 'beverley':
      condition = samples[col2].str.contains('Beverley_Brook')
    else:
      raise ValueError('Unknown river: %s' % rn)
    return DataSamples(samples[condition], data_cols=self.chems, dataset=self)
  def _add_column_campaign(self, df):
    divider_dates = ['2020-05-01', '2021-05-01']

    # Convert 'Sampling_Date' column to datetime if it's not already
    df['date'] = pd.to_datetime(df['date'])

    # Create a new column 'Group' based on the divider dates
    df['campaign'] = np.where(df['date'] < pd.to_datetime(divider_dates[0]), 1,
                          np.where(df['date'] < pd.to_datetime(divider_dates[1]), 2, 3))
    return df
  def _get_extent(self, pad=1000):
    x1 = min(self.df['x']) - pad
    x2 = max(self.df['x']) + pad
    y1 = min(self.df['y']) - pad
    y2 = max(self.df['y']) + pad    
    extent = [x1, x2, y1, y2]
    return extent
@logged
class DataSamples(Geospatial):
  def __init__(self, df: pd.DataFrame, data_cols=[], dataset=None):
    self.dataset = dataset
    self.df = df
    self.latlon = self._create_samples_dict_latlon()
    self.xy = self._create_samples_dict_xy()
    self.list = self._create_samples_list()
    self.dict = self._create_samples_dict()
    
    self.data_cols = data_cols
  def cluster(self, mode='dbscan', **kwargs):
    """
    _summary_

    Parameters
    ----------
    mode : str, optional
        _description_, by default 'dbscan'

    Returns
    -------
    pd.DataFrame
        Original df with 'cluster' column added to it.

    Raises
    ------
    ValueError
        _description_
    """
    if mode == 'dbscan':
      df = self._cluster_dbscan(**kwargs)
    else:
      raise ValueError('Unknown mode: %s' % mode)
    self.clusters = DataClusters(df, data_cols=self.data_cols, data_samples=self)
    return self.clusters
  def create_sample_graph(self, fname_d8): # NOTE 
    fname = 'tmp.csv'
    self._export_input_for_get_sample_graph(fname)
    sample_graph, _ = snu.get_sample_graphs(fname_d8, fname)
    os.remove(fname)
    self.graph = SampleGraph(sample_graph)
    if self.graph.no_nodes != len(self.df):
      self.__log.error('No. of graph nodes (%s) != no. of samples (%s). Try to decrease the threshold.' % \
        (self.graph.no_nodes, len(self.df)))     
  def get_catchment_areas(self, drainage):
    # in km3
    catch = Catchment(self.graph.graph, drainage.extent)
    d = catch.subcatchments
    d = {k: d[k] for k in sorted(d)}
    dx, dy = drainage.mg.dx, drainage.mg.dy
    areas = {k: v.astype(int).sum() * dx * dy * 1e-6 for k,v in d.items()}
    return areas
  def get_chem_dict(self, chem):
    return sort(self.df.set_index('id')[chem].to_dict())
  def get_ids(self):
    self.ids = list(sort(self.dict).keys()) # 
    return self.ids
  def get_sampling_dates(self):
    timestamps = set(self.df['date'])
    self.dates = [timestamp.to_pydatetime() for timestamp in timestamps]
    return self.dates
  def plot(self, ax=None, xcol='x', ycol='y', column=None, annotate=None, norm='log', cbar=False, cmap=None, \
    marker='o', vmin=None, vmax=None, edgecolor='k', zorder=10, s=100,\
      annotate_kw={}, s_factor=1, alpha=1., **kwargs):
    samples = self.df
    ax = get_ax(ax)

    x = samples[xcol].values
    y = samples[ycol].values

    if column is None:
      c = kwargs.get('c', 'r')
    elif column == 'id':
      ids = sorted(samples['id'])
      d = dict(list(zip(ids, range(len(ids)))))
      c = samples['id'].apply(lambda x: d[x])
    else:
      c = samples[column].values
      if s == 'val':
        s = c * s_factor

    if norm == 'log':
      norm = LogNorm(vmin=vmin, vmax=vmax)
      sc = ax.scatter(x, y, c=c, marker=marker, norm=norm, edgecolor=edgecolor, zorder=zorder, s=s, cmap=cmap, alpha=alpha)
    else:
      sc = ax.scatter(x, y, c=c, marker=marker, vmin=vmin, vmax=vmax, edgecolor=edgecolor, zorder=zorder, s=s, cmap=cmap, alpha=alpha)
    
    self.sc = sc
    if cbar:
      cb = plt.colorbar(sc, ax=ax)
      cb.set_label("Concentration (ng/l)")
    if annotate is not None:
      annotate_kw['zorder'] = annotate_kw.get('zorder', zorder + 1)
      self._annotate(annotate, ax=ax, **annotate_kw)
    
    ax.set_xlabel('Easting (km)')
    ax.set_ylabel('Northing (km)')
    AxesFormatter()._adjust(ax=ax, **kwargs)
    return ax
  def process(self, drainage, thresh:float, use_letters=False, relocate_dict={}):
    """
    Create an acyclic graph of samples after snapping them to the drainage
    given a certain upstream-area threshold. 
    Given the graph (topology), rename the samples (assign them new IDs)
    based on their order from upstream to downstream. 

    Parameters
    ----------
    drainage : _type_
        _description_
    thresh : float [m2]
        e.g. 1e7
    use_letters : bool, optional
        _description_, by default False
    relocate_dict : dict, optional
        _description_, by default {}
    """
    fname_d8 = drainage.fname
    self.snap_to_drainage(drainage, thresh=thresh)
    self.create_sample_graph(fname_d8)
    self.rename_based_on_topology(use_letters=use_letters)
    for sample_id, new_coords in relocate_dict.items():
      self.relocate_sample(sample_id=sample_id, new_coords=new_coords)
    self.snap_to_drainage(drainage, thresh=thresh)
    self.create_sample_graph(fname_d8)      
  def relocate_sample(self, sample_id, new_coords):
    x, y = new_coords
    df = self.df
    assert sample_id in df['id'].to_list()
    if x is not None:
      df.loc[df['id'] == sample_id, 'x'] = x
    if y is not None:
      df.loc[df['id'] == sample_id, 'y'] = y
    self.df = df
    # update all the attributes
    self.latlon = self._create_samples_dict_latlon()
    self.xy = self._create_samples_dict_xy()
    self.list = self._create_samples_list()
    self.dict = self._create_samples_dict()
  def rename_based_on_topology(self, use_letters=False, rjust=2):
    letters = [letter for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']
    df = self.df
    df.loc[:, 'id'] = df['id'].astype(str) # cause sample_name of the graph are str
    # print(df['id'])
    new_id = {}
    for i, (sample_name, my_data) in enumerate(snu.nx_topological_sort_with_data(self.graph.graph)):
      new_id[sample_name] = letters[i] if use_letters else str(int(i)+1).rjust(rjust, '0')        
    # update a column based on the mapping defined by the new_id dict
    self.df['id'] = self.df['id'].map(new_id)
    # update all the attributes
    self.latlon = self._create_samples_dict_latlon()
    self.xy = self._create_samples_dict_xy()
    self.list = self._create_samples_list()
    self.dict = self._create_samples_dict()
  def snap_to_drainage(self, drainage, thresh=1e7):  
    coords = np.array(self.df[['x', 'y']].values.tolist())
    snappd = ac.autosampler.snap_to_drainage(drainage.mg, coords, thresh)
    assert len(coords) == len(snappd)
    # self.__log.info('len coords vs. snappd', len(coords), len(snappd))
    self.df['x'] , self.df['y']  = ac.toolkit.model_xy_to_geographic_coords(
        (snappd[:, 0], snappd[:, 1]), drainage.mg)
  def zoom(self, id, ax, pad):
    self.dict[id].set_lims_around_sample(ax=ax, pad=pad)
    return ax
  def _add_cluster_id_to_df(self, df, cluster_dict):
    def find_cluster_id(id_value):
      for cluster_id, id_list in cluster_dict.items():
        if id_value in id_list:
          return cluster_id
      return None
    return df.assign(cluster=df['id'].apply(find_cluster_id))
  def _annotate(self, column, ax=None, xcol='x', ycol='y', \
    textcoords="offset points", xytext=(5,5), ha='center', va='center', **kwargs):
    if ax is None:
      ax = plt.gca()
    for i, row in self.df.iterrows():
      ax.annotate(row[column], (row[xcol], row[ycol]),\
         textcoords=textcoords, xytext=xytext, ha=ha,  **kwargs)
    return ax
  def _cluster_dbscan(self, max_dist, min_clusters=2):
    points_dict = self.xy
    points = np.array(list(points_dict.values()))
    dbscan = DBSCAN(eps=max_dist, min_samples=min_clusters)
    labels = dbscan.fit_predict(points)
    # Output formatting
    clusters = {}
    for idx, cluster_id in enumerate(labels):
        if cluster_id not in clusters:
            clusters[cluster_id] = []
        clusters[cluster_id].append(list(points_dict.keys())[idx])
    self.__log.info('Output: %s' % clusters)
    df = self._add_cluster_id_to_df(self.df, clusters)
    return df
  def _cluster_between_tributaries_wandle(dat):
    y1 = 172000
    y2 = 166000
    g1 = dat[dat['y_bng'] >= y1]
    g2 = dat[(dat['y_bng'] < y1) & (dat['y_bng'] > y2)]
    g3 = dat[dat['y_bng'] <= y2]
    groups = [g1, g2, g3]
    return groups
  def _create_samples_dict(self):
    samples = {}
    ids, xs, ys = self.df['id'], self.df['x'], self.df['y']
    lats, lons = self.df['lat'], self.df['lon']
    for (id, x, y, lat, lon) in zip(ids, xs, ys, lats, lons):
      samples[id] = DataSample(id, x, y, lat, lon, data_samples=self)
    return samples 
  def _create_samples_dict_latlon(self):
    latlon = {}
    ids, lats, lons = self.df['id'], self.df['lat'], self.df['lon']
    for (id, lat, lon) in zip(ids, lats, lons):
      latlon[id] =  (lat, lon)
    return latlon  
  def _create_samples_dict_xy(self):
    xy = {}
    ids, xs, ys = self.df['id'], self.df['x'], self.df['y']
    for (id, x, y) in zip(ids, xs, ys):
      xy[id] =  (x, y)
    return xy   
  def _create_samples_list(self):
    samples = []
    ids, xs, ys = self.df['id'], self.df['x'], self.df['y']
    lats, lons = self.df['lat'], self.df['lon']
    for (id, x, y, lat, lon) in zip(ids, xs, ys, lats, lons):
      samples.append(DataSample(id, x, y, lat, lon))
    return samples
  def _export_input_for_get_sample_graph(self, fname='tmp.csv'): # NOTE
    df = self.df.copy()[['id', 'x', 'y', *self.data_cols]]
    # yes, renaming is necessary
    df.rename(inplace=True, columns={
      'id': 'Sample.Code',
      'x': 'x_coordinate',
      'y': 'y_coordinate',
    })
    self.__log.debug('len(df)=%s' % len(df))
    df.to_csv(fname, sep=' ', index=False)
@logged
class DataSample(Geospatial):
  """
  Single data point measured in the field.
  """
  def __init__(self, id, x, y, lat, lon, data=None, data_samples=None):
    self.data_samples = data_samples
    self.id = id
    self.x = x
    self.y = y
    self.lat = lat
    self.lon = lon
    self.data = data
    self.df = self._create_df()
  def _create_df(self):
    """
    Return a DataFrame consisting of all attributes set
    in `__init__` 

    Returns
    -------
    pd.DataFrame
        _description_
    """
    d = {key: value for key, value in self.__dict__.items() if key != 'data'}
    if self.data is not None:
      d.update(self.data)
    return pd.DataFrame([d])
  def plot_xy(self, **kwargs):
    return plt.scatter([self.x], [self.y], **kwargs)
  def set_lims_around_sample(self, ax, pad):
    x, y = self.x, self.y
    x1, x2 = x - pad, x + pad
    y1, y2 = y - pad, y + pad
    ax.set_xlim(x1,x2)
    ax.set_ylim(y1,y2)
    return ax
@logged
class DataClusters(Geospatial):
  def __init__(self, df: pd.DataFrame, data_cols=[], data_samples=None):
    """
    Takes a dataframe containing 'cluster' column created by
    `DataSamples.cluster` method and then creates a dictionary 
    of `DataCluster` objects based on that.

    Parameters
    ----------
    df : pd.DataFrame
        _description_
    data_cols : list, optional
        _description_, by default []
    """
    self.data_samples = data_samples
    self.detection_threshold = self.data_samples.dataset.detection_threshold
    self.fillna = self.detection_threshold
    self.df = df
    self.data_cols = data_cols
    self.dict = self._create_clusters_dict()
  def export_statistic(self, statistic='mean'):
    dfs = []
    for key, cluster in self.dict.items():
      dfs.append(getattr(cluster, statistic).df)
    df = pd.concat(dfs, ignore_index=True)
    return DataSamples(df, data_cols=self.data_cols)
  def get_boxplot_df(self, chem):
    """
    Each cluster has a different no. of data samples 
    which makes prepping a single df a bit tricky.
    We concatenate dfs and fill them with NaNs so that 
    the outcome is box-diagonal.


    Parameters
    ----------
    chem : _type_
        _description_
    """
    di = {}
    for k, v in self.dict.items():
      df = v.df.drop([i for i in v.df.columns if i != chem], axis=1)
      df = df.rename(columns={chem: k})
      di[k] = np.log10(df)
    df = pd.concat(di.values(), axis=1)
    df.reset_index(drop=True, inplace=True)
    self.df_boxplot = df        
    return self.df_boxplot
  def get_df_for_boxplot(self, chem):
    """
    Each cluster has a different no. of data samples 
    which makes prepping a single df a bit tricky.
    We concatenate dfs and fill them with NaNs so that 
    the outcome is box-diagonal.


    Parameters
    ----------
    chem : _type_
        _description_
    """
    di = {}
    for k, v in self.dict.items():
      df = v.df.drop([i for i in v.df.columns if i != chem], axis=1)
      df = df.rename(columns={chem: k})
      di[k] = df
    df = pd.concat(di.values(), axis=1)
    df.reset_index(drop=True, inplace=True)
    self.df_boxplot = df        
    return self.df_boxplot
  def get_variance_ratios(self):
    self.vb, self.vw, self.vw_vb_ratio = {}, {}, {}
    for chem in self.data_cols:
      groups = self.get_groups_for_anova_etc(chem)
      # calculate between-group variance
      means = []
      for group in groups:
        means.append(np.mean(group))
      vb = np.var(means)
      # calculate mean within-group variance
      vw = []
      for group in groups:
        vw.append(np.var(group))
      vw = np.mean(vw)
      self.vb[chem], self.vw[chem], self.vw_vb_ratio[chem] = vb, vw, vw / vb
    return self.vb, self.vw, self.vw_vb_ratio
  def get_groups_for_anova_etc(self, chem):
    fillna = self.detection_threshold
    groups = []
    for k, cl in self.dict.items():
      groups.append(np.log10(cl.df[chem].fillna(fillna)))
    return groups
  def perform_anova(self):
    """
    Perform one-way ANOVA.

    Parameters
    ----------
    fillna : float, optional
        Small value (under detection threshold) to replace NaNs with, by default 1e-6

    Returns
    -------
    list of 2 dicts
        F_ratios and p_values
    """
    from scipy.stats import f_oneway
    F_ratios = {}
    p_values = {}
    for chem in self.data_cols:
      groups = self.get_groups_for_anova_etc(chem)
      F_ratios[chem], p_values[chem] = f_oneway(*groups)
    self.F_ratios = F_ratios
    self.p_values = p_values
    return self.F_ratios, self.p_values
  def plot(self):
    clrs = colors(len(self.dict))
    for key, val in self.dict.items():
      c = next(clrs)
  def plot_anova(self, chems: list, ax=None, vert=0, signif_level=0.05, style='bar',\
    yticklabels=False, c='k'):
    if ax is None:
      fig, ax = plt.subplots()      
    # select only contaminants of emerging concern (CECs)
    x, y = chems, []
    d = self.p_values
    for chem in chems:
      y.append(d[chem])
    # plot 
    if vert:
      x, y = y, x
      plt.plot(x, np.full(len(y), signif_level), 'r--', label='p=%s' % signif_level)
      plt.ylabel('p value')
    else:
      plt.vlines(x=signif_level, ymin=-1, ymax=len(chems), colors='k', linestyles='dotted', linewidth=1)
      plt.xlabel('p value')
   
    if style == 'bar':
      if vert:
        ax.bar(x,y)
      else:
        ax.barh(x,y)
    ax.plot(y, x, 'd', markerfacecolor='w', markeredgecolor=c)
    if not yticklabels:
      ax.set_yticklabels([])   
    return ax
  def plot_box(self, chem, c='b', **kwargs):
    # if not hasattr(self, 'df_boxplot'): # it would have to be for all chems
    self.get_boxplot_df(chem)
    ax = boxplot(self.df_boxplot, c, **kwargs)
    ax.set_xlabel('log10(c)')
    ax.set_ylabel('Subcatchment')
    ymin, ymax = ax.get_ylim()
    ax.vlines(x=np.log10(self.fillna), ymin=ymin, ymax=ymax, colors='k', linestyles='dotted', linewidth=2)
    
    return ax
  def plot_boxplot(self, chem, **kwargs):
    self.get_df_for_boxplot(chem)
    # print(self.df_boxplot)
    ax = plot_boxplot(self.df_boxplot, **kwargs)
    ax.set_xlabel('Concentration (ng/l)')
    ax.set_ylabel('Subcatchment')
    ymin, ymax = ax.get_ylim()
    ax.vlines(x=self.fillna, ymin=ymin, ymax=ymax, colors='k', linestyles='dotted', linewidth=2)
    ax.grid('none')
    return ax
  def plot_xy(self, ax=None, **kwargs):
    if ax is None:
      fig, ax = plt.subplots()  
    for key, cluster in self.dict.items():
      cluster.plot_xy(ax=ax, **kwargs)
    return ax
  def plot_variance(self, chems: list, ax=None, vert=0, cw='r', cb='Grey',\
    legend=False, recalculate=False, edgecolor='k'):
    if ax is None:
      fig, ax = plt.subplots()      
    if not (hasattr(self, 'vb') or recalculate):
      _ = self.get_variance_ratios() 
    vwn, vbn = [], []
    for chem in chems:
      vw = self.vw[chem]
      vb = self.vb[chem]
      N = vw + vb # normalisation
      vwn.append(100*vw/N)
      vbn.append(100*vb/N)  
    col_vw = 'vw'
    col_vb = 'vb'
    d = {'chem': chems, col_vw: vwn, col_vb: vbn}
    df = pd.DataFrame(d)
    colors = {'vw': cw, 'vb': cb}
    ax = df.plot(x='chem', ax=ax, kind='barh', stacked=True, legend=legend,\
      color=[colors[col] for col in df.columns[1:]], edgecolor=edgecolor,
    )
    ax.set_xlim(0,100)
    # plt.grid(axis='x', linestyle='--')
    ax.vlines(x=50, ymin=-1, ymax=len(chems), colors='k', linestyles='dashed', linewidth=1)
    ax.set_xlabel('% of variance')
    ax.set_ylabel('')
    return ax
  def update_ids(self, df: pd.DataFrame):
    assert ('cluster' in df.columns) and ('id' in df.columns)
    id_mapping = {k: v for (k,v) in zip(df['cluster'], df['id'])}
    self.df['cluster'] = self.df['cluster'].map(id_mapping)
    self.dict = self._create_clusters_dict()
  def _create_clusters_dict(self):
    cluster_dict = {}
    for cluster_id, cluster_df in self.df.groupby('cluster'):
      cluster_dict[cluster_id] = DataCluster(cluster_df, data_cols=self.data_cols, data_clusters=self)
    return cluster_dict
@logged
class DataCluster(Geospatial):
  """
  Aggregate of >1 data points measured in the field.
  """
  def __init__(self, df, data_cols=[], data_clusters=None):
    self.data_clusters = data_clusters
    self.detection_threshold = data_clusters.data_samples.dataset.detection_threshold
    self.df = df
    self.data_cols = data_cols
    self.xm = self.df['x'].mean()
    self.ym = self.df['y'].mean()
    cluster_id = self.create_id(mode='from_clustering')
    self.id = cluster_id
    self.nsamples = len(self.df)

    self.mean = DataSample(self.id, self.xm, self.ym, lat=None, lon=None, \
      data=dict(**self._get_statistic(func=np.mean), **{'cluster': cluster_id}))
    self.std = DataSample(self.id, self.xm, self.ym, lat=None, lon=None, \
      data=dict(**self._get_statistic(func=np.std), **{'cluster': cluster_id}))
  def _get_statistic(self, func=np.mean):
    """
    Calculate a statistic (e.g. mean or std) for each of the data columns
    stored in `self.data_cols` and return the result as a dictionary
    mapping the column name onto the statistic value.

    Parameters
    ----------
    func : _type_, optional
        converts a vector of numbers into a scalar, by default np.mean

    Returns
    -------
    dict
        {column: func(column), ...}
    """
    # "all of this is in the log space!"
    values = {} 
    for col in self.data_cols:
      values[col] = 10 ** func(np.log10(self.df[col])) 
    return values
  def plot(self, chem, vmin, vmax, ax=None, aspect='auto', s=50, pad=1):
    if ax is None:
      fig, ax = plt.subplots()     
    cid = self.id
    xm, ym = self.xm, self.ym
    
    #plt.annotate(cid, c='w', xy=np.array((xm,ym)), xytext=np.array((xm+100,ym+100))), #fontsize=12,  ha='center')
    # width, height = .3, .3
    # lef = offsets[cid][0]
    # bot = offsets[cid][1]
    # axins = ax.inset_axes([lef, bot, width, height], transform=ax.transAxes)
    df = self.df
    df_nan = df[df[chem].isna()]      
    df_notnan = df[df[chem].notna()] 
    x1, x2 = df['x'].min(), df['x'].max()
    y1, y2 = df['y'].min(), df['y'].max()
    
    norm = LogNorm(vmin=vmin, vmax=vmax)
    ax.scatter(df_notnan['x'], df_notnan['y'], s=s, edgecolor='k',
                  c=df_notnan[chem], norm=norm)
    ax.scatter(df_nan['x'], df_nan['y'], s=s, edgecolor='grey', facecolor='none')
    for i, row in df.iterrows():
      date_str = row['date'].strftime('%d.%m.%y')
      ax.annotate(date_str, xy=(row['x'], row['y']), xytext=(row['x'], row['y']))
    # def format_axins():
    #     pad = pads[cid]
    #     axins.set_title(cid)
    #     axins.set_xlim(xm-pad, xm+pad)
    #     axins.set_ylim(ym-pad, ym+pad)
    #     axins.set_xticks([])
    #     axins.set_yticks([])
    #     axins.set_xticklabels([])
    #     #plt.xticks(rotation=90)
    #     axins.set_yticklabels([])
    #     axins.set_aspect('equal')
    # format_axins()
    # ax.indicate_inset_zoom(axins, edgecolor="k", alpha=.1)  
    ax.set_xlim(x1-pad, x2+pad)
    ax.set_ylim(y1-pad, y2+pad)
    ax.set_aspect(aspect)
    return ax
  def plot_hist(self, ax=None, col='date', **kwargs):
    if ax is None:
      fig, ax = plt.subplots()    
    return ax.hist(self.df[col], **kwargs)
  def plot_conc(self, chem, dates, ax=None, vmin=1, vmax=1000, cbar=False, rot=0,\
    xtick_labels=True, ytick_labels=True, ylim=(-1.5,3.0), ylabel=True, decimate_xticklabels=None):
    if ax is None:
      fig, ax = plt.subplots() 
    df = self.df
    df['date'] = pd.to_datetime(df['date']).dt.strftime('%Y-%m-%d')
    xticks = dates
    
    x = [xticks.index(xi) for xi in df['date']]
    y = np.log10(df[chem])
    
    norm = LogNorm(vmin=vmin, vmax=vmax) # note, c is not log10-ed below
    self.sc = ax.scatter(x, y, c=df[chem], s=100, edgecolor='k', norm=norm)
    ax.plot(dates, np.full(len(dates), np.log10(self.detection_threshold)), 'r--')

    xticks = range(len(dates))
    labels = xticks# dates
    
    
    if xtick_labels:
      plt.xlabel('Date ID')
      plt.xticks(xticks, labels, rotation=rot)
    else:
      plt.xticks(xticks, [''] * len(xticks))
    
    if decimate_xticklabels is not None:
      ax = remove_every_nth_tick_label(ax=ax, n=decimate_xticklabels, axis='x')

    # if ytick_labels:
    #   plt.ylabel('log10(concentration) (ng/l)')
    # else:
    #   plt.ylabel()
    plt.tick_params(axis='x', which='minor', bottom=False)
    if ylabel:
      plt.ylabel('log10(concentration) (log10[ng/l])')
    plt.ylim(ylim)
    yticks = ax.get_yticks()
    if ytick_labels:
      plt.yticks(yticks, [int(i) for i in yticks])
    else:
      plt.yticks(yticks, [''] * len(yticks))
    plt.yticks(rotation=90, va='center')
    # plt.tick_params(axis='y', left=False, right=True, labelleft=False, labelright=True)
    ax.yaxis.set_label_position("left")
    ax.yaxis.tick_right()
    return ax
  def plot_pie(self, ax=None, col='date', cmap='Paired', **kwargs):
    if ax is None:
      fig, ax = plt.subplots()    
    df = self.df
    if col == 'date':
      df['date'] = pd.to_datetime(df['date'])  # Convert to datetime
      # Exclude time and calculate date counts
      counts = df['date'].dt.date.value_counts(normalize=False)
    else:
      counts = df[col].value_counts(normalize=False)
      # raise NotImplementedError('for col: %s' % col)
    total = counts.sum()
    colors = cm.get_cmap(cmap).colors
    return ax.pie(counts, colors=colors, labels=counts.index, \
      autopct=lambda p: '{:.0f}'.format(p * total / 100))
  def plot_zoom(self, ax=None):
    if ax is None:
      fig, ax = plt.subplots()    
    return ax    
  def plot_xy(self, col='date', mode='hist', ax=None, width=.05, height=.05, **kwargs):
    if ax is None:
      fig, ax = plt.subplots()
    if mode == 'hist':
      plot_fun = self.plot_hist
    elif mode == 'pie':
      plot_fun = self.plot_pie
    elif mode == 'zoom':
      plot_fun = self.plot_zoom
    else:
      raise ValueError('Unknown mode: %s' % mode)
    # The subtracted terms are essential for centering the plots
    lef = convert_to_ax_coords(self.xm, 'x', ax) - width/2
    bot = convert_to_ax_coords(self.ym, 'y', ax) - height/2
    inset_ax = ax.inset_axes([lef, bot, width, height], transform=ax.transAxes)
    return plot_fun(inset_ax, col=col, **kwargs)
  def create_id(self, mode='xy', **kwargs):
    if mode == 'xy':
      cluster_id = self._create_id_based_on_mean_xy(**kwargs)
    elif mode == 'from_clustering':
      cluster_id = list(set(self.df['cluster']))
      assert len(cluster_id) == 1
      cluster_id = cluster_id[0]
    else:
      raise ValueError('Unknown mode: %s' % mode)
    self.__log.debug('Created id %s' % cluster_id)
    return cluster_id
  def _create_id_based_on_mean_xy(self, acc=0.1):
    xm = int(round2accuracy(self.xm, acc))
    ym = int(round2accuracy(self.ym, acc))
    cluster_id = 'x%sy%s' % (xm, ym)
    return cluster_id

@logged
class Drainage(Geospatial):
  def __init__(self, domain):
    self.domain = domain
    self.extract(self.domain.extent)
  def extract(self, extent, fname_inp='../../alex/thames-sewage/input_dir/thames_d8.nc', fname_out=None,\
    pad_width=10, pad_value=0):
    self.__log.info(f'Extracting from {fname_inp}')
    if fname_out is None:
      root, extension = os.path.splitext(fname_inp)
      fname_out = '%s_%s_%s_%s_%s%s' % (root, *extent, extension)
    
    self.fname = fname_out

    d8 = extract_d8_subdomain(*extent, fname_inp, fname_out, pad_width=pad_width, pad_value=pad_value)
    
    self.extent = [float(i) for i in [d8.x.min(), d8.x.max(), d8.y.min(), d8.y.max()]]  
    self.mg = ac.toolkit.load_d8(fname_out)
    return self.mg
  def plot(self, ax=None, cmap='Greys', c='turquoise', cbar=True,\
    min_log10area_to_show=6, alpha_max=0.8, norm='log', **kwargs):
    ax = get_ax(ax)   
    grid = self.mg
    A = grid.at_node["drainage_area"].reshape(grid.shape)
    extent = self.get_extent()
    alpha = np.where(np.log10(A)>min_log10area_to_show, alpha_max, 0)
    
    if c is not None:
      # Create a custom colormap with a single color
      cmap = plt.cm.colors.ListedColormap([c]) 
    else:
      kwargs['norm'] = LogNorm(vmin=np.min(A), vmax=np.max(A))
    im = ax.imshow(A, extent=extent, cmap=cmap, alpha=alpha, **kwargs)
    
    if cbar and c is None:
      colorbar_flexible(im, ax, "Upstream area (m sq.)")
      
    AxesFormatter()._set_xylabels()
    return ax
  def get_extent(self):
    x1, y1 = self.mg.xy_of_lower_left
    dx, dy = self.mg.dx, self.mg.dy
    ny, nx = self.mg.shape # note the order is reversed in mg
    x2 = x1 + nx * dx
    y2 = y1 + ny * dy
    self.extent = [x1, x2, y1, y2]
    return self.extent

@logged
class Experiment(Geospatial):
  def __init__(self, experiment_id='thames-icl'):
    if experiment_id == 'thames-icl':
      self.dataset = DatasetICL()
      self.topography = Topography().extract(self.dataset.extent)  
  def plot(self, ax=None):
    ax = get_ax(ax)
    ax = self.topography.plot(ax=ax)
    ax = self.dataset.plot(ax=ax, norm=None)
    return ax
@logged
class FlowRate(Geospatial):
  def __init__(self, data_samples, means, drainage, river='Wandle'):
    self.data_samples = data_samples
    self.drainage = drainage
    di = Gauges().get_flow_rates(data_samples, means, drainage, river='Wandle')
    mean, std = {}, {}
    for k, v in di.items():
      mean[k] = v[0]
      std[k] = v[1]
    self.mean = ModelUnmix(means.graph.graph, {}, mean, drainage.extent)
    self.std = ModelUnmix(means.graph.graph, {}, std, drainage.extent)
  def get_data(self, date_start='2020-10-18', date_end='2020-11-18', \
    station_id=39003):
    gaug = Gauges()
    _ = gaug.read_stations_data([station_id])
    self.data = gaug.get_station_data_between_dates(station_id, \
      date_start, date_end)
    return self.data
  def plot(self, ax=None, date_start='2020-10-18', \
    date_end='2020-11-18', station_id=39003, c='k', **kwargs):
    ax = get_ax(ax)
    data = self.get_data(date_start, date_end, station_id)
    ax.plot(data.keys(), data.values(), '-o', c=c, **kwargs)
    ax.set_ylabel('Flow rate (m$^3$/s)')
    return ax
@logged
class Flux(Geospatial):
  def calculate(self, flowrate, flowratestd, mfinal, mfinalstd, means, \
    drainage_extent, unit='kg/day'):
    self.unit = unit
    self.__log.warning('Unit: %s' % unit)
    ng_per_l_to_kg_per_m3 = 1e-6
    if unit == 'kg/year':
      factor = 3600 * 24 * 365 * ng_per_l_to_kg_per_m3
    elif unit == 'kg/month':
      factor =  3600 * 24 * 30 * ng_per_l_to_kg_per_m3
    elif unit == 'kg/day':
      factor =  3600 * 24 * ng_per_l_to_kg_per_m3
    else:
      raise ValueError('Uknown unit: %s' % unit)
    mean, sigm = {}, {}
    for k, _ in flowrate.dict.items():
      f_mean = flowrate.dict[k]
      f_sigm = flowratestd.dict[k]
      c_mean = mfinal.dict[k] 
      c_sigm = mfinalstd.dict[k]
      mean[k] = f_mean * c_mean * factor
      sigm[k] = f_mean * c_sigm * factor + c_mean * f_sigm * factor
    self.mean = ModelUnmix(means.graph.graph, {}, mean, drainage_extent)
    self.std = ModelUnmix(means.graph.graph, {}, sigm, drainage_extent)   
    return self.mean, self.std

@logged
class Gauges(Geospatial):
  def get_catchment_area(self, station_id):
    # in km2
    df = self.df
    return float(df[df.id==station_id]['catchment-area'].iloc[0])
  def get_flow_rates(self, data_samples, means, drainage, river='Wandle'):
    self.read_metadata(river=river)
    sampling_dates = data_samples.get_sampling_dates()
    sample_areas = means.get_catchment_areas(drainage)
    if river == 'Wandle':
      station_id = 39003
    else:
      raise NotImplementedError(f'River {river}')
    df = self.read_stations_data([station_id])
    self.__log.info('Gauge catchment area: %s km2' % self.get_catchment_area(station_id))
    flow_rates = self.get_mean_and_std_flow_rates_from_basins(sample_areas, station_id, sampling_dates)
    return flow_rates 
  def get_mean_and_std_flow_rate(self, station_id, sampling_dates):
    di = self.get_station_data_from_sampling_dates(station_id, sampling_dates)
    vals = np.array(list(di.values()))
    mean, sigm = np.mean(vals), np.std(vals)
    return mean, sigm
  def get_mean_and_std_flow_rates_from_basins(self, basin_areas, \
    station_id, sampling_dates):
    di = {}
    for k, v in basin_areas.items():
      di[k] = self.get_mean_and_std_flow_rate_from_basin(v, station_id, sampling_dates)
    return sort(di)
  def get_mean_and_std_flow_rate_from_basin(self, basin_area, \
    station_id, sampling_dates):
    r_mean, r_sigm = self.get_mean_and_std_runoff(station_id, sampling_dates)
    return r_mean * basin_area, r_sigm * basin_area
  def get_mean_and_std_runoff(self, station_id, sampling_dates):
    area = self.get_catchment_area(station_id)
    flow_mean, flow_sigm = self.get_mean_and_std_flow_rate(station_id, sampling_dates)
    return flow_mean / area, flow_sigm / area
  def get_station_data_all_dates(self, station_id):
    df = self.station_data
    data = df[df.station_id == station_id]['time_series'].iloc[0]
    dates = [datetime.strptime(date, '%Y-%m-%d') for date in data[::2]]
    values = data[1::2]
    di = {k: v for (k,v) in zip(dates, values)}
    return di
  def get_station_data_from_sampling_dates(self, station_id, sampling_dates):
    di = self.get_station_data_all_dates(station_id)
    return {k: di[k] for k in sampling_dates}
  def get_station_data_between_dates(self, station_id, date_start, date_end):
    di = self.get_station_data_all_dates(station_id)

    date_start = datetime.strptime(date_start, '%Y-%m-%d')
    date_end = datetime.strptime(date_end, '%Y-%m-%d')
    delta = timedelta(days=1)
    
    newdi = {}
    current_date = date_start
    while current_date <= date_end:
      newdi[current_date] = di[current_date]
      current_date += delta  # Move to the next day

    return newdi
  def plot(self, ax=None, marker='>', c='tab:blue', annotate=False,\
    edgecolor='k', **kwargs):
    if ax is None:
      fig, ax = plt.subplots(figsize=[12,8])
    ax.scatter(self.df[self.xcol], self.df[self.ycol], marker=marker,\
      c=c,
      edgecolor=edgecolor, **kwargs)
    if annotate:
      self._annotate(c=c, **kwargs)
    return ax
  def read_metadata(self, fname='../../data/nrfa-station-metadata-2023-07-11.csv',\
     river='wandle'):
    """
    Read previously downloaded metadata of all NRFA gauging stations.

    Parameters
    ----------
    fname : str, optional
        _description_, by default '../../data/nrfa-station-metadata-2023-07-11.csv'

    Returns
    -------
    _type_
        _description_
    """
    self.xcol = 'easting'
    self.ycol = 'northing'
    gauges = pd.read_csv(fname)
    gauges = gauges[pd.isna(gauges['closed'])] # exclude closed ones
    river = 'Wandle' if river == 'wandle' else 'Wandle'
    self.df = gauges.loc[gauges['name'].str.contains(river)]
    return self.df
  def read_stations_data(self, ids, \
    base_url="https://nrfaapps.ceh.ac.uk/nrfa/ws"):
    station_metadata = []
    # Loop through the station IDs and fetch metadata
    for sid in ids:
        query = "station=%s&data-type=gdf&format=json-object" % sid
        stations_info_url = "{BASE}/time-series?{QUERY}".format(BASE=base_url,
                                                                QUERY=query)  
        # Send request and read response
        response = urllib.request.urlopen(stations_info_url).read()

        # Decode from JSON to Python dictionary
        response = json.loads(response)
        
        # Extract relevant metadata and time-series data
        station_info = response['station']
        data_type_info = response['data-type']
        station_metadata.append({
            "station_id": station_info['id'],
            "name": station_info['name'],
            "easting": station_info.get('easting', None),
            "northing": station_info.get('northing', None),
            "latitude": station_info.get('latitude', None),
            "longitude": station_info.get('longitude', None),
            "data_type_id": data_type_info['id'],
            "data_type_name": data_type_info['name'],
            "data_type_parameter": data_type_info['parameter'],
            "data_type_units": data_type_info['units'],
            "data_type_measurement_type": data_type_info['measurement-type'],
            "data_type_period": data_type_info['period'],
            "time_series": response['data-stream'],
        })

    self.station_data = pd.DataFrame(station_metadata)
    return self.station_data
  def _annotate(self, c, column='name', offset=None, **kwargs):
    if offset is None:
      offset = 500
      offset = [2*offset, offset]
    for i, row in self.df.iterrows():
      xy = np.array([row[self.xcol], row[self.ycol]])
      xytext = xy + offset
      plt.annotate(row[column], c=c, xy=xy, xytext=xytext, \
        # arrowprops=dict(facecolor='white', edgecolor='black', arrowstyle='->'), \
        bbox=dict(boxstyle='round,pad=0.3', edgecolor='none', facecolor='white', alpha=.5), \
        fontsize=12,  ha='center')

@logged
class Model(Geospatial):
  # Generic
  def create(self, mode='homog', *args, **kwargs):
    if mode == 'homog':
      out = self.create_homog(*args, **kwargs)
    elif mode == 'random':
      out = self.create_random(*args, **kwargs)      
    else:
      raise ValueError('Unknown mode: %s' % mode)
    return out
class ModelDrainageSingleParam(Model):
  def __init__(self, ids):
    self.ids = ids
  def create_homog(self, value) -> dict:
    self.dict = {k: value for k in self.ids}
    
    return self.dict
  def create_random(self, vmin=1, vmax=1000) -> dict:
    pass    
class ModelAreas(ModelDrainageSingleParam):
  pass
class ModelConcs(ModelDrainageSingleParam):
  def create_random(self, vmin=1, vmax=1000):
    """
    Create a mdoel uniformly distributed in the log-space
    between two values.

    Parameters
    ----------
    vmin : int, optional
        Before taking log, by default 1
    vmax : int, optional
        Before taking log, by default 1000

    Returns
    -------
    dict
        Mapping self.ids -> values
    """
    self.dict = {k: 10**random.uniform(np.log10(vmin), np.log10(vmax)) for k in self.ids}
    return self.dict
class ModelExRates(ModelDrainageSingleParam):
  def create_homog(self, value=1):
    return super().create_homog(value)
class ModelDrainageMultiParam: # TODO: add areas?
  def __init__(self, concs: ModelConcs, exrates: ModelExRates):
    self.concs = concs
    self.exrates = exrates
@logged
class ModelUnmix(Model):
  def __init__(self, network, exrates: dict, concentrations: dict, \
    drainage_extent):
    self.exrates = exrates
    self.dict = sort(concentrations)
    self.extent = drainage_extent
    self.catch = Catchment(network, drainage_extent)    
    self.subcatchments = self.catch.subcatchments
    self.ids = sorted(self.subcatchments.keys())
    self.upst_conc_map = snu.get_upstream_concentration_map(self.subcatchments, \
      self.dict)    
  def create_homog(self, value) -> dict:
    self.dict = {k: value for k in self.ids}
    self.upst_conc_map = snu.get_upstream_concentration_map(self.subcatchments, \
      self.dict) 
    return self.dict    
  def create_spike(self, id_peak, val_peak, val_bckg) -> dict:
    self.dict = self.create_homog(val_bckg)
    self.dict[id_peak] = val_peak
    self.upst_conc_map = snu.get_upstream_concentration_map(self.subcatchments, \
      self.dict)    
    return self.dict      
  def plot_new(self, cmap='viridis', ax=None, norm='log', catch=True,\
     vmin=None, \
    vmax=None, cbar=True, cbar_label="Concentration (ng/l)", extent=None, \
    detection_threshold=None, detection_threshold_c='w', **kwargs):  
    pass
  def plot(self, cmap='viridis', ax=None, norm='log', catch=True, vmin=None, \
    vmax=None, cbar=True, cbar_label="Concentration (ng/l)", extent=None, \
    detection_threshold=None, detection_threshold_c='w', **kwargs):
    ax = get_ax(ax)
    upstream_map = self.upst_conc_map
    cmap = plt.get_cmap(cmap)
    A = np.where(upstream_map==0, np.nan, upstream_map) # replace 0s with nans, so they're not plotted
    kwargs['cmap'] = cmap
    # vmin = kwargs.get('vmin', np.nanmin(A))
    # vmax = kwargs.get('vmax', np.nanmax(A))
    if norm == 'log':
      # kwargs['norm'] = LogNorm(vmin=np.nanmin(A), vmax=np.nanmax(A))
      kwargs['norm'] = LogNorm(vmin=vmin, vmax=vmax)
    else:
      kwargs['vmin'] = vmin
      kwargs['vmax'] = vmax

    im = ax.imshow(A, extent=self.extent, **kwargs)
    if cbar:
      cb = plt.colorbar(im, ax=ax)
      cb.set_label(cbar_label)
      cb.ax.minorticks_off()
      if detection_threshold is not None:
        cb.ax.axhline(detection_threshold, color=detection_threshold_c, linestyle='-', linewidth=2)


    if catch:
      ax = self.catch.plot(ax=ax)
    
    # axins = ax.inset_axes([0.7, 0.8, .3, .2])
    # axins.imshow()
    # fig.colorbar(means.sc, ax=ax, cax=cax, location='bottom', orientation='horizontal',
    #             label='Concentration (ng/l)')
    AxesFormatter()._convert_ticks_m2km(ax=ax)
    # remove_every_nth_tick_label(n=2, axis='x', ax=ax)
    yt_kw = dict(rotation=90, va='center')
    plt.yticks(**yt_kw)
    # plt.ylim(151500, 177000)
    ax.set_aspect('equal')    
    if extent is not None:
      ax.set_xlim(*extent[:2])
      ax.set_ylim(*extent[2:])
    return ax
  def plot_dict(self, ax=None, vert=True, c='k'):
    ax = get_ax(ax)
    d = sort(self.dict)
    x = np.log10(np.array(list(d.values())))
    y = range(1, len(x) + 1)
    if not vert:
      x, y = y, x
    ax.plot(x, y, c, linestyle='--')
    return ax
  def plot_inset(self, ax, extent, ticks=False, labels=False, \
    width=0.4, add_to_y=0.03, **kwargs):
    x1, x2, y1, y2 = extent
    aspect = (y2 - y1) / (x2 - x1)
    height = width * aspect 
    x0 = 1 - width
    y0 = 1 - height + add_to_y
    axins = ax.inset_axes([x0, y0, width, height])
    self.plot(ax=axins, cbar=0, **kwargs)
    axins.set_xlim(x1,x2)
    axins.set_ylim(y1,y2)
    if not ticks:
      axins.set_xticks([])
      axins.set_yticks([])
    if not labels:
      axins.set_xlabel(None)
      axins.set_ylabel(None)
    ax.indicate_inset_zoom(axins, edgecolor="black")
    return axins
  def read(self, fname: str):
    self.dict = read(fname)
    self.upst_conc_map = snu.get_upstream_concentration_map(self.subcatchments, self.dict) 
    return self.dict
  @classmethod
  def read_for_chems(cls, path, graph, drain, chems):
    m, mstd = {}, {}
    for chem in chems:
      m[chem] = cls(graph, {}, {}, drain.extent).read(f'{path}/{chem}-ModelFinalMean.json')
      mstd[chem] = cls(graph, {}, {}, drain.extent).read(f'{path}/{chem}-ModelFinalStd.json')
    return m, mstd
  def save(self, fname: str):
    pass
@logged
class ModelFinal(ModelUnmix):
  def read(self, path, chem, io='riverpy'):
    if io == 'riverpy':
      suffix = '-ModelFinalMean.json'
    else:
      raise ValueError(f'io:{io}')
    
    fname = f'{path}/{chem}{suffix}'
    return super().read(fname)
@logged
class ModelStd(ModelUnmix):
  def read(self, path, chem, io='riverpy'):
    if io == 'riverpy':
      suffix = '-ModelFinalStd.json'
    else:
      raise ValueError(f'io:{io}')
    
    fname = f'{path}/{chem}{suffix}'
    return super().read(fname)

@logged
class Rainfall:
  pass
@logged
class RainfallNW3:
  def href2date(self, href):
    year, rest = href.split('/wxhistday.php?year=')[1].split('&month=')
    month, day = rest.split('&day=')
    date = f'{year}-{int(month):02d}-{int(day):02d}'
    return date
  def get_all_data(self, years=np.arange(2019,2024)):
    data = {}
    for year in years:
      di = self.scrape(year)
      data = dict(data, **di)
    self.data = data
    return self.data
  def get_data(self, date_start, date_end):
    # assert date_start[:2] == date_end[:2] # make sure it is the same year, as we scrape one year at a time
    date_start = datetime.strptime(date_start, '%Y-%m-%d')
    date_end = datetime.strptime(date_end, '%Y-%m-%d')
    delta = timedelta(days=1)
    data = {}
    current_date = date_start
    while current_date <= date_end:
      # date = current_date.strftime("%Y-%m-%d")
      data[current_date] = self.data[current_date.strftime("%Y-%m-%d")]
      current_date += delta  # Move to the next day
    return data
  def plot(self, date_start, date_end, ax=None, c='k', **kwargs):
    ax = get_ax(ax)
    di = self.get_data(date_start, date_end)
    ax.plot(di.keys(), di.values(), '-o', c=c, **kwargs)
    ax.set_ylabel('Rainfall (mm)')
    return ax
  def scrape(self, year):
    url = f"http://nw3weather.co.uk/wxdataday.php?year={year}&vartype=rain"
    self.__log.info(f'Connecting to {url}')
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')
    table = soup.find('table', class_='table1')
    mm = {}
    for row in table.find_all('tr')[1:]:
      for cell in row.find_all('td'):
        anchor = cell.find('a')
        if anchor:
          href = anchor.get('href')
          date = self.href2date(href)
          mm[date] = float(cell.text.strip())
    return sort(mm)
@logged
class RainfallEA(Geospatial):
  def read(self):
    response = requests.get("http://environment.data.gov.uk/flood-monitoring/id/stations")
    data = response.json()
    stations = data['items']
    ids, labels, lats, lons = [], [], [], []
    for station in stations:
      if type(station.get('lat', np.nan)) != float:
        self.__log.debug('Skipping station with >1 lat/lon: ', station)
        continue
      ids.append(station['@id'].split('/')[-1])
      labels.append(station.get('label', 'N/A'))
      lats.append(station.get('lat', np.nan))
      lons.append(station.get('long', np.nan))
    df = pd.DataFrame({'id': ids, 'location': labels, 'lat': lats, 'lon': lons})
    df = df.dropna()
    transformer = create_coord_transform('latlon', 'bng')
    df['x'], df['y'] = transformer.transform(df['lon'].values, df['lat'].values) 
    return df
  def fetch_rainfall_data(self, station_id, start_date, end_date):
    base_url = "http://environment.data.gov.uk/flood-monitoring/id/stations"
    url = f"{base_url}/{station_id}/measures?parameter=rainfall&startdate={start_date}&enddate={end_date}"
    response = requests.get(url)
    data = response.json()
    return data
@logged
class Rivers(Geospatial):
  def plot(self, ax=None, c='turquoise', lw=1, zorder=1):
    if ax is None:
      fig, ax = plt.subplots(figsize=[5,10])         
    return self.gdf.plot(ax=ax, color=c, lw=lw, zorder=zorder)
@logged
class RiversWfd(Rivers):
  def read(self, *args, **kwargs):
    gdf = read('../../data/WFD_London.shp')
    gdf['WB_NAME'] = gdf['WB_NAME'].str.replace('Lea', 'Lee')
    self.gdf = gdf
    return self.gdf
@logged
class RiversOsd(Rivers):
  def read(self, river='wandle', extent=None):
    if river == 'wandle':
      self.gdf = read('../data/rivers_osd_521500_540000_151500_179500.shp')
    else:
      from shapely.geometry import box
      rivers = read('../data/rivers_osd.shp')
      x1, x2, y1, y2 = extent
      bbox = box(x1, y1, x2, y2) # sic!
      self.gdf = rivers[rivers.intersects(bbox)]
    return self.gdf
@logged
class RiversAll(Rivers):
  def __init__(self, *args, **kwargs):
    self.wfd = RiversWfd()
    _ = self.wfd.read(*args, **kwargs)
    self.osd = RiversOsd()
    _ = self.osd.read(*args, **kwargs)
    self.all = [self.wfd, self.osd]
  def plot(self, ax=None, **kwargs):
    ax = get_ax(ax)
    for rivers in self.all:
      rivers.plot(ax=ax, **kwargs)
    return ax

@logged
class SampleGraph(Geospatial):
  def __init__(self, graph):
    """
    Consider adding a method 'create' based on `create_sample_graph` 
    and `_export_input_for_get_sample_graph` of `DataSamples`.

    Parameters
    ----------
    graph : networkx.classes.digraph.DiGraph
        _description_
    """
    self.graph = graph
    self.count_nodes()
  def count_nodes(self):
    for i, (sample_name, my_data) in enumerate(self.graph.nodes(data=True)):
      pass
    n = i + 1
    self.no_nodes = n
    return self.no_nodes
  def get_topology(self) -> dict:
    di = {}
    for i, (id, dat) in enumerate(self.graph.nodes(data=True)):
      di[id] = dat['data'].downstream_node
    return di
  def get_topology_arrow_coords(self):
    di = self.get_topology()
    coords = {}
    for ida, idb in di.items():
      a = self.graph.nodes[ida]['data']
      if idb == '##ROOT##':
        coords[ida] = 'ROOT'
      else:
        b = self.graph.nodes[idb]['data']
        coords[ida] = (a.x, a.y), (b.x, b.y)
    self.arrow_coords = coords
    return self.arrow_coords
  def plot_spatial(self, ax=None, width=100, c='w', edgecolor='k', **kwargs):
    ax = get_ax(ax)
    if not hasattr(self, 'arrow_coords'):
      self.get_topology_arrow_coords()
    for start_node, coords in self.arrow_coords.items():
      if coords == 'ROOT':
        continue
      a, b = coords
      ax = plot_arrow_from_a_to_b(a, b, ax=ax, width=width, c=c, edgecolor=edgecolor, **kwargs)
    AxesFormatter()._set_xylabels()  
    return ax
  def plot(self, ax=None):
    if ax is None:
      fig, ax = plt.subplots(figsize=[5,10])      
    snu.plot_network(self.graph) #TODO transpose?
@logged
class Topography(Geospatial):
  def __init__(self, domain):
    self.domain = domain
    self.extract(self.domain.extent)
  def extract(self, extent, fname_inp='../../data/LIDAR_Composite_10m_DTM_2022.tif', fname_out=None):
    self.__log.info(f'Extracting from {fname_inp}')
    xmin, xmax, ymin, ymax = extent
    self.extent = extent
    if fname_out is None:
      root, extension = os.path.splitext(fname_inp)
      fname_out = '%s_%s_%s_%s_%s%s' % (root, *extent, extension)    
    if os.path.exists(fname_out):
      os.remove(fname_out)
    cmd = ['gdalwarp', '-te', 
           str(xmin), str(ymin), str(xmax), str(ymax), fname_inp, fname_out]
    completed_process = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = completed_process.stdout
    stderr = completed_process.stderr
    self.__log.debug(stdout)
    if len(stderr) > 0:
      self.__log.error(stderr)
    # topo = Topography()
    # topo.extent = extent
    # topo.domain = Domain(*extent)
    topo = self.read(fname_out)
    return topo
  def extract_around(self, x, y, pad):
    extent = [x-pad, x+pad, y-pad, y+pad]
    return Topography(Domain(*extent))
  def read(self, fname):
    self.xarr = read(fname)
    self.xarr = self.xarr.assign_coords(xkm=(self.xarr['x'] / 1000.0))
    self.xarr = self.xarr.assign_coords(ykm=(self.xarr['y'] / 1000.0))    
    self.__log.debug('Set self.xarr attribute')
    return self.xarr
  def plot(self, ax=None, cmap='Greys_r', aspect='equal', azdeg=45, altdeg=45, \
    unit='m',\
    xlabel='Easting (km)', ylabel='Northing (km)', label='Altitude (m a.s.l.)', **kwargs):

    ax = get_ax(ax)
    
    if unit == 'km':
      extent=[self.xarr.xkm.min(), self.xarr.xkm.max(), self.xarr.ykm.min(), self.xarr.ykm.max()]
    else:
      extent=[self.xarr.x.min(), self.xarr.x.max(), self.xarr.y.min(), self.xarr.y.max()]
    A = self.xarr.values.squeeze()

    self.shader = Shade()
    self.shader.plot(A.T, extent=extent, cmap=cmap, aspect=aspect, azdeg=azdeg, altdeg=altdeg,\
      label=label, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # plt.xticks(rotation=45)
    # print(ax.get_xlim(), ax.get_ylim())
    AxesFormatter()._convert_ticks_m2km(ax=ax)
    ax = remove_every_nth_tick_label(n=2, ax=ax, axis='x', start=0, remove_ticks=True)
    ax.set_xlim(self.extent[:2])
    ax.set_ylim(self.extent[2:])
    # ax.invert_yaxis()
    return ax
  def plot_all(self, domain, catch, drainage, wfd, osd, ds, cso, stw,\
    chem, ax=None, cmin=1, cmax=1000):
    if ax is None:
      fig, ax = plt.subplots(figsize=[12,8])
    ax = self.plot(ax=ax, cbar=0)
    if catch is not None:
      ax = catch.plot(ax=ax, c='k', ls='--', lw=1)
    ax = drainage.plot(ax=ax, cbar=0)
    ax = wfd.plot(ax=ax)
    ax = osd.plot(ax=ax)
    if ds is not None:
      ax = ds.plot(ax, column=chem, vmin=cmin, vmax=cmax, s=400, edgecolor='k', \
        cbar=0, annotate='id', annotate_kw={'fontsize':12, 'color': 'w'}, zorder=9)
    cso.plot(ax=ax, zorder=8, cbar=0, marker='X', s=100, c='k', alpha=0.6)
    stw.plot(ax=ax, zorder=7, s=400)
    ax.set_xlim(domain.extent[:2])
    ax.set_ylim(domain.extent[2:])    
    return ax
@logged
class TopographyWandle(Topography):
  def __init__(self):
    self.domain = DomainWandle()
    self.extract(self.domain.extent)
@logged
class WTWs(Geospatial):
  def __init__(self, *args, **kwargs):
    self.read(*args, **kwargs)
  def plot(self, ax=None, s=200, offset=400, zorder=None):
    if ax is None:
      fig, ax = plt.subplots(figsize=[12,8])
        
    df = self.df
    c_wtw = 'red'
    ax.scatter(df.x, df.y, marker='d', s=s, edgecolor='k', facecolor=c_wtw, zorder=zorder)
    xcol, ycol = 'x', 'y'
    zorder = zorder +1 if zorder is not None else zorder
    ax = self._annotate(ax=ax, c='w', df=df, xcol='x', ycol='y',\
      column='WWTP_NAME', offset=offset, zorder=zorder), 
    # column = 'WWTP_NAME'

    # offset = 400
    # for i, row in df.iterrows():
    #   ax.annotate(row[column], (row[xcol]+offset, row[ycol]+offset), c='w', 
    #                 fontsize=12)    
    return ax
  def read(self, extent, catch=None, restrict_to_catch=True):
    x1, x2, y1, y2 = extent
    # load wtws 
    fname = '../../data/HydroWASTE_v10/HydroWASTE_v10_exported.csv'
    df = pd.read_csv(fname, sep=',', encoding='utf-8') 
    df = df.loc[df.COUNTRY=='United Kingdom']
    self.__log.info('There are %s WTWs in the UK' % len(df))
    # bedd is actually more accurate in hydrowaste!
    # df.loc[df['WWTP_NAME'] == 'LONDON Beddington STW', 'LON_WWTP'] = -0.1507
    df.loc[df['WWTP_NAME'] == 'LONDON Mogden STW', 'LAT_WWTP'] = 51.4646354003
    df.loc[df['WWTP_NAME'] == 'LONDON Mogden STW', 'LON_WWTP'] = -0.34031514583
    df.loc[df['WWTP_NAME'] == 'LONDON HOGSMILL Valley STW', 'LAT_WWTP'] = 51.40138889
    df.loc[df['WWTP_NAME'] == 'LONDON HOGSMILL Valley STW', 'LON_WWTP'] = -0.28416666
    # Rename for plots
    df.loc[df['WWTP_NAME'] == 'LONDON Beddington STW', 'WWTP_NAME'] = 'Beddington STW'
    df.loc[df['WWTP_NAME'] == 'LONDON HOGSMILL Valley STW', 'WWTP_NAME'] = 'Hogsmill STW'
    df.loc[df['WWTP_NAME'] == 'LONDON Mogden STW', 'WWTP_NAME'] = 'Mogden STW'
    df.loc[df['WWTP_NAME'] == 'MERSTHAM STW', 'WWTP_NAME'] = 'Merstham STW'

    transformer = create_coord_transform('latlon', 'bng')
    df['x'], df['y'] = transformer.transform(df['LON_WWTP'].values, df['LAT_WWTP'].values) 
    df = df[
        (df['x'] >= x1) &
        (df['x'] <= x2) &
        (df['y'] >= y1) &
        (df['y'] <= y2)
    ]
    self.df = df
    if restrict_to_catch:
      assert catch is not None
      self.df = self._restrict_to_catchment(self.df, catch, xcol='x', ycol='y')
    return self.df
  def _annotate(self, ax, xcol, ycol, c, df, column, offset=400, \
    shorten=2, zorder=None):
    ax = get_ax(ax)
    for i, row in df.iterrows():
      txt = row[column][:shorten] if shorten is not None else row[column]
      ax.annotate(txt, (row[xcol]+offset, row[ycol]+offset), c=c, 
                    fontsize=12, zorder=zorder)  
    return ax  
