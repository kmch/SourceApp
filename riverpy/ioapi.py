from iogeo.api import read, save

class File:
  def init(self, core, path, extn, prefix='', suffix=''):
    self.core = core
    self.path = path
    self.extn = extn
    self.name = f'{prefix}{self.core}{suffix}.{self.extn}' 
    self.fname = f'{self.path}/{self.name}'
  def read(self, **kwargs):
    return read(self.fname, **kwargs)
  def save(self, data, **kwargs):
    return save(self.fname, data, **kwargs)

class FileDataSynFromModelFinalMean(File):
  def __init__(self, *args, **kwargs):
    core = 'DataSynFromModelFinalMean'
    self.init(core, *args, **kwargs)

class FileModelFinalMean(File):
  def __init__(self, *args, **kwargs):
    core = 'ModelFinalMean'
    self.init(core, *args, **kwargs)

class FileModelFinalStd(File):
  def __init__(self, *args, **kwargs):
    core = 'ModelFinalStd'
    self.init(core, *args, **kwargs)
