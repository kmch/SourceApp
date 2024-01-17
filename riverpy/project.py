from autologging import logged, traced
from nsh.generic import Dir
from nsh.utils import add_to_path
from riverpy.frame import *
from riverpy.input import *
from riverpy.solver import *

@logged
class Project:
  def __init__(self, name, path='./', problem='forward', discrete=True, \
    solver='funmix', ioapi='default', **kwargs):
    # the order matters
    self.name = name
    self._init_path(path)
    self._create_directory()
    self._init_kwargs(**kwargs)
    self._init_problem(problem, discrete)
    self._init_ioapi(ioapi)
    self._init_input()
    self._init_solver(solver)
    self._init_output()
  def _init_input(self):
    self.inp = Throughput(proj=self, inp_or_out='inp')
    self.problem._init_input(self.inp)
  def _init_ioapi(self, ioapi):
    pass
  def _init_kwargs(self, **kwargs):
    # Save the kwargs in case we want duplicate the project
    # Note, dict prevents deleting the path
    self.kws = dict(kwargs) 
  def _init_path(self, path='./'):
    self.path = add_to_path(path, self.name)
    self.__log.info(f'Project path set to {self.path}')
  def _init_problem(self, problem, discrete):
    if problem == 'forward':
      if discrete:
        problem = ProblemForwardDiscrete()
      else:
        problem = ProblemForwardContinuous()
    elif problem == 'inverse':
      if discrete:
        problem = ProblemInverseDiscrete()
      else:
        raise NotImplementedError(problem)
    else:
      raise ValueError(f'Problem: {problem}')
    self.problem = problem
  def _init_output(self):
    self.out = Throughput(proj=self, inp_or_out='out')
  def _init_solver(self, solver):    
    self.problem._init_solver(solver)
  def _create_directory(self):
    self.dir = Dir(self.path)
    self.dir.create()

@logged
class Throughput:
  def __init__(self, proj, inp_or_out:str, **kwargs):
    self.proj = proj
    self._init_path(inp_or_out)
    self._create_directory()
  def _create_directory(self):
    self.dir = Dir(self.path)
    self.dir.create()
  def _init_objects(self):
    pass
  def _init_path(self, inp_or_out:str):
    dirname = '/inp' if inp_or_out == 'inp' else '/out'
    self.path = add_to_path(self.proj.path, dirname)

class Problem:
  pass


class ProblemForward(Problem):
  def init_input(self):
    self.acq = 1
    self.mtrue = 1
    self.solver = 1
  def init_output(self):
    self.dsyn = 1
class ProblemForwardDiscrete(ProblemForward):
  def _init_input(self, inp):
    inp.acq = Acq()
    inp.mtrue = Model()
  def _init_solver(self, solver):
    self.solver = Solver()    
  def init_output(self, out):
    out.dsyn = Data()
class ProblemForwardContinuous(ProblemForward):
  pass

class ProblemInverse(Problem):
  def init_input(self):
    self.acq = 1
    self.dobs = 1
    self.mstart = 1
    self.solver = 1
  def init_output(self):
    self.mfinal = 1
    self.dsyn = 1
class ProblemInverseDiscrete(ProblemInverse):
  pass

class ProblemOED(Problem):
  """
  Find optimal experiment design, i.e. sample locations
  that subdivide the catchment in an optimal way.

  """
  def init_input(self):
    self.domain = 1
    self.most_downstream_point = 1
    self.topography


class InputForwardDiscrete(Throughput):
  pass
class InputForwardDiscrete(Throughput):
  pass
