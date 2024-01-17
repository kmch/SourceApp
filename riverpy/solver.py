from autologging import logged, traced
import sample_network_unmix as snu

class SolverForward:
  pass
class SolverForwardFunmixer(SolverForward):
  def __init__(self, sample_network, areas, upst_conc_map, exrates):
    self.sample_network = sample_network
    self.areas = areas
    self.upst_conc_map = upst_conc_map
    self.exrates = exrates  
  def run(self):
    dsyn_start, _ = snu.mix_downstream(self.sample_network, self.areas, self.upst_conc_map, self.exrates)
    return dsyn_start
class SolverForwardRiverPy(SolverForward):
  def __init__(self):
    pass    

class SolverInverse:
  pass
class SolverInverseFunmixer(SolverInverse): # TODO: IT'S NOT USED ACTUALLY
  def __init__(self, sample_network, areas, upst_conc_map, exrates):
    self.sample_network = sample_network
    self.areas = areas
    self.upst_conc_map = upst_conc_map
    self.exrates = exrates    
  def run(self, sample_network, element_data, solver='ecos'):
    problem = snu.SampleNetworkUnmixer(sample_network=sample_network, 
                                   continuous=False, use_regularization=True)
    regularizer_strength = self._get_optimal_regulariser_strength()
    element_pred_down, element_pred_upstream = problem.solve(element_data, \
      solver=solver, regularization_strength=regularizer_strength)
    
    area_dict = snu.get_unique_upstream_areas(sample_network)
    upstream_map = snu.get_upstream_concentration_map(area_dict, element_pred_upstream)     
    return upstream_map 
  def _get_optimal_regulariser_strength(self): #TODO
    snu.plot_sweep_of_regularizer_strength(problem, element_data, -5, -1, 11)
    regularizer_strength = None #FIXME
    return regularizer_strength    

