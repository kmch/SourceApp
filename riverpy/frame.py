"""
Inversion framework. The canonical sequence of steps is the following:
0. Initialise synthetic and inversion projects
1. Prepare input for synthetic calculation.
2. Run the synthetic calculation.
3. Plot the output of the synthetic calculation.
4. Prepare input for inversion.
5. Run the inversion.
6. Plot the output of the inversion.
"""
#FIXME: move to plot!
import cmocean as cm 
import colorcet as cc
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
# --------------------------------------
#TODO from numba import jit # over chemicals or for MC
import numpy as np

from plotea.mpl2d import *
from riverpy.input import *
from riverpy.solver import SolverForwardFunmixer #, SolverInverseFunmixer #FIXME
from riverpy.utils import sort, convert_xy_to_ij
import sample_network_unmix as snu #FIXME
from scipy.ndimage import gaussian_filter #FIXME

@logged
class InputPrep:
  @classmethod
  def run(cls, dataset='ICL', river='wandle', max_dist=100, detection_threshold=3):
    domain, topo, drain = cls.prep_domain_etc(river)
    dobs, dobs_unc, clusters, acq = cls.prep_data(drain, dataset, river, max_dist, detection_threshold)
    return domain, topo, drain, dobs, dobs_unc, clusters, acq
  @classmethod
  def prep_data(cls, drainage, dataset='ICL', river='wandle', max_dist=100, detection_threshold=3):
    dataset = DatasetICL(detection_threshold=detection_threshold) if dataset == 'ICL' else None
    data_samples = dataset.select_subset(river)
    clusters = data_samples.cluster(max_dist=max_dist)
    means = clusters.export_statistic('mean')
    sigms = clusters.export_statistic('std')
    if river == 'wandle':
      thresh = 1e7
      use_letters = True
      relocate_dict = {'F': [None, 170500]}
    kwargs = dict(thresh=thresh, relocate_dict=relocate_dict, use_letters=use_letters)
    means.process(drainage, **kwargs)
    sigms.process(drainage, **kwargs)
    clusters.update_ids(means.df)
    acq = Acquisition(means)
    return means, sigms, clusters, acq
  @classmethod
  def prep_domain_etc(cls, river='wandle'):
    if river == 'wandle':
      domain = DomainWandle()
    elif river == 'hogsmill':
      domain = DomainHogsmill()
    elif river == 'thames':
      domain = DomainThames()
    else:
      raise NotImplementedError('River: %s' % river)
    topo = Topography(domain)
    drain = Drainage(domain)
    return domain, topo, drain

@logged
class Forward:
  def __init__(self, mtrue, acquisition, solver='faster-unmixer'):
    self.mtrue = mtrue
    self.acquisition = acquisition
    if solver == 'faster-unmixer':
      sample_network = acquisition.network
      areas = mtrue.subcatchments
      exrates = mtrue.exrates
      upst_conc_map = mtrue.upst_conc_map

      self.solver = SolverForwardFunmixer(sample_network, areas, upst_conc_map, exrates)
  def plot_data(self, vmin, vmax, ax=None, **kwargs):
    if ax is None:
      fig, ax = plt.subplots()

    ax = self.mtrue.plot( ax=ax, vmin=vmin, vmax=vmax, **kwargs)
    if hasattr(self, 'dsyn'):
      kwargs['cbar'] = False
      ax = self.dsyn.plot(ax=ax, vmin=vmin, vmax=vmax, **kwargs)
    return ax
  def run(self, **kwargs):
    data_dict = self.solver.run(**kwargs)
    self.dsyn = Data(self.acquisition, data_dict)
    return self.dsyn

@logged
class Inversion:
  def __init__(self, acquisition, data, data_uncertainty=None, extent=None, code='faster-unmixer', exrates=None,\
    use_regularization=True, detection_threshold=None, **kwargs):
    self.detection_threshold = detection_threshold
    self.acquisition = acquisition
    self.data = data
    self.data_uncertainty = data_uncertainty
    self.exrates = exrates
    self.extent = extent
    self.code = code
    if code == 'faster-unmixer':
      self.problem = snu.SampleNetworkUnmixer(sample_network=acquisition.network,
                                   continuous=False, use_regularization=use_regularization)
    else:
      raise NotImplementedError('code: %s' %code)
  def get_fluxes(self, flow_rates: dict, unit='kg/year'):
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
    for k, (f_mean, f_sigm) in flow_rates.items():
      c_mean = self.mfinal.dict[k] 
      c_sigm = self.mfinalstd.dict[k]
      mean[k] = f_mean * c_mean * factor
      sigm[k] = f_mean * c_sigm * factor + c_mean * f_sigm * factor
    self.flux_mean = ModelUnmix(self.data.graph.graph, {}, mean, self.extent)
    self.flux_sigm = ModelUnmix(self.data.graph.graph, {}, sigm, self.extent)   
    return self.flux_mean, self.flux_sigm
  def get_misfit_vs_regul(self, chem, min_log=-5, max_log=-1, n=100, plot=True, **kwargs):
    vals = np.logspace(min_log, max_log, n)
    dobs = self._select_data(chem, self.data.df)
    misfits = []
    for val in vals:
      mfinal, dsyn = self.run(chem, 10**(val), **kwargs)
      misfit = self._calc_misfit(dsyn.dict, dobs)
      # Total L1 misfit
      misfit = np.sum(np.abs(np.array(list(misfit.values()))))
      misfits.append(misfit)
    if plot:
      plt.semilogx(vals, misfits, '.-')
    return vals, misfits
  def get_misfit_vs_rough(self, chem, min_log=-5, max_log=-1, n=11, plot=True):
    dobs = self._select_data(chem, self.data.df)
    vals, rough, misf = snu.plot_sweep_of_regularizer_strength(self.problem, \
      dobs, min_log, max_log, n, plot=plot)
    self.lcurve = [vals, rough, misf]
    return self.lcurve
  def plot(self, *args, **kwargs):
    return self.plot_model(*args, **kwargs)
  def plot_flux_arrows(self, flow_rates, angles=None, ax=None, unit='kg/day', \
    multiply_flux1=1e3, multiply_flux2=1e2):
    ax = get_ax(ax)
    fmean, fsigm = self.get_fluxes(flow_rates, unit=unit)
    for bid, flux in fmean.dict.items():
      ds = self.data.dict[bid]
      x, y = ds.x, ds.y
      width = np.log10(flux * multiply_flux1) * multiply_flux2
      # angle = random.uniform(-45, 45)
      angle = 10 if angles is None else angles[bid]
      ax = plot_arrow(x, y, cvalue=flux, ax=ax, angle=angle, width=width, pad=.1)
      # plt.scatter(x, y, s=100, c='r')
    return ax
  def plot_pseudo_lcurve(self, chem, min_log=-5, max_log=-1, n=11):
    dobs = self._select_data(chem, self.data.df)
    return snu.plot_sweep_of_regularizer_strength(self.problem, dobs, min_log, max_log, n)
  def plot_misfit(self, dobs, ax=None, log=True, annotate=True, cbar=True,\
    c='id', **kwargs):
    ax = get_ax(ax)
    dsyn = self.dsyn.dict
    dobs = dobs if type(dobs) == dict else dobs.dict
    
    syn = np.array(list(sort(dsyn).values()))
    obs = np.array(list(sort(dobs).values()))
    if c == 'misfit':
      misf = self._calc_misfit(dsyn, dobs)
      val = np.array(list(sort(misf).values()))
      cmap = 'Greys'
      vmin = 1
      vmax = 2
    elif c == 'id':
      #print(self.data.dict)
      val = np.arange(len(self.data.dict))
      cmap = 'rainbow'
      vmin, vmax = None, None
      cbar = False
    else:
      raise ValueError(f'c: {c}')
    x1 = np.min([min(syn), min(obs)])
    x2 = np.max([max(syn), max(obs)])
    pad = 1
    xlim = kwargs.get('xlim', (x1-pad, x2+pad))
    if log:
      syn, obs = np.log10(syn), np.log10(obs)
      xlim = np.log10(xlim) if xlim[0] > 0 else np.array((xlim[0], np.log10(xlim[1])))
      ylim = xlim
    sc = ax.scatter(syn, obs, c=val, cmap=cmap, edgecolor='k', s=200, \
      vmin=vmin, vmax=vmax)
    x = np.linspace(*xlim, 100)
    y = x
    ax.plot(x, y, 'k--') #, label='1:1 line')
    # if annotate:
    #   offset = .05 if log else 50
    #   for k, v in self.dsyn.dict.items():
        # y = self.dsyn.dict[k]
        # y = np.log10(y) if log else y
        # x = np.log10(v) if log else v
        # ax.annotate(k, (x + 2 * offset, y + offset), c='k', fontsize=12)

    if cbar:
      mode = kwargs.get('mode', 'hitcount')
      if mode == 'hitcount':
        cax = ax.inset_axes([0.4, 0.15, .5, .05])
        plt.colorbar(sc, ax=ax, cax=cax, location='bottom', orientation='horizontal',
                label='max(dsyn/dobs, dobs/dsyn)')
      else:
        plt.colorbar(sc, ax=ax,
                label='max(dsyn/dobs, dobs/dsyn)')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    # ax.yaxis.set_ticks_position('right')
    # ax.yaxis.set_label_position('right')
    plt.yticks(rotation=90, ha='center', va='center')
    ax.set_xlabel('log10(predicted concentration)')
    ax.set_ylabel('log10(observed concentration)')
    ax.set_aspect(kwargs.get('aspect', 'equal'))
    return ax
  def plot_misfit_map(self, vmin, vmax, ax=None, **kwargs):
    if ax is None:
      fig, ax = plt.subplots()
    dobs = self._select_data(self.chem, self.data.df)
    dsyn = self.dsyn.dict
    misf = {}
    for k, syn in dsyn.items():
      syn = np.log10(syn)
      obs = np.log10(dobs[k])
      misf[k] = np.max([obs / syn, syn / obs])
    misf = Data(self.acquisition, misf)
    # ax = self.mfinal.plot(ax=ax, vmin=vmin, vmax=vmax, **kwargs)
    ax = misf.plot(  ax=ax, vmin=vmin, vmax=vmax, **kwargs)
    return ax
  def plot_misfit_vs_rough(self, ax=None, s=200, offsetx=.7, offsety=.02, \
    cbar=True, cmap='plasma', annotate=False):
    ax = get_ax(ax)
    if hasattr(self, 'lcurve'):
      vals, rough, misf = self.lcurve
    else:
      raise ValueError('Run get_misfit_vs_rough first')
    # ax.text(.3, .9, 'Picked λ: %s' %\
    #    round(np.log10(self.regul_strength),3), fontsize=12, transform=ax.transAxes)
    
    # kw_scatt = dict(c=vals, norm=LogNorm(), edgecolor='k', s=s, zorder=2)
    kw_scatt = dict(c=vals, norm=LogNorm(), edgecolor='k', s=s, zorder=2)
    sc = ax.scatter(rough, misf, **kw_scatt, cmap=cmap)
    ax.plot(rough, misf, '--', c='k', zorder=1)
    if cbar:
      cax = ax.inset_axes([0.4, 0.9, .5, .05])
      cmap = plt.get_cmap(cmap)
      cb = plt.colorbar(sc, ax=ax, cax=cax, location='bottom',\
        orientation='horizontal', label='λ') #r'$\lambda$')
      cb.ax.minorticks_off()
    
    # axins = ax.inset_axes([0.5, 0.5, 0.5, 0.5], transform=ax.transAxes)
    # axins.scatter(rough, misf, **kw_scatt)

    # axins.set_xticks([])
    # axins.set_yticks([])
    # x1, x2 = ax.get_xlim()
    # y1, y2 = ax.get_ylim()
    # x2 = x1 + (x2-x1) / 3
    # y2 = y1 + (y2-y1) / 3
    # axins.set_xlim(x1,x2)
    # axins.set_ylim(y1,y2)
    # for i, val in enumerate(vals):
    #   x, y = rough[i], misf[i]
    #   if x > x1 and x < x2 and y > y1 and y < y2:
    #     axins.text(rough[i], misf[i], str(round(np.log10(vals[i]), 2)))    
    # ax.indicate_inset_zoom(axins, edgecolor="black")
    if annotate:
      for i, val in enumerate(vals):
        x, y = rough[i], misf[i]
        ax.text(rough[i]+offsetx, misf[i]+offsety, \
          str(round(np.log10(vals[i]), 2)), zorder=2, fontsize=10)  
    ax.set_xlabel("Model roughness")
    ax.set_ylabel("Data misfit")
    plt.yticks(rotation=90, ha='center', va='center')
    return ax
  def plot_model(self, vmin, vmax, ax=None, **kwargs):
    if ax is None:
      fig, ax = plt.subplots()
    ax = self.mfinal.plot(ax=ax, vmin=vmin, vmax=vmax, **kwargs)
    ax = self.dsyn.plot(  ax=ax, vmin=vmin, vmax=vmax, **kwargs)
    return ax
  def plot_recovery(self, mtrue, vmin, vmax, ax=None, figsize=(12,6), **kwargs):
    fig = plt.figure(figsize=figsize)
    ndigits = 2
    gs = gridspec.GridSpec(1,2, width_ratios=[1,1])
    ax1 = fig.add_subplot(gs[0])
    ylim1 = vmin if vmin is None else np.log10(vmin)
    ylim2 = vmax if vmax is None else np.log10(vmax)
    ax1 = self.plot_true_vs_final(mtrue, ax=ax1, ylim=(ylim1,ylim2), **kwargs)
    ax2 = fig.add_subplot(gs[1])
    ax2 = self.plot_model(vmin, vmax, ax=ax2)
    return ax1, ax2
  def plot_regularisation(self):
    pass
  def plot_true_vs_final(self, mtrue, ax=None, ylim=(-0.1, None)):
    if ax is None:
      fig, ax = plt.subplots()

    d = self.dobs
    # mfc='none'
    ax.plot(d.keys(), np.log10(list(d.values())), '.-', ms=11, label='observed data')
    d = self.dsyn.dict
    ax.plot(d.keys(), np.log10(list(d.values())), '.-', ms=15, label='synthetic data', alpha=.4)
    # d = self.dsyn_calc.dict # just checked that indeed dsyn_calc = dsyn
    # ax.plot(d.keys(), np.log10(list(d.values())), '.-', ms=15, label='synthetic data (calc)', alpha=.4)
    d = mtrue.dict
    ax.plot(d.keys(), np.log10(list(d.values())),'.-', c='g', ms=11, label='true model')
    d = self.mfinal.dict
    ax.plot(d.keys(), np.log10(list(d.values())),'.-', c='r', ms=15, label='final model', alpha=.4)
    ax.set_ylim(ylim)
    ax.set_ylabel('log10(concentration)')
    ax.legend()
    return ax
  def plot_old(self, output='model', ax=None, norm='log', **kwargs):
    if ax is None:
      fig, ax = plt.subplots()

    if output == 'model':
      output = self.final_model
      # log = True
      cmap = kwargs.get('cmap', 'viridis')
    elif output == 'uncertainty':
      output = self.final_model_std
      # log = False
      cmap = kwargs.get('cmap', 'magma')
    else:
      raise ValueError('Output: %s' % output)

    area_dict = snu.get_unique_upstream_areas(self.data.graph.graph)
    upstream_map = snu.get_upstream_concentration_map(area_dict, output)

    cmap = plt.get_cmap(cmap)
    A = np.where(upstream_map==0, np.nan, upstream_map)
    kwargs['cmap'] = cmap
    if norm == 'log':
      kwargs['norm'] = LogNorm(vmin=np.nanmin(A), vmax=np.nanmax(A))

    im = ax.imshow(A, extent=self.model.extent, **kwargs)
    cb = plt.colorbar(im, ax=ax)
    # cb.set_label("log10(C) of %s" % chem)
    # plt.xlim(x1,x2)
    # plt.ylim(y1,y2)
    return ax
  def run(self, chem, regularization_strength, solver='ecos', verbose=False, \
    **kwargs):
    self.regul_strength = regularization_strength
    self.chem = chem
    self.dobs = self._select_data(chem, self.data.df)
    dsyn, mfinal = self.problem.solve(self.dobs, export_rates=self.exrates, \
      regularization_strength=regularization_strength, solver=solver, verbose=verbose)
    self.dsyn = Data(self.acquisition, dsyn)
    self.mfinal = ModelUnmix(self.acquisition.network, self.exrates, mfinal, self.extent)
    # Checked that indeed dsyn_calc = dsyn (ALL GOOD)
    # self.dsyn_calc = Forward(self.mfinal, self.acquisition, solver=self.code).run()
    return self.mfinal, self.dsyn
  def run_mc(self, chem, regularization_strength, solver='ecos', n=1000):
    """
    WARNING
    This function is deprecated. Please use InversionMC().run() instead!

    Parameters
    ----------
    chem : _type_
        _description_
    regularization_strength : _type_
        _description_
    solver : str, optional
        _description_, by default 'ecos'
    n : int, optional
        _description_, by default 1000

    Returns
    -------
    _type_
        _description_
    """
    n = int(n)
    means = self._select_data(chem, self.data.df)
    stds = self._select_data(chem, self.data_uncertainty.df)

    result = {key: [] for key in means.keys()}

    for i in range(n):
      data = {}
      # print(dict(sorted(means.items())))
      # print(dict(sorted(stds.items())))
      for id, val in means.items():
        mean = np.log(means[id]) # NOT log10, as np.random.lognormal assumes natural one
        std = np.log(stds[id]) # try 0.00000000000001 for checks
        # print('%s: %s +/- %s' % (id, mean, std))
        data[id] = np.random.lognormal(mean, std, size=1)[0]
      # print(print(dict(sorted(data.items()))))

      dsyn_final, final_model = self.problem.solve(data, solver=solver,\
        regularization_strength=regularization_strength)
      for id, _ in result.items():
        result[id].append(final_model[id])

    result_mean = {}
    result_std = {}
    for id, all_iters in result.items():
      result_mean[id] = np.mean(np.log10(np.array(all_iters)))
      result_std[id] = np.std(np.log10(np.array(all_iters)))

    self.final_model = result_mean
    self.final_model_std = result_std
    # self.dsyn_final = downstream_means
    # self.dsyn_std = downstream_uncertainties

    return result_mean, result_std
  def run_mc_alex(self, chem, regularization_strength, solver='ecos', relative_error=10, num_repeats=50):
      # % relative error on observations
    element_pred_down_mc, element_pred_up_mc = self.problem.solve_montecarlo(\
      self._select_data(chem), relative_error=relative_error, num_repeats=num_repeats,\
        regularization_strength=regularization_strength, solver=solver)

    downstream_means, downstream_uncertainties = {}, {}
    for sample, values in element_pred_down_mc.items():
        downstream_uncertainties[sample] = np.std(values)
        downstream_means[sample] = np.mean(values)
    upstream_means, upstream_uncertainties = {}, {}
    for sample, values in element_pred_up_mc.items():
        upstream_uncertainties[sample] = np.std(values)
        upstream_means[sample] = np.mean(values)

    self.final_model = upstream_means
    self.final_model_std = upstream_uncertainties
    self.dsyn_final = downstream_means
    self.dsyn_std = downstream_uncertainties

    return self.final_model, self.dsyn_final
  def save(self, fname):
    path, name = os.path.split(fname)
    core, ext = os.path.splitext(name)
    dicts = [self.mfinal.dict, self.dsyn.dict]
    sufxs = ['mfinalmean', 'dsyn']
    if hasattr(self, 'mfinalstd'):
      dicts.append(self.mfinalstd.dict)
      sufxs.append('mfinalstd')
    for suffix, data in zip(sufxs, dicts):
      with open(path + '/' + core + '_' + suffix + ext, 'w') as file:
        json.dump(data, file)
  def _calc_lcurve(self, dobs):
    return
  def _calc_misfit(self, dobs: dict, dsyn: dict) -> dict:
    """
    As in Barnes & Lipp (2023).
    """
    misfit = {}
    for k in dsyn.keys():
      syn = dsyn[k] if dsyn[k] > 0 else self.detection_threshold
      obs = dobs[k] if dobs[k] > 0 else self.detection_threshold
      misfit[k] = np.max([obs / syn, syn / obs])
    return misfit
  def _calc_regul_term(self, dobs, dsyn):
    return
  def _select_data(self, chem, df):
    return df[['id', chem]].set_index('id')[chem].to_dict()
@logged
class InversionSyn(Inversion):
  def plot_dobs(self):
    return self.dobs.plot(**kwargs)
  def plot_dsyn(self):
    return self.dsyn.plot(**kwargs)
  def plot_dsyn_vs_dobs(self):
    pass
  def plot_dsyn_vs_dobs_boxplot(self):
    pass
  def plot_dsyn_vs_dobs_crossplot(self):
    pass
  def plot_dsyn_vs_dobs_diff(self):
    pass
  def plot_dsyn_vs_dobs_side_by_side(self):
    pass  
  def plot_mtrue(self, mtrue):
    pass  
  def plot_mfinal(self, **kwargs):
    return self.mfinal.plot(**kwargs)
  def plot_mfinal_vs_mtrue(self, mtrue, mode='crossplot', **kwargs):
    if mode == 'boxplot':
      out = self.plot_mfinal_vs_mtrue_boxplot(mtrue, **kwargs)
    elif mode == 'crossplot':
      out = self.plot_mfinal_vs_mtrue_crossplot(mtrue, **kwargs)
    else:
      raise ValueError('Unknown mode: %s' % s)
    return out
  def plot_mfinal_vs_mtrue_boxplot(self, mtrue, ax=None):
    pass
  def plot_mfinal_vs_mtrue_crossplot(self, mtrue, ax=None, annotate=True, c='turquoise', **kwargs):
    ax = get_ax(ax)
    mtrue = mtrue.dict
    mfinl = self.mfinal.dict
    ids, x, y = [], [], []
    for basen_id, val in sort(mtrue).items():
      ids.append(basen_id)
      x.append(mfinl[basen_id])
      y.append(mtrue[basen_id])
    if annotate:
      kwargs['annotations'] = ids
    ax = plot_crossplot(x, y, ax=ax, logx=True, logy=True, c=c, **kwargs)
    ax.set_xlabel('log10(recovered model)')
    ax.set_ylabel('log10(true model)') 
    return ax
  def plot_mfinal_vs_mtrue_diff(self, mtrue, ax=None):
    pass
  def plot_mfinal_vs_mtrue_side_by_side(self, mtrue, ax=None):
    pass
@logged
class InversionMC(Inversion):
  def __init__(self, acquisition, data, data_uncertainty, \
    extent=None, code='faster-unmixer',\
    exrates=None, use_regularization=True, **kwargs):
    super().__init__(acquisition, data, data_uncertainty, extent, code, exrates, use_regularization, **kwargs)
  def generate_data_lognormal(self, chem, std_is_relative: bool, n=1000):
    self.chem = chem
    n = int(n)
    means = self._select_data(chem, self.data.df)
    stds = self._select_data(chem, self.data_uncertainty.df)
    result = {key: [] for key in means.keys()}
    self.dobs_ensemble = []
    for i in range(n):
      data = {}
      for id, val in means.items():
        mean = np.log(means[id]) # NOT log10, as np.random.lognormal assumes natural one
        if std_is_relative:
          std = stds[id] * mean
        else:
          std = np.log(stds[id]) # try 0.00000000000001 for checks
        data[id] = np.random.lognormal(mean, std, size=1)[0]
      self.dobs_ensemble.append(sort(data))
    return self.dobs_ensemble
  def _add_subfigure_label(self, ax, label, pos='tl', fontsize=20):
    if pos == 'tl':
      pos = [0.03, 0.92]
    elif pos == 'tr':
      pos = [0.92, 0.92]
    elif pos == 'br':
      pos = [0.92, 0.08]
    elif pos == 'bl':
      pos = [0.03, 0.08]
    else:
      pass
    ax.text(*pos, label, transform=ax.transAxes,
            fontsize=fontsize, va='top')
    return ax
  def plot_all(self, vmin, vmax, figsize=(12,10), alpha=.5, mode='multicolor',\
     field_data=None, mtrue=None, vmin_std=None, vmax_std=None, title=False,\
      annotate_lcurve=False):
    
    x0 = 0.03
    
    fig = plt.figure(figsize=figsize)
    x1, x2, y1, y2 = self.extent
    aspect = (y2 - y1) / (x2 - x1)
    numbers = np.array(list(self._select_data(self.chem, self.data.df).values()))
    geom_mean = np.prod(numbers) ** (1 / len(numbers))
    gs = gridspec.GridSpec(2,3, width_ratios=[1,1,1])
    if title:
      fig.text(0.45, 0.9, self.chem, fontsize=14)
    # fig.suptitle(self.chem)
    # gs.tight_layout(fig, rect=[0, 0.01, 1, 0.99])
    ax1 = fig.add_subplot(gs[0,0]) # ----------------------------------------------
    
    c1, c2 = 'dodgerblue', 'red'
    self._plot_boxplot(self.df_dobs, c1, ax=ax1, alpha=alpha, connect_means=0)
    self._plot_boxplot(self.df_dsyn, c2, ax=ax1, alpha=alpha, connect_means=0)
    if field_data is not None:
      c3 = 'k'
      self._plot_boxplot(field_data, c3, ax=ax1, alpha=alpha, connect_means=0)
      ax1.text(0.05, 0.7, "field", c=c3, transform=ax1.transAxes) 

    plt.xlim(np.log10(vmin), np.log10(vmax))
    # ax1.text(0.05, 0.85, "observed", c=c1, transform=ax1.transAxes) #, fontsize=12)
    # ax1.text(0.05, 0.78, "predicted", c=c2, transform=ax1.transAxes) #, fontsize=12)
    ax1.text(x0, 0.88, "observed", c=c1, transform=ax1.transAxes) #, fontsize=12)
    ax1.text(x0, 0.78, "predicted", c=c2, transform=ax1.transAxes) #, fontsize=12)
    # ax1.text(0.05, 0.95, 'a', transform=ax1.transAxes,
    #         fontsize=12, fontweight='bold', va='top')
    self._add_subfigure_label(ax1, 'a', pos='bl')
    # ax1.set_xlabel('log10(c)')
    # ax1.set_ylabel('Subcatchment')
    plt.grid(None)
    ymin, ymax = ax1.get_ylim()
    ax1.vlines(x=np.log10(self.detection_threshold), ymin=ymin, ymax=ymax, colors='k', linestyles='dashed', linewidth=2)
    ax1.vlines(x=np.log10(geom_mean), ymin=ymin, ymax=ymax, colors='k', linestyles='dotted', linewidth=2)

    ax2 = fig.add_subplot(gs[0,1]) # ----------------------------------------------
    _ = self.plot_misfit(self._select_data(self.chem, self.data.df), ax=ax2, \
      xlim=[vmin, vmax], aspect='auto', mode=mode)
    ax2.set_aspect(aspect)
    self._add_subfigure_label(ax2, 'b', pos='bl')

    ax3 = fig.add_subplot(gs[0,2]) # ----------------------------------------------
    ax3 = self.plot_misfit_vs_rough(ax=ax3, annotate=annotate_lcurve)
    xlim3 = ax3.get_xlim()
    ylim3 = ax3.get_ylim()
    ratio  = (ylim3[1] - ylim3[0]) / (xlim3[1] - xlim3[0])
    ax3.set_aspect(aspect/ratio)
    self._add_subfigure_label(ax3, 'c', pos='bl')
    ax4 = fig.add_subplot(gs[1,0]) # ----------------------------------------------
    c = 'Grey'
    ax4 = self._plot_boxplot(self.df_mfin, c, alpha=alpha, vert=0, ax=ax4)
    ax4.text(x0, .88, "recovered", c=c, transform=ax4.transAxes)
    if mtrue is not None:
      c ='k'
      mtrue.plot_dict(ax=ax4, c=c)
      ax4.text(0.05, 0.83, "true", c=c, transform=ax4.transAxes)
    

    ax4.set_xlim(np.log10(vmin), np.log10(vmax))
    ax4.set_xlabel('log10(concentration)')
    ax4.set_ylabel('Subcatchment')
    plt.grid(None)
    ax4.vlines(x=np.log10(self.detection_threshold), ymin=ymin, ymax=ymax, colors='k', linestyles='dashed', linewidth=2)
    # ax4.vlines(x=np.log10(geom_mean), ymin=ymin, ymax=ymax, colors='k', linestyles='dotted', linewidth=2)
    self._add_subfigure_label(ax4, 'd', pos='br')

    ax5 = fig.add_subplot(gs[1,1]) # ----------------------------------------------
    # ax5.axis('off')
    inset_kws = dict(extent=[524e3, 527e3, 173.5e3, 175.5e3], width=0.5, \
      ticks=False, labels=False)
    # ax5 = self.plot_model(ax=ax5, vmin=vmin, vmax=vmax, cmap='viridis', cbar=0)
    # ax5 = self.mfinal.plot_inset(ax5, vmin=vmin, vmax=vmax, **inset_kws)
    # ax5 = self.plot_model(ax=ax5, vmin=vmin, vmax=vmax, cmap='viridis', cbar=0)
    cmap = 'viridis'
    cmap_m = cmap
    self.mfinal.plot(ax=ax5, vmin=vmin, vmax=vmax, cmap=cmap, cbar=0)
    self.dsyn.plot(ax=ax5, vmin=vmin, vmax=vmax, cmap=cmap, cbar=0,\
      annotate=None, annotate_kw=dict(c='w'))
    axi = self.mfinal.plot_inset(ax=ax5, vmin=vmin, vmax=vmax, cmap=cmap,\
      **inset_kws)
    self.dsyn.plot(ax=axi, vmin=vmin, vmax=vmax, cmap=cmap, cbar=0, \
      annotate=None, annotate_kw=dict(c='k'))
    AxesFormatter()._adjust(ax=axi, labels=0, ticklabels=0, ticks=0)
    if True:
      cax = ax5.inset_axes([0.79, 0.3, .05, .4])
      cmap = plt.get_cmap(cmap)
      norm = LogNorm(vmin, vmax)
      cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),\
        ax=ax5, cax=cax, location='right',\
        orientation='vertical', label='Concentration (ng/l)')
      cb.ax.minorticks_off()
      if self.detection_threshold is not None:
        cb.ax.axhline(self.detection_threshold, color='w', linestyle='-', linewidth=2)
      # plt.yticks(rotation=90, ha='center', va='center') # not working
      # cax.yaxis.set_major_formatter(ScalarFormatter())
      
    ax5.set_aspect('equal')
    ax5.set_xlim(x1, x2)
    ax5.set_ylim(y1, y2)
    # ax5.text(0.05, .9, self.chem, c='k', fontsize=12, transform=ax5.transAxes)
    # yticks, yticklabels = plt.yticks()
    # plt.ylabel(None)
    # plt.yticks(yticks, []*len(yticks))
    self._add_subfigure_label(ax5, 'e', pos='br')    

    ax6 = fig.add_subplot(gs[1,2]) # ----------------------------------------------
    # ax6.axis('off')
    cmap = 'Reds'
    # ax6 = self.plot_std(ax=ax6, vmin=vmin, vmax=vmax, vmin_std=vmin_std, vmax_std=vmax_std, \
    #   cmap=cmap, cbar=0, cbar_label='Concentration uncertainty (ng/l)')
    self.mfinalstd.plot(ax=ax6, vmin=vmin_std, vmax=vmax_std, cmap=cmap, cbar=0)
    self.dsyn.plot(ax=ax6, vmin=vmin, vmax=vmax, cmap=cmap_m, cbar=0,\
      annotate=None, annotate_kw=dict(c='k'))
    axi = self.mfinalstd.plot_inset(ax=ax6, vmin=vmin_std, vmax=vmax_std, cmap=cmap,\
      **inset_kws)
    self.dsyn.plot(ax=axi, vmin=vmin, vmax=vmax, cmap=cmap_m, cbar=0,\
      annotate=None, annotate_kw=dict(c='k'))
    AxesFormatter()._adjust(ax=axi, labels=0, ticklabels=0, ticks=0)
    if True:
      cax = ax6.inset_axes([0.79, 0.3, .05, .4])
      cmap = plt.get_cmap(cmap)
      norm = LogNorm(vmin_std, vmax_std)
      cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),\
        ax=ax5, cax=cax, location='right',\
        orientation='vertical', label='Concentration uncertainty (ng/l)')
      cb.ax.minorticks_off()
      # plt.yticks(rotation=90, ha='center', va='center') # not working
      # cax.yaxis.set_major_formatter(ScalarFormatter())  


    # ax6.set_xlabel(None)
    ax6.set_ylabel(None)
    ax6.set_aspect('equal')
    # self.mfinalstd.plot_inset(ax6, cmap=cmap, vmin=vmin_std, vmax=vmax_std, **inset_kws)
    ax6.set_xlim(x1, x2)
    ax6.set_ylim(y1, y2)
    self._add_subfigure_label(ax6, 'f', pos='br') 
    return None
  def plot_boxplots(self, figsize=(12,6), mtrue=None, df_dfield=None):
    df_dobs, df_dsyn, df_mfin = self.df_dobs, self.df_dsyn, self.df_mfin
    df_dobs, df_dsyn, df_mfin = self._parse_output()
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1,2, width_ratios=[1,1])
    ax1 = fig.add_subplot(gs[0])
    self._plot_boxplot(df_dobs, 'g', ax=ax1)
    self._plot_boxplot(df_dsyn, 'k', ax=ax1)
    if df_dfield is not None:
      self._plot_boxplot(df_dfield, 'r', ax=ax1)
    ax2 = fig.add_subplot(gs[1])
    self._plot_boxplot(df_mfin, 'k', vert=1, ax=ax2)
    if mtrue is not None:
      df_mtru = np.log10(pd.DataFrame({k: [v] for k, v in mtrue.dict.items()}))
      self._plot_boxplot(df_mtru, 'g', vert=1, ax=ax2)
  def plot(self, vmin, vmax, vmin_std, vmax_std, basin_id, ax=None, \
    figsize=(12,6), cmap_model='viridis', cmap_std='Reds', **kwargs):
    fig = plt.figure(figsize=figsize)
    ndigits = 2
    fig.suptitle('%s concentration in basin %s:\n %s +/- %s ng/l' % (self.chem, basin_id, \
      round(self.mfinal.dict[basin_id],ndigits), round(self.mfinalstd.dict[basin_id], ndigits)))
    gs = gridspec.GridSpec(1,2, width_ratios=[1,1])
    ax1 = fig.add_subplot(gs[0])
    ax1 = self.plot_model(ax=ax1, vmin=vmin, vmax=vmax, cmap=cmap_model)
    ax2 = fig.add_subplot(gs[1])
    ax2 = self.plot_std(ax=ax2, vmin=vmin, vmax=vmax, vmin_std=vmin_std, vmax_std=vmax_std, cmap=cmap_std) #'RdYlGn_r')'Reds'
    return ax1, ax2
  def plot_misfit(self, dobs, ax=None, log=True, mode='multicolor', **kwargs):
    if ax is None:
      fig, ax = plt.subplots()
     
    # ax.text(.05, .9, 'No. of MC iterations: 10^%s' % int(np.log10(len(self.df_dsyn))), \
    #   fontsize=12, transform=ax.transAxes)
    df1, df2 = self.df_dsyn, self.df_dobs

    if mode == 'multicolor':
      # nbinx, nbiny = 100, 100
      # hm, extent = self._get_hit_count(df1, df2, nbinx, nbiny)
      # im = ax.imshow(hm.T, extent=extent, origin="lower", cmap='Greys', vmin=0, vmax=5)
      # plt.colorbar(im, label="Hit count with %sx%s bins" % (nbinx, nbiny))
      clrs = colors(len(df1.columns), cmap='rainbow')
      for basen_id in df1.columns:
        ax.scatter(df1[basen_id], df2[basen_id], s=30, marker='.', alpha=.8, edgecolor='none', \
          facecolor=next(clrs), label=basen_id)
      plt.legend(loc='lower right', frameon=1, fontsize=10, alignment='left',\
        labelspacing=0.2, markerscale=3)
      
      
    elif mode == 'hitcount':
      x, y = [], []
      for basen_id in df1.columns:
        x += list(df1[basen_id])
        y += list(df2[basen_id])
      
      nbinx, nbiny = 100, 100
      heatmap, xedges, yedges = np.histogram2d(x, y, bins=(nbinx,nbiny))
      heatmap_smooth = gaussian_filter(heatmap, sigma=1)
      extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
      plt.imshow(heatmap_smooth.T, extent=extent, origin="lower", cmap=cc.cm.fire_r,\
        vmin=0, vmax=10)
      plt.colorbar(label="Hit count with %sx%s bins" % (nbinx, nbiny))
    else:
      raise ValueError('Mode: %s' % mode)
    kwargs['mode'] = mode
    ax = super().plot_misfit(dobs, ax, log, **kwargs)
  def plot_final_misfit(self, vmin, vmax, vmin_std, vmax_std, basin_id, \
    ax=None, figsize=(12,6), cmap_model='viridis', cmap_std='Reds',\
    alpha=.5, **kwargs):
    fig = plt.figure(figsize=figsize)
    # ndigits = 2
    # fig.suptitle('%s concentration in basin %s:\n %s +/- %s ng/l' % (self.chem, basin_id, \
    #   round(self.mfinal.dict[basin_id],ndigits), round(self.mfinalstd.dict[basin_id], ndigits)))
    gs = gridspec.GridSpec(1,3, width_ratios=[.5,1,1])
    ax1 = fig.add_subplot(gs[0])
    df_dobs, df_dsyn, df_mfin = self._parse_output()
    self._plot_boxplot(df_dobs, 'dodgerblue', ax=ax1, alpha=alpha)
    self._plot_boxplot(df_dsyn, 'red', ax=ax1, alpha=alpha)
    plt.xlim(np.log10(vmin), np.log10(vmax))
    ax1.set_xlabel('log10(c)')
    ax1.set_ylabel('Subcatchment')
    ax2 = fig.add_subplot(gs[1])
    xlim = [vmin, vmax]
    _ = self.plot_misfit(self._select_data(self.chem, self.data.df), ax=ax2, xlim=xlim, **kwargs)
    # plt.ylim(np.log10(vmin), np.log10(vmax))
    ax3 = fig.add_subplot(gs[2])
    ax3 = self.plot_misfit_vs_rough(ax=ax3)
    return ax1, ax2, ax3
  def plot_final_model(self, vmin, vmax, vmin_std, vmax_std, basin_id, \
    ax=None, figsize=(12,6), \
    cmap_model='viridis', cmap_std='Reds', alpha=1., c='turquoise', topo=None, **kwargs):
    fig = plt.figure(figsize=figsize)
    ndigits = 2
    # fig.suptitle('%s concentration in basin %s:\n %s +/- %s ng/l' % (self.chem, basin_id, \
      # round(self.mfinal.dict[basin_id],ndigits), round(self.mfinalstd.dict[basin_id], ndigits)))
    gs = gridspec.GridSpec(1,3, width_ratios=[.5,1,1])
    ax3 = fig.add_subplot(gs[0])
    df_dobs, df_dsyn, df_mfin = self._parse_output()
    ax3 = self._plot_boxplot(df_mfin, c, alpha=alpha, vert=0, ax=ax3)
    ax3.set_xlim(np.log10(vmin), np.log10(vmax))
    ax3.set_xlabel('log10(c)')
    ax3.set_ylabel('Subcatchment')
    ax1 = fig.add_subplot(gs[1])
    if topo is not None:
      ax1 = topo.plot(ax=ax1)
    ax1 = self.plot_model(ax=ax1, vmin=vmin, vmax=vmax, cmap=cmap_model)
    ax2 = fig.add_subplot(gs[2])
    ax2 = self.plot_std(ax=ax2, vmin=vmin, vmax=vmax, vmin_std=vmin_std, vmax_std=vmax_std, cmap=cmap_std) #'RdYlGn_r')'Reds'
    # plt.yticks([])
    yticks, yticklabels = plt.yticks()
    plt.ylabel(None)
    plt.yticks(yticks, []*len(yticks))
    return ax1, ax2, ax3
  def plot_fluxes(self, flow_rates, ax=None, figsize=(12,6), unit='kg/day', \
    vmin=0, vmax=5):
    fmean, fsigm = self.get_fluxes(flow_rates, unit=unit) 
    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(1,2, width_ratios=[1,1])
    ax = fig.add_subplot(gs[0,0]) 
    kwargs = dict(norm='linear', cmap='Reds', alpha=.5, extent=self.extent,
                  cbar_label='Mass flux (kg/day)', vmin=vmin, vmax=vmax)
    # ax = topo.plot_all(domain, None, drainage, wfd, osd, None, cso, stw, chem, ax=ax)
    ax = fmean.plot(ax=ax, **kwargs)
    ax = fig.add_subplot(gs[0,1]) 
    # ax = topo.plot_all(domain, None, drainage, wfd, osd, None, cso, stw, chem, ax=ax)
    ax = fsigm.plot(ax=ax, **dict(kwargs, cbar_label='Mass flux uncertainty (kg/day)'))
    return ax
  def plot_std(self, vmin, vmax, vmin_std, vmax_std, ax=None, **kwargs):
    if ax is None:
      fig, ax = plt.subplots()
    ax = self.mfinalstd.plot(ax=ax, vmin=vmin_std, vmax=vmax_std, **kwargs)
    ax = self.dsyn.plot(ax=ax, vmin=vmin, vmax=vmax, **kwargs)
    return ax
  def run(self, chem, regularization_strength, std_is_relative: bool, solver='ecos', n=1000, verbose=False):
    self.regul_strength = regularization_strength
    self.chem = chem
    n = int(n)
    means = self._select_data(chem, self.data.df)
    stds = self._select_data(chem, self.data_uncertainty.df)
    result = {key: [] for key in means.keys()}

    self.dobs_ensemble = []
    self.dsyn_ensemble = []
    self.mfin_ensemble = []
    for i in range(n):
      data = {}
      for id, val in means.items():
        # mean = np.log(means[id]) # NOT log10, as np.random.lognormal assumes natural one
        # std = np.log(stds[id]) # try 0.00000000000001 for checks
        # print('%s: %s +/- %s' % (id, mean, std))
        # data[id] = np.random.lognormal(mean, std, size=1)[0]
        # data[id] = np.random.normal(means[id], stds[id], size=1)[0]

        # if np.isnan(means[id]):
        #   means[id] = self.detection_threshold
        # if np.isnan(stds[id]):
        #   stds[id] = self.detection_threshold

        # ---- Method 1
        mean = np.log(means[id]) # NOT log10, as np.random.lognormal assumes natural one
        if std_is_relative:
          std = stds[id] * mean
        else:
          std = np.log(stds[id]) # try 0.00000000000001 for checks
        data[id] = np.random.lognormal(mean, std, size=1)[0]


        # ---- Method 2
        # dat = np.log(means[id])
        # if std_is_relative:
        #   std = stds[id] * dat
        # else:
        #   std = stds[id]
        # noise = np.random.normal(0, std, size=1)[0]
        # data[id] = np.exp(dat + noise)


        # ---------------
        if data[id] <= 0:
          data[id] = self.detection_threshold
      self.dobs_ensemble.append(data)
      dsyn_final, final_model = self.problem.solve(data, export_rates=self.exrates, \
      regularization_strength=regularization_strength, solver=solver, verbose=verbose)

      self.dsyn_ensemble.append(dsyn_final)
      self.mfin_ensemble.append(final_model)
      for id, _ in result.items():
        result[id].append(final_model[id])

    result_mean = {}
    result_std = {}
    for id, all_iters in result.items():
      result_mean[id] = 10**np.mean(np.log10(np.array(all_iters)))
      result_std[id] = 10**np.std(np.log10(np.array(all_iters)))

    self.mfinal =    ModelUnmix(self.acquisition.network, self.exrates, result_mean, self.extent)
    self.mfinalstd = ModelUnmix(self.acquisition.network, self.exrates, result_std, self.extent)
    # calculate synthetics from the final model
    self.dsyn = Forward(self.mfinal, self.acquisition, solver=self.code).run()
    # self.dsyn_final = downstream_means
    # self.dsyn_std = downstream_uncertainties
    self._parse_output()
    return self.mfinal, self.mfinalstd, self.dsyn
  def _get_hit_count(self, df1, df2, nbinx=100, nbiny=100):
    """
    Note: the data is binned using the tight data range.

    Parameters
    ----------
    df1 : _type_
        _description_
    df2 : _type_
        _description_
    nbinx : int, optional
        _description_, by default 100
    nbiny : int, optional
        _description_, by default 100

    Returns
    -------
    _type_
        _description_
    """
    x, y = [], []
    for basen_id in df1.columns:
      x += list(df1[basen_id])
      y += list(df2[basen_id])
    
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=(nbinx,nbiny))
    heatmap_smooth = gaussian_filter(heatmap, sigma=1)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    return heatmap_smooth, extent
  def _parse_output(self):
    inv = self
    df_dobs, df_dsyn, df_mfin = {}, {}, {}
    for k in sorted(inv.mfinal.dict.keys()):
      df_dobs[k] = []
      df_dsyn[k] = []
      df_mfin[k] = []
      for obs, syn, fin in zip(inv.dobs_ensemble, inv.dsyn_ensemble, inv.mfin_ensemble):
        df_dobs[k].append(obs[k])
        df_dsyn[k].append(syn[k])
        df_mfin[k].append(fin[k])
    self.df_dobs = np.log10(pd.DataFrame(df_dobs))
    self.df_dsyn = np.log10(pd.DataFrame(df_dsyn))
    self.df_mfin = np.log10(pd.DataFrame(df_mfin))
    return self.df_dobs, self.df_dsyn, self.df_mfin
  def _plot_boxplot(self, df, c, vert=0, ax=None, alpha=1., connect_means=True):
    if ax is None:
      fig, ax = plt.subplots()
    rot = 0 # if vert else 90
    c_med = 'k' if c == 'Grey' else c
    c_box = 'k' if c == 'Grey' else c
    df.boxplot(ax=ax, rot=rot, showmeans=1, vert=vert, patch_artist=True, # to fill in boxes
                boxprops={'facecolor': c, 'edgecolor': 'k', 'alpha': alpha},
           whiskerprops={'color': c,},
            capprops={'color': c,},
            flierprops={'markerfacecolor': c, 'markeredgecolor': 'none',
             'marker': '.'}, # 'markersize': 7},
                 meanprops={'markerfacecolor': c, 'markeredgecolor': 'k'},
                 medianprops={'color': c_med, "zorder":100},
          )
    if connect_means:
      means = df.mean()
      x, y = range(1, len(means) + 1), means
      if not vert:
        x, y = y, x
      ax.plot(x, y, c=c)
    return ax


@logged
class EnsembleInversion:
  def __init__(self, name, inv_cls, *args, **kwargs):
    self.name = name
    self.inv_cls = inv_cls
    self.inv_args = args
    self.inv_kwargs = kwargs
  def run(self, vals, *args, **kwargs):
    self.keys = vals
    self.inv = {}
    for val in vals:
      self.inv[val] = self.inv_cls(*self.inv_args, **self.inv_kwargs)
      self.inv[val].run(*args, **kwargs)
    return self.inv
  def plot(self, mtrue, ax=None):
    pass
  def save(self, path, plots=True, results=True, **kwargs):
    if plots:
      self.save_plots(path, **kwargs)
    if results:
      self.save_results(path)
@logged
class RegularisationTuning(EnsembleInversion):
  def run(self, vals, *args, **kwargs):
    self.keys = vals
    self.inv = {}
    for val in vals:
      self.inv[val] = self.inv_cls(*self.inv_args, **self.inv_kwargs)
      self.inv[val].run(*args, regularization_strength=10.**(-val), **kwargs)
    return self.inv
  def plot(self, mode='model', **kwargs):
    if mode == 'all':
      # out = self.(, **kwargs)
      pass
    elif mode == 'model':
      out = self.plot_model_crossplot(**kwargs)
    else:
      raise ValueError('Unknown mode: %s' % mode)
    return out    
  def plot_model_crossplot(self, mtrue, ax=None):
    ax = get_ax(ax)
    clrs = colors(len(self.inv))
    for key, inv in self.inv.items():
      ax = inv.plot_mfinal_vs_mtrue(mtrue, ax=ax, c=next(clrs), cbar=0, annotate=0, alpha=.5)
    return ax
@logged
class MultichemInversion(EnsembleInversion):
  def read_results(self, chems, path, extn='json'):
    mfinalmean, mfinalstd = {}, {}
    for chem in chems:
      mfinalmean[chem] = FileModelFinalMean(path, extn, prefix=f'{chem}-').read()
      mfinalstd[chem] = FileModelFinalStd(path, extn, prefix=f'{chem}-').read()
    return mfinalmean, mfinalstd
  def run(self, chems, regularisation:dict, save=True, path_results=None, **kwargs):
    self.keys = chems
    self.inv = {}
    for chem in chems:
      self.inv[chem] = self.inv_cls(*self.inv_args, **self.inv_kwargs)
      elbow = regularisation[chem]
      self.inv[chem].get_misfit_vs_rough(chem, max_log=0, n=10, plot=False)
      self.inv[chem].run(chem, regularization_strength=10.**(-elbow), **kwargs)
    if save:
      self.save(path_results, plots=True, results=True, extn='png')
    return self.inv
  def plot(self,  vmin=.5, vmax=1e3, vmin_std=1, vmax_std=100, alpha=.5, figsize=(14,10), \
    mode='multicolor'):
    for chem, inv in self.inv.items():
      inv.plot_all(vmin=vmin, vmax=vmax, vmin_std=vmin_std, vmax_std=vmax_std, alpha=alpha, \
        figsize=figsize, mode=mode)
  def save_plots(self, path, extn='png', show=False, dpi=300, vmin=.5, vmax=1e3, vmin_std=1, vmax_std=100, alpha=.5, \
    figsize=(14,10), mode='multicolor'):
    for chem, inv in self.inv.items():
      kws = dict(path=path, extn='json', prefix=f'{chem}-')    
      inv.plot_all(vmin=vmin, vmax=vmax, vmin_std=vmin_std, vmax_std=vmax_std, alpha=alpha, \
        figsize=figsize, mode=mode)
      plt.savefig(f'{path}/{chem}.{extn}', dpi=dpi)
      if not show:
        plt.close()
  def save_results(self, path):
    for chem, inv in self.inv.items():
      kws = dict(path=path, extn='json', prefix=f'{chem}-')
      FileModelFinalMean(**kws).save(inv.mfinal.dict)
      FileModelFinalStd(**kws).save(inv.mfinalstd.dict)
      FileDataSynFromModelFinalMean(**kws).save(inv.dsyn.dict)

@logged
class OED:
  """
  Optimal experiment design (OED).

  """
  def __init__(self):
    pass
  def run(self, target_area_m2, extent, most_downstream_xy, \
    path_to_elev_file='../../alex/thames-sewage/input_dir/thames_elev.nc'):
    fname_d8 = 'tmp.nc'
    self.extent = extent
    drain = self.extract_drainage_from_elev(extent, path_to_elev_file, fname_out=fname_d8)
    catch = self.whole_catch(most_downstream_xy, extent)
    ids, locx, locy, node_map = self.optimal_sample_xy(drain.mg, target_area_m2)
    ds, subs = self.constrain(catch, ids, locx, locy, node_map)
    self.ds = ds
    self.subs = subs
    return ds, subs
  def extract_drainage_from_elev(self, extent, fname_inp='../../alex/thames-sewage/input_dir/thames_elev.nc',\
    fname_out='tmp.nc'):
    self.__log.info('Extracting drainage from elevation.')
    d8 = extract_d8_subdomain(*extent, fname_inp, fname_out, pad_width=10, pad_value=0)
    mg = ac.toolkit.load_topo(fname_out)
    # ac.toolkit.viz_drainage_area(mg)
    # self.mg = mg
    drain = Drainage()
    # drain.extract(extent, fname_inp, fname_out, pad_width=10, pad_value=0)
    drain.fname = fname_out
    drain.mg = mg
    drain.extent = extent
    return drain
  def whole_catch(self, most_downstream_xy, extent):
   self.__log.info(f'Delineating the catchment using most-downstream point at {most_downstream_xy}')
   ds = DataSamples(pd.DataFrame({
       'id': ['1'], 
       'x': most_downstream_xy[0], 
       'y': most_downstream_xy[1],
       'lat': None,
       'lon': None,
   }))
   
   drainage = Drainage().extract(extent)
   ds.process(drainage, thresh=1e6, relocate_dict={}, use_letters=0)
   catch = Catchment(ds.graph.graph, drainage.extent)
   #ax = catch.plot()
   #ax = ds.plot(ax=ax, norm=None)
   return catch
  def optimal_sample_xy(self, mg, subcatch_area_m2):
    sample_nodes_catchments = ac.autosampler.get_sample_nodes_by_area(mg, subcatch_area_m2)
    localities, node_map = ac.autosampler.process_output_dict(sample_nodes_catchments, mg)
    # ac.autosampler.save_autosampler_results(localities, node_map, mg, ".")
    # ac.autosampler.viz_sample_site_results(localities, node_map, mg)
    ids, locx, locy = [str(int(i)) for i in localities[:,0]], localities[:,1],  localities[:,2]
    return ids, locx, locy, node_map
  def constrain(self, catch, ids, locx, locy, node_map):
    """
    Retain only the samples included in the catchment.

    Parameters
    ----------
    catch : _type_
        _description_
    ids : _type_
        _description_
    locx : _type_
        _description_
    locy : _type_
        _description_
    node_map : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    self.__log.info('Excluding samples outside of the catchment')
    d = {int(k): (x,y) for k, x, y in zip(ids, locx, locy)}
    dsubs = {k: node_map == k for k in d.keys()}
    subs = {}
    basin = np.logical_or.reduce(list(catch.subcatchments.values()))
    x_, y_ = [], []
    i_, j_ = [], []
    count = 0
    for id, x, y in zip(ids, locx, locy):
      i, j = convert_xy_to_ij(x, y, catch.extent, 50, 50)
      if basin[j, i]: # yes, reversed order111
        count += 1
        x_.append(x)
        y_.append(y)
        subs[str(count)] = dsubs[int(id)]
    ds2 = DataSamples(pd.DataFrame({
      'id': [str(i) for i in range(len(x_))], 
      'x': x_, 
      'y': y_,
      'lat': None,
      'lon': None,
      }))
    return ds2, CatchmentOptimal(subs, catch.extent)
  def plot(self, topo=None, ax=None):
    if topo is None:
      topo = Topography()
      topo.extract(domain.extent)
    
  def save(self, fname):
    # fname = '../data/oed-wandle-area%skm2.csv' % "{:.2f}".format(subcatch_area_m2 / 1e6).rjust(5,'0')
    self.ds.df.to_csv(fname, index=0)
@logged
class OEDthames(OED):
  def whole_catch(self, most_downstream_xy, extent):
    return None
  def constrain(self, catch, ids, locx, locy, node_map):
    self.__log.info('Excluding samples outside of the catchment')
    d = {int(k): (x,y) for k, x, y in zip(ids, locx, locy)}
    dsubs = {k: node_map == k for k in d.keys()}
    subs = {}
    x_, y_ = [], []
    i_, j_ = [], []    
    count = 0
    for id, x, y in zip(ids, locx, locy):
      if True:
        count += 1
        x_.append(x)
        y_.append(y)
        subs[str(count)] = dsubs[int(id)]
    ds2 = DataSamples(pd.DataFrame({
      'id': [str(i) for i in range(len(x_))], 
      'x': x_, 
      'y': y_,
      'lat': None,
      'lon': None,
      }))
    return ds2, CatchmentOptimal(subs, self.extent)

@traced
@logged
class PCA:
  def __init__(self, path_results:str, chems:list):
    self.path_results = path_results
    self.chems = chems
  def _check_data_matrix(self, df):
    m, n = df.shape 
    if m <= n:
      raise ValueError(f'No. of data points m={m} smaller than no. of dimensions n={n}')    
  def read_inversion_results(self, chems, path, extn='json'):
    dicts = []
    for chem in chems:
      di = FileModelFinalMean(path, extn, prefix=f'{chem}-').read()
      dicts.append(di)
    df = pd.DataFrame(dicts).T
    df = df.rename(columns={k: v for (k,v) in zip(list(df.columns), chems)})
    self._check_data_matrix(df)
    return df
  def transform_to_PC_basis(self, pca):
    dfpc = pd.DataFrame(dict(pca.clr_df.apply(pca.transform, axis=1))).T
    dfpc.columns = ['PC%s' % (i+1) for i in range(len(pca.clr_df.columns))] 
    self.dfpc = dfpc
    return self.dfpc   
  def run(self, plot=True):
    df = self.read_inversion_results(self.chems, self.path_results)
    pca = CompositionalPCA(df)
    pca.fit_pca()    
    self.transform_to_PC_basis(pca)
    return pca
  def plot(self, means, drainage_extent, pc='PC1', ax=None, title=None, show_title=True, cmap='PuOr', norm='linear',\
    vmin=-10, vmax=10, **kwargs):
    ax = get_ax(ax)
    kwargs['cmap'] = cmap
    kwargs['norm'] = norm
    kwargs['vmin'] = vmin
    kwargs['vmax'] = vmax
    title = pc if title is None else title
    if show_title:
      ax.set_title(title)
    m = ModelUnmix(means.graph.graph, {}, dict(self.dfpc[pc]), drainage_extent)
    ax = m.plot(ax=ax, cbar_label=f'{pc} score', **kwargs)
    ax = m.plot_inset(ax=ax, extent=[524e3, 527e3, 174e3, 175.5e3], ticks=0, \
      width=0.5, **kwargs)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')
    return ax

