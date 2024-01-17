
def within_group_variance(chems, groups, ddof):
  from math import isnan
  good_chems = []
  pvals = []
  within = {}
  for chem in chems:
      means, varss = [], []
      for g in groups:
          var = np.log10(g[chem]).var(skipna=True, ddof=ddof)
          mean = np.log10(g[chem]).mean(skipna=True) # btw. skipna=True by default
          means.append(mean)
          varss.append(var)
      if not (any(isnan(x) for x in means) or any(isnan(x) for x in varss)):
          within[chem] = means, varss
  return within

def between_group_variance(within_group_variance: dict, ddof: int):
  between = {}
  for key in within_group_variance.keys():
    means, varss = within_group_variance[key]
    between[key] = np.var(means, ddof=ddof)
  return between

def variance_ratios_old(samples, campaign, rivers, river, acc, ddof):

  chems = get_list_of_chemicals(samples)
  if campaign == 'all':
    df = samples
  else:
    df = samples[samples['campaign'] == campaign]
  riv, df = select_subset(df, rivers, river)
  groups = aggregate_by_proximity(df, acc=acc)

  within = within_group_variance(chems, groups, ddof)
  between = between_group_variance(within, ddof)
  ratios = {}
  for key in within.keys():
    ratios[key] = between[key] / np.array(within[key][1])

  fig, ax = plt.subplots(figsize=[13,5])
  ax.plot(ratios.keys(), ratios.values(), '.')
  F1 = 1
  ax.plot(ratios.keys(), np.full(len(ratios),F1), 'r--', label='ratio=%s' % F1)
  plt.ylabel('between/within-group variance ratio')
  # plt.legend()
  _ = plt.xticks(rotation=90)    
  plt.savefig('%s_campaign_%s_variance_ratios.png' % (river, campaign))
  return within, between, ratios
def variance_ratios(samples, acc, ddof):
  chems = get_list_of_chemicals(samples)
  groups = aggregate_by_proximity(samples, acc=acc)

  within = within_group_variance(chems, groups, ddof)
  between = between_group_variance(within, ddof)
  ratios = {}
  for key in within.keys():
    ratios[key] = between[key] / np.array(within[key][1])

  fig, ax = plt.subplots(figsize=[13,5])
  ax.plot(ratios.keys(), ratios.values(), '.')
  F1 = 1
  ax.plot(ratios.keys(), np.full(len(ratios),F1), 'r--', label='ratio=%s' % F1)
  plt.ylabel('between/within-group variance ratio')
  # plt.legend()
  _ = plt.xticks(rotation=90)    
  # plt.savefig('%s_campaign_%s_variance_ratios.png' % (river, campaign))
  return within, between, ratios

def plot_p_values(ax, chems, pvals, signif_level, *args, **kwargs):
  ax.plot(chems, pvals, *args, **kwargs)
  ax.plot(chems, np.full(len(pvals),signif_level), 'r--', label='p=%s' % signif_level)
  plt.legend()
  plt.ylabel('p value')
  _ = plt.xticks(rotation=90) 
  return ax

def anova_old(samples, campaign, rivers, river, signif_level=0.05, acc=100):
  from scipy.stats import f_oneway

  chems = get_list_of_chemicals(samples)
  if campaign == 'all':
    df = samples
  else:
    df = samples[samples['campaign'] == campaign]
  riv, df = select_subset(df, rivers, river)
  groups = aggregate_by_proximity(df, acc=acc)

  good_chems = []
  pvals = []
  Fratios = []
  for chem in chems:
    try:
      F, p = f_oneway(*[np.log10(group[chem].dropna()) for group in groups])
    except TypeError:
      continue
    if p is not np.nan:
      pvals.append(p)
      Fratios.append(F)
      good_chems.append(chem)

  fig, ax = plt.subplots(figsize=[13,5])
  plot_p_values(ax, good_chems, pvals, signif_level, '.')
  plt.savefig('%s_campaign_%s_pvalues.png' % (river, campaign))
  fig, ax = plt.subplots(figsize=[13,5])
  ax.plot(good_chems, Fratios, '.')
  F1 = 1
  ax.plot(good_chems, np.full(len(Fratios),F1), 'r--', label='F=%s' % F1)
  plt.ylabel('F ratio')
  plt.legend()
  _ = plt.xticks(rotation=90)
  plt.savefig('%s_campaign_%s_Fratios.png' % (river, campaign))

  return good_chems, Fratios, pvals
def anova(df, acc=100, signif_level=0.05, dropna=True):
  from scipy.stats import f_oneway
  chems = get_list_of_chemicals(df)
  groups = aggregate_by_proximity(df, acc=acc)
  good_chems = []
  pvals = []
  Fratios = []
  for chem in chems:
    try:
      if dropna:
        F, p = f_oneway(*[np.log10(group[chem].dropna()) for group in groups])
      else:
        F, p = f_oneway(*[np.log10(group[chem].fillna(1e-9)) for group in groups])
    except TypeError:
      print('ANOVA failed for %s' % chem)
      continue
    if p is not np.nan:
      pvals.append(p)
      Fratios.append(F)
      good_chems.append(chem)

  # fig, ax = plt.subplots(figsize=[13,5])
  # plot_p_values(ax, good_chems, pvals, signif_level, '.')
  # # plt.savefig('%s_campaign_%s_pvalues.png' % (river, campaign))
  # fig, ax = plt.subplots(figsize=[13,5])
  # ax.plot(good_chems, Fratios, '.')
  # F1 = 1
  # ax.plot(good_chems, np.full(len(Fratios),F1), 'r--', label='F=%s' % F1)
  # plt.ylabel('F ratio')
  # plt.legend()
  # _ = plt.xticks(rotation=90)
  # # plt.savefig('%s_campaign_%s_Fratios.png' % (river, campaign))

  return good_chems, Fratios, pvals
