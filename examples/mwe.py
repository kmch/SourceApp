from jupconfig.loggers import *
from riverpy import *
from plotea import *
log_lvl(20) # change verbosity of the logs (diagnostic messages), the higher the number, the less verbose
domain, topo, drain, dobs, dobs_unc, clusters, acq = InputPrep.run(dataset='ICL', river='wandle', 
                                                                   max_dist=100)
dataset = clusters.data_samples.dataset
catch = Catchment(dobs.graph.graph, drain.extent)
# dobs.plot()
wf = MultichemInversion('multichem01', 
                        InversionMC, # type of inversion (InversionMC - Monte Carlo inversion)
                        acq, # acquisition geometry object
                        detection_threshold=3, # ng/l
                        data=dobs, data_uncertainty=dobs_unc, 
                        extent=drain.extent, # spatial extent of apportionment
                        cp_norm=2, # p in the L^p, i.e. the vector norm used in the objective funciton
)
lambda_dict = {k: 1.7 for k in dataset.cecs}
for k in ['Benzoylecgonine', 'Cocaine', 'Salicylic acid']:
    lambda_dict[k] = 2.2
log_lvl(50)
wf.run(chems=dataset.cecs, # select contaminants of emerging concern analysed in the paper
       regularisation=lambda_dict, 
       n=1e2, # number of Monte Carlo runs
       std_is_relative=False, # define whether the input data_uncertainty is relative or absolute
       save=False, # save results as json files?
)
wf.save('./', plots=True, results=False, extn='png')
