function PC = SolvePC_TDPC(opts)

arguments
    opts.dist = "Normal";
    opts.rv_mean = 1;
    opts.rv_level = 0.1;
    opts.nSample = 10^3;
    opts.rv_name = "l1";
    opts.gauss_point = 50;
    opts.D = 1;
    opts.o = 4;
    opts.h = 1e-2;
    opts.endTime = 30;
    opts.th1_0 = 0;
    opts.stopcrt_t = 1e-2;
end
opts.rv_std = opts.rv_mean*opts.rv_level;

%% gauss quadrature
[loc,wts] = quadrature(opts.dist,opts.gauss_point);
%% basis quadrature & sample
basis = GenBasis(loc,opts.D,opts.o,opts.dist);
[XiSample,BasisSample] = BasisSampling(opts.nSample,opts.D,opts.o,opts.dist);
%% Galerkin method
Galerkin = PC_Pre(opts,basis,loc,wts,opts.rv_mean,opts.rv_std);
%% Slider crank model
[Prop,init] = SinglePendulum("h",opts.h,"endTime",opts.endTime,"th1_0",opts.th1_0,"L",opts.rv_mean);

PC = SolveStochasticDynamics_TDPC(Prop,init,Galerkin,XiSample,BasisSample,opts);
PC.XiSample = XiSample;
