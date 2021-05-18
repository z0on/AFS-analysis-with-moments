for PY in *_ga.py; do
cat $PY | perl -pe "s/# uses genetic algorithm from GADMA for optimization/# does NOT use GADMA/" | perl -pe "s/import gadma/#import gadma/" | perl -pe "s/^if len\(sys\.argv\)==9/\'\'\'\nif len\(sys\.argv\)==9/" | perl -pe "s/poptg=result.x/poptg=result\.x\n\'\'\'\nif len(sys\.argv)==9\:\n    params = np\.loadtxt(sys\.argv\[8\]\, delimiter=\" \"\, unpack=False\)\nelse\:\n        params=\[1\]\*\(len\(upper_bound\)-1\)\n    params\.append\(0\.01\)\n\nparams = moments\.Misc\.perturb_params(params\, fold=2\, upper_bound=upper_bound\, lower_bound=lower_bound\)\n\npoptg = moments\.Inference\.optimize_log\(params\, data\, func\,\n                                   lower_bound=lower_bound\,\n                                   upper_bound=upper_bound\,\n                                   verbose=False\, maxiter=30\)/" >${PY/_ga.py/}.py;
done

if len(sys.argv)==9:
    params = np.loadtxt(sys.argv[8], delimiter=" ", unpack=False)
else:
    params=[[1]*(len(upper_bound)-1),0.01]
