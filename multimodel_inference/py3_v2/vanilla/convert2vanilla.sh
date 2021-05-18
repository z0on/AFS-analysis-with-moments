for PY in *_ga.py; do
cat $PY | perl -pe "s/# uses genetic algorithm from GADMA for optimization/# does NOT use GADMA/" | perl -pe "s/import gadma/#import gadma/" | perl -pe "s/^par_labels/\'\'\'\npar_labels/" | perl -pe "s/poptg=result.x/poptg=result\.x\n\'\'\'\npoptg = moments\.Inference\.optimize_log\(params\, data\, func\,\n                                   lower_bound=lower_bound\,\n                                   upper_bound=upper_bound\,\n                                   verbose=False\, maxiter=30\)/" >${PY/_ga.py/}.py;
done

