

#F="sc_ga.py"

for F in [aIs]*_ga.py; do
cat $F |\
perl -pe 's/moments\.Spectrum\.from_file\(infile\)/moments\.Spectrum\.from_file\(infile\)\nfs=fs.fold\(\)/' |\
perl -pe 's/\,p_misid//' |\
perl -pe 's/..f_misid.+$/\)/' |\
perl -pe 's/return.+misid.+(fs2).+/return fs2/' |\
perl -pe 's/return.+misid.+(fs).+/return fs/' |\
perl -pe 's/^(.+er_bound.+)\,.+]$/$1\]/' > fold_$F;
done
chmod +x *.py

cat allmodels_unfolded | perl -pe 's/^(.+)/fold_$1/' > allmodels_folded

grep " = params" fold_*py | perl -pe 's/ =.+//' | perl -pe 's/fold_(\S+)_.*(\s.+)/$1$2/'>folded_params
grep " = params" [aIs]*py | perl -pe 's/ =.+//' | perl -pe 's/^(\S+)_.*(\s.+)/$1$2/'>unfolded_params
