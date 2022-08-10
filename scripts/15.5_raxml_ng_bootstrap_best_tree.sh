ml gcc/7.1.0
export PATH=/data/home/mpx469/software/raxml-ng-pthreads/raxml-ng/bin/:$PATH

# see best scoring trees
# grep "Final LogLikelihood:" raxml_ng/tree_search_*.raxml.log | sort -k 3 | head

# create symbolic link to best tree
ln -sf `grep "Final LogLikelihood:" ${PWD}/raxml_ng/tree_search_*.raxml.log | sort -k 3 -n -r | head -n 1 | cut -f 1 -d ":" | sed 's/log/bestTree/g'` ${PWD}/raxml_ng/best.tre

raxml-ng \
   --support \
   --tree raxml_ng/best.tre \
   --bs-trees raxml_ng/bootstrap_trees \
   --bs-metric fbp \
   --prefix raxml_ng/best.tre.fbp

raxml-ng \
   --support \
   --tree raxml_ng/best.tre \
   --bs-trees raxml_ng/bootstrap_trees \
   --bs-metric tbe \
   --prefix raxml_ng/best.tre.tbe

