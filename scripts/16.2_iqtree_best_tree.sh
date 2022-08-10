
# get tree likelihood scores
grep -e "BEST SCORE FOUND" iqtree/*log | sed -e 's/:BEST SCORE FOUND : /\t/g' -e 's/.log/.treefile/g' > iqtree/tree_likelihood.txt

# identify tree with the best score
BEST_TREE=$(sort -k 2 -n -r iqtree/tree_likelihood.txt | head -n 1 | cut -f 1)

# create symbolic link
ln -sf ${PWD}/$BEST_TREE ${PWD}/iqtree/best.tre.iqtree

