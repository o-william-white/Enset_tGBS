source /data/home/mpx469/software/python-virtualenv/ete3/bin/activate

mkdir -p best_trees

cp raxml_ng/best.tre.fbp.raxml.support best_trees
cp raxml_ng/best.tre.tbe.raxml.support best_trees
cp iqtree/best.tre.iqtree best_trees

python additional_scripts/root_and_ladderise_tree.py \
  best_trees/best.tre.fbp.raxml.support \
  best_trees/best.tre.fbp.raxml.support.rooted.newick \
  pop82,pop160,pop162 

python additional_scripts/root_and_ladderise_tree.py \
  best_trees/best.tre.tbe.raxml.support \
  best_trees/best.tre.tbe.raxml.support.rooted.newick \
  pop82,pop160,pop162 
  
python additional_scripts/root_and_ladderise_tree.py \
  best_trees/best.tre.iqtree \
  best_trees/best.tre.iqtree.support.rooted.newick \
  pop82,pop160,pop162 

