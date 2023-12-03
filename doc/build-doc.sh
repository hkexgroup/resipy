for f in ../jupyter-notebook/nb_*
do
fname=${f/\.\.\/jupyter-notebook/gallery}
cp $f $fname
done

cp -r "../jupyter-notebook/img/" "gallery/img"

sphinx-build -M html . _build
