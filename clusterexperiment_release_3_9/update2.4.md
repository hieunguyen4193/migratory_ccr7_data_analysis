# New to Version 2.3.0

## Change to `ClusterExperiment` Object

**Important** If you have objects created since 2.0.0 but with a version < 2.3.0 (i.e. including 2.0.0), you should run `updateObject` to update the class definition because there have been changes to the class definition since that version:

```
ceObj<-updateObject(ceObj)
```

*Warning* This command will, however, loose information saved about the last `mergeClusters` call that you have made if your object is from version < 2.1.4. You may want to save that information and manually update the slots. If you do so, make sure you call `validObject` to make sure that you have done so correctly (in particular, you will have to have a value for the slot `merge_demethod`, see `?ClusterExperiment` which is a new slot). For example,

```
ceObjNew<-updateObject(ceObj)
ceObjNew@merge_index<-ceObj@merge_index
<etc>
```

If you have objects from before 2.0.0 (when the class was called 'clusterExperiment'), you should construct a new object using the `ClusterExperiment` function. For example,

```
ceObjNew<-ClusterExperiment(
	se=as(ceObj,"SingleCellExperiment"), 
	clusterMatrix=ceObj@clusterMatrix, 
	<etc>
	)
```

See `?ClusterExperiment` for the names of the slots. 


## Other important changes

There have also been a number of changes and enhancements to the package. These are the most important (a complete list is detailed in the [NEWS](https://github.com/epurdom/clusterExperiment/blob/master/NEWS) file of the package -- all releases since May 1, 2018)

* We have changed the function `combineMany` to `makeConsensus`. This has resulted in changes to the names of the arguments of `RSEC`
	- `combineProportion` -> `consensusProportion` in `RSEC`
	- `combineMinSize` -> `consensusMinSize` in `RSEC`
* Add functionality to `getBestFeatures` to allow `edgeR` for DE, as well as weights used with `edgeR` for compatability with weights to handle zero-inflation. As part of this change  `isCount` argument has been replaced with more fine-grained `DEMethod` argument in `getBestFeatures`, `mergeClusters`; and the argument `mergeDEMethod` in `RSEC` is now available.
* We have changed the argument `sampleData` in various plotting commands to `colData` to better indicate that the argument is to identify columns in `colData` that should also be plotted. Furthermore `plotDendrogram` now takes the argument `colData` for plotting of information in `colData` with the dendrogram.
* We have changed the names of arguments related to unassigned (`-1` or `-2` assignments) to more consistently use the term "unassigned", as well as adding the function `assignUnassigned`:
	- argument `removeNegative` -> `removeUnassigned` in `getBestFeatures` 
	- argument `ignoreUnassignedVar` -> `filterIgnoresUnassigned` in `mergeClusters` (and other functions) for clarity.
	- function `removeUnclustered` -> `removeUnassigned`
* New plotting functions:
	- `plotTableClusters`
	- `plotFeatureScatter`
* Allow the arguments `subsample` and `sequential` to `RSEC` to allow for opting out of those options for large datasets (but default is `TRUE` unlike `clusterMany`)
* The argument `whichAssay` is added to most functions to allow the user to select the assay on which the operations will be performed.
* We've changed how we store the cluster hierarchies so that we now use the `phylo4d` class of `phylobase` package (previously we stored them as a `dendrogram` class). This makes it easier to store information about the dendrograms and manipulate them. There are various helper functions related to this change. See `?clusterDendrogram`. 
* We now store the coClustering matrix in the `coClustering` slot as a `sparseMatrix` class from the package `Matrix`. This will reduce the size of the object in memory. 






