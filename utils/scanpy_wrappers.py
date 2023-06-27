from biofile_handling import SampleDict, GeneListFile, CellAnnotFile
from string_functions import prefixify
from string_functions import make_gene_list
import scanpy as sc
import pandas as pd
from itertools import chain

class ScanpyMetaObject():
    """Object that collects standardized collections of scanpy functions.
    
    You can access the specific scanpy data in this object using `ScanpyMetaObject.adata`.  
    Use this variable like you would use `adata` in any scanpy operation.  

    Args:
        matrix (GxcFile | ExcFile): BioFile object of the matrix file.
        sampledict (SampleDict): a SampleDict object from the BioFileDocket.
    """
    def __init__(self, matrix, sampledict: SampleDict):
        
        self.matrix = matrix
        self.sampledict = sampledict
        self.matrixname = self.matrix.filename
        self.species = self.sampledict.species
        self.directory = self.sampledict.directory
        self.conditions = self.sampledict.conditions
    
    @property
    def species_prefix(self):
        if '_' not in self.species:
            return self.species
        else:
            return prefixify(self.species)
    
    @property
    def matrixpath(self):
        return self.directory + self.matrixname
    
    @property
    def datatype(self):
        return str(type(self.matrix)).split('.')[1].split('File')[0]
    
    @property
    def resultsfile(self):
        return self.matrixname + '_scresults.h5ad'
    
    @property
    def resultspath(self):
        return self.directory + self.resultsfile
    
    def read(self, delimiter = '\t', cache = True, transpose = True, filter_set = set()):
        """Reads the contained matrix into a scanpy.adata object, transposing if needed.
        
        Args:
            delimiter (str): a delimeter for the data, such as ',' or '\t'. Defaults to '\t'.
            cache (bool): whether or not to create a cache file for the matrix. Defaults to True.
            transpose (bool): whether or not to transpose the data matrix after loading. Defaults to True.
            filter_set (set): a set of gene ids to keep from the data, removing those that don't match.  
                Transposes matrix after filtering if transpose == True.
        """
        self.adata = sc.read(self.matrixpath, cache = cache, delimiter = delimiter)
        
        if len(filter_set) != 0:
            self.adata = self.adata[self.adata.obs.index.isin(filter_set)].copy()
        
        if transpose:
            self.adata = self.adata.transpose()
    
    def violin(self, x = 'total_counts', y = 'n_genes_by_counts', plot = True):
        """Runs sc.pp.calculate_qc_metrics, sc.pl.violin, and sc.pl.scatter.
        
        Args:
            x (str): label of the feature to be plotted on the x axis of the violin & scatter plot.  
                Defaults to 'total_counts'.
            y (str): label of the feature to be plotted on the y axis of the violin & scatter plot.
                Defaults to 'n_genes_by_counts'.
        """
        sc.pp.calculate_qc_metrics(self.adata, percent_top=None, log1p=False, inplace=True)
        
        if plot:
            sc.pl.violin(self.adata, [y, x],
             jitter=0.4, multi_panel=True)
            sc.pl.scatter(self.adata, x=x, y=y)
    
    def cellgene_filter(self, min_genes=100, min_cells=20):
        """Filters data by minimum number of genes expressed per cell and minimum expressing cells per gene.
        
        Args:
            min_genes (int): minimum number of genes expressed per cell as a cutoff.
            min_cells (int): minimum number of cells expressed per gene as a cutoff.
        """
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
    

    def normalize(self, max_n_genes_by_counts = 7000, target_sum = 1e4):
        """Normalizes data using a max_n_genes_by_counts cutoff and a target sum.

        Modifies underlying `adata` using cutoff and normalization.
        
        Args:
            max_n_genes_by_counts (int): removes cells that exceed this number of counds.
            target_sum (float): target number of counts per cell, using scientific notation.
        """
        self.adata = self.adata[self.adata.obs.n_genes_by_counts < max_n_genes_by_counts, :]
        sc.pp.normalize_total(self.adata, target_sum = target_sum)
        sc.pp.log1p(self.adata)
    
    def variable_filter(self, min_mean = 0.0125, max_mean = 3, min_disp = 0.1, max_disp = 10, plot = True):
        """Filters genes to use for dimensionality reduction using max/min of mean and max/min of dispersion, plotting optionally.
        
        Modifies underlying `adata` and saves original data under `adata.raw`.

        Args:
            min_mean (float | int): minimum mean expression of genes.
            max_mean (float | int): maximum mean expression of genes.
            min_disp (float | int): minimum dispersion of gene expression.
            max_disp (float | int): maximum dispersion of gene expression.
            plot (bool): whether or not to generate a scatter plot showing cutoffs.
        """

        sc.pp.highly_variable_genes(self.adata, min_mean = min_mean, max_mean=max_mean, min_disp=min_disp, max_disp=max_disp)
        if plot:
            sc.pl.highly_variable_genes(self.adata)
        
        self.adata.raw = self.adata
        self.adata = self.adata[:, self.adata.var.highly_variable]
    
    def regress_scale(self, how = ['total_counts'], max_value = 10):
        """Runs regression of specific features and scales data to a maximum value.
        
        Args:
            how (list): feature to regress out. Defaults to 'total_counts'. 
            max_value (float | int): maximum value to scale data.
        """

        sc.pp.regress_out(self.adata, how)
        sc.pp.scale(self.adata, max_value = max_value)
    
    def map_cellannots(self, cellannot):
        """Adds an additional .obs feature for cell annotation.

        Imports from a file that has two columns: `cell_barcode` and `cell_type`.  
        The `cell_barcode` field should be an exact match for cells in the data.  
        Cells without a matching barcode are given the `celltype` label `'Unlabeled'`.
        
        Args:
            cellannot (CellAnnotFile): CellAnnotFile object, usually in the same BioFileDocket.
        """

        cell_ids = pd.DataFrame({'cell_barcode':self.adata.obs.index})
        cell_annots = pd.read_csv(cellannot.path, sep = '\t')
        cell_ids = cell_ids.merge(cell_annots, on = 'cell_barcode', how = 'left')
        cell_ids.fillna('Unlabeled', inplace = True)
        self.adata.obs['celltype'] = cell_ids['celltype'].values
    
    def map_cellannots_multispecies(self, msd):
        """Adds an additional .obs feature for cell annotation from a MultiSpeciesBioFileDocket.

        Imports from a file that has two columns: `cell_barcode` and `cell_type`.  
        The `cell_barcode` field should be an exact match for cells in the data.  
        Cells without a matching barcode are given the `celltype` label `'Unlabeled'`.
        
        Args:
            msd (MultiSpeciesBioFileDocket): the docket containing information about files from all species in the dataset.  
                The function looks for the BioFile object at the key `cellannot` for each species.  
        """

        cell_ids = pd.DataFrame({'cell_barcode':self.adata.obs.index})
        
        cell_annots = pd.DataFrame()
        
        for i, pre in enumerate(msd.species_BioFileDockets):
            dummy = pd.read_csv(msd.species_BioFileDockets[pre].cellannot.path, sep = '\t')
            dummy['cell_barcode'] = pre + '_' + dummy['cell_barcode']
            dummy['celltype'] = pre + '_' + dummy['celltype']
            
            if i == 0:
                cell_annots = dummy
            else:
                cell_annots = pd.concat([cell_annots, dummy])
            
        cell_ids = cell_ids.merge(cell_annots, on = 'cell_barcode', how = 'left')
        cell_ids.fillna('Unlabeled', inplace = True)
        self.adata.obs['celltype'] = cell_ids['celltype'].values
    
    def pca_basic(self, svd_solver='arpack', color = [], plot = True, **kwargs):
        """Runs a simple PCA, displaying the first two components and variance plot.
        
        Args:
            svd_solver (str): SVD sovling function, passed to `sc.tl.pca()`.
            color (list): list of coloring schemes for points in the data; creates one plot per color scheme.
            plot (bool): whether or not to make a plot.
        """
        sc.tl.pca(self.adata, svd_solver=svd_solver)
        
        if not plot:
            return None
        
        if color != []:
            sc.pl.pca(self.adata, color = color, save = self.species_prefix + self.datatype + '_pca.pdf', **kwargs)
        else:
            sc.pl.pca(self.adata, save = self.species_prefix + self.datatype + '_pca.pdf', **kwargs)
            sc.pl.pca_variance_ratio(self.adata, log=True)
    
    def umap_leiden(self, n_neighbors=50, n_pcs = 40, legend_loc='on data', save = True, plot = True, save_file = '', **kwargs):
        """Runs Leiden clustering, followed by UMAP, on the data.
        
        Args:
            n_neighbors (int): number of neighbors to pass to `sc.pp.neighbors`.
            n_pcs (int): number of Principal Components to use from PCA, passed to `sc.pp.neighbors`.
            legend_loc (str): where the legend should go (e.g. `'on data'`.)
            save (bool): whether to save the plot.
            plot (bool): whether to make the plot.
        """
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata)
        
        if not plot:
            return None
        
        if save_file == '':
            save_file = '_'.join([self.species_prefix, self.datatype, 'leiden.pdf'])
        
        if save:
            sc.pl.umap(self.adata, color=['leiden'], legend_loc=legend_loc, save = save_file, **kwargs)
        else:
            sc.pl.umap(self.adata, color=['leiden'], legend_loc=legend_loc, **kwargs)
    
    def rank_genes(self, on = 'leiden', method = 't-test', plot = True, n_genes = 25, sharey = False):
        """Runs `sc.tl.rank_genes_groups` based on passed parameters.

        Parameters can be changed to alter the ranking scheme.

        Args:
            on (str): clustering feature to use. Defaults to 'leiden'.
            method (str): how to compare groups for gene ranking. Defaults to 't-test'.
            plot (bool): whether or not to plot the data.
            n_genes (int): number of genes to plot.
            sharey (bool): pass to `sc.pl.rank_genes_groups`.
        """

        sc.tl.rank_genes_groups(self.adata, on, method=method)
        sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=sharey)
        self.adata.write(self.resultsfile)
    
    def get_top_genes(self, datatype: str, top_number = 200, tofile = True):
        """Get the top_number genes per cluster across all clusters without repetition.  

        Also get the top marker gene per cluster.  
        Save these values as object attributes with a special key to the parent object.

        Args:
            datatype (str): descriptor for the data type, to be included in the outfile name.  
                Note: This is not checked for corretness!
            top_number (int): number of top genes to pull out from each cluster.
            tofile (bool): whether to save a file. (Not implemented yet.)
        
        Returns:
            (tuple[str, str]): keys of the top gene list files.
        """
        from biofile_handling import GeneListFile
        
        self.marker_genes_df = pd.DataFrame(self.adata.uns['rank_genes_groups']['names'])
        top_genelist = list(set(chain.from_iterable([i for i in self.marker_genes_df.iloc[0:top_number].values.tolist()])))
        top_marker_genes = list(set(self.marker_genes_df.iloc[0].values))
        
        genelisttype = '_'.join(['top', str(top_number) + 'DE', datatype, 'list'])
        markergenetype = '_'.join(['top', 'marker', datatype, 'list'])
        setattr(self, genelisttype, top_genelist)
        setattr(self, markergenetype, top_marker_genes)
        
        return genelisttype, markergenetype
    
    def export_top_genes(self, key: str, filename_overwrite = ''):
        """Exports the top genes list to a text file, assigning the object as an attribute of the parent object.
        
        Args:
            key (str): a key for the gene list, ending in `'_list'`.
        
        Raises:
            Exception: when the key does not end in `'_list'`.
        """
        if '_list' not in key:
            raise Exception('key must be an attribute ending in "_list"')
        
        top_type = key.split('_')[1]
        datatype = key.split('_')[2]
        
        top_genelist_filename = '_'.join([self.species_prefix, self.datatype, 'top', top_type, datatype, 'ids.txt'])
        
        if filename_overwrite != '':
            top_genelist_filename = filename_overwrite
        
        top_genelist_file = GeneListFile(filename = top_genelist_filename, 
                                         sampledict = self.sampledict, 
                                         sources = self.matrix, 
                                         genes = getattr(self, key), 
                                         identifier = datatype,
                                         autoname = False)
        setattr(self, key + '_file', top_genelist_file)
    
    # Using an IdmmFile object, convert IDs from one type to another
    # Usually you'll be converting from gene_name to the embedding (e.g. Orthogroup)
    # Prints a table of the mappings
    # Returns a list of ids for the ones actually present in the scanpy object
    def map_gene_to_id(self, idmm, gene_list: list, from_id: str, to_id: str, check_ids = True):
        """Maps IDs from one type into another using an IdmmFile object.
        
        Args:
            idmm (IdmmFile): the IdmmFile object you're converting IDs between.
            gene_list (list of str): list of genes in the dataset you want to convert.
            from_id (str): starting feature column in the idmm.
            to_id (str): ending feature in the idmm.
            check_ids (bool): whether to check if the ids in gene_list are present in the parent object.  
                If true, only returns mappings that are actually found in the parent data.
        """

        import pandas as pd
        idmm = pd.read_csv(idmm.path, sep = '\t')

        id_table = idmm[idmm[from_id].isin(gene_list)][[from_id, to_id]].drop_duplicates()
        display(id_table)
        
        if check_ids:
            ids = [i for i in id_table[to_id].values if i in list(self.adata.var.index)]
        else:
            ids = [i for i in id_table[to_id].values]
        
        return ids

def diagonalize_df(df):
    """Sorts values of a 2d dataframe along their diagonal for prettier plotting."""
    max_order = [list(df.loc[i]).index(max(df.loc[i])) for i in df.index]
    reordered = df.copy(deep = True)
    reordered['max_col'] = max_order
    reordered = reordered.sort_values(axis = 'index', by = 'max_col')
    reordered.drop('max_col', axis = 1, inplace = True)
    
    return reordered