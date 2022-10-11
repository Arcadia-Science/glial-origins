from biofile_handling import SampleDict
from biofile_handling import GeneListFile
from string_functions import prefixify
from string_functions import make_gene_list
import scanpy as sc
import pandas as pd
from itertools import chain

# A new class to contain the metadata and functions associated with a scanpy analysis
# To directly access the normal scanpy object within this class, you can just use ScanpyMetaObject.adata
# Use this in place of the 'adata' variable you would see in most scanpy examples
class ScanpyMetaObject():
    def __init__(self, matrix, sampledict: SampleDict):
        
        self.matrix = matrix
        self.sampledict = sampledict
        self.matrixname = self.matrix.filename
        self.species = self.sampledict.species
        self.directory = self.sampledict.directory
        self.conditions = self.sampledict.conditions
    
    @property
    def species_prefix(self):
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
    
    # read in a file using scanpy read, transposing it if needed
    def read(self, delimiter = '\t', cache = True, transpose = True):
        self.adata = sc.read(self.matrixpath, cache = cache, delimiter = delimiter)
        
        if transpose:
            self.adata = self.adata.transpose()
    
    # generate a violin and scatter plot for the x and y variables passed
    # Usually used for total_counts vs. n_genes_by_counts
    def violin(self, x = 'total_counts', y = 'n_genes_by_counts'):
        sc.pp.calculate_qc_metrics(self.adata, percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(self.adata, [y, x],
             jitter=0.4, multi_panel=True)
        sc.pl.scatter(self.adata, x=x, y=y)
    
    # Filter out the data by minimum genes per cell and minimum cells per gene
    def cellgene_filter(self, min_genes=100, min_cells=20):
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
    
    # Normalize the data by removing cells with too many counts
    # Then scale counts per cell to the same value
    def normalize(self, max_n_genes_by_counts = 7000, target_sum = 1e4):
        self.adata = self.adata[self.adata.obs.n_genes_by_counts < max_n_genes_by_counts, :]
        sc.pp.normalize_total(self.adata, target_sum = target_sum)
        sc.pp.log1p(self.adata)
    
    # Filter highly variable genes using mean and dispersion, saving a raw copy if needed
    def variable_filter(self, min_mean = 0.0125, max_mean = 3, min_disp = 0.1, max_disp = 10):
        sc.pp.highly_variable_genes(self.adata, min_mean = min_mean, max_mean=max_mean, min_disp=min_disp, max_disp=max_disp)
        sc.pl.highly_variable_genes(self.adata)
        
        self.adata.raw = self.adata
        self.adata = self.adata[:, self.adata.var.highly_variable]
    
    # Regress out total counts from the analysis and scale
    def regress_scale(self, how = ['total_counts'], max_value = 10):
        sc.pp.regress_out(self.adata, how)
        sc.pp.scale(self.adata, max_value = max_value)
    
    # Perform PCA, displaying the first two PCs and the PCA variance ratio
    def pca_basic(self, svd_solver='arpack'):
        sc.tl.pca(self.adata, svd_solver=svd_solver)
        sc.pl.pca(self.adata, save = self.species_prefix + self.datatype + '_pca.pdf')
        sc.pl.pca_variance_ratio(self.adata, log=True)
    
    # Perform UMAP and then Leiden clustering using parameters passed in
    # Save the plot by default
    def umap_leiden(self, n_neighbors=50, n_pcs = 40, legend_loc='on data', save = True):
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata)
        if save:
            sc.pl.umap(self.adata, color=['leiden'], legend_loc=legend_loc, save = '_'.join([self.species_prefix, self.datatype, 'leiden.pdf']))
        else:
            sc.pl.umap(self.adata, color=['leiden'], legend_loc=legend_loc)
    
    # Rank the genes by their importance for leiden clusters using a t-test
    # Parameters can be changed to alter the ranking scheme
    def rank_genes(self, on = 'leiden', method = 't-test', plot = True, n_genes = 25, sharey = False):
        sc.tl.rank_genes_groups(self.adata, on, method=method)
        sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=sharey)
        self.adata.write(self.resultsfile)
    
    # Get the top_number genes per cluster across all clusters without repetition
    # Also get the top marker gene per cluster
    # Save these values as object attributes with a special key
    # Return these keys, which can be used with getattr() to access the actual list of names
    # The datatype string is only used to name the output file; it's not cross-checked to make sure it's correct
    def get_top_genes(self, datatype: str, top_number = 200, tofile = True, export = True):
        from biofile_handling import GeneListFile
        
        self.marker_genes_df = pd.DataFrame(self.adata.uns['rank_genes_groups']['names'])
        top_genelist = list(set(chain.from_iterable([i for i in self.marker_genes_df.iloc[0:top_number].values.tolist()])))
        top_marker_genes = list(set(self.marker_genes_df.iloc[0].values))
        
        genelisttype = '_'.join(['top', str(top_number) + 'DE', datatype, 'list'])
        markergenetype = '_'.join(['top', 'marker', datatype, 'list'])
        setattr(self, genelisttype, top_genelist)
        setattr(self, markergenetype, top_marker_genes)
        
        return genelisttype, markergenetype
    
    # Export the top genes list to a text file
    # Generates a GeneListFile object and assigns that as an object attribute
    def export_top_genes(self, key):
        if '_list' not in key:
            raise Exception('key must be an attribute ending in "_list"')
        
        top_type = key.split('_')[1]
        datatype = key.split('_')[2]
        
        top_genelist_filename = '_'.join([self.species_prefix, self.datatype, 'top', top_type, datatype, 'ids.txt'])
        top_genelist_file = GeneListFile(top_genelist_filename, self.sampledict, self.matrix, getattr(self, key), datatype)
        setattr(self, key + '_file', top_genelist_file)
    
    # Using an IdmmFile object, convert IDs from one type to another
    # Usually you'll be converting from gene_name to the embedding (e.g. Orthogroup)
    # Prints a table of the mappings
    # Returns a list of ids for the ones actually present in the scanpy object
    def map_gene_to_id(self, idmm, gene_list: list, from_id: str, to_id: str):
        import pandas as pd
        idmm = pd.read_csv(idmm.path, sep = '\t')

        id_table = idmm[idmm[from_id].isin(gene_list)][[from_id, to_id]].drop_duplicates()
        display(id_table)
        
        ids = [i for i in id_table[to_id].values if i in list(self.adata.var.index)]
        
        return ids
