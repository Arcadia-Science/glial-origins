from biofile_handling import SampleDict
import scanpy as sc

class ScanpyMetaObject():
    def __init__(self, matrix, sampledict: SampleDict):
        from string_functions import prefixify
        
        self.matrix = matrix
        self.sampledict = sampledict
        self.matrixname = self.matrix.filename
        self.species = self.sampledict.species
        self.directory = self.sampledict.directory
        self.conditions = self.sampledict.conditions
        self.species_prefix = prefixify(self.species)
        self.matrixpath = self.directory + self.matrixname
        self.datatype = str(type(self.matrix)).split('.')[1].split('File')[0]
        
        self.resultsfile = self.matrixname + '_scresults.h5ad'
        self.resultspath = self.directory + self.resultsfile
    
    def read(self, delimiter = '\t', cache = True, transpose = True):
        self.adata = sc.read(self.matrixpath, cache = cache, delimiter = delimiter)
        
        if transpose:
            self.adata = self.adata.transpose()
            
    def violin(self):
        sc.pp.calculate_qc_metrics(self.adata, percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(self.adata, ['n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True)
        sc.pl.scatter(self.adata, x='total_counts', y='n_genes_by_counts')
        
    def cellgene_filter(self, min_genes=100, min_cells=20):
        sc.pp.filter_cells(self.adata, min_genes=100)
        sc.pp.filter_genes(self.adata, min_cells=20)
    
    def normalize(self, max_n_genes_by_counts = 7000, target_sum = 1e4):
        self.adata = self.adata[self.adata.obs.n_genes_by_counts < max_n_genes_by_counts, :]
        sc.pp.normalize_total(self.adata, target_sum = target_sum)
        sc.pp.log1p(self.adata)
    
    def variable_filter(self, min_mean = 0.0125, max_mean = 3, min_disp = 0.1, max_disp = 10):
        sc.pp.highly_variable_genes(self.adata, min_mean = min_mean, max_mean=max_mean, min_disp=min_disp, max_disp=max_disp)
        sc.pl.highly_variable_genes(self.adata)
        
        self.adata.raw = self.adata
        self.adata = self.adata[:, self.adata.var.highly_variable]
        
    def regress_scale(self, how = ['total_counts'], max_value = 10):
        sc.pp.regress_out(self.adata, how)
        sc.pp.scale(self.adata, max_value = max_value)
    
    def pca_basic(self, svd_solver='arpack'):
        sc.tl.pca(self.adata, svd_solver=svd_solver)
        sc.pl.pca(self.adata, save = self.species_prefix + self.datatype + '_pca.pdf')
        sc.pl.pca_variance_ratio(self.adata, log=True)
    
    def umap_leiden(self, n_neighbors=50, n_pcs = 40, legend_loc='on data'):
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs = n_pcs)
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata)
        sc.pl.umap(self.adata, color=['leiden'], legend_loc=legend_loc, save = '_'.join([self.species_prefix, self.datatype, 'leiden.pdf']))