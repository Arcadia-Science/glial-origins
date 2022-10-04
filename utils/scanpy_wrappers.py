from biofile_handling import SampleDict
from biofile_handling import GeneListFile
from string_functions import prefixify
from string_functions import make_gene_list
import scanpy as sc
import pandas as pd
from itertools import chain

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
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
    
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
    
    def umap_leiden(self, n_neighbors=50, n_pcs = 40, legend_loc='on data', save = True):
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(self.adata)
        sc.tl.leiden(self.adata)
        if save:
            sc.pl.umap(self.adata, color=['leiden'], legend_loc=legend_loc, save = '_'.join([self.species_prefix, self.datatype, 'leiden.pdf']))
        else:
            sc.pl.umap(self.adata, color=['leiden'], legend_loc=legend_loc)
    
    def rank_genes(self, on = 'leiden', method = 't-test', plot = True, n_genes = 25, sharey = False):
        sc.tl.rank_genes_groups(self.adata, on, method=method)
        sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=sharey)
        self.adata.write(self.resultsfile)
        
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
    
    def export_top_genes(self, key):
        if '_list' not in key:
            raise Exception('key must be an attribute ending in "_list"')
        
        top_type = key.split('_')[1]
        datatype = key.split('_')[2]
        
        top_genelist_filename = '_'.join([self.species_prefix, self.datatype, 'top', top_type, datatype, 'ids.txt'])
        top_genelist_file = GeneListFile(top_genelist_filename, self.sampledict, self.matrix, getattr(self, key), datatype)
        setattr(self, key + '_file', top_genelist_file)
    
    def map_gene_to_id(self, idmm, gene_list: list, from_id: str, to_id: str):
        import pandas as pd
        idmm = pd.read_csv(idmm.path, sep = '\t')

        id_table = idmm[idmm[from_id].isin(gene_list)][[from_id, to_id]].drop_duplicates()
        display(id_table)
        
        ids = [i for i in id_table[to_id].values if i in list(self.adata.var.index)]
        
        return ids