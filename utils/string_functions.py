def convert_fields_to_dict_gtf(string: str):
    """ Takes a `string` from the GTF additional fields (column 8) and parses into a `dict`.  
    Tolerates strings with and without a terminal ';' separator.

    Args:
        string (str): String in the GTF additional fields format, as below.  
            `'feature1 "entry1"; feature2 "entry2";'`
    
    Returns:
        (dict): Unpacked string into a dict with key:value format, as below.  
            `{'feature1': 'entry1', 'feature2': 'entry2'}`
    
    Examples:
        >>> convert_fields_to_dict_gtf('gene_name "Hh"; gene_id "123"')
        {'gene_name': 'Hh', 'gene_id': '123'}
    """
    # Additional fields are demarcated with ';', splits them by ';'
    # Additional fields string ends in a terminal ';', causing problems later; this handles that
    entries = string.split(";")
    if entries[-1] == '':
        entries = entries[:-1]
    # Additional fields are irregularly padded with whitespace and should be cleaned by this
    output = {entry.split('"')[0].replace(' ', ''): entry.split('"')[1] for entry in entries}
    return output

def convert_fields_to_dict_gff(string: str):
    """Takes a `string` from the GFF additional fields (column 8) and parses into a `dict`.

    Args:
        string (str): String in the GFF additional fields format, as below.  
            `'feature1=entry1;feature2=entry2'`
    
    Returns:
        (dict): Unpacked string into a dict with key:value format, as below.  
            `{'feature1': 'entry1', 'feature2': 'entry2'}`
    
    Examples:
        >>> convert_fields_to_dict_gff('gene_name=Hh;gene_id=123')
        {'gene_name': 'Hh', 'gene_id': '123'}
    """

    # Additional fields are demarcated with ';', splits them by ';'
    entries = string.split(";")
    entries = [entry for entry in entries if '=' in entry]
    # Additional fields are irregularly padded with whitespace and should be cleaned by this
    output = {entry.split('=')[0]: entry.split('=')[1] for entry in entries}
    return output

def convert_dbxref_to_dict(string: str):
    """Takes a `string` from the Dbxref additional fields and parses into `dict`.
    
    Args:
        string (str): String in the Dbxref format, as below.  
            `'feature1:entry1,feature2:entry2'`
    
    Returns:
        (dict): Unpacked string into a dict with key:value format, as below.  
            `{'feature1': 'entry1', 'feature2': 'entry2'}`

    Examples:
        >>> convert_dbxref_to_dict('Ensembl:123104,UniProtKB:Q7D56')
        {'Ensembl': '123104', 'UniProtKB': 'Q7D56'}
    """
    if string == None:
        return None
    # Additional fields are demarcated with ',', splits them by ','
    entries = string.split(",")
    # Key-value pairs are specified with a ':', split into dictionary based on that sep
    output = {entry.split(':')[0]: entry.split(':')[-1] for entry in entries}
    return output

def convert_dict_to_fields_gtf(dictionary: dict):
    """Takes a `dict` and converts it to a `string` for the additional fields section of GTF format.
    
    Args:
        dictionary (dict): a dict with key:value format, as below.  
            `'feature1 "entry1"; feature2 "entry2"'`
    
    Returns:
        (str): String in GTF additional fields format, as below:  
            `'feature1=entry1;feature2=entry2'`

    Examples:
        >>> convert_dict_to_fields_gtf({'gene_name': 'Hh', 'gene_id': '123'})
        'gene_name "Hh"; gene_id "123"'
    """
    return '; '.join([key + ' "' + value + '"' for key,value in dictionary.items()])

# Function takes a dictionary and returns an additional-fields string for GTF
# For example, a dictionary {'gene_name': 'Hh', 'gene_id': '123'} becomes
# 'gene_name=Hh;gene_id=123'
def convert_dict_to_fields_gff(dictionary: dict):
    """Takes a `dict` and converts it to a `string` for the additional fields section of GFF format.

    Args:
        dictionary (dict): a dict with key:value format, as below.  
            `'feature1 "entry1"; feature2 "entry2"'`
    
    Returns:
        (str): String in GFF additional fields format, as below:  
            `'feature1=entry1;feature2=entry2'`

    Examples:
        >>> convert_dict_to_fields_gff({'gene_name': 'Hh', 'gene_id': '123'})
        'gene_name=Hh;gene_id=123'
    """
    return ';'.join([str(key) + '=' + str(value) for key,value in dictionary.items()])

def make_gene_list(lst: list, filename: str):
    """Makes a gene list file using a `list` of genes and a `string` filename or path.  
    File will display a list of genes, one per line.

    Args:
        lst (list): list of gene ids to be saved.
        filename (str): name or path to destination file.

    Examples:
        >>> make_gene_list(['Hh', 'bcd', 'nos'], 'fly_genes.txt')
        >>> with open('fly_genes.txt', 'r') as file:
        >>>     for line in file:
        >>>         print(line)
        'Hh'
        'bcd'
        'nos'
    """

    # Creates a file with name `filename`
    with open(filename, 'w') as fp:
        for item in lst:
        # write each item on a new line
            fp.write("%s\n" % item)
        print('Wrote', len(lst), 'gene ids to', filename)

# 
# Takes the first letter of "Genus" -> G
# Takes the first three letters of "species" -> spe
# Concatenates -> Gspe
def prefixify(species:str):
    """Converts a  `string` of format "Genus_species" to "Gspe" species prefix format.  

    Note:
        If the `species` string lacks an underscore, will return the entire original string.

    Args:
        species (str): name of species as `Genus_species`
    
    Returns:
        (str): species prefix, e.g. `Gspe`

    Examples:
        >>> prefixify('Genus_species')
        'Gspe'
        >>> prefixify('MmusDrerXlae')
        'MmusDrerXlae'
    """
    try:
        output = species.split('_')[0][0] + species.split('_')[1][0:3]
        return output
    except IndexError:
        return species