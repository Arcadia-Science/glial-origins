# Function takes a string from the GTF additional fields (column 8) and parses into dict
# For example, a string of format 'gene_name "Hh"; gene_id "123"' becomes
# {'gene_name': 'Hh', 'gene_id': '123'}
def convert_fields_to_dict_gtf(string):
    # Additional fields are demarcated with ';', splits them by ';'
    # Additional fields string ends in a terminal ';', causing problems later; this handles that
    entries = string.split(";")
    if entries[-1] == '':
        entries = entries[:-1]
    # Additional fields are irregularly padded with whitespace and should be cleaned by this
    output = {entry.split('"')[0].replace(' ', ''): entry.split('"')[1] for entry in entries}
    return output

# Function takes a string from the GFF additional fields (column 8) and parses into dict
# For example, a string of format 'gene_name=Hh;gene_id=123' becomes
# {'gene_name': 'Hh', 'gene_id': '123'}
def convert_fields_to_dict_gff(string):
    # Additional fields are demarcated with ';', splits them by ';'
    entries = string.split(";")
    entries = [entry for entry in entries if '=' in entry]
    # Additional fields are irregularly padded with whitespace and should be cleaned by this
    output = {entry.split('=')[0]: entry.split('=')[1] for entry in entries}
    return output

# Function takes a string from the Dbxref additional fields and parses into dict
# For example, a string of format 'Ensembl:123104,UniProtKB:Q7D56' becomes
# {'Ensembl': '123104', 'UniProtKB': 'Q7D56'}
def convert_dbxref_to_dict(string):
    if string == None:
        return None
    # Additional fields are demarcated with ',', splits them by ','
    entries = string.split(",")
    # Key-value pairs are specified with a ':', split into dictionary based on that sep
    output = {entry.split(':')[0]: entry.split(':')[-1] for entry in entries}
    return output

# Function takes a dictionary and returns an additional-fields string for GTF
# For example, a dictionary {'gene_name': 'Hh', 'gene_id': '123'} becomes
# 'gene_name "Hh"; gene_id "123"'
def convert_dict_to_fields_gtf(dictionary):
    return '; '.join([key + ' "' + value + '"' for key,value in dictionary.items()])

# Function takes a dictionary and returns an additional-fields string for GTF
# For example, a dictionary {'gene_name': 'Hh', 'gene_id': '123'} becomes
# 'gene_name=Hh;gene_id=123'
def convert_dict_to_fields_gff(dictionary):
    return ';'.join([str(key) + '=' + str(value) for key,value in dictionary.items()])

# Makes a gene list .txt file from a python list
# Each item of the list is written to a new line
# For example, a list ['Hh', 'bcd', 'nos'] becomes a text file:
# Hh
# bcd
# nos
def make_gene_list(lst, filename):
    # Creates a file with name `filename`
    with open(filename, 'w') as fp:
        for item in lst:
        # write each item on a new line
            fp.write("%s\n" % item)
        print('Wrote', len(lst), 'gene ids to', filename)
    return True

# Converts a string of format "Genus_species" to "Gspe" species prefix format
# Takes the first letter of "Genus" -> G
# Takes the first three letters of "species" -> spe
# Concatenates -> Gspe
def prefixify(species):
    try:
        output = species.split('_')[0][0] + species.split('_')[1][0:3]
        return output
    except IndexError:
        return species