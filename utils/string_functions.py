# Function takes a string from the Ensembl additional fields and parses into dict
def convert_fields_to_dict_gtf(string):
    # Additional fields are demarcated with ';', splits them by ';'
    # Additional fields string ends in a terminal ';', causing problems later; this handles that
    entries = string.split(";")
    if entries[-1] == '':
        entries = entries[:-1]
    # Additional fields are irregularly padded with whitespace and shoudl be cleaned by this
    output = {entry.split('"')[0].replace(' ', ''): entry.split('"')[1] for entry in entries}
    return output

# Function takes a string from the Ensembl additional fields and parses into dict
def convert_fields_to_dict_gff(string):
    # Additional fields are demarcated with ';', splits them by ';'
    entries = string.split(";")
    # Additional fields are irregularly padded with whitespace and should be cleaned by this
    output = {entry.split('=')[0]: entry.split('=')[1] for entry in entries}
    return output

# Function takes a string from the Ensembl additional fields and parses into dict
def convert_dbxref_to_dict(string):
    # Additional fields are demarcated with ';', splits them by ';'
    entries = string.split(",")
    # Additional fields are irregularly padded with whitespace and should be cleaned by this
    output = {entry.split(':')[0]: entry.split(':')[1] for entry in entries}
    return output

def convert_dict_to_fields_gtf(dictionary):
    return '; '.join([key + ' "' + value + '"' for key,value in dictionary.items()])

def make_gene_list(lst, filename):
    # Stolen from Jase's mouse code, writes to file
    with open(filename, 'w') as fp:
        for item in lst:
        # write each item on a new line
            fp.write("%s\n" % item)
        print('Done')
    return True

def prefixify(species):
    return species.split('_')[0][0] + species.split('_')[1][0:3]