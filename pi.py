from Bio import Entrez
from Bio import Medline
import pandas as pd
import re
import pubmed_parser as pp

Entrez.email = '***@gmail.com'
api_key = '***'


def getMEDLINE(term):
    handle_test = Entrez.esearch(db='pubmed', term=term+'[ti]', retmax = 3000)
    result_test = Entrez.read(handle_test)
    handle_test.close()
    pmid = result_test['IdList']
    recs = []
    for i in pmid:
        handle = Entrez.efetch(db='pubmed', id=i, rettype='medline', retmode='text')
        r = handle.read()
        recs.append(r)
    return recs


#For a given input term, get all previous indexing MeSH terms as a list from MEDLINE
def getPIs(term):
    recs = getMEDLINE(term)
    mesh = 'MH  -'
    PI_recs = [i for i in recs if mesh in i] #check only articles with MeSH for MeSH terms
    PI_recs = [i.split('\n') for i in recs]
    combined = [i for sublist in PI_recs for i in sublist]
    mesh = [i for i in combined if i.startswith('MH  - ')]
    mesh = [i.split('/') for i in mesh]
    PIs = []
    for i in mesh:
        for w in i:
            if w.startswith('MH  - *'):
              w = (w[7:])
              PIs.append(w)
            if w.startswith('MH  - '):
              w = (w[6:])
              PIs.append(w)
            else:
              pass
    return PIs


#Create a dictionary from previous indexing MeSH terms with term as key, total count of term occurrences as value
def PI_dict(term):
    dict = {}
    PIs = getPIs(term)
    for i in PIs:
        if i in dict:
            dict[i] += 1
        else:
            dict[i] = 1      
    return dict


#Convert previous indexing counts to percentages of previous indexing terms total
#Create new dictionary with top previous indexing terms 1-10, sorted highest to lowest by percentage of previous indexing
#Break ties in percentages by sorting alphabetically by previous indexing term
def PI_data(term):
    PIs = PI_dict(term)
    PIs = dict(sorted(PIs.items(), reverse=True, key=lambda t: t[::-1]))
    total = sum(PIs.values())

    for key, value in PIs.items():
        pct = value * 100 / total
        PIs[key] = round(pct, 2)

    PI_k = sorted(PIs, key=lambda x: (-PIs[x], x))
    PI_v = list(PIs.values())

    PI_sorted_dict = {'Prioritized New Concept':term,
                      'PI 1':PI_k[0], 'PI 1 %':PI_v[0], 'PI 2':PI_k[1], 'PI 2 %':PI_v[1], 'PI 3':PI_k[2], 'PI 3 %':PI_v[2],
                      'PI 4':PI_k[3],'PI 4 %':PI_v[3], 'PI 5':PI_k[4], 'PI 5 %':PI_v[4], 'PI 6':PI_k[5], 'PI 6 %':PI_v[5], 
                      'PI 7':PI_k[6], 'PI 7 %': PI_v[6], 'PI 8':PI_k[7], 'PI 8 %':PI_v[7], 'PI 9':PI_k[8], 'PI 9 %':PI_v[8], 
                      'PI 10':PI_k[9], 'PI 10 %':PI_v[9], 
                      'Total PIs':total}
    return PI_sorted_dict 


#Create a dataframe composed of dictionaries of previous indexing terms and percentages from a given list
def PI_report(terms):
    dict_list = []
    for i in terms:
        dict = PI_data(i)
        dict_list.append(dict)
    df = pd.DataFrame(dict_list)
    return df


df_source = pd.read_csv('SP and SDOH with Postings.csv')

newterms = df_source['Prioritized New Concept'].tolist()

df_report = PI_report(newterms)

df = pd.concat([df_source, df_report], axis=1)

df.to_csv('SP and SDOH with Postings and PIs.csv', index=False)



