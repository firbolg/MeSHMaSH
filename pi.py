from Bio import Entrez
from Bio import Medline
import pandas as pd
import time
import re
import numpy as np
from ratelimiter import RateLimiter
import pubmed_parser as pp

Entrez.email = 'levijdolan@gmail.com'
api_key = 'c34c74b04d76cdf36ea6835f39e191a3a108'


@RateLimiter(max_calls=10, period=1)
def getPMIDs (term):
    handle_test = Entrez.esearch(db='pubmed', term=term+'[ti]', retmax = 3000)
    result_test = Entrez.read(handle_test)
    handle_test.close()
    pmid = result_test['IdList']
    return pmid


@RateLimiter(max_calls=10, period=1)
def returnPMRecords (ids):
    tester = []
    rate_limiter = RateLimiter(max_calls=3, period=1)

    for x in ids :
        with rate_limiter:
            tester.append(pp.parse_xml_web(x, save_xml=False))
    return tester


#For a given input term, get all previous indexing MeSH terms as a list
@RateLimiter(max_calls=10, period=1)
def getPIs(term):
    word_lst = []
    ids = getPMIDs(term)
    for i in ids:
        rec = returnPMRecords(i)
        mesh_lst = [i['keywords'] for i in rec]
        for i in mesh_lst:
            words = re.split(':|;', i)
            words = words[1::2]
            word_lst = word_lst + words
        return word_lst


#Create a dictionary from previous indexing MeSH terms with term as key, total count of term occurrences as value
@RateLimiter(max_calls=10, period=1)
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
#Create new dictionary with top previous indexing terms 1-15, sorted highest to lowest by percentage of previous indexing
#Break ties in percentages by sorting alphabetically by previous indexing term
@RateLimiter(max_calls=10, period=1)
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
                      'PI 10':PI_k[9], 'PI 10 %':PI_v[9], 'PI 11':PI_k[10], 'PI 11 %':PI_v[10], 'PI 12':PI_k[11], 'PI 12 %':PI_v[11], 
                      'PI 13':PI_k[12], 'PI 13 %':PI_v[12], 'PI 14':PI_k[13], 'PI 14 %':PI_v[13], 'PI 15':PI_k[14], 'PI 15 %':PI_v[14],
                      'Total PIs':total}
    return PI_sorted_dict 


#Create a dataframe composed of dictionaries of previous indexing terms and percentages from a given list
@RateLimiter(max_calls=10, period=1)
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

df.to_csv('SP and SDOH.csv', index=False)

