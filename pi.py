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
def getPMIDs (term) :
    handle_test = Entrez.esearch(db='pubmed', term=term+'[ti]', retmax = 1000)
    result_test = Entrez.read(handle_test)
    handle_test.close()
    pmid = result_test['IdList']
    return pmid

@RateLimiter(max_calls=10, period=1)
def returnPMRecords (ids) :
    tester = []
    rate_limiter = RateLimiter(max_calls=3, period=1)

    for x in ids :
        with rate_limiter:
            tester.append(pp.parse_xml_web(x, save_xml=False))
    return tester

newterm = 'Federal Poverty Level' #TODO insert list of proposed terms

ids = getPMIDs(newterm)

#TODO create dict with all MeSH as keys, total counts as values
#take a PMID, get MeSH
#check empty dict for each MeSH term
#if no key, create key
#value +1

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

def PI_dict(term):
    dict = {}
    PIs = getPIs(term)
    for i in PIs:
        if i in dict:
            dict[i] += 1
        else:
            dict[i] = 1
    return dict

#print(PI_dict(newterm))

#TODO create dataframe with all desired outputs
#create dataframe with cols: #1 word/%/#2 word/%/#3 word/%/#4 word/%/#5 word/%...#10 word/%/total count of MeSH occurrences 
#create new row, but newterm in first column
#look through dict, count all values, send this total to total count column
#for each key, take value, divide by total count
#express this as percentage and replace value for each key with this percentage
#sort from greatest to smallest by value
#for first 25 keys, take key, check if next dataframe column empty
#if empty, put key there
#put value in next column
#take Dan's csv, to df, then concat dataframes and output as new csv

def PI_report(term):
    column_names = ['Proposed Term','1','%','2','%','3','%','4','%','5','%','6','%','7','%','8','%','9','%','10','%','11','%','12','%','13','%','14','%','15','%','Total PIs']
    df = pd.DataFrame(columns = column_names)

    PIs = PI_dict(term)
    PIs = dict(sorted(PIs.items(), key=lambda item: item[1], reverse=True))
    #PIs = dict(sorted(PIs.items(), reverse=True, key=lambda t: t[::-1]))
    total = sum(PIs.values())
    for key, value in PIs.items():
        pct = value * 100 / total
        PIs[key] = round(pct, 2)
    return PIs, total

print(PI_report(newterm))



#TODO create visualizations of PIs
#for each row of previous dataframe:
#use plotly to make a pie chart with top 5 #word column values as labels, % column values as pie slices, whole chart labelled with index/newterm
#create x,y axis of x=total Title/Abstract of newterm  y=total MeSh terms
#put pie charts on this graph as bubbles; when you hover over it, expand 
#underneath, print larger pie charts alphabetically with x,y info underneath
#output as .html

