import os
import sys
import argparse
from pymongo import MongoClient

parser = argparse.ArgumentParser(description='push genic values to the database')
parser.add_argument('-bed', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='bed file with 7 columns (6 regular + 1 value column)')
parser.add_argument('-config', required=False, default= "mongoconfig.json", help='mongo config file (json) with "key" and the link')
parser.add_argument('-genome', required= True, help='Genome version')
parser.add_argument('-title', required=True, type=str, help='title')

args = parser.parse_args()


def genome2organism(genome):
    genome2organismDict = {
        'hg19': 'human',
        'tair10': 'Arabidopsis',
        'mm10': 'mouse'
    }
    return genome2organismDict[genome]

if os.path.isfile(args['config']):
    connection = MongoClient(json.load(open(args['config']))['key'])
    db = connection.data
    collection = db['gene']

    for line in args.bed:
        ll = line.split('\t')
        chromosome = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        name = ll[3]
        score = ll[4]
        strand = ll[5]
        value = float(ll[6])
        row = {
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'name': name,
            'score': score,
            'strand': strand,
            'value': value,
            'genome': args.genome,
            'organism': genome2organism(args.genome)
        }

        cursor = collection.find({'$and':[{"genome":row['genome']}, {"name": row['name']}]})

        if(cursor.count() == 0): #insert entry if doesn't exist in db
            
            iterData = dict.copy(row)
            for key in iterData.keys():
                if(iterData[key] is None):
                    row.pop(key,None)

            collection.insert(row)
        elif(cursor.count() > 1): #there should only be 1 document for each genome/gene pair
            raise ValueError("Duplicate document found. More than one found for Gene " + row['name'] + " for " + args.genome)

        for oldDoc in cursor: 
            for key in row.keys():
        
                oldDoc[key] = row[key]
                uqId = oldDoc['_id']
                if(row[key] is None):
                    oldDoc.pop(key, None) #get rid of null fields
                else:
                    #update using document's unique ID
                    collection.update({"_id": uqId}, {"$set":{key:oldDoc[key]}})
else:
    print("No such file " + args['config'] + ". Skipping this script.")