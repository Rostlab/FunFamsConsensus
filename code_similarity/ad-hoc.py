import pandas as pd

path = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\bindPredict_performance.txt'
df = pd.read_csv(path, sep=' ', header=0, engine='python')
#print(df.head())

uniprot_ids = df['id']
print(len(uniprot_ids))

path2 = r'C:\Users\Linus\LRZ Sync+Share\UniversitätMünchen\Bioinformatik\6. Semester\Bachelorarbeit\evaluation_full_new.tsv'
data = pd.read_csv(path2, sep='\t', header=0, engine='python')
#print(data.head())

data = data[data['protein'].isin(uniprot_ids)]
print(data.shape)
print(data.mean())

