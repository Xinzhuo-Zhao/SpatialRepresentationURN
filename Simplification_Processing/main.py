from Graph_structure import *
import warnings
warnings.simplefilter("ignore")
#%%
fileaddress = "./TestFile/TestFile.shp"
savename = 'testings'
data = gpd.read_file(fileaddress).to_crs('epsg:3857')
type = data['highway'].value_counts().to_dict()
desired_type = [
    'cycleway',
    'living_street',
    'primary',
    'raceway',
    'residential',
    'rest_area',
    'road',
    'secondary',
    'service',
    'tertiary',
    'track',
    'trunk',
    'unclassified',
    'tertiary_link',
    'motorway'
]

data = data[data['highway'].apply(lambda x: True if x in desired_type else False)]
data = data[['OBJECTID', 'geometry']]
data['nodes'] = data.apply(lambda x: [], axis=1)
data = data[data['geometry'].notna()]
data = data.explode()
data = data.set_index('OBJECTID').to_dict('index')


geo_keys = [i for i in data]
for o in tqdm(range(len(geo_keys))):
    for d in range(o + 1, len(geo_keys)):
        line1_idx = o
        line2_idx = d
        line1 = data[geo_keys[o]]['geometry']
        line2 = data[geo_keys[d]]['geometry']
        if line1.intersects(line2):
            intersectionPT = line1.intersection(line2)
            if intersectionPT.type == 'MultiPoint':
                for i in intersectionPT.geoms:
                    data[geo_keys[o]]['nodes'].append(i)
                    data[geo_keys[d]]['nodes'].append(i)
            elif intersectionPT.type == 'Point':
                data[geo_keys[o]]['nodes'].append(intersectionPT)
                data[geo_keys[d]]['nodes'].append(intersectionPT)

separatLi, _ = Crossing_Checking(data, geo_keys)

linkstg = []
for i in range(len(separatLi)):
    linkstg.append(separatLi[i])


urn_graph = GraphBuilder(linkstg,dict(), dict(), list(),savename)
urn_graph.run()
urn_graph.save()
