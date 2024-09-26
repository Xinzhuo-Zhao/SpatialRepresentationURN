from Graph_structure import *
import warnings
warnings.simplefilter("ignore")
#%%
def Crossing_Checking(data, geo_key):
    '''
    This function Generate
    :param data: dict road_idx: "node": node_geos(list), "geometry": road_geos(shapely)
    :param geo_key:
    :return: separateLi(list)
    :return: separateMap(dict) 原有的geo_keys向separateLi的idx映射
    '''
    separatLi = list()
    separatDict = dict()

    for li in tqdm(range(len(geo_key))):
        rd_idx = geo_key[li]
        rd_geo = data[rd_idx]['geometry']
        cross_pt = data[rd_idx]['nodes']
        '''
        所有点入列
        '''
        if len(cross_pt) == 0:
            continue
        pt_que = Queue()
        rd_list = list()
        rd_list.append(rd_geo)
        for pt in cross_pt:
            pt_que.Que_in(pt)
        '''
        按照类型处理每一个点跟线的情况，
        until all pts are used
        '''
        while not pt_que.Que_isEmpty():
            current_pt = pt_que.Que_out()
            l_idx = 0
            while l_idx < len(rd_list):
                current_li = rd_list[l_idx]
                dist_to_curve = check_pt_on_curve(current_pt, current_li)
                #   开始检测currentli和currentpt之间的关系
                if dist_to_curve > 1:
                    #   如果距离大于1m，说明跟这个线段是无关的，直接找下一条线
                    l_idx += 1
                    continue
                if check_endpt(current_pt, current_li):
                    #   检查p是不是endpt
                    #   是：  对线不做任何处理，直接算下一个点
                    #   不是: 对线进行split
                    break
                else:
                    #   是切割点且找到了需要被切割的线
                    #   但是点不一定直接在线上
                    if dist_to_curve == 0:
                        #   正好在线上
                        break_pt = current_pt
                        spt_list = split_polyline(break_pt, current_li)
                        rd_list.pop(l_idx)
                        rd_list = rd_list + spt_list
                        break
                    else:
                        #   不在线上
                        break_pt = nearest_pt(current_pt, current_li)
                        spt_list = split_polyline(break_pt, current_li)
                        #   有可能找不到最近点

                        rd_list.pop(l_idx)
                        rd_list = rd_list + spt_list
                        break
                # l_idx += 1
        '''
        点已经完全分割了线
        线存储在
        转移到separatLi里
        并保留映射关系
        '''
        for i in rd_list:
            idx = len(separatLi)
            separatLi.append(i)
            separatDict[idx] = li

    return separatLi, separatDict
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
# 清理出想要的线型
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
