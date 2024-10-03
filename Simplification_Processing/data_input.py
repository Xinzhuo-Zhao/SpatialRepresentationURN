import geopandas as gpd
import pickle
import os


def save_testingFile(geo_list, filename):
    filelist = os.listdir()
    if 'result' not in filelist:
        os.makedirs('result')
    geos = list()
    for i in geo_list:
        geos.append(i)
    p = gpd.GeoDataFrame(geos, columns=['geometry'], crs='epsg:3857').reset_index()
    filename = filename + ".geojson"
    path = os.path.join(r"./result",filename)
    p.to_file(path)


def load_testingFile(filename):
    filename = filename + ".geojson"
    load_data = gpd.read_file(filename, crs='epsg:3857')
    return load_data.geometry.to_list()


def save_idx(idxs, storage, filename):
    filelist = os.listdir()
    if 'result' not in filelist:
        os.makedirs('result')
    geos = {'fid': list(), 'geometry': list()}
    if idxs.__class__ == dict:
        for i in idxs:
            if idxs[i] is None:
                continue
            geos['fid'].append(i)
            geos['geometry'].append(storage[i])
    elif idxs == list or set:
        for i in idxs:
            geos['fid'].append(i)
            geos['geometry'].append(storage[i])
    else:
        raise TypeError("Check the type of input value 'idxs'")
    p = gpd.GeoDataFrame.from_dict(geos, crs='epsg:3857')
    filename = filename + ".geojson"
    path = os.path.join(r"./result", filename)
    p.to_file(path)


def save_data_structure(data, filename='data.pkl'):
    """
    save the connectivity configurations
    :param data: 要保存的数据结构。
    :param filename: 保存的文件名，默认为'data.pkl'。
    """
    filelist = os.listdir()
    if 'result' not in filelist:
        os.makedirs('result')
    path = os.path.join(r"./result", filename)
    with open(path, 'wb') as file:
        pickle.dump(data, file)


def load_data_structure(filename='data.pkl'):
    """
    从本地文件中读取数据结构。
    :param filename: 要读取的文件名，默认为'data.pkl'。
    :return: 读取的数据结构。
    """
    try:
        with open(filename, 'rb') as file:
            data = pickle.load(file)
        return data
    except FileNotFoundError:
        print(f"The file {filename} does not exist.")
        return None
