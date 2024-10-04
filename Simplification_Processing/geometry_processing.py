from shapely.ops import split, nearest_points
import shapely
from tqdm import tqdm
from numpy.linalg import norm
from data_structure import Queue
from Angle import *
from math import pi


def check_endpt(pt, li):
    end1, end2 = endPts(li)
    if end1.distance(pt) < 1 or end2.distance(pt) < 1:
        return True
    else:
        return False


def endPts(s):
    end1 = Point(s.xy[0][0], s.xy[1][0])
    end2 = Point(s.xy[0][-1], s.xy[1][-1])
    return end1, end2


def compute_box(data):
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    return q1, q3, iqr, lower, upper


def normalize_vector(vector):
    length = np.linalg.norm(vector, 1)
    if length == 0:
        return vector
    return vector / length


def AddPoint(line_list):
    li1 = line_list[0]
    li2 = line_list[1]
    pts = []
    for i in range(len(li1.xy[0])):
        pts.append((li1.xy[0][i], li1.xy[1][i]))
    pts.append((li2.xy[0][-1], li2.xy[1][-1]))
    return LineString(pts)


def cosine_similarity(vec1, vec2):
    cos_sim = vec1.dot(vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    return cos_sim


def midPt(s):
    return Point(((s.xy[0][-1] + s.xy[0][0]) / 2, (s.xy[1][-1] + s.xy[1][0]) / 2))



def LineStringVec(pts):
    if pts.__class__ == np.ndarray:
        pass
    elif pts.__class__ == shapely.geometry.linestring.LineString:
        pts = np.array(pts.xy)
    vec = np.diff(pts).T
    vec_length = norm(vec, ord=2, axis=1)
    angle = np.array([cosine_similarity(vec[i], vec[i + 1]) for i in range(len(vec) - 1)])
    return vec, vec_length, angle


def check_pt_on_curve(pt, inputli):
    return inputli.distance(pt)


def split_polyline(p, l):
    sp_list = split(l, p)
    li_list = list()
    #   Check whether the l is split successfully
    if len(sp_list.geoms) < 2:
        #   not split successfully
        bk = nearest_pt(p, l)
        sp_list = split(l, bk)
    if len(sp_list.geoms) < 2:
        #   cannot use nearest pt to achieve the split mission
        #   manually split using intermediate pts
        shortest_distance = 1
        close_seg = None
        for i in range(len(l.xy[0]) - 1):
            curli = LineString([(l.xy[0][i], l.xy[1][i]), (l.xy[0][i + 1], l.xy[1][i + 1])])
            if curli.distance(p) < shortest_distance:
                close_seg = curli
                shortest_distance = curli.distance(p)
        if close_seg == None:
            raise ValueError("未知问题1:找不到符合标准的线段")
            return l
        vector_x = (close_seg.xy[0][0] - close_seg.xy[0][-1]) / close_seg.length
        vector_y = (close_seg.xy[1][0] - close_seg.xy[1][-1]) / close_seg.length
        verticalVec = (vector_y * 5, -vector_x * 5)
        crossing = LineString(
            [(p.x + verticalVec[0], p.y + verticalVec[1]), (p.x - verticalVec[0], p.y - verticalVec[1])])
        sp_list = split(l, crossing)

    if len(sp_list.geoms) < 2:
        raise ValueError("未知问题2:仍然无法切割线段")
        return l
    for i in sp_list.geoms:
        li_list.append(i)
    return li_list


def nearest_pt(p, l):
    brk, _ = nearest_points(l, p)
    return brk


def get_closest_nd(nd, target_nd):
    '''
    this function find the closest nd of nd among target_nd
    :param nd:
    :param target_nd: candidating select nodes
    :return: closest_idx: the idx of closest nd in target_nd
    :return: target_nd[closest_idx]: the geometry of the closest pt
    '''
    closest_dist = None
    closest_idx = None
    for n_idx, nd_geom in enumerate(target_nd):
        if len(target_nd) > 3:
            if n_idx == 0 or n_idx == (len(target_nd) - 1):  # ignore the start and end node
                continue
        if closest_idx is None or closest_dist is None:
            closest_dist = nd.distance(nd_geom)
            closest_idx = n_idx
        new_dist = nd.distance(nd_geom)
        if new_dist < closest_dist:
            closest_idx = n_idx
            closest_dist = new_dist
    return closest_idx, target_nd[closest_idx]


def get_mid_nd(n1, n2):
    n1 = np.array(n1.coords[0])
    n2 = np.array(n2.coords[0])
    return Point((n1 + n2) / 2)


def separate_lines(s):
    li_list = []
    pt_list = Points_in_Line(s)
    for i in range(len(pt_list) - 1):
        li_list.append(LineString([(pt_list[i].x, pt_list[i].y), (pt_list[i + 1].x, pt_list[i + 1].y)]))
    return li_list


def Points_in_Line(s):
    pt_list = []
    for i in range(len(s.xy[0])):
        pt_list.append(Point(s.xy[0][i], s.xy[1][i]))
    return pt_list


def LineSimplication(s, l_id):
    '''
    line_Simplication only simplify the intermediate pts of segment geometry. Reduce the description and storage burden
    of geometric information
    :param l_id:
    :param s:   LineString
    :return res_lines:  store all lines shapely
    :return res_node:  store all nodes shapely
    :return node_map:  relations between links and nodes(nodes to links)
    :return bool: False is for unlinked lines
    '''
    li_list = separate_lines(s)
    li_list = [i for i in li_list if i.length > 0]
    res_lines = list()
    node_map = dict()
    res_node = list()

    times = 0
    if len(li_list) < 2:
        return None, None, None, False

    res_node.append(Point(s.xy[0][0], s.xy[1][0]))
    node_map[0] = {0}
    while len(li_list) > 1:
        li1 = li_list[0]
        li2 = li_list[1]
        angle = AngleCal(li1, li2)

        if angle < 0.8:

            nd_idx = len(res_node)
            li_idx = len(res_lines)
            current_crossing = Point(li1.xy[0][-1], li1.xy[1][-1])
            node_map[nd_idx] = {li_idx, li_idx + 1}
            res_node.append(current_crossing)
            res_lines.append(li_list[0])
            li_list.pop(0)
        elif (li1.length + li2.length) < 50:
            li_list[0] = LineString([(li1.xy[0][0], li1.xy[1][0]), (li2.xy[0][1], li2.xy[1][1])])
            li_list.pop(1)
        else:
            li_list[0] = AddPoint([li1, li2])
            li_list.pop(1)
    if len(li_list) > 1:
        raise ValueError('wrong li_list length:{}'.format(len(li_list)))
    for i in li_list:
        nd_idx = len(res_node)
        li_idx = len(res_lines)
        current_crossing = Point(i.xy[0][-1], i.xy[1][-1])
        node_map[nd_idx] = {li_idx}
        res_node.append(current_crossing)
        res_lines.append(i)
    return res_lines, res_node, node_map, True


def nodes_split_curve(nodesList, curve):
    '''
    This function Generate
    :param data: dict road_idx: "node": node_geos(list), "geometry": road_geos(shapely)
    :param geo_key:
    :return: separateLi(list)
    :return: separateMap(dict) Original geo_keys to separateLi idx mapping
    '''
    separatLi = list()

    rd_geo = curve
    cross_pt = nodesList
    pt_que = Queue()
    rd_list = list()
    rd_list.append(rd_geo)
    for pt in cross_pt:
        pt_que.Que_in(pt)
    '''
    use pts to split the link accordingly
    until all pts are used
    '''
    while not pt_que.Que_isEmpty():
        current_pt = pt_que.Que_out()
        l_idx = 0
        while l_idx < len(rd_list):
            current_li = rd_list[l_idx]
            dist_to_curve = check_pt_on_curve(current_pt, current_li)
            #   check the intersection releationship between the currentli and currentpt
            if dist_to_curve > 1:
                #   if distance larger than tolerance(1m), than it means this split pt is not correct.
                l_idx += 1
                continue
            if check_endpt(current_pt, current_li):
                #   check currentpt is endpt or not
                #   True：  pass
                #   False: currentpt should split the current_li
                break
            else:
                #   current_pt should split the currentli
                #   current_pt is direct on current li or not
                if dist_to_curve == 0:
                    #   current_pt is on the currentli
                    break_pt = current_pt
                    spt_list = split_polyline(break_pt, current_li)
                    rd_list.pop(l_idx)
                    rd_list = rd_list + spt_list
                    break
                else:
                    #   current_pt not on currentli
                    break_pt = nearest_pt(current_pt, current_li)
                    spt_list = split_polyline(break_pt, current_li)

                    rd_list.pop(l_idx)
                    rd_list = rd_list + spt_list
                    break
    for i in rd_list:
        idx = len(separatLi)
        separatLi.append(i)

    return separatLi


def CombinePassLinks(l1, l2):
    l1_list = Points_in_Line(l1)
    l2_list = Points_in_Line(l2)
    if l1_list[0].distance(l2_list[0]) < l1_list[-1].distance(l2_list[0]):
        if l1_list[0].distance(l2_list[0]) < 1:
            l1_list.reverse()
            return LineString(l1_list + l2_list)
        else:
            return LineString(l2_list + l1_list)
    else:
        if l1_list[-1].distance(l2_list[0]) < 1:
            return LineString(l1_list + l2_list)
        else:
            l2_list.reverse()
            return LineString(l1_list + l2_list)


def LineModification(s, end_identify=3):
    '''
    Based on  LineSimplification
    :param s:
    :return: a LineString
    '''
    total_length = s.length
    pts = np.array(s.xy)
    if pts.shape[1] < 4:
        return s
    vec, vec_length, angle = LineStringVec(pts)
    if pts.shape[1] == 4:
        #   四个点三段线
        if vec_length[0] < vec_length[1] / end_identify and vec_length[2] < vec_length[1] / end_identify:
            pts = np.delete(pts, [1, 2], axis=1)
            return LineString(pts.T)
        else:
            return s
    #   先检测波动情况
    #   波动的情况意味着两个角度的cosine值的正负值是不一样的
    num_angle = angle.shape[0]
    angle_idx = 0
    while angle_idx < num_angle:
        if angle[angle_idx] < 0.8:
            #   条件1：说明对应点没有起到形变作用，中间点应该取消
            #   条件2：说明是回折点，应当被清除
            #   构成这个角的中间点应当被remove
            #   重新计算一遍pts、vec和angle
            mid_pt_idx = angle_idx + 1
            pts = np.delete(pts, mid_pt_idx, axis=1)
            vec, vec_length, angle = LineStringVec(pts)
            num_angle = angle.shape[0]
        elif vec_length[angle_idx] + vec_length[angle_idx + 1] < 50:
            mid_pt_idx = angle_idx + 1
            pts = np.delete(pts, mid_pt_idx, axis=1)
            vec, vec_length, angle = LineStringVec(pts)
            num_angle = angle.shape[0]
        else:
            angle_idx += 1
    return LineString(pts.T)


def roundness(s):
    '''
    roundness or circularity
    :param s: s should be the list of point
    :return:
    '''
    geom = Polygon(s)
    area = geom.area
    c = geom.length
    return 4 * pi * area / np.square(c)


def straighten_links(pt_list):
    #   ensuring the start and end pt position
    if len(pt_list) < 3:
        return pt_list
    pt_list = pt_list.copy()
    # block = set()
    res_list = [pt_list[0]]
    idx = 1
    while idx < len(pt_list) - 1:
        prev = res_list[-1]
        current = pt_list[idx]
        next_pt = pt_list[idx + 1]
        if current.distance(prev) == 0:
            idx = idx + 1
            continue
        if current.distance(next_pt) == 0:
            res_list.append(current)
            idx = idx + 2
            continue
        angle = AngleCal_nd(prev, current, next_pt)
        if angle > 0:
            res_list.append(current)
        idx = idx + 1
    # pt_list = [j for i,j in enumerate(pt_list) if i not in block]
    res_list.append(pt_list[-1])
    return res_list


def Crossing_Checking(data, geo_key):
    '''
    This function Generate
    :param data: dict road_idx: "node": node_geos(list), "geometry": road_geos(shapely)
    :param geo_key:
    :return: separateLi(list)
    :return: separateMap(dict) he original geo_keys are mapped to separateLi idx
    '''
    separatLi = list()
    separatDict = dict()

    for li in tqdm(range(len(geo_key))):
        rd_idx = geo_key[li]
        rd_geo = data[rd_idx]['geometry']
        cross_pt = data[rd_idx]['nodes']

        if len(cross_pt) == 0:
            continue
        pt_que = Queue()
        rd_list = list()
        rd_list.append(rd_geo)
        for pt in cross_pt:
            pt_que.Que_in(pt)
        '''
        Handling each point-following-line case by type.
        until all pts are used
        '''
        while not pt_que.Que_isEmpty():
            current_pt = pt_que.Que_out()
            l_idx = 0
            while l_idx < len(rd_list):
                current_li = rd_list[l_idx]
                dist_to_curve = check_pt_on_curve(current_pt, current_li)
                #    Start detecting the relationship between currentli and currentpt
                if dist_to_curve > 1:
                    #   If the distance is greater than 1m, it means that it is not related to this line segment,
                    #   and you can find the next line directly.
                    l_idx += 1
                    continue
                if check_endpt(current_pt, current_li):
                    #   check whether the pt is the endpt
                    #   True：  do not split
                    #   False: pt should split
                    break
                else:
                    #   find split pt and corresponding link
                    #   but the split pt may not locate on the the link geometry
                    if dist_to_curve == 0:
                        #   the split geometry is on the link geometry
                        break_pt = current_pt
                        spt_list = split_polyline(break_pt, current_li)
                        rd_list.pop(l_idx)
                        rd_list = rd_list + spt_list
                        break
                    else:
                        #   not on the link geometry
                        break_pt = nearest_pt(current_pt, current_li)
                        spt_list = split_polyline(break_pt, current_li)

                        rd_list.pop(l_idx)
                        rd_list = rd_list + spt_list
                        break
        '''
        The point has completely split the line
        The line is stored in separatLi
        the mapping relationship is preserved in separatDict
        '''
        for i in rd_list:
            idx = len(separatLi)
            separatLi.append(i)
            separatDict[idx] = li

    return separatLi, separatDict