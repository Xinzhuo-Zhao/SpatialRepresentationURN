from tqdm import tqdm
from geometry_processing import *
from Angle import *
from copy import deepcopy
from data_structure import *
from data_input import *

class GraphBuilder:
    def __init__(self, link_storage, links, node, node_storage, filename):
        self.link_storage = link_storage
        self.links = links
        self.node = node
        self.node_storage = node_storage
        self.crossing_management = dict()
        self.merging_dict = dict()
        self.merging_storage = list()
        self.prev_merging_dict = dict()
        self.prev_merging_stg= list()
        self.filename = filename
    
    def save(self):
        save_testingFile(self.link_storage, self.filename+ 'links')
        save_testingFile(self.node_storage, self.filename + 'nodes')
        save_data_structure(self.links, self.filename + "Graph_links.pkl")
        save_data_structure(self.node, self.filename + "Graph_node.pkl")
    
    def run(self):
        self.graph_connections()
        self.graph_establishment()
        
        times = len(self.link_storage)
        for i in tqdm(range(times)):
            if self.link_storage[i] is None:
                continue
            self.run_linesimplification(i)
        
        for node_id, node_set in tqdm(self.node.items()):
            self.run_Combine_PassLinks(node_id, node_set)
        
        for link_id in tqdm(range(len(self.link_storage))):
            if self.link_storage[link_id] is None:
                continue
            if len(self.links[link_id]) == 0:
                self.link_storage[link_id] = None
                if len(self.links[link_id]) > 0:
                    kept = self.links[link_id].pop()
                    for nd in self.links[link_id]:
                        self.node[kept] = self.node[kept].union(self.node[nd])
                        self.node.pop(nd)
                        self.node_storage[nd] = None
                    self.node[kept].discard(link_id)
                self.links.pop(link_id)
        
        for l, n in tqdm(self.links.items()):
            if not self.link_storage[l]:
                continue
            self.run_remove_anomalous_shortlinks(l, n)
        self.check_graph_validation()

        for i in tqdm(range(len(self.link_storage))):
            if self.link_storage[i] is not None:
                self.link_storage[i] = LineModification(self.link_storage[i])

        for val, items in tqdm(self.node.items()):
            self.run_IdentifyParrallelCrossing(val, items)

        times = 0
        while times < len(self.merging_storage):
            if self.merging_storage[times] is None:
                times += 1
                continue
            self.run_find_close_cycles_on_graph(times)
            times += 1

        self.prev_merging_dict = deepcopy(self.merging_dict)
        self.prev_merging_stg = deepcopy(self.merging_storage)
        self.run_cycle_simplify()

        self.run_finishing()
        
    def run_finishing(self, looptimes=0):
        print("Currently in {} times".format(looptimes))
        print("run_Combine_PassLinks in {} time".format(looptimes))
        for n_idx, link_set in tqdm(self.node.items()):
            self.run_Combine_PassLinks(n_idx, link_set)

        self.check_graph_validation()

        print("run_remove_anomalous_shortlinks in {} time".format(looptimes))
        for l, n in tqdm(self.links.items()):
            if self.link_storage[l] is None:
                continue
            self.run_remove_anomalous_shortlinks(l, n)
        self.check_graph_validation()
        print("run_IdentifyParrallelCrossing in {} time".format(looptimes))
        self.merging_dict = dict()
        self.merging_storage = list()
        for n_idx, link_set in tqdm(self.node.items()):
            self.run_IdentifyParrallelCrossing(n_idx, link_set)

        print("run_find_close_cycles_on_graph in {} time".format(looptimes))
        times = 0
        while times < len(self.merging_storage):
            if self.merging_storage[times] is None:
                times += 1
                continue
            self.run_find_close_cycles_on_graph(times)
            times += 1

        update_status = self.Identify_MergingDifference(self.prev_merging_dict, 
                                                        self.prev_merging_stg, 
                                                        self.merging_dict, 
                                                        self.merging_storage)
        if update_status:
            self.prev_merging_dict = deepcopy(self.merging_dict)
            self.prev_merging_stg = deepcopy(self.merging_storage)
            print("run_cycle_simplify in {} time".format(looptimes))
            self.run_cycle_simplify()
            self.run_finishing( looptimes + 1)
        else:
            return
    
    def run_graph_establishment(self):
        for l, nd_set in self.crossing_management.items():
            spt_list = nodes_split_curve([self.node_storage[i] for i in nd_set], self.link_storage[l])
            add_idx = list()
            for i in spt_list:
                current_index = len(self.link_storage)
                add_idx.append(current_index)
                self.link_storage.append(i)
                self.links[current_index] = set()
            self.link_storage[l] = None
            for relate_nd in self.links[l]:
                for relate_li in add_idx:
                    if check_endpt(self.node_storage[relate_nd], self.link_storage[relate_li]):
                        self.node[relate_nd].add(relate_li)
                        self.links[relate_li].add(relate_nd)
                self.node[relate_nd].remove(l)
                if l in self.node[relate_nd]:
                    raise ValueError("未知错误：未找到继承关系线")
            self.links.pop(l)
    
    def run_linesimplification(self, link_id):
        spt_list, spt_node, nd2li, status = LineSimplication(self.link_storage[link_id], link_id)
        if not status:
            return
        if len(spt_list) == 1:
            self.link_storage[link_id] = spt_list[0]
            return
        related_nd = self.links.pop(link_id)
        self.link_storage[link_id] = None
        dup_record = dict()
        #   检测有没有被纪录过
        for nd in related_nd:
            self.node[nd].remove(link_id)
        for nd_idx in range(len(spt_node)):
            li_set = nd2li[nd_idx]
            current_node = spt_node[nd_idx]
            current_nidx = len(self.node_storage)
            if nd_idx == 0 or nd_idx == len(spt_node) - 1:
                #   说明这是开始点和终止点
                #   说明原图中已经有对应点了
                if len(li_set) != 1:
                    raise ValueError('未知错误：端点有多条线')
                li = list(li_set)[0]
                for nd in related_nd:
                    #   nd is the index of self.node_storage
                    if current_node.distance(self.node_storage[nd]) < 1:
                        #   说明currentnode和原node是同一个点
                        if li not in dup_record:
                            #   li第一次被遍历到
                            current_lidx = len(self.link_storage)
                            dup_record[li] = current_lidx
                            self.links[current_lidx] = {nd}
                            self.node[nd].add(current_lidx)
                            self.link_storage.append(spt_list[li])
                        else:
                            #   li已经被遍历过了
                            current_lidx = dup_record[li]
                            self.links[current_lidx].add(nd)
                            self.node[nd].add(current_lidx)
                        break
            else:
                #   说明是中间点
                #   说明是新点
                self.node[current_nidx] = set()
                self.node_storage.append(current_node)
                for li in li_set:
                    #   li is the index of spt_list
                    if li not in dup_record:
                        #   li第一次被遍历到
                        current_lidx = len(self.link_storage)
                        dup_record[li] = current_lidx
                        self.links[current_lidx] = {current_nidx}
                        self.node[current_nidx].add(current_lidx)
                        self.link_storage.append(spt_list[li])
                    else:
                        #   li已经被遍历过了
                        current_lidx = dup_record[li]
                        self.links[current_lidx].add(current_nidx)
                        self.node[current_nidx].add(current_lidx)

    def run_Combine_PassLinks(self, n_idx, n_set):
        if n_set is None:
            return
        if len(n_set) == 2:
            li = list(n_set)
            angle = abs(AngleCal(self.link_storage[li[0]], self.link_storage[li[1]]))
            if angle > 0.6:
                new_link = CombinePassLinks(self.link_storage[li[0]], self.link_storage[li[1]])
                idx = len(self.link_storage)
                new_node = set()
                self.link_storage.append(new_link)
                for l in n_set:
                    self.links[l].discard(n_idx)
                    new_node = new_node.union(self.links[l])
                    self.link_storage[l] = None
                    for n in self.links[l]:
                        self.node[n].discard(l)
                        self.node[n].add(idx)
                    self.links[l] = None
                    # self.links.pop(l)
                self.links[idx] = new_node
                self.node[n_idx] = None
                # self.node.pop(n_idx)
                self.node_storage[n_idx] = None

    def run_remove_anomalous_shortlinks(self, l_idx, n_set, length_limits=50):
        if self.link_storage[l_idx].length < length_limits:
            # meanlength = 0
            # count = 0
            depth1, depth_end1 = self.Line_forward_Line(l_idx)
            depth2, _ = self.Line_forward_Line(depth1, depth_end1)
            related_li = depth2 + depth1
            length_record = list()
            for li_geo_idx in related_li:
                length_record.append(self.link_storage[li_geo_idx].length)
            # if l_idx == 344749:
            #     print(related_li)
            #     print(length_record)
            if len(length_record) == 0:
                #   经过修改后l_idx成孤立线段了：
                # del_links([l_idx])
                # #   同时你n_set也都是无效的
                # for n_id in n_set:
                #     if len(self.node[n_id]) == 0:
                #         self.node.pop(n_id)
                #         self.node_storage[n_id] = None
                return
            q1, q3 = np.percentile(length_record, [33, 75])
            meanlength = np.mean(length_record)
            midlength = np.median(length_record)
            # for i in n:
            #     for j in self.node[i]:
            #         if j != l:
            #             meanlength += self.link_storage[j].length
            #             count += 1
            # meanlength = meanlength / count
            if (meanlength / self.link_storage[l_idx].length > 4) and (self.link_storage[l_idx].length < q1):
                # 过短的奇异线条删除：
                # 包括线条的storage替换成None，线条目录删除，对应node修改
                # node修改：
                # 该线条连接的node，坍缩成一个点，同时返回修改对应的线条
                # 新的点继承原来点的链接线条
                # 相关线条修改:
                # 从被修改的点返回修改被改动的线
                # 1连接的点idx被修改成新的
                # 2对线的断点进行延伸
                pt_idx = len(self.node_storage)
                # 建立新节点为中心点
                new_pt = midPt(self.link_storage[l_idx])
                self.node_storage.append(new_pt)
                new_node = set()
                # 继承原有线段两个端点的连接线段关系
                for i in n_set:
                    self.node[i].discard(l_idx)
                    new_node = new_node.union(self.node[i])
                    # self.node[i] = None
                    self.node.pop(i)
                    self.node_storage[i] = None
                self.node[pt_idx] = new_node
                # 对连接的所有线段进行修改
                for fixline in new_node:
                    # 几何关系修改
                    l1 = self.link_storage[fixline]
                    if l1.length == 0:
                        self.link_storage[fixline] = None
                        if fixline in self.links:
                            for n in self.links[fixline]:
                                self.node[n].discard(fixline)
                            self.links.pop(fixline)
                        continue
                    if l1 is None:
                        continue
                    l1_list = Points_in_Line(l1)
                    startPt = l1_list[0]
                    endPt = l1_list[-1]
                    if startPt.distance(new_pt) > endPt.distance(new_pt):
                        BefVec = np.array(l1_list[-1].coords[0]) - np.array(l1_list[-2].coords[0])
                        AftVec = np.array(new_pt.coords[0]) - np.array(l1_list[-1].coords[0])
                        ChangedVec = BefVec - AftVec
                        BefLen = norm(BefVec, ord=2, axis=0)
                        ChangedLen = norm(ChangedVec, ord=2, axis=0)
                        while BefLen == 0:  # -1位置存在重合点
                            l1_list.pop()
                            BefVec = np.array(l1_list[-1].coords[0]) - np.array(l1_list[-2].coords[0])
                            AftVec = np.array(new_pt.coords[0]) - np.array(l1_list[-1].coords[0])
                            ChangedVec = BefVec - AftVec
                            BefLen = norm(BefVec, ord=2, axis=0)
                            ChangedLen = norm(ChangedVec, ord=2, axis=0)
                        if ChangedLen < 20 and BefLen < 20:
                            l1_list[-1] = new_pt
                        else:
                            ept_vec = np.array(l1_list[-1].coords[0]) - BefVec / BefLen * 10
                            l1_list[-1] = Point(ept_vec[0], ept_vec[1])
                            l1_list.append(new_pt)
                    else:
                        BefVec = np.array(l1_list[0].coords[0]) - np.array(l1_list[1].coords[0])
                        AftVec = np.array(new_pt.coords[0]) - np.array(l1_list[1].coords[0])
                        ChangedVec = BefVec - AftVec
                        BefLen = norm(BefVec, ord=2, axis=0)
                        ChangedLen = norm(ChangedVec, ord=2, axis=0)
                        while BefLen == 0:  # -1位置存在重合点
                            l1_list.pop(0)
                            BefVec = np.array(l1_list[0].coords[0]) - np.array(l1_list[1].coords[0])
                            AftVec = np.array(new_pt.coords[0]) - np.array(l1_list[1].coords[0])
                            ChangedVec = BefVec - AftVec
                            BefLen = norm(BefVec, ord=2, axis=0)
                            ChangedLen = norm(ChangedVec, ord=2, axis=0)
                        if ChangedLen < 20 and BefLen < 20:
                            l1_list[0] = new_pt
                        else:
                            ept_vec = np.array(l1_list[0].coords[0]) - BefVec / BefLen * 10
                            l1_list[0] = Point(ept_vec[0], ept_vec[1])
                            l1_list.insert(0, new_pt)
                    self.link_storage[fixline] = LineString(l1_list)
                    # 连接关系修改
                    self.links[fixline] = self.links[fixline] - n_set
                    self.links[fixline].add(pt_idx)
                self.links[l_idx] = None
                # self.links.pop(l)
                self.link_storage[l_idx] = None

    def run_IdentifyParrallelCrossing(self, n_idx, n_set):
        if n_set is None:
            return
        loop = list(n_set)
        for l1 in range(len(loop)):
            for l2 in range(l1 + 1, len(loop)):
                # 保证角度计算是从交点出发
                t1 = [Point(self.link_storage[loop[l1]].xy[0][0], self.link_storage[loop[l1]].xy[1][0]),
                      Point(self.link_storage[loop[l1]].xy[0][-1], self.link_storage[loop[l1]].xy[1][-1])]
                t2 = [Point(self.link_storage[loop[l2]].xy[0][0], self.link_storage[loop[l2]].xy[1][0]),
                      Point(self.link_storage[loop[l2]].xy[0][-1], self.link_storage[loop[l2]].xy[1][-1])]
                if t1[0].distance(self.node_storage[n_idx]) > 1:
                    t1.reverse()
                if t2[0].distance(self.node_storage[n_idx]) > 1:
                    t2.reverse()
                angle = AngleCal(LineString(t1), LineString(t2))
                if angle > 0.9:
                    if (loop[l1] in self.merging_dict.keys()) and (loop[l2] in self.merging_dict.keys()):
                        idx = len(self.merging_storage)
                        self.merging_storage.append(
                            self.merging_storage[self.merging_dict[loop[l1]]].union(
                                self.merging_storage[self.merging_dict[loop[l2]]]))
                        self.merging_storage[self.merging_dict[loop[l1]]] = None
                        self.merging_storage[self.merging_dict[loop[l2]]] = None
                        for i in self.merging_storage[-1]:
                            self.merging_dict[i] = idx
                    elif loop[l1] in self.merging_dict.keys():
                        self.merging_dict[loop[l2]] = self.merging_dict[loop[l1]]
                        self.merging_storage[self.merging_dict[loop[l1]]].add(loop[l2])
                    elif loop[l2] in self.merging_dict.keys():
                        self.merging_dict[loop[l1]] = self.merging_dict[loop[l2]]
                        self.merging_storage[self.merging_dict[loop[l2]]].add(loop[l1])
                    else:
                        self.merging_dict[loop[l1]] = len(self.merging_storage)
                        self.merging_dict[loop[l2]] = len(self.merging_storage)
                        self.merging_storage.append({loop[l1], loop[l2]})

    def run_find_close_cycles_on_graph(self, storage_id):
        #   先检查li_set
        li_set = self.merging_storage[storage_id].copy()
        pts_status, brk_pts = self.check_breaknd(li_set)
        if len(brk_pts) == 0:
            #   不一定闭合，但是没有可以往外延伸的点了
            #   例如位于边缘的link也没有闭合，但是也不具备bkps
            return
        while brk_pts:
            n = brk_pts.pop()
            if pts_status[n]['status'] != 'breaknd':
                continue

            path_list, pts_status, brk_pts = self.explore_roads(n, pts_status, brk_end=brk_pts)
            #   status表示是否寻找到其他merging_dict的组，如果有则更新status 并更新需要找的pts_status和brk_pts
            status = 0

            if len(path_list):  # 如果有符合要求的路径被找到了
                #   开始找路径上的所有线段并对merging_dict和merging_storage进行更新
                for path in path_list:
                    #   获得节点之间的link
                    newlinks = self.nodeseq_to_linkseq(path)
                    li_set = li_set.union(newlinks)
                    # print(li_set)
                    #   获取所有path相关的merging_storage内容
                    for i in newlinks:
                        if i in self.merging_dict:
                            #   找到之前被记录的其他组
                            status = 1
                            li_set = self.merging_storage[self.merging_dict[i]].union(li_set)
                    #   更新node_storage
                    self.merging_storage[storage_id] = li_set
                    for i in li_set:
                        #   对于新更新进来的线段
                        if i not in self.merging_dict:
                            #   说明i是新找到的link
                            self.merging_dict[i] = storage_id
                        elif self.merging_dict[i] != storage_id:
                            #   说明找到的是其他组内的link
                            self.merging_storage[self.merging_dict[i]] = None
                            self.merging_dict[i] = storage_id
            if status:
                pts_status, brk_pts = self.check_breaknd(li_set, pts_status, brk_pts)

    def run_cycle_simplify(self):
        for idx, testing_g in tqdm(enumerate(self.merging_storage)):
            if testing_g is None:
                continue
            for i in testing_g:
                if i not in self.links:
                    self.merging_storage[idx] = None
                    break
            if self.merging_storage[idx] is None:
                continue
            sub_l, sub_n = self.pre_cycle_modification(testing_g)
            sub_node_graph = self.node_graph_construction(sub_n, sub_l)
            node_manage = self.nodetype_in_sub(sub_n, self.node)
            if len(sub_n) < 2:  # 不构成回路
                continue
            n_idx = self.select_best_initial_junction(sub_n)
            cycle_list = self.find_cycles(n_idx, sub_node_graph, sub_n)
            cycle_list = self.cycle_filter(cycle_list)
            while len(cycle_list):
                min_cycle = self.minimal_cycle(cycle_list)
                sub_n, sub_l = self.cycle_simplify(min_cycle, sub_n, sub_l, node_manage)
                sub_node_graph = self.node_graph_construction(sub_n, sub_l)
                node_manage = self.nodetype_in_sub(sub_n, self.node)
                n_idx = self.select_best_initial_junction(sub_n)
                cycle_list = self.find_cycles(n_idx, sub_node_graph, sub_n)
                cycle_list = self.cycle_filter(cycle_list)
            self.merging_storage[idx] = None

    def graph_establishment(self):
        for l, nd_set in self.crossing_management.items():
            spt_list = nodes_split_curve([self.node_storage[i] for i in nd_set], self.link_storage[l])
            add_idx = list()
            for i in spt_list:
                current_index = len(self.link_storage)
                add_idx.append(current_index)
                self.link_storage.append(i)
                self.links[current_index] = set()
            self.link_storage[l] = None
            for relate_nd in self.links[l]:
                for relate_li in add_idx:
                    if check_endpt(self.node_storage[relate_nd], self.link_storage[relate_li]):
                        self.node[relate_nd].add(relate_li)
                        self.links[relate_li].add(relate_nd)
                self.node[relate_nd].remove(l)
                if l in self.node[relate_nd]:
                    raise ValueError("未知错误：未找到继承关系线")
            self.links.pop(l)

    def Identify_MergingDifference(self, prev_merging_dict, prev_merging_stg, merging_d, merging_stg):
        for i in merging_d.keys():
            if i in prev_merging_dict:
                #   说明之前这个点已经被辨识过
                if prev_merging_stg[prev_merging_dict[i]] != merging_stg[merging_d[i]]:
                    #   说明需要进一步合并
                    return True
            else:
                return True
        return False
    
    def graph_connections(self):
        crossing_management = dict()
        pt_idx = [0, 0]
        for o in tqdm(range(len(self.link_storage))):
            for d in range(o + 1, len(self.link_storage)):
                line1 = self.link_storage[o]
                line2 = self.link_storage[d]

                if o not in self.links.keys():
                    self.links[o] = set()
                if d not in self.links.keys():
                    self.links[d] = set()
                if pt_idx[0] not in self.node.keys():
                    self.node[pt_idx[0]] = set()

                if line1.intersects(line2):
                    intersectionGEOM = line1.intersection(line2)
                    if intersectionGEOM.type == 'Point':
                        crossing_management, pt_idx = self.intersection_of_pt(intersectionGEOM, crossing_management, pt_idx,
                                                                         o, d)
                    elif intersectionGEOM.type == 'MultiPoint':
                        crossing_management, pt_idx = self.intersection_of_multipt(intersectionGEOM, crossing_management,
                                                                              pt_idx,
                                                                              o, d)
                    elif intersectionGEOM.type == "LineString":
                        crossing_management, pt_idx = self.intersection_of_line(intersectionGEOM, crossing_management,
                                                                           pt_idx,
                                                                           o, d)
                    elif intersectionGEOM.type == "MultiLineString":
                        crossing_management, pt_idx = self.intersection_of_multili(intersectionGEOM, crossing_management,
                                                                              pt_idx,
                                                                              o, d)
                    elif intersectionGEOM.type == "GeometryCollection":
                        crossing_management, pt_idx = self.intersection_of_geom_collection(intersectionGEOM,
                                                                                      crossing_management, pt_idx, o, d)
                    else:
                        raise TypeError("出现预期外的交叉类型：{}".format(intersectionGEOM.type))
        self.crossing_management = crossing_management

    def intersection_of_geom_collection(self, geom_collection, crossing_management, pt_idx, o, d):
        for geom in geom_collection.geoms:
            if geom.type == 'Point':
                crossing_management, pt_idx = self.intersection_of_pt(geom, crossing_management, pt_idx, o, d)
            elif geom.type == 'MultiPoint':
                crossing_management, pt_idx = self.intersection_of_multipt(geom, crossing_management, pt_idx, o, d)
            elif geom.type == "LineString":
                crossing_management, pt_idx = self.intersection_of_line(geom, crossing_management, pt_idx, o, d)
            elif geom.type == "MultiLineString":
                crossing_management, pt_idx = self.intersection_of_multili(geom, crossing_management, pt_idx, o, d)
            else:
                raise TypeError("出现预期外的交叉类型：{}".format(geom.type))
        return crossing_management, pt_idx

    def intersection_of_multili(self, intersectionLis, crossing_management, pt_idx, o, d):
        for n in intersectionLis.geoms:
            crossing_management, pt_idx = self.intersection_of_line(n, crossing_management, pt_idx, o, d)
        return crossing_management, pt_idx

    def intersection_of_line(self, intersectionLi, crossing_management, pt_idx, o, d):
        nd1_geom = Point(intersectionLi.xy[0][0], intersectionLi.xy[1][0])
        nd2_geom = Point(intersectionLi.xy[0][1], intersectionLi.xy[1][1])
        crossing_management, pt_idx = self.intersection_of_multipt([nd1_geom, nd2_geom], crossing_management, pt_idx, o, d)
        return crossing_management, pt_idx

    def intersection_of_multipt(self, intersectionPTS, crossing_management, pt_idx, o, d):
        if intersectionPTS.__class__ == list:
            for n in intersectionPTS:
                crossing_management, pt_idx = self.intersection_of_pt(n, crossing_management, pt_idx, o, d)
            return crossing_management, pt_idx

        for n in intersectionPTS.geoms:
            crossing_management, pt_idx = self.intersection_of_pt(n, crossing_management, pt_idx, o, d)
        return crossing_management, pt_idx

    def intersection_of_pt(self, intersectionPT, crossing_management, pt_idx, o, d):
        new_pt = 1  # 默认找到的新的点位是(new_pt = 1)
        for i in range(len(self.node_storage)):
            if self.node_storage[i].distance(intersectionPT) < 1:
                intersectionPT = self.node_storage[i]
                pt_idx[1] = i
                new_pt = 0
                break
        #   检测相交点是在线端点还是中间
        if new_pt == 0:
            self.node[pt_idx[1]].add(o)
            self.node[pt_idx[1]].add(d)
            self.links[o].add(pt_idx[1])
            self.links[d].add(pt_idx[1])
            if not check_endpt(intersectionPT, self.link_storage[o]):
                if o in crossing_management.keys():
                    crossing_management[o].add(pt_idx[1])
                else:
                    crossing_management[o] = {pt_idx[1]}
            if not check_endpt(intersectionPT, self.link_storage[d]):
                if d in crossing_management.keys():
                    crossing_management[d].add(pt_idx[1])
                else:
                    crossing_management[d] = {pt_idx[1]}
            # new_pt = 1
        else:
            if pt_idx[0] not in self.node.keys():
                self.node[pt_idx[0]] = set()
            self.node[pt_idx[0]].add(o)
            self.node[pt_idx[0]].add(d)
            self.links[o].add(pt_idx[0])
            self.links[d].add(pt_idx[0])
            self.node_storage.append(intersectionPT)
            if not check_endpt(intersectionPT, self.link_storage[o]):
                if o in crossing_management.keys():
                    crossing_management[o].add(pt_idx[0])
                else:
                    crossing_management[o] = {pt_idx[0]}
            if not check_endpt(intersectionPT, self.link_storage[d]):
                if d in crossing_management.keys():
                    crossing_management[d].add(pt_idx[0])
                else:
                    crossing_management[d] = {pt_idx[0]}
            pt_idx[0] += 1
        return crossing_management, pt_idx

    def explore_roads(self, startnd, cross_status, brk_end, blockednd=None, timelimits=15):
        '''
        一个能增量的广度优先算法
        '''

        if blockednd is None:
            blockednd = set()

        brk_end = deepcopy(brk_end)
        cross_status = deepcopy(cross_status)

        end_positions = set(cross_status.keys()).difference(set(brk_end))
        # 定义一个队列，用于存储遍历的节点
        queue = Queue()
        # 将起点加入队列中
        queue.Que_in(startnd)

        # 定义一个字典，用于记录每个节点的前一个节点，最终用于查找最短路径
        # 并记录到startli的总距离
        prev = {startnd: {'prev': None, 'distance': 0, 'step': 0}}
        neighbor, _, _ = self.Node_forward_Node(startnd, self.node, self.links, end_positions)
        # print(neighbor)
        currentnd = None
        _, prevnd = self.Link_to_Node(list(cross_status[startnd]['self.links'])[0], startnd, self.links)
        # print(prevnd)
        end_targets = set()

        while not queue.Que_isEmpty():
            # print("Current self.node id is : {}".format(currentnd))
            # print("Current self.node id's Neighbor is : {}".format(neighbor))
            if currentnd is None:
                currentnd = queue.Que_out()
                blockednd.add(currentnd)
            elif prev[currentnd]['step'] < timelimits:
                currentnd = queue.Que_out()
                neighbor, _, _ = self.Node_forward_Node(currentnd, self.node, self.links, blockednd)
            else:
                #   如果现有的li已经超过timelimits的情况，不需要继续延伸了
                currentnd = queue.Que_out()
                neighbor, _, _ = self.Node_forward_Node(currentnd, self.node, self.links, blockednd)
                continue

            for nextnd in neighbor:
                inter_li, _ = self.Link_between_Nodes(currentnd, nextnd)
                new_distance = self.link_storage[inter_li].length + prev[currentnd]['distance']
                # angle = AngleCal_nd(self.node_storage[prevnd], self.node_storage[currentnd], self.node_storage[nextnd])
                # print(angle)
                # if abs(angle) > 0.6:
                if (nextnd not in blockednd):  # nextli是第一次找到
                    #   保证不会反向传播
                    prev[nextnd] = dict()
                    prev[nextnd]['prev'] = currentnd
                    prev[nextnd]['distance'] = new_distance
                    prev[nextnd]['step'] = prev[currentnd]['step'] + 1
                    blockednd.add(nextnd)
                    if nextnd in brk_end:
                        end_targets.add(nextnd)
                        # print("新找到的点是：{}".format(nextnd))
                    else:
                        queue.Que_in(nextnd)

                else:  # nextli不是第一次找到
                    if new_distance < prev[nextnd]['distance']:
                        #   说明范围内找到了距离更短的点
                        prev[nextnd]['prev'] = currentnd
                        prev[nextnd]['distance'] = new_distance
                        prev[nextnd]['step'] = prev[currentnd]['step'] + 1

        routine = list()
        if len(end_targets) == 0:
            return routine, cross_status, brk_end
        for nd in end_targets:
            path = list()
            currentnd = nd
            while currentnd is not None:
                path.append(currentnd)
                currentnd = prev[currentnd]['prev']
            if len(path) == 0:
                continue
            tem_path = self.find_whole_path_within_seq(path, cross_status)
            if len(tem_path) < 3:  # 说明不成环
                continue
            path_roundness = roundness(tem_path)
            if path_roundness < 0.25:
                routine.append(path)
                cross_status[startnd]['status'] = "normalnd"
                cross_status[nd]['status'] = "normalnd"
                brk_end.remove(nd)
        if cross_status[startnd]['status'] == 'breaknd':
            cross_status[startnd]['status'] = 'passnd'
        if startnd in brk_end:
            brk_end.remove(startnd)
        return routine, cross_status, brk_end

    def find_whole_path_within_seq(self, seq, cross):
        '''
        用于找到seq中所有空间铆钉点，并形成一个闭环回路
        :param seq:
        :param cross:
        :return:
        '''
        seq = seq.copy()
        start = seq[0]
        end = seq[-1]
        whole_path_node = list()
        #   补充每两个节点之间的线段中的位置铆钉点，首尾点与node重合因此不需要。
        prevnd = start

        # print("cross:{}".format(cross))
        # print("end:{}".format(end))

        #   先找到start 和 end中间的点
        for nd, link_record in cross.items():
            if nd == start or nd == end:
                continue
            if link_record['self.links'].intersection(cross[start]['self.links']) and (nd not in seq):
                seq = [nd] + seq
            elif link_record['self.links'].intersection(cross[end]['self.links']) and (nd not in seq):
                seq.append(nd)
        if seq[0] != start and seq[-1] != end:
            #   收尾点都更新的情况,最多还存在一个junction
            if not cross[seq[0]]["self.links"].intersection(cross[seq[0]]["self.links"]):
                #   说明start 和 end 不能形成闭环
                nd1, _, _ = self.Node_forward_Node(seq[0], self.node, self.links, endpt=set(seq))
                nd2, _, _ = self.Node_forward_Node(seq[-1], self.node, self.links, endpt=set(seq))
                missing_nd = set(nd1).intersection(set(nd2))
                if len(missing_nd) == 1:
                    seq.append(missing_nd.pop())
                else:
                    raise ValueError("Unable to close between {} and {}".format(seq[0], seq[-1]))

        for currentnd in seq:
            if currentnd == prevnd:
                #   从第二位置遍历开始
                whole_path_node.append(self.node_storage[prevnd])
                continue
            interlink: set = self.node[prevnd].intersection(self.node[currentnd])
            if len(interlink) > 1:
                #   多条相交线段中，返回最短那条
                shortest_link = interlink.pop()
                shortest_length = self.link_storage[shortest_link].length
                while interlink:
                    current_link = interlink.pop()
                    current_length = self.link_storage[current_link].length
                    if current_length < shortest_length:
                        shortest_length = current_length
                        shortest_link = current_link
                interlink = shortest_link
            elif len(interlink) == 1:
                interlink = interlink.pop()
            else:
                raise TypeError("Check the link relationship between {} and {} and the path is {}".format(prevnd,
                                                                                                          currentnd,
                                                                                                          seq))
            interlink_pts = Points_in_Line(self.link_storage[interlink])
            if len(interlink_pts) == 2:
                whole_path_node.append(self.node_storage[currentnd])
            elif len(interlink_pts) > 2:
                #   需要补足中间的铆钉点
                #   先确认中间点的顺序
                if interlink_pts[0].distance(self.node_storage[prevnd]) > 1:
                    interlink_pts.reverse()
                for idx, nd in enumerate(interlink_pts):
                    if idx == 0:
                        continue
                    whole_path_node.append(nd)
            else:
                raise TypeError("Check the geom of link {}".format(interlink))

            prevnd = currentnd
        return whole_path_node

    def check_breaknd(self, li_set, crossing_record=None, breakpt_record=None):
        if crossing_record is None:
            crossing_record = dict()
        if breakpt_record is None:
            breakpt_record = list()
        for li in li_set:
            for nd in self.links[li]:
                num_crossing = 0
                if nd not in crossing_record:
                    crossing_nd_set = set()
                    for ndtolink in self.node[nd]:
                        if ndtolink in li_set:
                            num_crossing += 1
                            crossing_nd_set.add(ndtolink)
                    crossing_record[nd] = {"self.links": crossing_nd_set}
                    if num_crossing < 2:
                        crossing_record[nd]['status'] = 'breaknd'
                        breakpt_record.append(nd)
                    else:
                        crossing_record[nd]['status'] = 'normalnd'
        return crossing_record, breakpt_record

    def pre_cycle_modification(self, bf_merge_set):
        sub_link, sub_node = self.subgraph_filtering(bf_merge_set)
        terminal_map = self.find_terminals(sub_link)
        while len(terminal_map):
            sub_link, sub_node = self.remove_terminals(sub_link, sub_node, terminal_map)
            terminal_map = self.find_terminals(sub_link)
        return sub_link, sub_node

    def subgraph_filtering(self, involving_links: set):
        sub_link = dict()
        sub_node = dict()
        #   初始化sub_link的数据结构
        for l in involving_links:
            sub_link[l] = set()
        for l in sub_link.keys():
            for n in self.links[l]:
                if len(self.node[n].intersection(involving_links)) > 1:
                    #   说明n这个节点有至少两条连接并不构成孤点
                    sub_link[l].add(n)
                    if n in sub_node:
                        sub_node[n].add(l)
                    else:
                        sub_node[n] = {l}
        return sub_link, sub_node

    def find_terminals(self, link_dict):
        terminal_dict = dict()
        for l_idx, n_set in link_dict.items():
            if len(n_set) == 1:
                terminal_dict[l_idx] = list(n_set)[0]
        return terminal_dict

    def node_graph_construction(self, n_links, l_node):
        node_set = set(n_links.keys())
        node_graph = dict()
        for n_idx, link_set in n_links.items():
            next_nd, _, _ = self.Node_forward_Node(n_idx, n_links, l_node)
            next_nd = set(next_nd).intersection(node_set)
            node_graph[n_idx] = next_nd
        return node_graph

    def remove_terminals(self, link_dict: dict, node_dict: dict, changing_map: dict):
        for l_idx, n_idx in changing_map.items():
            link_dict, node_dict = self.update_remove_from_link(l_idx, link_dict, node_dict)
            if n_idx in node_dict:
                #   terminal is not updated, Thus check n_idx is removed
                if len(node_dict[n_idx]) == 1:
                    #   terminal Node should be removed
                    link_dict, node_dict = self.update_remove_from_node(n_idx, link_dict, node_dict)
            #   因为运行这个的前提是terminal node只跟一条线相连，所以此时n_idx已经为空了
            # node_dict.pop(n_idx)
        return link_dict, node_dict

    def select_best_initial_junction(self, node_graph):
        best_nid = None
        num_links = None
        for nid, link_set in node_graph.items():
            if best_nid is None:
                best_nid = nid
                num_links = len(link_set)
                continue
            if len(link_set) > num_links:
                best_nid = nid
                num_links = len(link_set)
        return best_nid

    def nodetype_in_sub(self, sub_n, par_n):
        '''
        Identify the status of each nodes in the sub_graph.
        fixed_nd means that they connect to self.node out of sub_graph
        eliminated_nd means that they do not connect to any self.node out of sub_graph
        :param sub_n: problematic nodes cliped from the whole graph
        :param par_n: whole self.node graph
        :return:      dict of self.node status
        '''
        #   区分调整点和消失点（转折点放进cycle中单段管理）
        #   调整点：表示与外界连线还有连接，需要谨慎空间变动的点； 用“fixed_nd”代表
        #   消失点：只存在子图内部的点；用“eliminated_nd”代表
        node_type = dict()
        for n_id, l_set in sub_n.items():
            #   消失点由于只跟子图内部link产生连接;因此l_par_set - l_set = empty
            #   调整点与子图中外面的点还有其他的link; 因此 l_par_set - l_set还有其他的元素
            l_par_set = par_n[n_id]
            difference = l_par_set.difference(l_set)
            if len(difference) > 0:
                #   n_id is the "fixed_nd"
                node_type[n_id] = "fixed_nd"
            else:
                node_type[n_id] = "eliminated_nd"
        return node_type

    def find_cycles(self, current_nd, node_graph, node_link):
        stack = list()
        #   nid入栈
        stack.append(current_nd)
        record = self.build_visited_graph(node_graph)
        cycle = list()

        while len(stack):
            top = stack[-1]
            if len(stack) == 1:
                prev_nd = None
            if len(stack) > 1:
                prev_nd = stack[-2]
            nextpt = record[top].keys()
            for i in nextpt:
                if (i not in stack or i == current_nd) and (record[top][i]['status'] == 0):
                    #   找到一个有效点
                    #   not in routine and it has never been passed
                    #   or reach the start nd and form a cycle
                    stack.append(i)
                    record[top][i]['status'] = 1
                    dist = list(self.node[top].union(self.node[i]))[0]
                    #   if there are more than 1 self.links between two self.node
                    #   choose either of them as distance since they are likely to be similar
                    if prev_nd is None:
                        record[top][i]["total_length"] = self.link_storage[dist].length
                    else:
                        record[top][i]["total_length"] = record[prev_nd][top]["total_length"] + self.link_storage[dist].length
                    break
            if current_nd == stack[-1] and len(stack) > 2:
                #   when current_nd is founded
                #   output the whole cycle
                if len(stack) == 3:
                    # 此时有两种可能：
                    # 情况一：在一条线上来回，
                    # 情况二：在两个节点中存在复数重复线段
                    # 情况一的情况不需要更新，只更新情况二
                    first_n = stack[0]
                    second_n = stack[1]
                    if len(node_link[first_n].intersection(node_link[second_n])) > 1:
                        cycle.append((stack.copy(), record[top][stack[-1]]["total_length"]))
                        # cycle[stack] = record[top][stack[-1]]["total_length"]
                        stack.pop()
                        continue
                    else:
                        stack.pop()
                        continue

                cycle.append((stack.copy(), record[top][stack[-1]]["total_length"]))
                # cycle[stack] = record[top][stack[-1]]["total_length"]
                stack.pop()
                continue
            last = stack[-1]
            if last == top:
                # no more nd can be found nor no cycle has been found(no pop operation)
                for nd_idx in record[top]:
                    # update all status of nd_idx
                    record[top][nd_idx]["status"] = 0
                    record[top][nd_idx]["total_length"] = 0
                stack.pop()
        return cycle

    def cycle_filter(self, cycle_dict):
        cycle_res = list()
        for data in cycle_dict:
            cycle, dist = data
            cycle_geom = self.find_whole_cycle_geom(cycle)
            cycle_roundness = roundness(cycle_geom)
            if cycle_roundness < 0.25:
                cycle_res.append((cycle, dist))
        return cycle_res

    def find_whole_cycle_geom(self, cycle):
        cycle = cycle.copy()
        whole_path_node = list()
        prevnd = cycle[0]
        for currentnd in cycle:
            if currentnd == prevnd:
                #   从第二位置遍历开始
                whole_path_node.append(self.node_storage[prevnd])
                continue
            interlink: set = self.node[prevnd].intersection(self.node[currentnd])
            if len(interlink) > 1:
                #   多条相交线段中，返回最短那条
                shortest_link = interlink.pop()
                shortest_length = self.link_storage[shortest_link].length
                while interlink:
                    current_link = interlink.pop()
                    current_length = self.link_storage[current_link].length
                    if current_length < shortest_length:
                        shortest_length = current_length
                        shortest_link = current_link
                interlink = shortest_link
            elif len(interlink) == 1:
                interlink = interlink.pop()
            else:
                raise TypeError("Check the link relationship between {} and {}".format(prevnd, currentnd))

            interlink_pts = Points_in_Line(self.link_storage[interlink])
            if len(interlink_pts) == 2:
                whole_path_node.append(self.node_storage[currentnd])
            elif len(interlink_pts) > 2:
                #   需要补足中间的铆钉点
                #   先确认中间点的顺序
                if interlink_pts[0].distance(self.node_storage[prevnd]) > 1:
                    interlink_pts.reverse()
                for idx, nd in enumerate(interlink_pts):
                    if idx == 0:
                        continue
                    whole_path_node.append(nd)
            else:
                raise TypeError("Check the geom of link {}".format(interlink))

            prevnd = currentnd
        return whole_path_node

    def minimal_cycle(self, cycle_dict: list):
        minimal_length = None
        minimal_cycle = None
        for data in cycle_dict:
            cycle, dist = data
            if minimal_cycle is None:
                minimal_length = dist
                minimal_cycle = cycle
            if dist < minimal_length:
                minimal_length = dist
                minimal_cycle = cycle
        return minimal_cycle

    def build_visited_graph(self, graph):
        '''
        this method is to manage the "visited" situation
        :param graph:
        :return: record is a dict_like object:{nid:{nds_linked_to_nid:{statue,"total_length"}}}
        '''
        record = dict()
        for n_id, next_nd_set in graph.items():
            record[n_id] = dict()
            for next in next_nd_set:
                record[n_id][next] = {"status": 0, "total_length": 0}
        return record

    def find_internal_nd(self, seq_o2d, seq_d2o, o_nd, d_nd):
        '''
        建立从O到D之间，共有左右两支节点链；输入时长度未知，但顺序都是从o 到 d; 这个方法是要建立位于从O到D的两个链路之间的中点链接；
        并记录可能的断点。
        this step follows the function of "separate_nd_seq"
        '''
        #   first, get all pts within links from both side
        #   Then return those pts along the sequence from o to d
        seq_head = [o_nd] + seq_o2d + [d_nd]
        seq_tail = [o_nd] + seq_d2o + [d_nd]

        #   status is the source of those pts, if the pt is mapped to a real node, then nid; if it is from a mid pt within
        #           links , then None
        #   record: status of links whether this link has been passed or not
        #   storage: new nodes geometry
        head_status, record, head_storage = self.get_inter_nd(seq_head)
        tail_status, record, tail_storage = self.get_inter_nd(seq_tail, record)
        if len(head_storage) > len(tail_storage):
            return head_storage, head_status, tail_storage, tail_status
        else:
            return tail_storage, tail_status, head_storage, head_status

    def get_inter_nd(self, seq_nd: list, stack: set = None, mode: str = "shortest"):
        '''

        :param seq_nd: （list）seq_o2d or seq_d2o
        :param stack: 记录那些已经被选过的link，避免重复的选择
        :return:status: dict{}
        '''
        if stack is None:
            stack = set()
        step = 0
        times = len(seq_nd) - 1  # check step and step+1 Thus total length need -1
        status = {0: seq_nd[0]}
        #   初始点入列
        seq_nd_storage = [self.node_storage[seq_nd[0]]]
        while step < times:
            nd1 = seq_nd[step]
            nd2 = seq_nd[step + 1]
            inter_link, stack = self.get_link_between_nds(nd1, nd2, stack, mode)
            #   get all pts within the inter_link
            inter_points = Points_in_Line(self.link_storage[inter_link])
            if self.node_storage[nd1].distance(inter_points[0]) > self.node_storage[nd2].distance(inter_points[0]):
                #   keep inter_points sequence from o to d
                inter_points.reverse()
            for i in range(1, len(inter_points) - 1):
                status[len(seq_nd_storage)] = 'mid_pt'
                seq_nd_storage.append(inter_points[i])
            status[len(seq_nd_storage)] = nd2
            seq_nd_storage.append(self.node_storage[nd2])
            step = step + 1
        return status, stack, seq_nd_storage

    def get_link_between_nds(self, nd1, nd2, record=None, mode: str = "shortest", sub_graph=None):
        '''
        这个方法返回两个节点之间的link，但是存在两种情况，
        1   存在复数link时：
            1.1 mode 1：“updating mode”：直接将待合并的线段直接更新合并成同一线段
            1.2 mode 2：“shortest mode”：直接返回最小的link值
            1.3 mode 3：“All mode”：返回全部相连线段
        2   record中存储需要排除的情况
        :param nd1:
        :param nd2:
        :param record:需要排除的link id
        :param mode: 默认返回最短的相连线段
                1.1 mode 1：“updating”：直接将待合并的线段直接更新合并成同一线段
                1.2 mode 2：“shortest”：直接返回最小的link值
                1.3 mode 3：“all”：返回全部相连线段(return list)
        :param sub_graph: default as None, None的时候不需要更新子图拓扑，非None则更新
                            subgraph[0]表示从节点向链接的映射
                            subgraph[1]表示从链接向节点饿映射
        :return:
        '''
        #   expecting to return a single line between nd1 and nd2
        #   case: 1 return the link not in Record
        #         2 return the shorter one
        if record is None:
            record = set()
        inter_link_sets = self.node[nd1].intersection(self.node[nd2])
        if len(inter_link_sets) > 1:
            #   one cycle contains only two nds
            if mode == "shortest":
                '''
                shortest mode: return the minimal length link
                '''
                shortest_link = None
                shortest_len = None
                for li in inter_link_sets:
                    if li not in record:
                        if shortest_link is None:
                            #   initialization
                            shortest_link = li
                            shortest_len = self.link_storage[li].length
                        elif self.link_storage[li].length < shortest_len:
                            shortest_link = li
                            shortest_len = self.link_storage[li].length
                if shortest_link is not None:
                    record.add(shortest_link)
                    return shortest_link, record
                else:  # 防御性语句
                    raise Exception("Not Found Non-Duplication Path, All possible links has been recorded")
            elif mode == "all":
                record = record.union(inter_link_sets)
                return inter_link_sets, record
            # elif mode == "updating":
            #     #   record 中现在没有输出需求，所以并没有更新record 中的信息，如果未来需要输出record
            #     #   则需要再updating中去除被
            #     if sub_graph is None:
            #         updated_link, new_link_id = merging_duplicate_links(nd1, nd2)
            #         for i in new_link_id:
            #             record.add(i)
            #         return updated_link, record
            #     else:
            #         sub_node = sub_graph[0]
            #         sub_link = sub_graph[1]
            #         updated_link, new_link_id, sub_node, sub_link = merging_duplicate_links(nd1, nd2,[sub_node, sub_link])
            #         for i in new_link_id:
            #             record.add(i)
            #         return updated_link, record, sub_node, sub_link

            else:
                raise Exception("Check mode. Not accept the {} mode as input".format(mode))
        elif len(inter_link_sets) == 1:  # Most Cases
            out_link = inter_link_sets.pop()
            record.add(out_link)
            return out_link, record
        else:  # 防御性语句
            raise Exception('No link between n1 and n2 can be found')

    def cycle_simplify(self, min_cycle, node_graph, link_graph, node_manage):
        '''
        修正minimal_cycle_list，并将其修正成单一线段。
        '''
        sub_n = deepcopy(node_graph)
        sub_l = deepcopy(link_graph)
        #   先检测每个node的属性
        #   找到转折点先
        manage = dict()
        endnd = list()  # 记录所有已经遍历过的位置
        node_angle = dict()  # minimal_circle_node_id : angle
        prev_id = min_cycle[-2]
        for idx in range(len(min_cycle) - 1):
            #   minimal_cycle是针对节点的有序列表
            #   前一个点和后一个点取交集分别跟当前点取交集，即是其相连的线段
            current_id = min_cycle[idx]
            aft_id = min_cycle[idx + 1]
            l1 = sub_n[current_id].intersection(sub_n[prev_id])
            l2 = sub_n[current_id].intersection(sub_n[aft_id])
            angle = AngleCal_turning(self.link_storage[l1.pop()],
                                     self.link_storage[l2.pop()],
                                     self.node_storage[current_id])
            node_angle[current_id] = angle
            if angle > 0.6:  # 说明是转折点
                endnd.append(current_id)
            manage[current_id] = node_manage[current_id]
            prev_id = current_id

        def quick_sort_endnd(s):
            if len(s) <= 1:
                return s
            else:
                pivot = s[0]
                less = [x for x in s[1:] if node_angle[x] <= node_angle[pivot]]
                greater = [x for x in s[1:] if node_angle[x] > node_angle[pivot]]
            return quick_sort_endnd(less) + [pivot] + quick_sort_endnd(greater)

        if len(endnd) != 2:
            end_nd_tem = quick_sort_endnd([i for i in min_cycle[1:]])
            endnd = [end_nd_tem[-1], end_nd_tem[-2]]

        if len(min_cycle) > 3:
            o_nd, d_nd, seq_o2d, seq_d2o = self.separate_nd_seq(min_cycle, endnd)
            mid_nodes = seq_o2d + seq_d2o
            l_storage, l_status, s_storage, s_status = self.find_internal_nd(seq_o2d, seq_d2o, o_nd, d_nd)
            new_link, link_status, upd_nd_stg, upd_nd_map, end_update = self.construct_mid_link(l_storage,
                                                                                           l_status,
                                                                                           s_storage,
                                                                                           s_status,
                                                                                           mid_nodes)

            #   开始更新新产生的节点与link里的信息
            #   先更新节点
            for geom_id, target_nid_set in upd_nd_map.items():
                new_node_geom = upd_nd_stg[geom_id]
                #   任意取出一个点作为target_nid，其他部分的nodeid信息全部迁移到当前点
                target_nid = target_nid_set.pop()
                if len(target_nid_set):  # 说明有其他点需要被合并
                    sub_n, sub_l = self.update_replacing_node(target_nid, target_nid_set, [sub_n, sub_l])
                    for i in target_nid_set:
                        if min_cycle[0] == i:
                            min_cycle.pop(0)
                            min_cycle.append(min_cycle[0])
                        if i in min_cycle:
                            min_cycle.remove(i)
                    for key, val in link_status.items():
                        link_status[key] = val.difference(target_nid_set)
                #   更新node_geometry和其相连的link端点位置
                self.update_change_node_geom(target_nid, new_node_geom)
            #   更新终点；看终点是否需要更新
            for end, replaced_set in end_update.items():
                sub_n, sub_l = self.update_replacing_node(end, replaced_set, [sub_n, sub_l])
                for i in replaced_set:
                    if min_cycle[0] == i:
                        min_cycle.pop(0)
                        min_cycle.append(min_cycle[0])
                    if i in min_cycle:
                        min_cycle.remove(i)
            #   再更新link
            current_step = 0
            next_step = 1
            while next_step < len(min_cycle):
                link_list = self.node[min_cycle[current_step]].intersection(self.node[min_cycle[next_step]])
                for li_idx in link_list:
                    #   删除拓扑信息
                    _, sub_n, sub_l = self.del_links(li_idx, [sub_n, sub_l])
                current_step = current_step + 1
                next_step = next_step + 1
            #   将新找到的线的拓扑和几何信息更新
            _, sub_n, sub_l = self.add_links(new_link, link_status, [sub_n, sub_l])
            return sub_n, sub_l
        else:  # 直接由两个点构成的link
            o_nd, d_nd = min_cycle[0], min_cycle[1]
            _, _, sub_n, sub_l = self.merging_duplicate_links(o_nd, d_nd, [sub_n, sub_l])
            return sub_n, sub_l

    def merging_duplicate_links(self, nd1, nd2, sub_graph=None):
        '''
        This function is to merge self.links when multiple self.links exist between nd1 and nd2
        And update the messages in graphs and subgraphs.
        :param nd1:
        :param nd2:
        :param sub_graph: default as None, None的时候不需要更新子图拓扑，非None则更新
                            subgraph[0]表示从节点向链接的映射
                            subgraph[1]表示从链接向节点饿映射
        :return:
        '''
        duplication_links = list(self.node[nd1].intersection(self.node[nd2]))
        if len(duplication_links) < 2:
            raise Exception("No duplication founded")
        times = len(duplication_links)
        step = 1
        l1 = duplication_links[0]
        l1_inter_points = Points_in_Line(self.link_storage[l1])
        new_link = list()
        while step < times:
            l2 = duplication_links[step]
            l2_inter_points = Points_in_Line(self.link_storage[l2])
            if l1_inter_points[0].distance(l2_inter_points[0]) > l1_inter_points[0].distance(l2_inter_points[-1]):
                #   保证 l1 和 l2 顺序是统一的
                #   由于只跟两个点有关, 不需要更新拓扑关系
                l2_inter_points.reverse()
            l1_inter_status = dict()
            l2_inter_status = dict()
            for idx, i in enumerate(l1_inter_points):
                l1_inter_status[idx] = "mid_pt"
            for idx, i in enumerate(l2_inter_points):
                l2_inter_status[idx] = "mid_pt"
            new_link = self.construct_mid_link(l1_inter_points, l1_inter_status, l2_inter_points, l2_inter_status, None)
            l1_inter_points = Points_in_Line(new_link)
            step = step + 1
        '''
        更新拓扑关系
        '''
        if sub_graph is None:
            self.del_links(duplication_links)
            newid = self.add_links([new_link], {0: {nd1, nd2}})
            return new_link, newid
        else:
            sub_node = sub_graph[0]
            sub_link = sub_graph[1]
            _, sub_node, sub_link = self.del_links(duplication_links, [sub_node, sub_link])
            newid, sub_node, sub_link = self.add_links([new_link], {0: {nd1, nd2}}, [sub_node, sub_link])
            return new_link, newid, sub_node, sub_link

    def construct_mid_link(self, l_storage, l_status, s_storage, s_status, junctions):
        '''
        find mid pts between l_storage and s_storage
        first, from the long storage to find their corresponding short storage
        reverse to check whether any nds from short storage missed
        add missed mid pt into the right position of new link
        :param l_storage: the longer pt list with geometry in it
        :param l_status: the source of the pts from self.node or "mid pts"
        :param s_storage: the shorter pt list with geometry
        :param s_status: the source of the pts from self.node or "mid pts"
        :return: new_links: storage the geometry of new self.links
        :return: links_status: each link self.links to which nodes
        '''

        #   storage the pts geometry between two self.node and splite at self.node
        tem_links = [{'status': l_status[0], 'geometry': l_storage[0]}]
        #   o pts and d pts are not considering within the self.links
        l_step = 1
        total_step = len(l_storage) - 1
        #   where last s was found in order to make sure no nd was missing
        s_step = 1
        '''
            找到合并后的geometry
        '''
        while l_step < total_step:
            current_nd = l_storage[l_step]
            #   s_idx is latest founded closest pt in s_storage
            #   s_idx to s_step are skipping ranges

            s_idx, another_nd = get_closest_nd(current_nd, s_storage)
            #   new_nd_geom is one pt on new_link
            new_nd_geom = get_mid_nd(current_nd, another_nd)
            '''
                check skips between s_idx and s_steps
                if there is a nd in s_storage, then split at nd, tem_link into new_link and tem_link initial at nd
            '''
            if s_idx - s_step > 1:
                #   there at least one nd
                for i in range(s_step + 1, s_idx):
                    # check i pt is from nodes or not
                    if s_status[i] == 'mid_pt':
                        #   mid_pt should be skipped
                        continue
                    #   i is from self.node
                    #   first update fixed_nd
                    #   tem_nd store the geometry of i
                    tem_nd = s_storage[i]
                    tem_l_idx, tem_another_nd = get_closest_nd(tem_nd, l_storage)
                    tem_new_nd = get_mid_nd(tem_nd, tem_another_nd)
                    tem_links.append({'status': 'mid_pt', 'geometry': tem_new_nd})  # todo 这是旧的

            else:
                #   When no skipped nds within the
                tem_links.append({'status': 'mid_pt', 'geometry': new_nd_geom})

            s_step = s_idx
            l_step = l_step + 1
        #   put endpt into the link and d must be self.node
        tem_links.append({'status': l_status[len(l_storage) - 1], 'geometry': l_storage[-1]})
        #   tem_link中已经获取了完整的geometry
        tem_link_geom = LineString([i['geometry'] for i in tem_links])
        if junctions is None:
            return tem_link_geom
        '''
            将中间的node找到新geometry上的最近点
        '''
        junctions = junctions.copy()
        new_node_storage = list()
        new_node_map = dict()
        for nid in junctions:
            old_geom = self.node_storage[nid]
            new_geom = nearest_pt(old_geom, tem_link_geom)
            junction_id = len(new_node_storage)
            #   判断junction之间是否过近导致合并
            for id, geom in enumerate(new_node_storage):
                if new_geom.distance(geom) < 20:
                    #   需要被合并
                    junction_id = id
            if junction_id != len(new_node_storage):
                #   说明junction之间需要被合并
                new_node_map[junction_id].add(nid)
            else:
                new_node_map[junction_id] = {nid}
                new_node_storage.append(new_geom)

        '''
            控制junction按照顺序插入
        '''
        current_links = tem_links
        end_replacing = dict()
        for junc_id, junc_geom in enumerate(new_node_storage):
            current = 0
            next = 1
            tem_links = deepcopy(current_links)
            current_links = [tem_links[0]]
            used_status = 0
            while next < len(tem_links):
                if used_status == 0:
                    current_geom = tem_links[current]["geometry"]
                    next_geom = tem_links[next]["geometry"]
                    '''
                        一下内容是为了避免现有nd——list中与junction几何完全重合
                    '''
                    if current_geom.distance(junc_geom) == 0:
                        if tem_links[current]['status'] == 'mid_pt':
                            #   currentpt直接被juncnd 取代
                            used_status = 1
                            current_links.append({'status': new_node_map[junc_id].copy(),
                                                  'geometry': junc_geom})
                            current = current + 1
                            next = next + 1
                            continue
                        else:  # 说明currentpt和juncnd都是nd，他们需要被合并.用junc替代currentpt
                            used_status = 1
                            current_links.append({'status': new_node_map[junc_id].copy(),
                                                  'geometry': junc_geom})
                            new_node_map[junc_id].add(tem_links[current]['status'])
                            current = current + 1
                            next = next + 1
                            continue
                    if next_geom.distance(junc_geom) == 0:
                        # nextpt如果是重合的
                        current_links.append(tem_links[next])
                        current = current + 1
                        next = next + 1
                        continue
                    '''
                        一下内容是为了避免现有nd——list中与junction
                    '''

                    angle = AngleCal_nd(current_geom, junc_geom, next_geom)
                    if angle > 0 or angle is None:
                        #   junc位于current和next之间
                        used_status = 1
                        current_links.append({'status': new_node_map[junc_id].copy(),
                                              'geometry': junc_geom})
                current_links.append(tem_links[next])
                current = current + 1
                next = next + 1
            if used_status == 0:
                #   说明之前没有找到合适的插入点
                #   一般发生在拐角处检测他是更靠近哪一个端点
                #   直接让端点信息取代原有信息
                if l_storage[0].distance(junc_geom) < l_storage[-1].distance(junc_geom):
                    if l_status[0] not in end_replacing:
                        end_replacing[l_status[0]] = new_node_map[junc_id]
                    else:
                        end_replacing[l_status[0]] = end_replacing[l_status[0]].union(new_node_map[junc_id])
                else:
                    if l_status[len(l_storage) - 1] not in end_replacing:
                        end_replacing[l_status[len(l_storage) - 1]] = new_node_map[junc_id]
                    else:
                        end_replacing[l_status[len(l_storage) - 1]] = end_replacing[l_status[len(l_storage) - 1]].union(
                            new_node_map[junc_id])

        '''
            将split link 分段输出
        '''
        #   links_status record the new_link idx in new_link and map to nds
        links_status = {0: {l_status[0]}}
        #   new_links store the geometry of new self.links
        new_links = list()
        tem_links = [current_links[0]['geometry']]
        for i in range(1, len(current_links) - 1):
            current_status = current_links[i]['status']
            current_geom = current_links[i]['geometry']
            if current_status == 'mid_pt':
                tem_links.append(current_geom)
            else:
                tem_links.append(current_geom)
                links_status[len(new_links)] = links_status[len(new_links)].union(current_status)
                new_links.append(LineString(tem_links))
                links_status[len(new_links)] = current_status
                tem_links = [current_geom]

        links_status[len(new_links)].add(l_status[len(l_storage) - 1])
        tem_links.append(current_links[-1]['geometry'])
        new_links.append(LineString(tem_links))

        return new_links, links_status, new_node_storage, new_node_map, end_replacing

    def separate_nd_seq(self, minimal_seq_list, endnd):
        #   将minimal_cycle给划分成两部分
        #   以endnd[0]为“O”, endnd[1]为“D”(Origin and Destination)
        #   O to D : 找到O之后的D之前的所有node
        #   剩下的部分是D to O
        #   最后为了保证方向一直，都从起点到终点作为， D to O进行reverse， 也成为O to D的顺序，构成环的另外一边
        o = None
        d = None
        seq_o2d = LinkedList()
        seq_d2o = LinkedList()
        length = None
        minimal_seq_list = minimal_seq_list.copy()
        minimal_seq_list.pop()
        for nd in minimal_seq_list:
            if nd in endnd:
                #   第一个遇到的是 O
                if o is None:
                    o = nd
                    #   避免minimal——Seq和之前endnd的顺序不一致
                    if nd != endnd[0]:
                        #   应该是o的点但在endnd顺序中错位了
                        endnd[1] = endnd[0]
                        endnd[0] = o
                    continue
                else:  # 永远是O先更新；then D更新，因此不用考虑endnd的顺序问题，直接赋值即可
                    d = nd
                    length = seq_d2o._length  # length是为了方便seq_d2o后置插入
                    continue
            if o is None and d is None:
                #   在找到o之前前置插入
                seq_d2o.insert(nd, 0)
            elif o is not None and d is None:
                #   OD 之间，seqo2d正向插入
                seq_o2d.append(nd)
            elif o is not None and d is not None:
                if length == None:
                    raise TypeError("Length should not be None, Please check the len(endnd):length = "
                                    "seq_d2o._length was not updated")
                seq_d2o.insert(nd, length)
        return o, d, seq_o2d.to_list(), seq_d2o.to_list()

    def update_remove_from_link(self, l_idx: int, l_dict: dict, n_dict: dict):
        if l_idx not in l_dict:
            return l_dict, n_dict
        if len(l_dict[l_idx]) == 0:
            #   Preventing：l_idx has been cleared
            return l_dict, n_dict
        for n_idx in l_dict[l_idx]:
            if l_idx in n_dict[n_idx]:
                n_dict[n_idx].remove(l_idx)
        if l_idx in l_dict:
            l_dict.pop(l_idx)
        return l_dict, n_dict

    def update_remove_from_node(self, n_idx: int, l_dict: dict, n_dict: dict):
        if n_idx not in n_dict:
            return l_dict, n_dict
        if len(n_dict[n_idx]) == 0:
            #   Preventing：n_idx has been cleared
            return l_dict, n_dict
        for l_idx in n_dict[n_idx]:
            if n_idx in l_dict[l_idx]:
                l_dict[l_idx].remove(n_idx)
        if n_idx in n_dict:
            n_dict.pop(n_idx)
        return l_dict, n_dict

    def straighten_links(self, pt_list):
        #   首先保证初始点和终结点
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
    #%%
    def Line_forward_Line(self, lid, endli=None):
        if endli is None:
            endli = set()
        if lid.__class__ == int:
            li_list = list(self.links[lid])
            endli.add(lid)
        elif lid.__class__ == list:
            li_list = set()
            for i in lid:
                li_list = li_list.union(self.links[i])
                endli.add(i)
        else:
            raise TypeError("Check input types: nid should be list and endpt should be set")
        result = set()
        for i in li_list:
            for j in self.node[i]:
                if j not in endli:
                    result.add(j)
        result = list(result)
        return result, endli

    def nodeseq_to_linkseq(self, seq):
        '''
        seq should be closed or ringed path
        :param seq:
        :return:
        '''
        prevnd = seq[0]
        linkseq = set()
        for currentnd in seq:
            if currentnd == prevnd:
                continue
            linkseq = linkseq.union(self.node[prevnd].intersection(self.node[currentnd]))
            prevnd = currentnd
        return set(linkseq)

    def Node_forward_Node(self, nid, node_graph, link_graph, endpt=None, endlk=None):
        if endpt is None:
            endpt = set()

        if endlk is None:
            endlk = set()
        if nid.__class__ == int:
            li_list = list(node_graph[nid])
            endpt.add(nid)
        elif nid.__class__ == list:
            li_list = set()
            for i in nid:
                li_list = li_list.union(node_graph[i])
                endpt.add(i)
        else:
            raise TypeError("Check input types: nid should be list and endpt should be set")
        res_node = set()
        for i in li_list:
            if i not in endlk:
                endlk.add(i)
                for j in link_graph[i]:
                    if j not in endpt:
                        res_node.add(j)
        res_node = list(res_node)

        return res_node, endpt, endlk

    def Link_between_Nodes(self, nd1, nd2, record=None, mode='shortest'):
        '''
        这个方法返回两个节点之间的link，但是存在两种情况，
        1   存在复数link时：
            1.2 mode 2：“shortest mode”：直接返回最小的link值
            1.3 mode 3：“All mode”：返回全部相连线段
        2   record中存储需要排除的情况
        :param nd1:
        :param nd2:
        :param record:需要排除的link id
        :param mode: 默认返回最短的相连线段
                1.1 mode 1：“updating”：直接将待合并的线段直接更新合并成同一线段
                1.2 mode 2：“shortest”：直接返回最小的link值
                1.3 mode 3：“all”：返回全部相连线段(return list)
        '''
        if record is None:
            record = set()
        inter_link_sets = self.node[nd1].intersection(self.node[nd2])
        if len(inter_link_sets) > 1:
            #   one cycle contains only two nds
            if mode == "shortest":
                '''
                shortest mode: return the minimal length link
                '''
                shortest_link = None
                shortest_len = None
                for li in inter_link_sets:
                    if li not in record:
                        if shortest_link is None:
                            #   initialization
                            shortest_link = li
                            shortest_len = self.link_storage[li].length
                        if self.link_storage[li].length < shortest_len:
                            shortest_link = li
                            shortest_len = self.link_storage[li].length
                if shortest_link is not None:
                    record.add(shortest_link)
                    return shortest_link, record
                else:  # 防御性语句
                    raise Exception("Not Found Non-Duplication Path, All possible self.links has been recorded. "
                                    "And the problem self.node ids are {} and {}.the intersect self.links are {}"
                                    .format(nd1, nd2, record))
            elif mode == "all":
                record = record.union(inter_link_sets)
                return inter_link_sets, record
        elif len(inter_link_sets) == 1:
            outli = inter_link_sets.pop()
            record.add(outli)
            return outli, record
        else:
            return None, record

    def Link_to_Node(self, lid, nid, link_garph, block=None):
        '''
        return another self.node id of current link
        :param lid:
        :param nid:
        :param link_garph:
        :return: mode 'multi_pt': return the list of self.node ids and the string 'multi_pt'
                 mode 'pt': find only one nid and return it and the string 'pt'
                 mode 'end': cannot find any nid and return None and the string 'end'
        '''
        if block is None:
            block = {nid}
        else:
            block.add(nid)
        other_node = list()
        for nd in link_garph[lid]:
            if nd not in block:
                other_node.append(nd)
        if len(other_node) > 1:
            return "multi_pt", other_node
        elif len(other_node) == 1:
            return "pt", other_node[0]
        else:
            return 'end', None

    def update_change_node_geom(self, nid: int, new_geom: Point):
        if nid not in self.node:
            raise IndexError("input nid: {} not exist in self.node".format(nid))
        if self.node_storage[nid] is None:
            raise ValueError("{} self.node is none".format(nid))
        #   优先更新geom本身
        # oringin_node = self.node_storage[nid]
        #   更新相连link的geom
        for li in self.node[nid]:
            pt_list = Points_in_Line(self.link_storage[li])
            #   检测新增的点应该是头还是尾
            if pt_list[0].distance(new_geom) < pt_list[-1].distance(new_geom):
                #   应当加在头
                pt_list[0] = new_geom
            else:
                #   应当加在尾
                pt_list[-1] = new_geom
            #   重新更新li的geom
            if len(pt_list) > 3:
                pt_list = straighten_links(pt_list)
            self.link_storage[li] = LineString(pt_list)
        self.node_storage[nid] = new_geom

    def update_replacing_node(self, target_nid, replaced_set, sub_graph=None):
        if replaced_set.__class__ is list:
            replaced_set = set(replaced_set)

        #   提取出对应的link id,
        for replaced_id in replaced_set:
            if replaced_id not in self.node:
                #   说明已经被更新过了
                continue
            link_set = self.node[replaced_id]
            for li_idx in link_set:
                self.links[li_idx].remove(replaced_id)
                self.links[li_idx].add(target_nid)
                self.node[target_nid].add(li_idx)
            self.node.pop(replaced_id)
            self.node_storage[replaced_id] = None
        if sub_graph is not None:
            node_graph = sub_graph[0]
            link_graph = sub_graph[1]
            for replaced_id in replaced_set:
                if replaced_id not in node_graph:
                    #   说明已经被更新过了
                    continue
                link_set = node_graph[replaced_id]
                for li_idx in link_set:
                    link_graph[li_idx].remove(replaced_id)
                    link_graph[li_idx].add(target_nid)
                    node_graph[target_nid].add(li_idx)
                node_graph.pop(replaced_id)
            return node_graph, link_graph

    def add_links(self, l_geom_list: list, mapping_relations, sub_graph=None):
        if l_geom_list.__class__ is int:
            l_geom_list = [l_geom_list]
        if l_geom_list.__class__ is not list:
            raise TypeError("plz check the the type of 'l_id_list'")

        newid = list()
        #   先把l_id_list单独的编码获得
        for tem_idx, link_geometry in enumerate(l_geom_list):
            #   存储并获得ID
            real_idx = len(self.link_storage)
            newid.append(real_idx)
            self.link_storage.append(link_geometry)
            #   更新link的拓扑信息
            self.links[real_idx] = mapping_relations[tem_idx]
            node_set = mapping_relations[tem_idx]
            #   更新节点拓扑信息
            for i in node_set:
                if real_idx not in self.node[i]:
                    self.node[i].add(real_idx)
            if sub_graph is not None:
                #   需要对子图信息进行更新
                sub_node = sub_graph[0]
                sub_link = sub_graph[1]
                #   更新子图link的拓扑信息
                sub_link[real_idx] = {i for i in node_set if i in sub_node}
                #   更新子图node的拓扑信息
                for i in sub_link[real_idx]:
                    if real_idx not in sub_node[i]:
                        sub_node[i].add(real_idx)

        if sub_graph is None:
            return newid
        else:
            return newid, sub_node, sub_link

    def del_links(self, l_id_list: list, sub_graph=None):
        '''

        :param l_id_list:
        :param sub_graph:   默认为None，非None时接收以长度为2的list形式的输入
                            subgraph[0]表示从节点向链接的映射
                            subgraph[1]表示从链接向节点饿映射
        :return:
        '''
        if l_id_list.__class__ is int:
            l_id_list = [l_id_list]
        if l_id_list.__class__ is not list:
            raise TypeError("plz check the the type of 'l_id_list'")
        del_relationships = dict()  # 保留被删除的links相关关系
        for li_id in l_id_list:
            node_set = self.links[li_id]
            del_relationships[li_id] = node_set
            if li_id in self.links:
                self.links.pop(li_id)
            for n_id in node_set:
                self.node[n_id].discard(li_id)
            if sub_graph is not None:
                sub_node = sub_graph[0]
                sub_link = sub_graph[1]
                if li_id in sub_link:
                    sub_link.pop(li_id)
                for n_id in node_set:
                    if n_id in sub_node:
                        sub_node[n_id].discard(li_id)
            self.link_storage[li_id] = None
        if sub_graph is None:
            return del_relationships
        else:
            #   子图的拓扑信息更新
            return del_relationships, sub_node, sub_link

#%%
    def check_graph_validation(self):
        status_n = self.check_node_validation()
        status_l = self.check_link_validation()
        if status_l and status_n:
            # 都通过了检测
            return
        else:
            self.check_graph_validation()

    def check_link_validation(self, status=1):
        for l_id in range(len(self.link_storage)):
            if self.link_storage[l_id] is None:
                continue
            if len(self.links[l_id]) == 0:
                self.link_storage[l_id] = None
                self.links.pop(l_id)
                status = 0  # 说明这次检查中有不合格的情况
                continue
            if self.link_storage[l_id].length == 0:
                status = 0  # 说明这次检查中有不合格的情况
                self.link_storage[l_id] = None
                keep_nid = self.links[l_id].pop()
                for n_id in self.links[l_id]:
                    self.node[keep_nid] = self.node[keep_nid].union(self.node[n_id])
                    for l in self.node[n_id]:
                        self.links[l].discard(n_id)
                        self.links[l].add(keep_nid)
                    if n_id in self.node:
                        self.node.pop(n_id)
                        self.node_storage[n_id] = None
                self.node[keep_nid].discard(l_id)
                if l_id in self.links:
                    self.links.pop(l_id)
        return status

    def check_node_validation(self, status=1):
        for n_idx, stg in enumerate(self.node_storage):
            if stg is None:
                continue
            if len(self.node[n_idx]) < 2:
                #   尽端点不存在
                for i in self.node[n_idx]:
                    self.links[i].discard(n_idx)
                self.node_storage[n_idx] = None
                self.node.pop(n_idx)
                status = 0  # 说明这次检查中有不合格的情况
        return status

