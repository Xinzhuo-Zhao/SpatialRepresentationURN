from tqdm import tqdm
from geometry_processing import *
from Angle import *
from copy import deepcopy
from data_structure import *
from data_input import *

class GraphBuilder:
    def __init__(self, link_storage, link, node, node_storage, filename):
        self.link_storage = link_storage
        self.link = link
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
        save_data_structure(self.link, self.filename + "Graph_links.pkl")
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
            if len(self.link[link_id]) == 0:
                self.link_storage[link_id] = None
                if len(self.link[link_id]) > 0:
                    kept = self.link[link_id].pop()
                    for nd in self.link[link_id]:
                        self.node[kept] = self.node[kept].union(self.node[nd])
                        self.node.pop(nd)
                        self.node_storage[nd] = None
                    self.node[kept].discard(link_id)
                self.link.pop(link_id)
        
        for l, n in tqdm(self.link.items()):
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
        for l, n in tqdm(self.link.items()):
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
                self.link[current_index] = set()
            self.link_storage[l] = None
            for relate_nd in self.link[l]:
                for relate_li in add_idx:
                    if check_endpt(self.node_storage[relate_nd], self.link_storage[relate_li]):
                        self.node[relate_nd].add(relate_li)
                        self.link[relate_li].add(relate_nd)
                self.node[relate_nd].remove(l)
                if l in self.node[relate_nd]:
                    raise ValueError("unfound related links")
            self.link.pop(l)
    
    def run_linesimplification(self, link_id):
        spt_list, spt_node, nd2li, status = LineSimplication(self.link_storage[link_id], link_id)
        if not status:
            return
        if len(spt_list) == 1:
            self.link_storage[link_id] = spt_list[0]
            return
        related_nd = self.link.pop(link_id)
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
                #   nd_idx is the start or end pt within a line
                #   There are already corresponding points in the original graph
                if len(li_set) != 1:
                    raise ValueError('未知错误：端点有多条线')
                li = list(li_set)[0]
                for nd in related_nd:
                    #   nd is the index of self.node_storage
                    if current_node.distance(self.node_storage[nd]) < 1:
                        #   the current node and the original node are the same point
                        if li not in dup_record:
                            #   li is traversed for the first time
                            current_lidx = len(self.link_storage)
                            dup_record[li] = current_lidx
                            self.link[current_lidx] = {nd}
                            self.node[nd].add(current_lidx)
                            self.link_storage.append(spt_list[li])
                        else:
                            #   li was traversed before
                            current_lidx = dup_record[li]
                            self.link[current_lidx].add(nd)
                            self.node[nd].add(current_lidx)
                        break
            else:
                #   this is the intermediate and new established pt
                self.node[current_nidx] = set()
                self.node_storage.append(current_node)
                for li in li_set:
                    #   li is the index of spt_list
                    if li not in dup_record:
                        #   li is traversed for the first time
                        current_lidx = len(self.link_storage)
                        dup_record[li] = current_lidx
                        self.link[current_lidx] = {current_nidx}
                        self.node[current_nidx].add(current_lidx)
                        self.link_storage.append(spt_list[li])
                    else:
                        #   li was traversed before
                        current_lidx = dup_record[li]
                        self.link[current_lidx].add(current_nidx)
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
                    self.link[l].discard(n_idx)
                    new_node = new_node.union(self.link[l])
                    self.link_storage[l] = None
                    for n in self.link[l]:
                        self.node[n].discard(l)
                        self.node[n].add(idx)
                    self.link[l] = None
                self.link[idx] = new_node
                self.node[n_idx] = None
                self.node_storage[n_idx] = None

    def run_remove_anomalous_shortlinks(self, l_idx, n_set, length_limits=50):
        if self.link_storage[l_idx].length < length_limits:
            depth1, depth_end1 = self.Line_forward_Line(l_idx)
            depth2, _ = self.Line_forward_Line(depth1, depth_end1)
            related_li = depth2 + depth1
            length_record = list()
            for li_geo_idx in related_li:
                length_record.append(self.link_storage[li_geo_idx].length)
            if len(length_record) == 0:
                return
            q1, q3 = np.percentile(length_record, [33, 75])
            meanlength = np.mean(length_record)
            if (meanlength / self.link_storage[l_idx].length > 4) and (self.link_storage[l_idx].length < q1):
                #   abnormal too short link identification and aggregation.
                pt_idx = len(self.node_storage)
                #   aggregating at their centroid
                new_pt = midPt(self.link_storage[l_idx])
                self.node_storage.append(new_pt)
                new_node = set()
                for i in n_set:
                    self.node[i].discard(l_idx)
                    new_node = new_node.union(self.node[i])
                    self.node.pop(i)
                    self.node_storage[i] = None
                self.node[pt_idx] = new_node
                # modify relevant links
                for fixline in new_node:
                    # geometric information update
                    l1 = self.link_storage[fixline]
                    if l1.length == 0:
                        self.link_storage[fixline] = None
                        if fixline in self.link:
                            for n in self.link[fixline]:
                                self.node[n].discard(fixline)
                            self.link.pop(fixline)
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
                    #   Graph connectivity update
                    self.link[fixline] = self.link[fixline] - n_set
                    self.link[fixline].add(pt_idx)
                self.link[l_idx] = None
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
            return
        while brk_pts:
            n = brk_pts.pop()
            if pts_status[n]['status'] != 'breaknd':
                continue

            path_list, pts_status, brk_pts = self.explore_roads(n, pts_status, brk_end=brk_pts)
            #   Status indicates whether other merging groups have been found.
            #   If so, update the status and update the pts_status and brk_pts that need to be found
            status = 0

            if len(path_list):  #   if any cycle is identified
                #   update cycles into  merging dict and merging storage
                for path in path_list:
                    newlinks = self.nodeseq_to_linkseq(path)
                    li_set = li_set.union(newlinks)
                    for i in newlinks:
                        if i in self.merging_dict:
                            #   this cycle is twisted with other cycles
                            status = 1
                            li_set = self.merging_storage[self.merging_dict[i]].union(li_set)
                    self.merging_storage[storage_id] = li_set
                    for i in li_set:
                        #   对于新更新进来的线段
                        if i not in self.merging_dict:
                            #   i is a new founded link in cycles
                            self.merging_dict[i] = storage_id
                        elif self.merging_dict[i] != storage_id:
                            #   i is also part of other cycles
                            self.merging_storage[self.merging_dict[i]] = None
                            self.merging_dict[i] = storage_id
            if status:
                pts_status, brk_pts = self.check_breaknd(li_set, pts_status, brk_pts)

    def run_cycle_simplify(self):
        for idx, testing_g in tqdm(enumerate(self.merging_storage)):
            if testing_g is None:
                continue
            for i in testing_g:
                if i not in self.link:
                    self.merging_storage[idx] = None
                    break
            if self.merging_storage[idx] is None:
                continue
            sub_l, sub_n = self.pre_cycle_modification(testing_g)
            sub_node_graph = self.node_graph_construction(sub_n, sub_l)
            node_manage = self.nodetype_in_sub(sub_n, self.node)
            if len(sub_n) < 2:  # no cycle found
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
                self.link[current_index] = set()
            self.link_storage[l] = None
            for relate_nd in self.link[l]:
                for relate_li in add_idx:
                    if check_endpt(self.node_storage[relate_nd], self.link_storage[relate_li]):
                        self.node[relate_nd].add(relate_li)
                        self.link[relate_li].add(relate_nd)
                self.node[relate_nd].remove(l)
                if l in self.node[relate_nd]:
                    raise ValueError("未知错误：未找到继承关系线")
            self.link.pop(l)

    def Identify_MergingDifference(self, prev_merging_dict, prev_merging_stg, merging_d, merging_stg):
        for i in merging_d.keys():
            if i in prev_merging_dict:
                if prev_merging_stg[prev_merging_dict[i]] != merging_stg[merging_d[i]]:
                    #   need more aggregation
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

                if o not in self.link.keys():
                    self.link[o] = set()
                if d not in self.link.keys():
                    self.link[d] = set()
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
        new_pt = 1  # default the intersection pt is a new pt(new_pt = 1)
        for i in range(len(self.node_storage)):
            if self.node_storage[i].distance(intersectionPT) < 1:
                intersectionPT = self.node_storage[i]
                pt_idx[1] = i
                new_pt = 0
                break
        #   check the new pt locate at end of link or mid of link
        if new_pt == 0:
            self.node[pt_idx[1]].add(o)
            self.node[pt_idx[1]].add(d)
            self.link[o].add(pt_idx[1])
            self.link[d].add(pt_idx[1])
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
        else:
            if pt_idx[0] not in self.node.keys():
                self.node[pt_idx[0]] = set()
            self.node[pt_idx[0]].add(o)
            self.node[pt_idx[0]].add(d)
            self.link[o].add(pt_idx[0])
            self.link[d].add(pt_idx[0])
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
        queue = Queue()
        queue.Que_in(startnd)

        # this dict is to record the previous li and used to find the minimal distance trace
        # record distance from current li to startli
        prev = {startnd: {'prev': None, 'distance': 0, 'step': 0}}
        neighbor, _, _ = self.Node_forward_Node(startnd, self.node, self.link, end_positions)
        currentnd = None
        _, prevnd = self.Link_to_Node(list(cross_status[startnd]['self.links'])[0], startnd, self.link)
        end_targets = set()

        while not queue.Que_isEmpty():
            if currentnd is None:
                currentnd = queue.Que_out()
                blockednd.add(currentnd)
            elif prev[currentnd]['step'] < timelimits:
                currentnd = queue.Que_out()
                neighbor, _, _ = self.Node_forward_Node(currentnd, self.node, self.link, blockednd)
            else:
                #   search times more than the the timelimits, the serch should be stopped
                currentnd = queue.Que_out()
                neighbor, _, _ = self.Node_forward_Node(currentnd, self.node, self.link, blockednd)
                continue

            for nextnd in neighbor:
                inter_li, _ = self.Link_between_Nodes(currentnd, nextnd)
                new_distance = self.link_storage[inter_li].length + prev[currentnd]['distance']
                if (nextnd not in blockednd):  # this nextli is first searched
                    #   ensuring this search process wont backward
                    prev[nextnd] = dict()
                    prev[nextnd]['prev'] = currentnd
                    prev[nextnd]['distance'] = new_distance
                    prev[nextnd]['step'] = prev[currentnd]['step'] + 1
                    blockednd.add(nextnd)
                    if nextnd in brk_end:
                        end_targets.add(nextnd)
                    else:
                        queue.Que_in(nextnd)

                else:  # this nextli is searched before
                    if new_distance < prev[nextnd]['distance']:
                        #   it means shorter path is found
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
            if len(tem_path) < 3:  # no cycle found
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
        find all intermediate pts to form a cycle
        :param seq:
        :param cross:
        :return:
        '''
        seq = seq.copy()
        start = seq[0]
        end = seq[-1]
        whole_path_node = list()
        #   start pt and end pts both overlap with node position
        prevnd = start

        #   find all intermediate pts from start pt to end pt
        for nd, link_record in cross.items():
            if nd == start or nd == end:
                continue
            if link_record['self.links'].intersection(cross[start]['self.links']) and (nd not in seq):
                seq = [nd] + seq
            elif link_record['self.links'].intersection(cross[end]['self.links']) and (nd not in seq):
                seq.append(nd)
        if seq[0] != start and seq[-1] != end:
            #   if start and end are both updated,in this case, there might exist a junction at most.
            if not cross[seq[0]]["self.links"].intersection(cross[seq[0]]["self.links"]):
                #   only start and end cannot form a cycle path
                nd1, _, _ = self.Node_forward_Node(seq[0], self.node, self.link, endpt=set(seq))
                nd2, _, _ = self.Node_forward_Node(seq[-1], self.node, self.link, endpt=set(seq))
                missing_nd = set(nd1).intersection(set(nd2))
                if len(missing_nd) == 1:
                    seq.append(missing_nd.pop())
                else:
                    raise ValueError("Unable to close between {} and {}".format(seq[0], seq[-1]))

        for currentnd in seq:
            if currentnd == prevnd:
                whole_path_node.append(self.node_storage[prevnd])
                continue
            interlink: set = self.node[prevnd].intersection(self.node[currentnd])
            if len(interlink) > 1:
                #   Among intersecting links, return the shortest one
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
                #   need add in intermediate pts
                #   check the sequence of intermediate pts
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
            for nd in self.link[li]:
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
        for l in involving_links:
            sub_link[l] = set()
        for l in sub_link.keys():
            for n in self.link[l]:
                if len(self.node[n].intersection(involving_links)) > 1:
                    #   This node n has at least two connections and does not constitute a solitary node
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
            #   Because the premise for running this is that the terminal node is only connected to one line,
            #   so at this point, n_idx is already empty
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
        #   Distinguish between adjustment points and vanishing points
        #   (turning points are placed in the cycle for single segment management)
        #   Adjustment point: Refers to a point that is still connected to the outside world and requires careful
        #   consideration of spatial changes; Represented by 'fixed_nd'
        #   Vanishing point: a point that only exists within the subgraph; Represented by 'eliminated_nd'
        node_type = dict()
        for n_id, l_set in sub_n.items():
            #   The vanishing point is only connected to the internal links of the subgraph;
            #   Therefore,l_par_set - l_set = empty
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
                    first_n = stack[0]
                    second_n = stack[1]
                    if len(node_link[first_n].intersection(node_link[second_n])) > 1:
                        cycle.append((stack.copy(), record[top][stack[-1]]["total_length"]))
                        stack.pop()
                        continue
                    else:
                        stack.pop()
                        continue

                cycle.append((stack.copy(), record[top][stack[-1]]["total_length"]))
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
                #   return the shortest link at this junction
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
                #   add intermediate pts
                #   ensuring the sequence of intermediate pts
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
        Establish two node chains from O to D, one on the left and one on the right; The length of the input is unknown,
        but the order is from o to d; This method is to establish a midpoint link between two links from O to D; And
        record possible breakpoints.
        this step follows the function of "separate_nd_seq"
        '''
        #   first, get all pts within links from both side
        #   Then return those pts along the sequence from o to d
        seq_head = [o_nd] + seq_o2d + [d_nd]
        seq_tail = [o_nd] + seq_d2o + [d_nd]

        #   status is the source of those pts, if the pt is mapped to a real node, then nid; if it is from a mid pt
        #   within links , then None
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
        get inter nodes
        :param seq_nd: （list）seq_o2d or seq_d2o
        :param stack: Record the links that have already been selected to avoid duplicate selections
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
        This method returns the link between two nodes, but there are two situations,
        When there is a complex link:
        1.1 mode 1: "Updating mode": Directly update and merge the line segments to be merged into the same line segment
        1.2 mode 2: "shortest mode": directly returns the smallest link value
        1.3 mode 3: "All mode": returns all connected line segments
        2 situations that need to be excluded when storing in the record
        :param nd1:
        :param nd2:
        : paramrecord: Link IDs that need to be excluded
        : param mode: default returns the shortest connected line segment
        1.1 Mode 1: "Updating": Directly update and merge the line segments to be merged into the same line segment
        1.2 mode 2: "shortest": directly returns the smallest link value
        1.3 mode 3: "all": Return all connected line segments (return list)
        : param_ subgraph: default as None. When None, there is no need to update the subgraph topology, otherwise it will be updated
        Subgraph [0] represents the mapping from nodes to links
        Subgraph [1] represents the mapping from links to nodes
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
                else:
                    raise Exception("Not Found Non-Duplication Path, All possible links has been recorded")
            elif mode == "all":
                record = record.union(inter_link_sets)
                return inter_link_sets, record
            else:
                raise Exception("Check mode. Not accept the {} mode as input".format(mode))
        elif len(inter_link_sets) == 1:  # Most Cases
            out_link = inter_link_sets.pop()
            record.add(out_link)
            return out_link, record
        else:
            raise Exception('No link between n1 and n2 can be found')

    def cycle_simplify(self, min_cycle, node_graph, link_graph, node_manage):
        '''
        Aggregate the minimal cycles into a distinct segment
        '''
        sub_n = deepcopy(node_graph)
        sub_l = deepcopy(link_graph)
        #   Find the turning pts
        manage = dict()
        endnd = list()  # record searched elements
        node_angle = dict()  # minimal_circle_node_id : angle
        prev_id = min_cycle[-2]
        for idx in range(len(min_cycle) - 1):
            current_id = min_cycle[idx]
            aft_id = min_cycle[idx + 1]
            l1 = sub_n[current_id].intersection(sub_n[prev_id])
            l2 = sub_n[current_id].intersection(sub_n[aft_id])
            angle = AngleCal_turning(self.link_storage[l1.pop()],
                                     self.link_storage[l2.pop()],
                                     self.node_storage[current_id])
            node_angle[current_id] = angle
            if angle > 0.6:  # turning pt
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

            #   update the connectivity information of links and nodes
            for geom_id, target_nid_set in upd_nd_map.items():
                new_node_geom = upd_nd_stg[geom_id]
                #   Take any point as the target_nid, and migrate all other parts of the nodeid information to the
                #   current point
                target_nid = target_nid_set.pop()
                if len(target_nid_set):  # Other node should be aggregated
                    sub_n, sub_l = self.update_replacing_node(target_nid, target_nid_set, [sub_n, sub_l])
                    for i in target_nid_set:
                        if min_cycle[0] == i:
                            min_cycle.pop(0)
                            min_cycle.append(min_cycle[0])
                        if i in min_cycle:
                            min_cycle.remove(i)
                    for key, val in link_status.items():
                        link_status[key] = val.difference(target_nid_set)
                #   Update node_geometry and its connected link endpoint positions
                self.update_change_node_geom(target_nid, new_node_geom)
            #   Update endpoint; Check if the endpoint needs to be updated
            for end, replaced_set in end_update.items():
                sub_n, sub_l = self.update_replacing_node(end, replaced_set, [sub_n, sub_l])
                for i in replaced_set:
                    if min_cycle[0] == i:
                        min_cycle.pop(0)
                        min_cycle.append(min_cycle[0])
                    if i in min_cycle:
                        min_cycle.remove(i)
            #   update link
            current_step = 0
            next_step = 1
            while next_step < len(min_cycle):
                link_list = self.node[min_cycle[current_step]].intersection(self.node[min_cycle[next_step]])
                for li_idx in link_list:
                    #   delete the connections
                    _, sub_n, sub_l = self.del_links(li_idx, [sub_n, sub_l])
                current_step = current_step + 1
                next_step = next_step + 1
            #   Update the topology and geometric information of newly found lines
            _, sub_n, sub_l = self.add_links(new_link, link_status, [sub_n, sub_l])
            return sub_n, sub_l
        else:  # start and end forming a link
            o_nd, d_nd = min_cycle[0], min_cycle[1]
            _, _, sub_n, sub_l = self.merging_duplicate_links(o_nd, d_nd, [sub_n, sub_l])
            return sub_n, sub_l

    def merging_duplicate_links(self, nd1, nd2, sub_graph=None):
        '''
        This function is to merge self.links when multiple self.links exist between nd1 and nd2
        And update the messages in graphs and subgraphs.
        :param nd1:
        :param nd2:
        :param sub_graph: default as None. If None, subgraph should not be updated.
                            subgraph[0]: node -> link
                            subgraph[1]: link -> node
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
        update connectivity configuration
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
                    tem_links.append({'status': 'mid_pt', 'geometry': tem_new_nd})

            else:
                #   When no skipped nds within the
                tem_links.append({'status': 'mid_pt', 'geometry': new_nd_geom})

            s_step = s_idx
            l_step = l_step + 1
        #   put endpt into the link and d must be self.node
        tem_links.append({'status': l_status[len(l_storage) - 1], 'geometry': l_storage[-1]})
        #   tem_link geometry
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
            #   identify whether the junctions are too close to each other, thus need aggregation.
            for id, geom in enumerate(new_node_storage):
                if new_geom.distance(geom) < 20:
                    #   should be aggregated
                    junction_id = id
            if junction_id != len(new_node_storage):
                #   junctions should be aggregated
                new_node_map[junction_id].add(nid)
            else:
                new_node_map[junction_id] = {nid}
                new_node_storage.append(new_geom)

        '''
            insert the junctions according to certain sequence
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
                        This part is to avoid the overlapping of geometry in nd_list and junctions
                    '''
                    if current_geom.distance(junc_geom) == 0:
                        if tem_links[current]['status'] == 'mid_pt':
                            #   currentpts are replaced by juncnd
                            used_status = 1
                            current_links.append({'status': new_node_map[junc_id].copy(),
                                                  'geometry': junc_geom})
                            current = current + 1
                            next = next + 1
                            continue
                        else:
                            '''
                                currentpt and juncnd are same, which means that junc should be replcaed by currentpt,
                                avoiding duplication
                            '''
                            used_status = 1
                            current_links.append({'status': new_node_map[junc_id].copy(),
                                                  'geometry': junc_geom})
                            new_node_map[junc_id].add(tem_links[current]['status'])
                            current = current + 1
                            next = next + 1
                            continue
                    if next_geom.distance(junc_geom) == 0:
                        # nextpt overlaps with juncpt
                        current_links.append(tem_links[next])
                        current = current + 1
                        next = next + 1
                        continue

                    angle = AngleCal_nd(current_geom, junc_geom, next_geom)
                    if angle > 0 or angle is None:
                        #   junc locates between current and next
                        used_status = 1
                        current_links.append({'status': new_node_map[junc_id].copy(),
                                              'geometry': junc_geom})
                current_links.append(tem_links[next])
                current = current + 1
                next = next + 1
            if used_status == 0:
                #   Explanation: No suitable insertion point was found before
                #   Usually occurs at the corner to detect which endpoint it is closer to
                #   Directly replace the original information with endpoint information
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
            Segment the split_link for output
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
        """
        divide the minimal cycle into two separate traces
        the O is the endnd[0] as start pt and the D is the end[1] as end pt of each trace.
        from O to D is to find the intermediate pts of one trace
        from D to O is the intermediate pts of another trace
        ensuring the direction align the principle of start to end, the 'from D to O' pt sequnce should be revert
        :param minimal_seq_list:  the sequence of pt to form a cycle
        :param endnd: the geometric turning of the cycle
        :return:
        """
        o = None
        d = None
        seq_o2d = LinkedList()
        seq_d2o = LinkedList()
        length = None
        minimal_seq_list = minimal_seq_list.copy()
        minimal_seq_list.pop()
        for nd in minimal_seq_list:
            if nd in endnd:
                if o is None:
                    o = nd
                    if nd != endnd[0]:
                        endnd[1] = endnd[0]
                        endnd[0] = o
                    continue
                else:
                    d = nd
                    length = seq_d2o._length
                    continue
            if o is None and d is None:
                seq_d2o.insert(nd, 0)
            elif o is not None and d is None:
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
        if len(pt_list) < 3:
            return pt_list
        pt_list = pt_list.copy()
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
        res_list.append(pt_list[-1])
        return res_list
    #%%
    def Line_forward_Line(self, lid, endli=None):
        if endli is None:
            endli = set()
        if lid.__class__ == int:
            li_list = list(self.link[lid])
            endli.add(lid)
        elif lid.__class__ == list:
            li_list = set()
            for i in lid:
                li_list = li_list.union(self.link[i])
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
        """
        This method returns the link between two nodes, but there are two cases:,
        When there are multiple links:
            1.2 mode 2: "shortest mode": directly return the smallest link value
            1.3 mode 3: "All mode": return all connected segments
        The conditions to be excluded are stored in the record 2
        :param nd1:
        :param nd2:
        :param record: link id to be excluded
        :param mode: return the shortest connected line segment by default
                1.1 mode 1: "updating": directly update and merge the segments to be merged into the same segment
                1.2 mode 2: "shortest": directly return the smallest link value
                1.3 mode 3: "all": return all connected segments (return list)
        """
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
                else:
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
        """
        return another self.node id of current link
        :param lid:
        :param nid:
        :param link_garph:
        :return: mode 'multi_pt': return the list of self.node ids and the string 'multi_pt'
                 mode 'pt': find only one nid and return it and the string 'pt'
                 mode 'end': cannot find any nid and return None and the string 'end'
        """
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
        #   update the nid geometry first
        #   then update the links connected to this node
        for li in self.node[nid]:
            pt_list = Points_in_Line(self.link_storage[li])
            #   Detect whether the newly added point is the head or the tail
            if pt_list[0].distance(new_geom) < pt_list[-1].distance(new_geom):
                #   at head
                pt_list[0] = new_geom
            else:
                #   at tail
                pt_list[-1] = new_geom
            #   update the geometry of li
            if len(pt_list) > 3:
                pt_list = straighten_links(pt_list)
            self.link_storage[li] = LineString(pt_list)
        self.node_storage[nid] = new_geom

    def update_replacing_node(self, target_nid, replaced_set, sub_graph=None):
        if replaced_set.__class__ is list:
            replaced_set = set(replaced_set)

        for replaced_id in replaced_set:
            if replaced_id not in self.node:
                #   was updated before
                continue
            link_set = self.node[replaced_id]
            for li_idx in link_set:
                self.link[li_idx].remove(replaced_id)
                self.link[li_idx].add(target_nid)
                self.node[target_nid].add(li_idx)
            self.node.pop(replaced_id)
            self.node_storage[replaced_id] = None
        if sub_graph is not None:
            node_graph = sub_graph[0]
            link_graph = sub_graph[1]
            for replaced_id in replaced_set:
                if replaced_id not in node_graph:
                    #   was updated before
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
            #   update connectivity of link
            self.link[real_idx] = mapping_relations[tem_idx]
            node_set = mapping_relations[tem_idx]
            #   update connectivity of node
            for i in node_set:
                if real_idx not in self.node[i]:
                    self.node[i].add(real_idx)
            if sub_graph is not None:
                #   update connectivity of subgraph
                sub_node = sub_graph[0]
                sub_link = sub_graph[1]
                #   update connectivity of link in subgraph
                sub_link[real_idx] = {i for i in node_set if i in sub_node}
                #   update connectivity of node in subgraph
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
        :param sub_graph:   subgraph[0] means nodes to links
                            subgraph[1] means links to nodes
        :return:
        '''
        if l_id_list.__class__ is int:
            l_id_list = [l_id_list]
        if l_id_list.__class__ is not list:
            raise TypeError("plz check the the type of 'l_id_list'")
        del_relationships = dict()  # 保留被删除的links相关关系
        for li_id in l_id_list:
            node_set = self.link[li_id]
            del_relationships[li_id] = node_set
            if li_id in self.link:
                self.link.pop(li_id)
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
            #   return the subgraphs
            return del_relationships, sub_node, sub_link

#%%
    def check_graph_validation(self):
        status_n = self.check_node_validation()
        status_l = self.check_link_validation()
        if status_l and status_n:
            # both link and node are all valid
            return
        else:
            self.check_graph_validation()

    def check_link_validation(self, status=1):
        for l_id in range(len(self.link_storage)):
            if self.link_storage[l_id] is None:
                continue
            if len(self.link[l_id]) == 0:
                self.link_storage[l_id] = None
                self.link.pop(l_id)
                status = 0  # there has invalid connectivity in this check
                continue
            if self.link_storage[l_id].length == 0:
                status = 0  # means there has invalid links in this check
                self.link_storage[l_id] = None
                keep_nid = self.link[l_id].pop()
                for n_id in self.link[l_id]:
                    self.node[keep_nid] = self.node[keep_nid].union(self.node[n_id])
                    for l in self.node[n_id]:
                        self.link[l].discard(n_id)
                        self.link[l].add(keep_nid)
                    if n_id in self.node:
                        self.node.pop(n_id)
                        self.node_storage[n_id] = None
                self.node[keep_nid].discard(l_id)
                if l_id in self.link:
                    self.link.pop(l_id)
        return status

    def check_node_validation(self, status=1):
        for n_idx, stg in enumerate(self.node_storage):
            if stg is None:
                continue
            if len(self.node[n_idx]) < 2:
                for i in self.node[n_idx]:
                    self.link[i].discard(n_idx)
                self.node_storage[n_idx] = None
                self.node.pop(n_idx)
                status = 0  # means nodes are all valid
        return status

