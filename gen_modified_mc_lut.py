#!/usr/bin/python
# -*- coding: UTF-8 -*-
import  math
from itertools import combinations, permutations

# references: 
# 1.https://www.cs.upc.edu/~virtual/SGI/docs/1.%20Theory/Unit%2010.%20Volume%20models.%20Marching%20Cubes/Marching%20Cubes.pdf
# 2.A modified look-up table for implicit disambiguation of marching cubes
'''
    vertex number:
         4--------------5
        /|             /|
       / |            / |
      7--|-----------6  |
      |  |           |  |
      |  0-----------|--1
      | /            | /
      |/             |/ 
      3--------------2
'''
vertices = [
    (-1.0, -1.0, -1.0),
    (1.0, -1.0, -1.0),
    (1.0, -1.0, 1.0),
    (-1.0, -1.0, 1.0),
    (-1.0, 1.0, -1.0),
    (1.0, 1.0, -1.0),
    (1.0, 1.0, 1.0),
    (-1.0, 1.0, 1.0)
]

'''
    edge number:
         -------4-------
        7|             5|
       / 8            / |
      ---|-----6-----/  9
     11  |          10  |
      |  -------0----|--|
      | 3            | 1
      |/             |/ 
      ---------2-----/
'''
vertex2edge = [
    (0, 3, 8),
    (0, 1, 9),
    (1, 2, 10),
    (2, 3, 11),
    (4, 7, 8),
    (4, 5, 9),
    (5, 6, 10),
    (6, 7, 11)
]

edge2vertex = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 0),
    (4, 5),
    (5, 6),
    (6, 7),
    (7, 4),
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7)
]


Face_Left   = 0
Face_Right  = 1
Face_Top    = 2
Face_Bottom = 3
Face_Front  = 4
Face_Back   = 5

facets = [
    (1, 2, 6, 5), #Face_Left
    (0, 3, 7, 4), #Face_Right
    (5, 6, 7, 4), #Face_Top
    (1, 2, 3, 0), #Face_Bottom
    (6, 2, 3, 7), #Face_Front
    (5, 1, 0, 4) #Face_Back
]

oppo_facets = {
    Face_Left : Face_Right,
    Face_Right : Face_Left,
    Face_Top : Face_Bottom,
    Face_Bottom : Face_Top,
    Face_Front : Face_Back,
    Face_Back : Face_Front
}

def has_bit(value, bit):
    return (value) & (1 << bit)

def indices_from_bit(value):
    indices = [i for i in range(8) if has_bit(value, i)]
    return indices

def indices_from_exclude(indices):
    exindices = [i for i in range(8) if i not in indices]
    return exindices

def indices_from_exclude_indices(indices, exindices):
    remains = []
    for i in range(len(indices)):
        if indices[i] not in exindices:
            remains.append(indices[i])
    return remains

def edge_gen_bits(value):
    ebits = list([0 for k in range(12)])
    for j in range(8):
        if (has_bit(value, j)):
            ebits[vertex2edge[j][0]] ^= 1
            ebits[vertex2edge[j][1]] ^= 1
            ebits[vertex2edge[j][2]] ^= 1
    return ebits

def edge_gen_bits_from_vertex_indices(indices):
    ebits = [0 for k in range(12)]
    for i in range(len(indices)):
        vertex_index = indices[i]
        edges = vertex2edge[vertex_index]
        ebits[edges[0]] ^= 1
        ebits[edges[1]] ^= 1
        ebits[edges[2]] ^= 1
    return ebits

def edge_cut_count(ebits):
    count = 0
    for i in range(12):
        if ebits[i]:
            count += 1
    return count

def edge_gen_number(ebits):
    value = 0
    for i in range(12):
        if ebits[i]:
            value |= (1 << i)
    return value

def edge_gen_edge_indices(ebits):
    edge_indices = [i for i in range(12) if ebits[i]]
    return edge_indices

def edge_gen_edge_indices_from_vertex_indices(indices):
    ebits = edge_gen_bits_from_vertex_indices(indices)
    return edge_gen_edge_indices(ebits)

def edge_query_edge_indices_from_vertex_indice(eindices, vindex):
    eindices_ = []
    for i in range(len(eindices)):
        eindex = eindices[i]
        if vindex in edge2vertex[eindex]:
            eindices_.append(eindex)
    return eindices_

def edge_is_one_cut(ebits, vindex):
    edges = vertex2edge[vindex]
    if ebits[edges[0]] + ebits[edges[1]] + ebits[edges[2]] == 1:
        if ebits[edges[0]]:
            return (True, edges[0])
        elif ebits[edges[1]]:
            return (True, edges[1])
        else:
            return (True, edges[2])
        return (False, None)
    return (False, None)

def edge_is_two_cut(ebits, vindex):
    edges = vertex2edge[vindex]
    if ebits[edges[0]] + ebits[edges[1]] + ebits[edges[2]] == 2:
        eindices = []
        for i in range(3):
            if ebits[edges[i]]:
                eindices.append(edges[i])
        return (True, eindices)
    return (False, None)

def edge_is_full_cut(ebits, vindex):
    edges = vertex2edge[vindex]
    if ebits[edges[0]] + ebits[edges[1]] + ebits[edges[2]] == 3:
        return (True, list(edges))
    return (False, None)

def edge_is_zero_cut(ebits, vindex):
    edges = vertex2edge[vindex]
    if ebits[edges[0]] + ebits[edges[1]] + ebits[edges[2]] == 0:
        return (True, [])
    return (False, None)

def subtract(P1, P2): # P1 - P2
    return (P1[0] - P2[0], P1[1] - P2[1], P1[2] - P2[2])

def length(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

def normalized(v):
    len = length(v)
    return (v[0] / len, v[1] / len, v[2] / len)

def dot(P, Q):
    return P[0]*Q[0] + P[1]*Q[1] + P[2]*Q[2]

def cross(P, Q):
    return (P[1]*Q[2] - P[2]*Q[1], P[2]*Q[0] - P[0]*Q[2], P[0]*Q[1] - Q[0]*P[1])

def gen_vector_from_edge_index(eindex):
    vindices =edge2vertex[eindex]
    v1 = vertices[vindices[0]]
    v2 = vertices[vindices[1]]
    return normalized(subtract(v1, v2))

def gen_midpoint_from_edge_index(eindex):
    vindices = edge2vertex[eindex]
    v1 = vertices[vindices[0]]
    v2 = vertices[vindices[1]]
    return ((v1[0] + v2[0]) / 2.0, (v1[1] + v2[1]) / 2.0, (v1[2] + v2[2]) / 2.0)

def is_same_face_from_vertex_index_2(vindex1, vindex2):
    for i in range(len(facets)):
        if (vindex1 in facets[i]) and (vindex2 in facets[i]):
            return True
    return False

def is_same_face_from_vertex_index_4(vindex1, vindex2, vindex3, vindex4):
    for i in range(len(facets)):
        if (vindex1 in facets[i]) and (vindex2 in facets[i]) and\
            (vindex3 in facets[i]) and (vindex4 in facets[i]):
            return True
    return False

'''
    确定三角面的顶点索引(假设三个顶点为P0, P1, P2, 环绕顺序)
    根据平面方程, 判断关联的顶点是否在三角面确定的平面(等值面)下, 即dot(N, P) + D < 0(其中, D = dot(N, P0))
'''
def is_below_triangle(P0, P1, P2, vertex):
    P0P1 = (P1[0] - P0[0], P1[1] - P0[1], P1[2] - P0[2])
    P0P2 = (P2[0] - P0[0], P2[1] - P0[1], P2[2] - P0[2])

    N1 = cross(P0P1, P0P2)
    D1 = -dot(N1, P0)
    if dot(N1, vertex) + D1 <= 0:
        return True

    return False

def is_regular_triangle(eindices, vertex):
    eindex0 = eindices[0]
    eindex1 = eindices[1]
    eindex2 = eindices[2]
    P0 = gen_midpoint_from_edge_index(eindex0)
    P1 = gen_midpoint_from_edge_index(eindex1)
    P2 = gen_midpoint_from_edge_index(eindex2)
    return is_below_triangle(P0, P1, P2, vertex)

def gen_triangle_indices__(eindex0, eindex1, eindex2, indices):
    count = 0
    vlen = len(indices)
    eindices = [eindex0, eindex1, eindex2]
    for i in range(vlen):
        if is_regular_triangle(eindices, vertices[indices[i]]):
            count += 1
    if count == vlen:
        return eindices

    count = 0
    eindices = [eindex0, eindex2, eindex1]
    for i in range(vlen):
        if is_regular_triangle(eindices, vertices[indices[i]]):
            count += 1
    if count == vlen:
        return eindices

    print("gen_triangle_indices!fatal error, (%d, %d, %d, #%d)"%(eindex0, eindex1, eindex2, vlen))
    return []

def gen_triangle_indices(eindices, indices):
    used_vindices = []
    for i in range(len(eindices)):
        eindex = eindices[i]
        edges = edge2vertex[eindex]
        contains0 = edges[0] in indices
        contains1 = edges[1] in indices
        if contains0 or contains1:
            if contains0 and edges[0] not in used_vindices:
                used_vindices.append(edges[0])
            elif contains1 and edges[1] not in used_vindices:
                used_vindices.append(edges[1])
    return gen_triangle_indices__(eindices[0], eindices[1], eindices[2], used_vindices)

def query_closest_edge(eindex, eindices):
    min_dist = None
    min_eindex = None
    P0 = gen_midpoint_from_edge_index(eindex)
    for i in range(len(eindices)):
        eindex_ = eindices[i]
        if eindex_ != eindex:
            P1 = gen_midpoint_from_edge_index(eindex_)
            dist = length(subtract(P0, P1))
            if min_dist == None or min_dist > dist:
                min_dist = dist
                min_eindex = eindex_
    if min_eindex != None:
        eindices.remove(min_eindex)
    else:
        print("query_closest_edge!fatal error, {eindex}: {eindices}".format(eindex=eindex, eindices=eindices))
    return min_eindex, eindices

def query_closest_vertex(vindex, vindices):
    min_dist = None
    min_vindex = None
    P0 = vertices[vindex]
    for i in range(len(vindices)):
        vindex_ = vindices[i]
        if vindex_ != vindex:
            P1 = vertices[vindex_]
            dist = length(subtract(P0, P1))
            if min_dist == None or min_dist > dist:
                min_dist = dist
                min_vindex = vindex_
    if min_vindex == None:
        print("query_closest_vertex!fatal error, {vindex}: {vindices}".format(vindex=vindex, vindices=vindices))
    return min_vindex

def query_farest_vertex(vindex, vindices):
    max_dist = None
    max_vindex = None
    P0 = vertices[vindex]
    for i in range(len(vindices)):
        vindex_ = vindices[i]
        if vindex_ != vindex:
            P1 = vertices[vindex_]
            dist = length(subtract(P0, P1))
            if max_dist == None or max_dist < dist:
                max_dist = dist
                max_vindex = vindex_
    if max_vindex == None:
        print("query_farest_vertex!fatal error, {vindex}: {vindices}".format(vindex=vindex, vindices=vindices))
    return max_vindex

def query_farest_edge_from_edge(eindex, eindices):
    max_dist = None
    max_eindex = None
    P0 = gen_midpoint_from_edge_index(eindex)
    for i in range(len(eindices)):
        eindex_ = eindices[i]
        if eindex != eindex_:
            P1 = gen_midpoint_from_edge_index(eindex_)
            dist = length(subtract(P0, P1))
            if max_dist == None or max_dist < dist:
                max_dist = dist
                max_eindex = eindex_
    if max_eindex != None:
        eindices.remove(max_eindex)
    else:
        print("query_farest_edge_from_edge!fatal error, {eindex}: {eindices}".format(eindex=eindex, eindices=eindices))
    return max_eindex, eindices

def query_farest_edge_from_vertex(vindex, eindices):
    max_dist = None
    max_eindex = None
    P0 = vertices[vindex]
    for i in range(len(eindices)):
        eindex = eindices[i]
        P1 = gen_midpoint_from_edge_index(eindex)
        dist = length(subtract(P0, P1))
        if max_dist == None or max_dist < dist:
            max_dist = dist
            max_eindex = eindex
    if max_eindex != None:
        eindices.remove(max_eindex)
    else:
        print("query_farest_edge_from_vertex!fatal error, {vindex}: {eindices}".format(vindex=vindex, eindices=eindices))
    return max_eindex, eindices

def query_cloest_edge_from_vertex(vindex, eindices):
    min_dist = None
    min_eindex = None
    P0 = vertices[vindex]
    for i in range(len(eindices)):
        eindex = eindices[i]
        P1 = gen_midpoint_from_edge_index(eindex)
        dist = length(subtract(P0, P1))
        if min_dist == None or min_dist > dist:
            min_dist = dist
            min_eindex = eindex
    if min_eindex != None:
        eindices.remove(min_eindex)
    else:
        print("query_cloest_edge_from_vertex!fatal error, {eindex}: {eindices}".format(eindex=min_eindex, eindices=eindices))
    return min_eindex, eindices

def triangle_indices_fill(triangles):
    count = len(triangles)
    return triangles + [-1 for i in range(16 - count)]

def gen_modified_mc_lut_case0(indices):
    return [-1 for i in range(16)]

def gen_modified_mc_lut_case1(indices):
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    triangles = []
    triangles.extend(gen_triangle_indices(eindices, indices))
    if len(triangles) / 3 != 1:
        print("gen_modified_mc_lut_case1!fatal error, triangle count=%d"%(len(triangles) / 3))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case2A(indices):
    vindex0 = indices[0]
    vindex1 = indices[1]
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    other_eindices = edge_query_edge_indices_from_vertex_indice(eindices, vindex0)
    other_eindices += edge_query_edge_indices_from_vertex_indice(eindices, vindex1)
    
    triangles = []
    min_eindex = other_eindices[0]
    other_eindices.pop(0)
    used_indices = [min_eindex]
    while len(other_eindices) > 0:
        min_eindex, other_eindices = query_closest_edge(min_eindex, other_eindices)
        used_indices.append(min_eindex)
        if len(used_indices) == 3:
            triangles.extend(gen_triangle_indices(used_indices, indices))
            used_indices.pop(1)

    if len(triangles) / 3 != 2:
        print("gen_modified_mc_lut_case2A!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(other_eindices):
        print("gen_modified_mc_lut_case2A!fatal error, other_eindices={other_eindices}".format(other_eindices=other_eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case2B(indices):
    vindex0 = indices[0]
    vindex1 = indices[1]
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    eindices0 = edge_query_edge_indices_from_vertex_indice(eindices, vindex0)
    eindices1 = edge_query_edge_indices_from_vertex_indice(eindices, vindex1)

    # 找到离vindex1距离最远的边
    max_eindex, eindices0 = query_farest_edge_from_vertex(vindex1, eindices0)

    triangles = []
    used_eindices = [max_eindex]
    min_eindex, eindices0 = query_closest_edge(max_eindex, eindices0)
    used_eindices.append(min_eindex)

    min_eindex, eindices1 = query_closest_edge(min_eindex, eindices1)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    used_eindices.pop(1)

    max_eindex1, eindices1 = query_farest_edge_from_vertex(vindex0, eindices1)
    used_eindices.append(max_eindex1)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    used_eindices.pop(1)

    used_eindices = [max_eindex, max_eindex1]
    min_eindex, eindices0 = query_closest_edge(max_eindex, eindices0)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    used_eindices.pop(0)

    min_eindex, eindices1 = query_closest_edge(min_eindex, eindices1)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))

    if len(triangles) / 3 != 4:
        print("gen_modified_mc_lut_case2B!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices0):
        print("gen_modified_mc_lut_case2B!fatal error, eindices0={eindices0}".format(eindices0=eindices0))
    if len(eindices1):
        print("gen_modified_mc_lut_case2B!fatal error, eindices1={eindices1}".format(eindices1=eindices1))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case2C(indices):
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    vindex0 = indices[0]
    vindex1 = indices[1]
    eindices0 = edge_query_edge_indices_from_vertex_indice(eindices, vindex0)
    eindices1 = edge_query_edge_indices_from_vertex_indice(eindices, vindex1)
    triangles = []
    triangles.extend(gen_triangle_indices(eindices0, indices))
    triangles.extend(gen_triangle_indices(eindices1, indices))
    if len(triangles) / 3 != 2:
        print("gen_modified_mc_lut_case2C!fatal error, triangle count=%d"%(len(triangles) / 3))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case3A(indices):
    # 找到one cut的顶点
    ebits = edge_gen_bits_from_vertex_indices(indices)
    one_vindex = None
    one_eindex = None
    for i in range(len(indices)):
        vindex = indices[i]
        isok, eindex = edge_is_one_cut(ebits, vindex)
        if isok:
            one_vindex = vindex
            one_eindex = eindex
            break

    vindex1 = None
    vindex2 = None
    for i in range(len(indices)):
        vindex = indices[i]
        if vindex != one_vindex:
            if None == vindex1:
                vindex1 = vindex
            else:
                vindex2 = vindex

    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    other_eindices = edge_query_edge_indices_from_vertex_indice(eindices, vindex1)
    other_eindices += edge_query_edge_indices_from_vertex_indice(eindices, vindex2)

    triangles = []
    min_eindex = one_eindex
    used_indices = [one_eindex]
    while len(other_eindices) > 0:
        min_eindex, other_eindices = query_closest_edge(min_eindex, other_eindices)
        used_indices.append(min_eindex)
        if len(used_indices) == 3:
            triangles.extend(gen_triangle_indices(used_indices, indices))
            used_indices.pop(1)
    
    if len(triangles) / 3 != 3:
        print("gen_modified_mc_lut_case3A!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(other_eindices):
        print("gen_modified_mc_lut_case3A!fatal error, other_eindices={other_eindices}".format(other_eindices=other_eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case3B(indices):
    # 找到full cut的顶点
    ebits = edge_gen_bits_from_vertex_indices(indices)
    vindex0 = None
    eindices0 = None
    for i in range(len(indices)):
        vindex = indices[i]
        isok, eindices = edge_is_full_cut(ebits, vindex)
        if isok:
            vindex0 = vindex
            eindices0 = eindices
            break

    # 找到离vindex0最近的顶点
    min_vindex = query_closest_vertex(vindex0, indices)

    # 找到最后一个剩余的顶点
    last_vindex = None
    for i in range(len(indices)):
        vindex = indices[i]
        if vindex != vindex0 and vindex != min_vindex:
            last_vindex = vindex
            break

    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    min_eindices = edge_query_edge_indices_from_vertex_indice(eindices, min_vindex)
    last_eindices = edge_query_edge_indices_from_vertex_indice(eindices, last_vindex)

    # 找到离min_vindex最远的边
    keep_min_eindex, eindices0 = query_farest_edge_from_vertex(min_vindex, eindices0)

    triangles = []
    used_eindices = [keep_min_eindex]

    min_eindex, eindices0 = query_closest_edge(keep_min_eindex, eindices0)
    used_eindices.append(min_eindex)
    last_min_eindex, last_eindices = query_closest_edge(min_eindex, last_eindices)
    used_eindices.append(last_min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    used_eindices.pop(0)

    min_eindex, min_eindices = query_closest_edge(last_min_eindex, min_eindices)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))

    used_eindices = [keep_min_eindex, last_min_eindex]
    min_eindex, last_eindices = query_closest_edge(last_min_eindex, last_eindices)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    used_eindices.remove(last_min_eindex)

    min_eindex, eindices0 = query_closest_edge(min_eindex, eindices0)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    used_eindices.remove(keep_min_eindex)

    min_eindex, min_eindices = query_closest_edge(min_eindex, min_eindices)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    
    if len(triangles) / 3 != 5:
        print("gen_modified_mc_lut_case3B!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices0):
        print("gen_modified_mc_lut_case3B!fatal error, eindices0={eindices0}".format(eindices0=eindices0))
    if len(min_eindices):
        print("gen_modified_mc_lut_case3B!fatal error, min_eindices={min_eindices}".format(min_eindices=min_eindices))
    if len(last_eindices):
        print("gen_modified_mc_lut_case3B!fatal error, last_eindices={last_eindices}".format(last_eindices=last_eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case3C(indices):
    complement_indices = indices_from_exclude(indices)
    ebits = edge_gen_bits_from_vertex_indices(indices)
    # 找到zero cut和full cut的顶点
    zero_vindex = None
    full_vindex = None
    full_eindices = None
    for i in range(len(complement_indices)):
        vindex = complement_indices[i]
        if None == full_vindex:
            isok, full_eindices = edge_is_full_cut(ebits, vindex)
            if isok:
                full_vindex = vindex
        if None == zero_vindex:
            isok, _ = edge_is_zero_cut(ebits, vindex)
            if isok:
                zero_vindex = vindex
        if full_vindex and zero_vindex:
            break
    
    # 从complement_indices找到离full cut最远的顶点
    farest_vindex = query_farest_vertex(full_vindex, complement_indices)
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)

    # 从eindices剔除full_eindices的边
    for i in range(len(full_eindices)):
        eindex = full_eindices[i]
        eindices.remove(eindex)

    # 随机从eindices里面找到离farest_vindex最近的边
    min_eindex, eindices = query_cloest_edge_from_vertex(farest_vindex, eindices)
    
    triangles = []
    triangles.extend(gen_triangle_indices(full_eindices, indices))

    used_eindices = [min_eindex]
    while len(eindices) > 0:
        min_eindex, eindices = query_closest_edge(min_eindex, eindices)
        used_eindices.append(min_eindex)
        if len(used_eindices) == 3:
            triangles.extend(gen_triangle_indices(used_eindices, indices))
            used_eindices.pop(1)

    if len(triangles) / 3 != 5:
        print("gen_modified_mc_lut_case3C!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices):
        print("gen_modified_mc_lut_case3C!fatal error, eindices={eindices}".format(eindices=eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case4A(indices):
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    min_eindex = eindices[0]
    eindices.pop(0)
    used_eindices = [min_eindex]
    triangles = []
    while len(eindices) > 0:
        min_eindex, eindices = query_closest_edge(min_eindex, eindices)
        used_eindices.append(min_eindex)
        if len(used_eindices) == 3:
            triangles.extend(gen_triangle_indices(used_eindices, indices))
            used_eindices.pop(1)

    if len(triangles) / 3 != 2:
        print("gen_modified_mc_lut_case4A!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices):
        print("gen_modified_mc_lut_case4A!fatal error, eindices={eindices}".format(eindices=eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case4B(indices):
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    min_eindex = eindices[0]
    eindices.pop(0)
    used_eindices = [min_eindex]
    triangles = []
    while len(eindices) > 0:
        min_eindex, eindices = query_closest_edge(min_eindex, eindices)
        used_eindices.append(min_eindex)
        if len(used_eindices) == 3:
            triangles.extend(gen_triangle_indices(used_eindices, indices))
            used_eindices.pop(1)

    if len(triangles) / 3 != 4:
        print("gen_modified_mc_lut_case4B!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices):
        print("gen_modified_mc_lut_case4B!fatal error, eindices={eindices}".format(eindices=eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case4C(indices):
    complement_indices = indices_from_exclude(indices)
    # 找到complement_indices里互相最近的两组顶点
    nearest_indices0 = []
    nearest_indices1 = []
    nearest_indices0.append(complement_indices[0])
    complement_indices.pop(0)
    min_dist = None
    min_vindex = None
    P0 = vertices[nearest_indices0[0]]
    for i in range(len(complement_indices)):
        vindex = complement_indices[i]
        P1 = vertices[vindex]
        dist = length(subtract(P0, P1))
        if min_dist == None or dist < min_dist:
            min_dist = dist
            min_vindex = vindex
    if min_vindex:
        nearest_indices0.append(min_vindex)
        complement_indices.remove(min_vindex)
        nearest_indices1.append(complement_indices[0])
        nearest_indices1.append(complement_indices[1])

    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    other_eindices0 = edge_query_edge_indices_from_vertex_indice(eindices, nearest_indices0[0])
    other_eindices0 += edge_query_edge_indices_from_vertex_indice(eindices, nearest_indices0[1])
    other_eindices1 = edge_query_edge_indices_from_vertex_indice(eindices, nearest_indices1[0])
    other_eindices1 += edge_query_edge_indices_from_vertex_indice(eindices, nearest_indices1[1])

    triangles = []
    min_eindex = other_eindices0[0]
    other_eindices0.pop(0)
    used_indices = [min_eindex]
    while len(other_eindices0) > 0:
        min_eindex, other_eindices0 = query_closest_edge(min_eindex, other_eindices0)
        used_indices.append(min_eindex)
        if len(used_indices) == 3:
            triangles.extend(gen_triangle_indices(used_indices, indices))
            used_indices.pop(1)

    min_eindex = other_eindices1[0]
    other_eindices1.pop(0)
    used_indices = [min_eindex]
    while len(other_eindices1) > 0:
        min_eindex, other_eindices1 = query_closest_edge(min_eindex, other_eindices1)
        used_indices.append(min_eindex)
        if len(used_indices) == 3:
            triangles.extend(gen_triangle_indices(used_indices, indices))
            used_indices.pop(1)

    if len(triangles) / 3 != 4:
        print("gen_modified_mc_lut_case4C!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(other_eindices0):
        print("gen_modified_mc_lut_case4C!fatal error, other_eindices0={other_eindices0}".format(other_eindices0=other_eindices0))
    if len(other_eindices1):
        print("gen_modified_mc_lut_case4C!fatal error, other_eindices1={other_eindices1}".format(other_eindices1=other_eindices1))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case4D(indices):
    # 找到one cut的顶点
    one_vindex = None
    one_eindex = None
    ebits = edge_gen_bits_from_vertex_indices(indices)
    for i in range(len(indices)):
        vindex = indices[i]
        isok, one_eindex = edge_is_one_cut(ebits, vindex)
        if isok:
            one_vindex = vindex
            one_eindex = one_eindex
            break

    # 剔除one_eindex
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    eindices.remove(one_eindex)

    triangles = []
    used_eindices0 = [one_eindex]
    min_eindex, eindices = query_closest_edge(one_eindex, eindices)
    used_eindices0.append(min_eindex)
    min_eindex, eindices = query_closest_edge(min_eindex, eindices)
    used_eindices0.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices0, indices))
    used_eindices0.pop(1)

    used_eindices1 = [one_eindex]
    min_eindex, eindices = query_closest_edge(one_eindex, eindices)
    used_eindices1.append(min_eindex)
    min_eindex, eindices = query_closest_edge(min_eindex, eindices)
    used_eindices1.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices1, indices))
    used_eindices1.pop(1)

    min_eindex, eindices = query_closest_edge(min_eindex, eindices)
    used_eindices1.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices1, indices))
    triangles.extend(gen_triangle_indices(used_eindices0 + [min_eindex], indices))

    if len(triangles) / 3 != 4:
        print("gen_modified_mc_lut_case4D!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices):
        print("gen_modified_mc_lut_case4D!fatal error, eindices={eindices}".format(eindices=eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case4E(indices):
    # 从indices找full cut的顶点
    full_vindex = None
    full_eindices = None
    ebits = edge_gen_bits_from_vertex_indices(indices)
    for i in range(len(indices)):
        vindex = indices[i]
        isok, full_eindices = edge_is_full_cut(ebits, vindex)
        if isok:
            full_vindex = vindex
            break
    
    # 找到离full_vindex最远的顶点
    farest_vindex = query_farest_vertex(full_vindex, indices)

    # 从eindces离剔除full_eindices, farest_eindex
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    for i in range(len(full_eindices)):
        eindices.remove(full_eindices[i])

    farest_eindex = edge_query_edge_indices_from_vertex_indice(eindices, farest_vindex)[0]
    eindices.remove(farest_eindex)

    triangles = []
    triangles.extend(gen_triangle_indices(full_eindices, indices))

    used_eindices = [farest_eindex]
    min_eindex, eindices = query_closest_edge(farest_eindex, eindices)
    used_eindices.append(min_eindex)
    min_eindex, eindices = query_closest_edge(farest_eindex, eindices)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    used_eindices.pop(0)

    min_eindex, eindices = query_closest_edge(min_eindex, eindices)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))
    used_eindices.pop(0)

    min_eindex, eindices = query_closest_edge(min_eindex, eindices)
    used_eindices.append(min_eindex)
    triangles.extend(gen_triangle_indices(used_eindices, indices))

    if len(triangles) / 3 != 4:
        print("gen_modified_mc_lut_case4E!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices):
        print("gen_modified_mc_lut_case4E!fatal error, eindices={eindices}".format(eindices=eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case4F(indices):
    # 找到4组full cut的顶点
    full_vindices = []
    full_eindices = []
    ebits = edge_gen_bits_from_vertex_indices(indices)
    complement_indices = indices_from_exclude(indices)
    for i in range(len(complement_indices)):
        vindex = complement_indices[i]
        isok, eindices = edge_is_full_cut(ebits, vindex)
        if isok:
            full_vindices.append(vindex)
            full_eindices.append(eindices)

    triangles = []
    for i in range(len(full_vindices)):
        triangles.extend(gen_triangle_indices(full_eindices[i], indices))

    if len(triangles) / 3 != 4:
        print("gen_modified_mc_lut_case4F!fatal error, triangle count=%d"%(len(triangles) / 3))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case5A(indices): # 和 case 3A 互补
    complement_indices = indices_from_exclude(indices)
    # 找到one cut的顶点
    ebits = edge_gen_bits_from_vertex_indices(indices)
    one_vindex = None
    one_eindex = None
    for i in range(len(complement_indices)):
        vindex = complement_indices[i]
        isok, eindex = edge_is_one_cut(ebits, vindex)
        if isok:
            one_vindex = vindex
            one_eindex = eindex
            break

    # 剔除one cut顶点的边
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    eindices.remove(one_eindex)

    triangles = []
    min_eindex = one_eindex
    used_indices = [one_eindex]
    while len(eindices) > 0:
        min_eindex, eindices = query_closest_edge(min_eindex, eindices)
        used_indices.append(min_eindex)
        if len(used_indices) == 3:
            triangles.extend(gen_triangle_indices(used_indices, indices))
            used_indices.pop(1)

    if len(triangles) / 3 != 3:
        print("gen_modified_mc_lut_case5A!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices):
        print("gen_modified_mc_lut_case5A!fatal error, eindices={eindices}".format(eindices=eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case5B(indices):
    complement_indices = indices_from_exclude(indices)
    # 从complement_indices找到full cut的顶点
    full_vindex = None
    full_eindices = None
    ebits = edge_gen_bits_from_vertex_indices(indices)
    for i in range(len(complement_indices)):
        vindex = complement_indices[i]
        isok, full_eindices = edge_is_full_cut(ebits, vindex)
        if isok:
            full_vindex = vindex
            break

    # 从eindices里面剔除full_eindices
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    for i in range(len(full_eindices)):
        eindices.remove(full_eindices[i])
    
    triangles = []
    triangles.extend(gen_triangle_indices(full_eindices, indices))

    min_eindex = eindices[0]
    eindices.pop(0)
    used_indices = [min_eindex]
    while len(eindices) > 0:
        min_eindex, eindices = query_closest_edge(min_eindex, eindices)
        used_indices.append(min_eindex)
        if len(used_indices) == 3:
            triangles.extend(gen_triangle_indices(used_indices, indices))
            used_indices.pop(1)

    if len(triangles) / 3 != 3:
        print("gen_modified_mc_lut_case5B!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(eindices):
        print("gen_modified_mc_lut_case5B!fatal error, eindices={eindices}".format(eindices=eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case5C(indices):
    # 从complement_indices里找到3个full cut顶点
    full_vindices = []
    full_eindices = []
    ebits = edge_gen_bits_from_vertex_indices(indices)
    complement_indices = indices_from_exclude(indices)
    for i in range(len(complement_indices)):
        vindex = complement_indices[i]
        isok, eindices = edge_is_full_cut(ebits, vindex)
        if isok:
            full_vindices.append(vindex)
            full_eindices.append(eindices)

    triangles = []
    for i in range(len(full_vindices)):
        triangles.extend(gen_triangle_indices(full_eindices[i], indices))

    if len(triangles) / 3 != 3:
        print("gen_modified_mc_lut_case5C!fatal error, triangle count=%d"%(len(triangles) / 3))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case6A(indices): # case 2A 互补
    complement_indices = indices_from_exclude(indices)
    vindex0 = complement_indices[0]
    vindex1 = complement_indices[1]
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    other_eindices = edge_query_edge_indices_from_vertex_indice(eindices, vindex0)
    other_eindices += edge_query_edge_indices_from_vertex_indice(eindices, vindex1)
    
    triangles = []
    min_eindex = other_eindices[0]
    other_eindices.pop(0)
    used_indices = [min_eindex]
    while len(other_eindices) > 0:
        min_eindex, other_eindices = query_closest_edge(min_eindex, other_eindices)
        used_indices.append(min_eindex)
        if len(used_indices) == 3:
            triangles.extend(gen_triangle_indices(used_indices, indices))
            used_indices.pop(1)

    if len(triangles) / 3 != 2:
        print("gen_modified_mc_lut_case6A!fatal error, triangle count=%d"%(len(triangles) / 3))
    if len(other_eindices):
        print("gen_modified_mc_lut_case6A!fatal error, other_eindices={other_eindices}".format(other_eindices=other_eindices))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case6B(indices):
    # 从complement_indices里找到2个full cut顶点
    full_vindices = []
    full_eindices = []
    ebits = edge_gen_bits_from_vertex_indices(indices)
    complement_indices = indices_from_exclude(indices)
    for i in range(len(complement_indices)):
        vindex = complement_indices[i]
        isok, eindices = edge_is_full_cut(ebits, vindex)
        if isok:
            full_vindices.append(vindex)
            full_eindices.append(eindices)

    triangles = []
    for i in range(len(full_vindices)):
        triangles.extend(gen_triangle_indices(full_eindices[i], indices))

    if len(triangles) / 3 != 2:
        print("gen_modified_mc_lut_case6B!fatal error, triangle count=%d"%(len(triangles) / 3))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case6C(indices): # 和 case 2C互补
    complement_indices = indices_from_exclude(indices)
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    vindex0 = complement_indices[0]
    vindex1 = complement_indices[1]
    eindices0 = edge_query_edge_indices_from_vertex_indice(eindices, vindex0)
    eindices1 = edge_query_edge_indices_from_vertex_indice(eindices, vindex1)
    triangles = []
    triangles.extend(gen_triangle_indices(eindices0, indices))
    triangles.extend(gen_triangle_indices(eindices1, indices))

    if len(triangles) / 3 != 2:
        print("gen_modified_mc_lut_case6C!fatal error, triangle count=%d"%(len(triangles) / 3))
    return triangle_indices_fill(triangles)

def gen_modified_mc_lut_case7(indices):
    eindices = edge_gen_edge_indices_from_vertex_indices(indices)
    triangles = []
    triangles.extend(gen_triangle_indices(eindices, indices))
    if len(triangles) / 3 != 1:
        print("gen_modified_mc_lut_case7!fatal error, triangle count=%d"%(len(triangles) / 3))
    return triangle_indices_fill(triangles)

def generate_edge_tables():
    edge_tables = []
    for i in range(256):
        ebits = edge_gen_bits(i)
        edge_tables.append(edge_gen_number(ebits))
    return edge_tables

def generate_triangle_tables():
    triangle_tables = []
    total_count_case0 = 0
    total_count_case1 = 0
    total_count_case2A = 0
    total_count_case2B = 0
    total_count_case2C = 0
    total_count_case3A = 0
    total_count_case3B = 0
    total_count_case3C = 0
    total_count_case4A = 0
    total_count_case4B = 0
    total_count_case4C = 0
    total_count_case4D = 0
    total_count_case4E = 0
    total_count_case4F = 0
    total_count_case5A = 0
    total_count_case5B = 0
    total_count_case5C = 0
    total_count_case6A = 0
    total_count_case6B = 0
    total_count_case6C = 0
    total_count_case7 = 0
    for i in range(256):
        indices = indices_from_bit(i)
        if i == 5:
            vi = 2;
        total_indices = len(indices)
        ebits = edge_gen_bits_from_vertex_indices(indices)
        cut_count = edge_cut_count(ebits)
        if 1 == total_indices: # case 1
            total_count_case1 += 1
            triangle_tables.append(gen_modified_mc_lut_case1(indices))
        elif 2 == total_indices: # case 2A/2B/2C
            if 4 == cut_count: # case 2A
                total_count_case2A += 1
                triangle_tables.append(gen_modified_mc_lut_case2A(indices))
            elif 6 == cut_count:
                if is_same_face_from_vertex_index_2(indices[0], indices[1]): # case 2B
                    total_count_case2B += 1
                    triangle_tables.append(gen_modified_mc_lut_case2B(indices))
                else: # case 2C
                    total_count_case2C += 1
                    triangle_tables.append(gen_modified_mc_lut_case2C(indices))
            else:
                print("generate_triangle_tables!case 2X fatal error, cut_count={cut_count}".format(cut_count=cut_count))
        elif 3 == total_indices: # case 3A/3B/3C
            if 5 == cut_count: # case 3A
                total_count_case3A += 1
                triangle_tables.append(gen_modified_mc_lut_case3A(indices))
            elif 7 == cut_count: # case 3B
                total_count_case3B += 1
                triangle_tables.append(gen_modified_mc_lut_case3B(indices))
            elif 9 == cut_count: # case 3C
                total_count_case3C += 1
                triangle_tables.append(gen_modified_mc_lut_case3C(indices))
            else:
                print("generate_triangle_tables!case 3X fatal error, cut_count={cut_count}".format(cut_count=cut_count))
        elif 4 == total_indices: # case 4A/4B/4C/4D/4F
            if 4 == cut_count: # case 4A
                total_count_case4A += 1
                triangle_tables.append(gen_modified_mc_lut_case4A(indices))
            elif 6 == cut_count: # case 4B/4D/4E
                two_cut_count = 0
                for i in range(len(indices)):
                    edges = vertex2edge[indices[i]]
                    if ebits[edges[0]] + ebits[edges[1]] + ebits[edges[2]] == 2:
                        two_cut_count += 1
                if 3 == two_cut_count: # case 4B
                    total_count_case4B += 1
                    triangle_tables.append(gen_modified_mc_lut_case4B(indices))
                elif 2 == two_cut_count: # case 4D
                    total_count_case4D += 1
                    triangle_tables.append(gen_modified_mc_lut_case4D(indices))
            elif 8 == cut_count: # case 4C/4E
                two_cut_count = 0
                for i in range(len(indices)):
                    edges = vertex2edge[indices[i]]
                    if ebits[edges[0]] + ebits[edges[1]] + ebits[edges[2]] == 2:
                        two_cut_count += 1
                if 4 == two_cut_count:
                    total_count_case4C += 1
                    triangle_tables.append(gen_modified_mc_lut_case4C(indices))
                elif 2 == two_cut_count: # case 4E
                    total_count_case4E += 1
                    triangle_tables.append(gen_modified_mc_lut_case4E(indices))
            elif 12 == cut_count: # case 4F
                total_count_case4F += 1
                triangle_tables.append(gen_modified_mc_lut_case4F(indices))
            else:
                print("generate_triangle_tables!case 4X fatal error, cut_count={cut_count}".format(cut_count=cut_count))
        elif 5 == total_indices: # case 5A/5B/5C
            if 5 == cut_count: # case 5A
                total_count_case5A += 1
                triangle_tables.append(gen_modified_mc_lut_case5A(indices))
            elif 7 == cut_count: # case 5B
                total_count_case5B += 1
                triangle_tables.append(gen_modified_mc_lut_case5B(indices))
            elif 9 == cut_count: # case 5C
                total_count_case5C += 1
                triangle_tables.append(gen_modified_mc_lut_case5C(indices))
            else:
                print("generate_triangle_tables!case 5X fatal error, cut_count={cut_count}".format(cut_count=cut_count))
        elif 6 == total_indices: # case 6A/6B/6C
            if 4 == cut_count: # case 6A
                total_count_case6A += 1
                triangle_tables.append(gen_modified_mc_lut_case6A(indices))
            elif 6 == cut_count: # case 6B/6C
                is_same_facet_4 = False
                for vi1, vi2, vi3, vi4 in combinations(indices, 4):
                    if is_same_face_from_vertex_index_4(vi1, vi2, vi3, vi4):
                        is_same_facet_4 = True
                        break
                if is_same_facet_4:
                    total_count_case6B += 1
                    triangle_tables.append(gen_modified_mc_lut_case6B(indices))
                else:
                    total_count_case6C += 1
                    triangle_tables.append(gen_modified_mc_lut_case6C(indices))
            else:
                print("generate_triangle_tables!case 3X fatal error, cut_count={cut_count}".format(cut_count=cut_count))
        elif 7 == total_indices: # case 7
            if 3 == cut_count:
                total_count_case7 += 1
                triangle_tables.append(gen_modified_mc_lut_case7(indices))
            else:
                print("generate_triangle_tables!case 7X fatal error, cut_count={cut_count}".format(cut_count=cut_count))
        else:
            total_count_case0 += 1
            triangle_tables.append(gen_modified_mc_lut_case0(indices))
    
    print("total_case0: %d"%total_count_case0)
    print("total_case1: %d"%total_count_case1)
    print("total_case2A: %d"%total_count_case2A)
    print("total_case2B: %d"%total_count_case2B)
    print("total_case2C: %d"%total_count_case2C)
    print("total_case3A: %d"%total_count_case3A)
    print("total_case3B: %d"%total_count_case3B)
    print("total_case3C: %d"%total_count_case3C)
    print("total_case4A: %d"%total_count_case4A)
    print("total_case4B: %d"%total_count_case4B)
    print("total_case4C: %d"%total_count_case4C)
    print("total_case4D: %d"%total_count_case4D)
    print("total_case4E: %d"%total_count_case4E)
    print("total_case4F: %d"%total_count_case4F)
    print("total_case5A: %d"%total_count_case5A)
    print("total_case5B: %d"%total_count_case5B)
    print("total_case5C: %d"%total_count_case5C)
    print("total_case6A: %d"%total_count_case6A)
    print("total_case6B: %d"%total_count_case6B)
    print("total_case6C: %d"%total_count_case6C)
    print("total_case7: %d"%total_count_case7)
    return triangle_tables

def modified_mc_lut_save_to_cxx(path, edge_tables, triangle_tables):
    macros_defines_start_templ = '''
#ifndef __MARCHING_CUBES_LUT__
#define __MARCHING_CUBES_LUT__
    '''
    macros_defines_ending_templ = '''
#endif // __MARCHING_CUBES_LUT__
    '''
    edge_table_formatted_templ = '''
int edgeTable[256] = {value_list}
    '''
    triangle_table_formatted_templ = '''
int triTable[256][16] = {value_list}
    '''

    edge_value_list = "{\n"
    for i in range(32):
        sub_edge_value = edge_tables[i*8:i*8+8]
        edge_value_list += ", ".join([hex(sub_edge_value[j]) for j in range(8)])
        edge_value_list += ",\n"
    edge_value_list += "};"

    triangle_value_list = "{\n"
    for i in range(256):
        triangle_value_list += "{"
        triangle_value_list += ", ".join([str(triangle_tables[i][j]) for j in range(16)])
        triangle_value_list += "},\n"
    triangle_value_list += "};"

    with open(path, mode="w+") as fp:
        fp.write(macros_defines_start_templ)
        fp.write(edge_table_formatted_templ.format(value_list=edge_value_list))
        fp.write(triangle_table_formatted_templ.format(value_list=triangle_value_list))
        fp.write(macros_defines_ending_templ)

def modified_mc_lut_save_to_csharp(path, edge_tables, triangle_tables):
    edge_table_formatted_templ = '''
static const int edges[256] = {value_list}
    '''
    triangle_table_formatted_templ = '''
static const int triangulation[256][16] = {value_list}
    '''

    edge_value_list = "{\n"
    for i in range(32):
        sub_edge_value = edge_tables[i*8:i*8+8]
        edge_value_list += ", ".join([hex(sub_edge_value[j]) for j in range(8)])
        edge_value_list += ",\n"
    edge_value_list += "};"

    triangle_value_list = "{\n"
    for i in range(256):
        triangle_value_list += "{"
        triangle_value_list += ", ".join([str(triangle_tables[i][j]) for j in range(16)])
        triangle_value_list += "},\n"
    triangle_value_list += "};"

    with open(path, mode="w+") as fp:
        fp.write(edge_table_formatted_templ.format(value_list=edge_value_list))
        fp.write(triangle_table_formatted_templ.format(value_list=triangle_value_list))

if __name__ == "__main__":
    edge_tables = generate_edge_tables()
    print(edge_tables)
    triangle_tables = generate_triangle_tables()
    print(len(triangle_tables))
    #for i in range(len(triangle_tables)):
        #print("{item}\n".format(item=triangle_tables[i]))
    modified_mc_lut_save_to_csharp("./mc_lut.cs", edge_tables, triangle_tables)