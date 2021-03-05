import numpy as np
import copy
from AnglePi import clip_angle, normalize_angle, difference_angle

def MorseArc(final_contour_points, contour_bounds, scope, step):
    # final_contour_points: list[list[list]] [contour 1:[ylist,zlist], contour 2:....]
    # contour_bounds: list[list] [contour 1: [ymin, ymax, zmin, zmax], contour 2:...]
    # scope: ymin, ymax
    # step: np.arange(ymin, ymax, step)
    ymin, ymax = scope
    y_interval = np.arange(ymin,ymax+step,step)
    sort_z = []
    for m in range(y_interval.shape[0]):
        y = y_interval[m]
        tz = []
        for i in range(len(final_contour_points)):
            if y < contour_bounds[i][1] and y >= contour_bounds[i][0]:
                ttz = []
                for j in range(len(final_contour_points[i][0])-1):
                    p1 = [final_contour_points[i][0][j],final_contour_points[i][1][j]]
                    p2 = [final_contour_points[i][0][j+1],final_contour_points[i][1][j+1]]
                    if abs(p2[1]-p1[1]) > np.pi:
                        if p1[1] < 0 :
                            p1[1] += np.pi*2
                        else:
                            p2[1] += np.pi*2
                    if y <= p1[0] and y > p2[0]:
                        ttz.append((y-p1[0])/(p2[0]-p1[0])*(p2[1]-p1[1])+p1[1])
                        if ttz != []:
                            if ttz[-1] > np.pi:
                                ttz[-1] -= np.pi*2
                    elif y > p1[0] and y <= p2[0]:
                        ttz.append((y-p1[0])/(p2[0]-p1[0])*(p2[1]-p1[1])+p1[1])
                        if ttz != []:
                            if ttz[-1] > np.pi:
                                ttz[-1] -= np.pi*2
                if ttz != []:
                    # assume all the angles sweeping a certain contour are less than 180 degree.
                    if abs(np.max(ttz) - np.min(ttz)) > np.pi:
                        tz.extend([np.max(ttz),np.min(ttz)])
                    else:
                        tz.extend([np.min(ttz),np.max(ttz)])
        if tz == []:
            continue   
        tz.sort()
        sort_z.append(tz)

    link = []
    stack_link = []
    sets = []
    for i in range(len(sort_z)):
        if sets == []:
            tz = []
            for j in range(int(len(sort_z[i])/2)):
                tz.append([[y_interval[i]],[sort_z[i][2*j]],[sort_z[i][2*j+1]]])
                stack_link.append(j)
            sets.extend(tz)
        else:  
            new_stack_link = []
            tz = []
            ttl = []
            # all arcs 
            counts = np.zeros(len(stack_link))
            for j in range(int(len(sort_z[i])/2)):
                tz.append([[y_interval[i]],[sort_z[i][2*j]],[sort_z[i][2*j+1]]])
                tl = []
                # each arc belongs to which stack_link
                p1 = sort_z[i][2*j+0]
                p2 = sort_z[i][2*j+1]
                if p2 < p1:
                    p2 += np.pi*2
                for k in range(len(stack_link)):
                    q1 = sets[stack_link[k]][1][-1]
                    q2 = sets[stack_link[k]][2][-1]
                    if q2 < q1 :
                        q2 += 2*np.pi
                    if np.ceil((q1-p1)/2/np.pi) <= np.floor((q2-p1)/2/np.pi):
                        counts[k] += 1
                        tl.append(k)
                    elif np.ceil((q1-p2)/2/np.pi) <= np.floor((q2-p2)/2/np.pi):
                        counts[k] += 1
                        tl.append(k)
                    elif np.ceil((p1-q1)/2/np.pi) <= np.floor((p2-q1)/2/np.pi):
                        counts[k] += 1
                        tl.append(k)
                    elif np.ceil((p1-q2)/2/np.pi) <= np.floor((p2-q2)/2/np.pi):
                        counts[k] += 1
                        tl.append(k)
                ttl.append(tl)
            for j in range(len(ttl)):
                tl = ttl[j]
                if tl == []:
                    new_stack_link.append(len(sets))
                    sets.append(tz[j])
                elif len(tl) > 1 :
                    new_stack_link.append(len(sets))
                    for k in range(len(tl)):
                        link.append([len(sets),stack_link[tl[k]]])
                    sets.append(tz[j])
                elif len(tl) == 1:
                    if counts[tl[0]] == 1:
                        for k in range(3):
                            sets[stack_link[tl[0]]][k].append(tz[j][k][0])
                        new_stack_link.append(stack_link[tl[0]])
                    elif counts[tl[0]] > 1:
                        new_stack_link.append(len(sets))
                        link.append([len(sets),stack_link[tl[0]]])
                        sets.append(tz[j])
                    else:
                        print('counts == 0:',counts)
            stack_link = new_stack_link
    
    return sets

def MorseD(final_contour_points, contour_bounds, scope, step):
    # final_contour_points: list[list[list]] [contour 1:[ylist,zlist], contour 2:....]
    # contour_bounds: list[list] [contour 1: [ymin, ymax, zmin, zmax], contour 2:...]
    # scope: ymin, ymax
    # step: np.arange(ymin, ymax, step)
    ymin, ymax = scope
    y_interval = np.arange(ymin,ymax,step)
    sort_z = []
    empty_index = []
    for m in range(y_interval.shape[0]):
        y = y_interval[m]
        tz = []
        for i in range(len(final_contour_points)):
            if y < contour_bounds[i][1] and y >= contour_bounds[i][0]:
                for j in range(len(final_contour_points[i][0])-1):
                    p1 = [final_contour_points[i][0][j],final_contour_points[i][1][j]]
                    p2 = [final_contour_points[i][0][j+1],final_contour_points[i][1][j+1]]
                    if y <= p1[0] and y > p2[0]:
                        tz.append((y-p1[0])/(p2[0]-p1[0])*(p2[1]-p1[1])+p1[1])
                    elif y > p1[0] and y <= p2[0]:
                        tz.append((y-p1[0])/(p2[0]-p1[0])*(p2[1]-p1[1])+p1[1])
        tz.sort()
        sort_z.append(tz)
        if tz == []:
            empty_index.append(m)
    for i in range(len(empty_index)):
        index = empty_index[-1-i]
        sort_z.pop(index)
        y_interval = np.delete(y_interval, index)
    nums = []
    link = []
    sets = []
    pre_num = 0
    for i in range(len(sort_z)):
        if len(sort_z[i]) % 2 == 0 and len(sort_z[i]) > 0:
            nums.append(len(sort_z[i]))
            if sets == []:
                tz = []
                for j in range(int(len(sort_z[i])/2)):
                    tz.append([[y_interval[i]],[sort_z[i][2*j]],[sort_z[i][2*j+1]]])
                sets.append(tz)
                link.append([])
            else:  
                if nums[-1] == pre_num:
                    for j in range(int(len(sort_z[i])/2)):
                        sets[-1][j][0].append(y_interval[i])
                        sets[-1][j][1].append(sort_z[i][2*j])
                        sets[-1][j][2].append(sort_z[i][2*j+1])
                else:
                    tz = []
                    tll = []
                    for j in range(int(len(sort_z[i])/2)):
                        tz.append([[y_interval[i]],[sort_z[i][2*j]],[sort_z[i][2*j+1]]])
                        tl = []
                        p1 = sort_z[i][2*j+0]
                        p2 = sort_z[i][2*j+1]
                        for k in range(len(sets[-1])):
                            q1 = sets[-1][k][1][-1]
                            q2 = sets[-1][k][2][-1]
                            if p1 <= q2 and p1 > q1:
                                if k not in tl:
                                    tl.append(k)
                            elif p2 <= q2 and p2 > q1:
                                if k not in tl:
                                    tl.append(k)
                            elif q1 <= p2 and q1 > p1:
                                if k not in tl:
                                    tl.append(k)
                            elif q2 <= p2 and q2 > p1:
                                if k not in tl:
                                    tl.append(k)
                        tll.append(tl)
                    link.append(tll)
                    sets.append(tz)
            pre_num = nums[-1]
#     if show == 1:
#         for i in range(len(final_contour)):
#             x, y = final_contour_points[i]
#             plt.plot(x,y,color = 'k')
#         plt.axis('equal')
#         for i in range(len(sets)):
#             for j in range(len(sets[i])):
#                 plt.plot(sets[i][j][0]+sets[i][j][0][::-1],sets[i][j][1]+sets[i][j][2][::-1])
#         plt.show()
    
    new_sets = []
    stack_link = []
    new_link = []
    for i in range(len(sets)):
        if new_sets == []:
            for j in range(len(sets[i])):
                new_sets.append(copy.deepcopy(sets[i][j]))
                stack_link.append(j)
        else:
            t_stack_link = [-1]*len(link[i])
            j = 0
            while 1:
                single_set_link = link[i][j]
                if len(link[i][j]) == 0:
                    t_stack_link[j] = len(new_sets)
                    new_sets.append(copy.deepcopy(sets[i][j]))
                elif len(link[i][j]) > 1:
                    t_stack_link[j]= len(new_sets)
                    new_sets.append(copy.deepcopy(sets[i][j]))
                    for index in link[i][j]:
                        new_link.append([t_stack_link[j], stack_link[index]])
                elif len(single_set_link) == 1:
                    tk = 0
                    tp = 0
                    for k in range(len(link[i])-j-1):
                        if len(link[i][j+k+1]) == 1 and single_set_link[0] == link[i][j+k+1][0]:
                            tk += 1
                            continue
                        else:
                            break
                    if tk == 0:
                        t_stack_link[j] = stack_link[single_set_link[0]]
                        for m in range(3):
                            new_sets[t_stack_link[j]][m].extend(copy.deepcopy(sets[i][j][m]))
                    else:
                        for p in range(tk+1):
                            tp = p
                            t_stack_link[j+p] = len(new_sets)
                            new_sets.append(copy.deepcopy(sets[i][j+p]))
                            new_link.append([t_stack_link[j+p],stack_link[single_set_link[0]]])
    #                 print(tp,tk)
                    j = j+tp
                j += 1
                if j == len(link[i]):
                    break
            stack_link = t_stack_link
#     if show == 2 :
#         plt.axis('equal')
#         for i in range(len(new_sets)):
#             p = plt.plot(new_sets[i][0]+new_sets[i][0][::-1],new_sets[i][1]+new_sets[i][2][::-1],linewidth = 1)
#             color = p[0].get_color()
#             plt.plot([new_sets[i][0][0],new_sets[i][0][0]],[new_sets[i][1][0],new_sets[i][2][0]],color = color,linewidth = 1)
#         plt.show()
    return new_sets, new_link

