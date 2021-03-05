import numpy as np
from AnglePi import clip_angle, normalize_angle, difference_angle

def lines_arc(new_sets,radius_of_tool):
    new_sets_y_bounds = []
    for i in range(len(new_sets)):
        new_sets_y_bounds.append([np.min(new_sets[i][0]),np.max(new_sets[i][0])])
    line_set = []
    for i in range(len(new_sets)):
        y1 = new_sets_y_bounds[i][0]
        y2 = new_sets_y_bounds[i][1]
        if y2 - y1 > radius_of_tool:  
            y_interval_t = np.arange(y1+radius_of_tool,y2,radius_of_tool*2)

            if y2 - y_interval_t[-1] > radius_of_tool:
                y_interval_t = np.hstack((y_interval_t,np.asarray([np.mean([y2, y_interval_t[-1]])])))
            line_y = []
            line_z = []
            for ty in y_interval_t:
                for j in range(len(new_sets[i][0])-1):
                    if ty >= new_sets[i][0][j] and ty < new_sets[i][0][j+1]:
                        line_y.extend([ty,ty])
                        zmin1 = new_sets[i][1][j]
                        zmin2 = new_sets[i][1][j+1]
                        zmax1 = new_sets[i][2][j]
                        zmax2 = new_sets[i][2][j+1]
                        if abs(zmin2-zmin1) > np.pi:
                            if zmin1 < 0:
                                zmin1 += np.pi*2
                            else:
                                zmin2 += np.pi*2
                        if abs(zmax2-zmax1) > np.pi:
                            if zmax1 < 0:
                                zmax1 += np.pi*2
                            else:
                                zmax2 += np.pi*2
                        zmin = zmin1+(zmin2-zmin1)*(ty-new_sets[i][0][j])/(new_sets[i][0][j+1]-new_sets[i][0][j])
                        zmax = zmax1+(zmax2-zmax1)*(ty-new_sets[i][0][j])/(new_sets[i][0][j+1]-new_sets[i][0][j])
                        if zmin > np.pi:
                            zmin -= np.pi*2
                        if zmax > np.pi:
                            zmax -= np.pi*2
                        if zmax > zmin:
                            if zmax - zmin > 2*radius_of_tool/ty:
                                line_z.extend([zmin+radius_of_tool/ty,zmax-radius_of_tool/ty])
                            else:
                                line_z.extend([(zmin+zmax)/2]*2)
                        else:
                            zmax += np.pi*2
                            if zmax - zmin > 2*radius_of_tool/ty:
                                line_z.extend([zmin+radius_of_tool/ty,zmax-radius_of_tool/ty])
                                if line_z[-1] > np.pi:
                                    line_z[-1] -= np.pi*2
                            else:
                                line_z.extend([(zmin+zmax)/2]*2)

                        break

        else:
            line_y = [(new_sets_y_bounds[i][1] + new_sets_y_bounds[i][1])/2]*2
            line_z = [np.mean(new_sets[i][1]),np.mean(new_sets[i][2])]
            if line_z[1]-line_z[0]> 2* radius_of_tool/ty:
                line_z[0] += radius_of_tool/ty
                line_z[1] -= radius_of_tool/ty
            else:
                mid = np.mean(line_z)
                line_z = [mid]*2
        line_set.append([line_y,line_z])
    return line_set

def path_arc(path, flag = 1):
    y,z = path
    path_y = []
    path_z = []
    l = 0
    clock = []
    if flag == 1:
        path_y.append(y[0])
        path_z.append(z[0])
        path_y.append(y[1])
        path_z.append(z[1])
        l += y[0]*difference_angle(z[0],z[1])
        clock.extend([1,0])
    elif flag == 2:
        path_y.append(y[1])
        path_z.append(z[1])
        path_y.append(y[0])
        path_z.append(z[0])
        l += y[0]*difference_angle(z[0],z[1])
        clock.extend([-1,0])
    elif flag == 3:
        path_y.append(y[-1])
        path_z.append(z[-1])
        path_y.append(y[-2])
        path_z.append(z[-2])
        l += y[0]*difference_angle(z[0],z[1])
        clock.extend([-1,0])
    elif flag == 4:
        path_y.append(y[-2])
        path_z.append(z[-2])
        path_y.append(y[-1])
        path_z.append(z[-1])
        l += y[0]*difference_angle(z[0],z[1])
        clock.extend([1,0])
        
    for i in range(int(len(y)/2)-1):
        if flag == 1 or flag == 2:
            index = 2*i + 2
        elif flag == 3 or flag == 4:
            index = -2*i-4
        p1 = np.asarray([path_y[-1]*np.cos(path_z[-1]),path_y[-1]*np.sin(path_z[-1])])
        p2 = np.asarray([y[index]*np.cos(z[index]),y[index]*np.sin(z[index])])
        p3 = np.asarray([y[index+1]*np.cos(z[index+1]),y[index+1]*np.sin(z[index+1])])
        l12 = np.linalg.norm(p1-p2)
        l13 = np.linalg.norm(p1-p3)
        if l12< l13:
            path_y.append(y[index])
            path_z.append(z[index])
            path_y.append(y[index+1])
            path_z.append(z[index+1])
            clock.extend([1,0])
            l += l12
        else:
            path_y.append(y[index+1])
            path_z.append(z[index+1])
            path_y.append(y[index])
            path_z.append(z[index])
            clock.extend([-1,0])
            l += l13
        l += y[index]*difference_angle(z[index],z[index+1])
    sf = [[path_y[0]*np.cos(path_z[0]),path_y[0]*np.sin(path_z[0])],[path_y[-1]*np.cos(path_z[-1]),path_y[-1]*np.sin(path_z[-1])]]
    new_path = [path_y,path_z]
    return [new_path, l, sf, clock ]
    
def lines_vertical(new_sets,radius_of_tool):
    new_sets_y_bounds = []
    for i in range(len(new_sets)):
        new_sets_y_bounds.append([np.min(new_sets[i][0]),np.max(new_sets[i][0])])

    line_set = []
    for i in range(len(new_sets)):
        y1 = new_sets_y_bounds[i][0]
        y2 = new_sets_y_bounds[i][1]
        if y2 - y1 > radius_of_tool:
            y_interval_t = np.arange(y1+radius_of_tool,y2,radius_of_tool*2)
            if y2 - y_interval_t[-1] > radius_of_tool:
                y_interval_t = np.hstack((y_interval_t,np.asarray([np.mean([y2, y_interval_t[-1]])])))
            line_y = []
            line_z = []
            for ty in y_interval_t:
                for j in range(len(new_sets[i][0])-1):
                    if ty >= new_sets[i][0][j] and ty < new_sets[i][0][j+1]:
                        line_y.extend([ty,ty])
                        zmin = new_sets[i][1][j]+(new_sets[i][1][j+1]-new_sets[i][1][j])*(ty-new_sets[i][0][j])/(new_sets[i][0][j+1]-new_sets[i][0][j])
                        zmax = new_sets[i][2][j]+(new_sets[i][2][j+1]-new_sets[i][2][j])*(ty-new_sets[i][0][j])/(new_sets[i][0][j+1]-new_sets[i][0][j])
                        if zmax - zmin > 2*radius_of_tool:
                            line_z.extend([zmin+radius_of_tool,zmax-radius_of_tool])
                        else:
                            line_z.extend([(zmin+zmax)/2]*2)
                        break

        else:
            line_y = [(new_sets_y_bounds[i][1] + new_sets_y_bounds[i][1])/2]*2
            line_z = [np.mean(new_sets[i][1]),np.mean(new_sets[i][2])]
            if line_z[1]-line_z[0]> 2* radius_of_tool:
                line_z[0] += radius_of_tool
                line_z[1] -= radius_of_tool
            else:
                mid = np.mean(line_z)
                line_z = [mid]*2
        line_set.append([line_y,line_z])
    return line_set

def path_vertical(path, flag = 1):
    y,z = path
    path_y = []
    path_z = []
    l = 0
    if flag == 1:
        path_y.append(y[0])
        path_z.append(z[0])
        path_y.append(y[1])
        path_z.append(z[1])
        l += np.linalg.norm(np.asarray([y[0]-y[1],z[0]-z[1]]))
    elif flag == 2:
        path_y.append(y[1])
        path_z.append(z[1])
        path_y.append(y[0])
        path_z.append(z[0])
        l += np.linalg.norm(np.asarray([y[0]-y[1],z[0]-z[1]]))
    elif flag == 3:
        path_y.append(y[-1])
        path_z.append(z[-1])
        path_y.append(y[-2])
        path_z.append(z[-2])
        l += np.linalg.norm(np.asarray([y[-1]-y[-2],z[-1]-z[-2]]))
    elif flag == 4:
        path_y.append(y[-2])
        path_z.append(z[-2])
        path_y.append(y[-1])
        path_z.append(z[-1])
        l += np.linalg.norm(np.asarray([y[-1]-y[-2],z[-1]-z[-2]]))
        
    for i in range(int(len(y)/2)-1):
        if flag == 1 or flag == 2:
            index = 2*i + 2
        elif flag == 3 or flag == 4:
            index = -2*i-4
        p1 = np.asarray([path_y[-1],path_z[-1]])
        p2 = np.asarray([y[index],z[index]])
        p3 = np.asarray([y[index+1],z[index+1]])
        l12 = np.linalg.norm(p1-p2)
        l13 = np.linalg.norm(p1-p3)
        if l12< l13:
            path_y.append(y[index])
            path_z.append(z[index])
            path_y.append(y[index+1])
            path_z.append(z[index+1])
            l += l12
        else:
            path_y.append(y[index+1])
            path_z.append(z[index+1])
            path_y.append(y[index])
            path_z.append(z[index])
            l += l13
        l += np.linalg.norm(p2-p3)
    sf = [[path_y[0],path_z[0]],[path_y[-1],path_z[-1]]]
    new_path = [path_y,path_z]
    return [new_path, l, sf ]