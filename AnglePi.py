import numpy as np

def clip_angle(angle, amin = 0, amax = np.pi):
    if angle< amin:
        angle = angle
    elif angle > amax:
        angle = amax
    return angle

def normalize_angle(theta):
    if theta > np.pi:
        while 1:
            theta -= np.pi*2
            if theta < np.pi:
                break
    if theta < -np.pi:
        while 1:
            theta += np.pi*2
            if theta > -np.pi:
                break
    return theta
def difference_angle(theta1, theta2):
    if theta2 > theta1 :
        return theta2 - theta1
    else:
        return theta2-theta1 + np.pi*2
