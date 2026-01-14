
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt

def read_path(fname):
    return np.loadtxt(fname, delimiter=",", dtype=np.complex_)

def discretize_path(path, to_pixel, steps=20):
    pixmap = defaultdict(list)
    for i, (z0, z1) in enumerate(zip(path[:-1], path[1:])):
        ts = np.linspace(0, 1, steps, endpoint=True)
        zs = z0 + (z1 - z0) * ts
        for z in zs:
            pixmap[to_pixel(z)].append(i)  # keep starting index
    return pixmap

def intersect(path1, path2, resolution=1000, steps=20):
    all_points = np.concatenate([path1, path2])
    minx, maxx = np.min(all_points.real), np.max(all_points.real)
    miny, maxy = np.min(all_points.imag), np.max(all_points.imag)

    def to_pixel(z):
        x = int((z.real - minx) / (maxx - minx) * (resolution - 1))
        y = int((z.imag - miny) / (maxy - miny) * (resolution - 1))
        return x, y

    pixmap1 = discretize_path(path1, to_pixel, steps)
    pixmap2 = discretize_path(path2, to_pixel, steps)

    intersections = []
    for pix in set(pixmap1.keys()) & set(pixmap2.keys()):
        for i in pixmap1[pix]:
            for j in pixmap2[pix]:
                intersections.append((pix, i, j))
    return intersections

def det(a, b):
    return a[0]*b[1] - a[1]*b[0]

def segment_intersection(p1, p2, q1, q2):
    p1, p2, q1, q2 = map(lambda z: np.array([z.real, z.imag]), (p1, p2, q1, q2))
    r, s = p2 - p1, q2 - q1
    denom = det(r, s)
    if denom == 0:
        return None
    t = det(q1 - p1, s) / denom
    u = det(q1 - p1, r) / denom
    if 0 <= t <= 1 and 0 <= u <= 1:
        return p1 + t * r
    return None

def find_intersection(path1, path2):
    intersections = []
    for i in range(len(path1) - 1):
        for j in range(len(path2) - 1):
            inter = segment_intersection(path1[i], path1[i+1], path2[j], path2[j+1])
            if inter is not None:
                return (i, j, inter[0] + 1j*inter[1])
    return None

def resolve_intersection(p1, p2, path_segment_i, path_segment_j, I, II):
    """Try to find precise intersection in the given segment ranges
       and decide orientation."""
    inter = find_intersection(
        p1[path_segment_i[0]:path_segment_i[1] + 1, 0],
        p2[path_segment_j[0]:path_segment_j[1] + 1, 0]
    )
    if inter is None:
        return None
    
    di, dj, r = inter
    i0 = path_segment_i[0] + di
    j0 = path_segment_j[0] + dj
    
    v1 = [(p1[i0+1,0] - p1[i0,0]).real, (p1[i0+1,0] - p1[i0,0]).imag]
    v2 = [(p2[j0+1,0] - p2[j0,0]).real, (p2[j0+1,0] - p2[j0,0]).imag]
    
    
    y_ori = (min(abs(np.exp(p1[i0,1]) - np.exp(p2[j0,1])), 
                 abs(np.exp(p1[i0,2]) - np.exp(p2[j0,2]))) 
           - min(abs(np.exp(p1[i0,1]) - np.exp(p2[j0,2])), 
                 abs(np.exp(p1[i0,2]) - np.exp(p2[j0,1]))))
    x_ori = det(v1, v2)

    min_sheet_diff = min(min(abs(np.exp(p1[i0,1]) - np.exp(p2[j0,2])), 
                 abs(np.exp(p1[i0,2]) - np.exp(p2[j0,1]))), 
           min(abs(np.exp(p1[i0,1]) - np.exp(p2[j0,1])), 
                 abs(np.exp(p1[i0,2]) - np.exp(p2[j0,2]))))
    if min_sheet_diff > 0.1:
        print(min_sheet_diff)        
        print(p1[i0])
        print(p2[j0])
        orientation = "Keine intersection"
    elif y_ori * x_ori > 0: 
        orientation = f"{I} nach {II}"
    else:
        orientation = f"{II} nach {I}"
    
    return (i0, j0, r, orientation)

# Example usage

cycle = {"A" : 0,
         "B" : 1,
         "C" : 2,
         "D" : 3,
         "E" : 4,
         "F" : 5,
         }
I = "C"
II = "F"

p1 = read_path(f"data/path_data_{cycle[I]}/path_data_12.csv")
p2 = read_path(f"data/path_data_{cycle[II]}/path_data_10.csv")
intersections = intersect(p1[:,0], p2[:,0], resolution=1000, steps=20)
print(len(intersections))

path_segment_i = [0,0]
path_segment_j = [0,0]
last_pix = (-100, -100)

for pix, i, j in intersections:
    d = max(abs(pix[0] - last_pix[0]), abs(pix[1] - last_pix[1]))
    if d <= 10:
        path_segment_i[1] = i + 1
        path_segment_j[1] = j + 1
        continue
    last_pix = pix
    
    res = resolve_intersection(p1, p2, path_segment_i, path_segment_j, I, II)
    if res is None:
        path_segment_i = [i - 2, i + 1]
        path_segment_j = [j - 2, j + 1]
        continue
    i0, j0, r, orientation = res
    print(f"Intersection at {r}, orientation: {orientation}")
    
    path_segment_i = [i - 2, i + 1]
    path_segment_j = [j - 2, j + 1]

# Final check for leftover segment
res = resolve_intersection(p1, p2, path_segment_i, path_segment_j, I, II)
if res:
    i0, j0, r, orientation = res
    print(f"Final intersection at {r}, orientation: {orientation}")

