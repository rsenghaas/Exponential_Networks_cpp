import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

# --- Step 1: exponent triples ---
def exponent_triples(F, variables):
    F = sp.expand(F)
    triples = []
    for m in F.as_ordered_terms():
        powers = m.as_powers_dict()
        i = powers.get(variables[0], 0)
        j = powers.get(variables[1], 0)
        k = powers.get(variables[2], 0)
        triples.append((int(i), int(j), int(k)))
    return triples

# --- Step 2: compute lower hull faces ---
def lower_hull_faces(triples, t):
    points = np.array([[i, j, k*t] for (i,j,k) in triples])
    hull = ConvexHull(points)
    lower_faces = []

    for simplex in hull.simplices:
        p1, p2, p3 = points[simplex]
        v1 = p2 - p1
        v2 = p3 - p1
        normal = np.cross(v1, v2)

        # lower hull: downward pointing normal
        if normal[2] < 0:
            lower_faces.append(simplex)
    return points, lower_faces

# --- Step 3: project lower faces and build dual graph ---
def dual_graph(points, faces):
    face_centers = []
    edges = set()

    for face in faces:
        pts = points[face][:,:2]  # project to (i,j)
        center = pts.mean(axis=0)
        face_centers.append(center)

    # naive edge: connect faces sharing at least 2 points
    for i, f1 in enumerate(faces):
        for j, f2 in enumerate(faces):
            if i >= j:
                continue
            if len(set(f1) & set(f2)) >= 2:
                edges.add((i,j))

    # ensure 2D array
    face_centers = np.array(face_centers, dtype=float)
    if face_centers.ndim == 1:
        face_centers = face_centers.reshape(1, -1)

    return face_centers, list(edges)

def external_edges(points, hull_vertices):
    # hull_vertices: coordinates of Newton polytope vertices
    # add rays along edges connecting boundary points

    edges = []
    N = len(hull_vertices)
    for i in range(N):
        p1 = hull_vertices[i]
        p2 = hull_vertices[(i+1)%N]
        edges.append((p1, p2))
    return edges

# --- Step 4: plot dual tropical curve ---
def plot_dual_graph(face_centers, edges):
    if len(face_centers) == 0:
        print("No lower faces → degenerate tropical curve.")
        return

    face_centers = np.array(face_centers, dtype=float)
    if face_centers.ndim == 1:
        face_centers = face_centers.reshape(1,-1)

    plt.figure(figsize=(6,6))
    plt.scatter(face_centers[:,0], face_centers[:,1], color='red')

    for i,j in edges:
        p1 = face_centers[i]
        p2 = face_centers[j]
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color='blue')

    plt.gca().set_aspect('equal')
    plt.xlabel('i')
    plt.ylabel('j')
    plt.title('Dual Graph / Tropical Curve')
    plt.show()

# --- Example usage ---
x,y,Q = sp.symbols('x y Q')

F = (1 + x + y + Q * x * y)

triples = exponent_triples(F, (x,y,Q))
t = 1.0
points, lower_faces = lower_hull_faces(triples, t)
face_centers, edges = dual_graph(points, lower_faces)
plot_dual_graph(face_centers, edges)
