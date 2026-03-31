import numpy as np
import matplotlib.pyplot as plt

# Adjacency matrix
A = np.array([
    [0, 0, 0, 2, 1, 1],
    [0, 0, 0, 1, 2, 1],
    [0, 0, 0, 0, 1, 2],
    [0, 1, 1, 0, 0, 0],
    [1, 2, 1, 0, 0, 0],
    [0, 1, 0, 0, 0, 0]
])

num_nodes = A.shape[0]

# 3+3 bipartite layout
group1 = [0, 1, 2]
group2 = [3, 4, 5]

y_spacing = 1.0
pos = {}
for idx, node in enumerate(group1):
    pos[node] = np.array([-1, y_spacing*(len(group1)-1-idx)])
for idx, node in enumerate(group2):
    pos[node] = np.array([1, y_spacing*(len(group2)-1-idx)])

fig, ax = plt.subplots(figsize=(6,4))

# Draw nodes
for node in group1:
    ax.scatter(*pos[node], s=1200, color='skyblue', zorder=2)
for node in group2:
    ax.scatter(*pos[node], s=1200, color='lightgreen', zorder=2)

# Draw labels
for i in range(num_nodes):
    ax.text(pos[i][0], pos[i][1], str(i+1), ha='center', va='center', fontsize=12, weight='bold', zorder=3)

# Draw arrows with **only perpendicular offset**
for i in range(num_nodes):
    for j in range(num_nodes):
        forward = A[i,j]
        backward = A[j,i]
        total = forward + backward
        if total == 0 or i >= j:  # handle each pair once
            continue
        start = pos[i] + i*0.01
        end = pos[j] + i *0.01
        vec = end - start
        length = np.linalg.norm(vec)
        if length == 0:
            continue
        u = vec / length
        perp = np.array([0, 1])  # perpendicular vector

        for k in range(total):
            offset = 0.08*(k - (total-1)/2)
            # Decide direction for forward/backward
            if k < forward:
                s = start + perp*offset
                e = end + perp*offset
            else:
                s = end + perp*offset
                e = start + perp*offset
            ax.annotate(
                '', xy=e, xytext=s,
                arrowprops=dict(arrowstyle='-|>', color='gray', lw=1.0, shrinkA=0, shrinkB=0)
            )

ax.set_xlim(-1.5,1.5)
ax.set_ylim(-0.5, y_spacing*(len(group1)-0.5))
ax.axis('off')
ax.set_title("Bipartite Quiver with Perpendicular Offset Only")
plt.show()
