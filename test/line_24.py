import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# 立方体の8つの頂点を定義
cube_vertices = np.array([[-0.5, -0.5, -0.5],
                          [0.5, -0.5, -0.5],
                          [0.5, 0.5, -0.5],
                          [-0.5, 0.5, -0.5],
                          [-0.5, -0.5, 0.5],
                          [0.5, -0.5, 0.5],
                          [0.5, 0.5, 0.5],
                          [-0.5, 0.5, 0.5]])

# 立方体のエッジを描画
edges = [
    [cube_vertices[0], cube_vertices[1]], [cube_vertices[1], cube_vertices[2]],
    [cube_vertices[2], cube_vertices[3]], [cube_vertices[3], cube_vertices[0]],
    [cube_vertices[4], cube_vertices[5]], [cube_vertices[5], cube_vertices[6]],
    [cube_vertices[6], cube_vertices[7]], [cube_vertices[7], cube_vertices[4]],
    [cube_vertices[0], cube_vertices[4]], [cube_vertices[1], cube_vertices[5]],
    [cube_vertices[2], cube_vertices[6]], [cube_vertices[3], cube_vertices[7]]
]
for edge in edges:
    ax.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], [edge[0][2], edge[1][2]], 'k-')

# 中心から立方体の面や頂点方向に線を描画
directions = [
    (1, 0.37796447, 0.37796447), (1, 0.37796447, -0.37796447), (1, -0.37796447, 0.37796447), (1, -0.37796447, -0.37796447),
    (-1, 0.37796447, 0.37796447), (-1, 0.37796447, -0.37796447), (-1, -0.37796447, 0.37796447), (-1, -0.37796447, -0.37796447),
    (0.37796447, 1, 0.37796447), (-0.37796447, 1, 0.37796447), (0.37796447, 1, -0.37796447), (-0.37796447, 1, -0.37796447),
    (0.37796447, -1, 0.37796447), (-0.37796447, -1, 0.37796447), (0.37796447, -1, -0.37796447), (-0.37796447, -1, -0.37796447),
    (0.37796447, 0.37796447, 1), (-0.37796447, 0.37796447, 1), (0.37796447, -0.37796447, 1), (-0.37796447, -0.37796447, 1),
    (0.37796447, 0.37796447, -1), (-0.37796447, 0.37796447, -1), (0.37796447, -0.37796447, -1), (-0.37796447, -0.37796447, -1)
    # 必要に応じて他の方向も指定してください
]

for d in directions:
    ax.plot([0, 0.5 * d[0]], [0, 0.5 * d[1]], [0, 0.5 * d[2]], 'b-')

# 軸とグリッドを非表示
ax.set_axis_off()

# グラフの表示

plt.savefig('24line.png')
