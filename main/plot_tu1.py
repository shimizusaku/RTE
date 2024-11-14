import matplotlib.pyplot as plt
import numpy as np

# load
f = np.load('/mnt/solar07a/c0287shimizu/data/d001/MultiRay/fortran.npz')

# point of tu=1
x_tu1  = f['x_tu1']
x_tu1_ = x_tu1[1:, 150]
x_ = np.arange(256)


# 四角形の頂点
x = [0, 255, 255, 0, 0]
y = [0, 0, 255, 255, 0]

# 四角形を描画
plt.plot(x, y, color='black')
plt.plot(x_, x_tu1_, color='red')

plt.axis('off')
plt.axis('equal')


# 描画を表示
save_dir = '/mnt/solar07a/c0287shimizu/fig/test'
file_name = 'tu1_point.png'
file_path = os.path.join(save_dir, file_name)

plt.tight_layout()
plt.savefig(file_path)
plt.close()