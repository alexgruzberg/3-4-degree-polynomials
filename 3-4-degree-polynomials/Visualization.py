from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

# Visualization of the roots of cubic
def visualizer(rotat_angle, corr, view_angle, a, b, c, d):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    # Range of X and Y
    X = np.arange(-corr, corr, 0.15)
    Y = np.arange(-corr, corr, 0.15)
    X, Y = np.meshgrid(X, Y)

    Z1 = np.empty_like(X)
    Z2 = np.empty_like(X)
    C1 = np.empty_like(X, dtype=object)
    C2 = np.empty_like(X, dtype=object)

    for i in range(len(X)):
        for j in range(len(X[0])):
            x = X[i][j]
            y = Y[i][j]
            z1 = (a*x+b)*(x*x-y*y)+x*(c-2*a*y*y)+d
            z2 = y*(a*(3*x*x-y*y)+2*b*x+c)
            Z1[i, j] = z1
            Z2[i, j] = z2

            C1[i, j] = plt.get_cmap("Greens")(0.5)
            C2[i, j] = plt.get_cmap("Blues")(z2)

    # Create a transparent bridge region
    X_bridge = np.vstack([X[-1, :], X[-1, :]])
    Y_bridge = np.vstack([Y[-1, :], Y[-1, :]])
    Z_bridge = np.vstack([Z1[-1, :], Z2[-1, :]])

    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_zlabel("Z axis")

    ax.set_zlim(-8 * corr * corr, 8 * corr * corr)
    color_bridge = np.empty_like(Z_bridge, dtype=object)
    color_bridge.fill((1, 1, 1, 0))

    # Join the two surfaces flipping one of them (using also the bridge)
    X_full = np.vstack([X, X_bridge, np.flipud(X)])
    Y_full = np.vstack([Y, Y_bridge, np.flipud(Y)])
    Z_full = np.vstack([Z1, Z_bridge, np.flipud(Z2)])
    color_full = np.vstack([C1, color_bridge, np.flipud(C2)])

    surf_full = ax.plot_surface(X_full, Y_full, Z_full, rstride=1, cstride=1, shade=True,
                                facecolors=color_full, linewidth=1,
                                antialiased=False)
    ax.view_init(view_angle, rotat_angle)
    fig.set_size_inches(8, 8)
    ax.grid(False)

    plt.draw()

a, b, c, d = map(float, input().split())
for angle in range(0, 4):
    visualizer(30 * angle, max(3, (abs(a) + abs(b) + abs(c) + abs(d)) // 3), 45 + 15 * angle, a, b, c, d)
    plt.pause(0.05)



