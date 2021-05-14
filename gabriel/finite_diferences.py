import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

class VectorGrid2:
    '''Class representing staggered grid, u and v component are stored at face
    centers, left and bottom respectively
    '''
    def __init__(self, resolution, grid_spacing, sample=None, data=None):
        self.resolution = resolution
        self.grid_spacing = grid_spacing
        self.sample = sample

        x_points = int(resolution[0]/float(grid_spacing[0]))
        y_points = int(resolution[1]/float(grid_spacing[1]))
        x = np.linspace(0, resolution[0] - grid_spacing[0], x_points)
        y = np.linspace(0, resolution[1] - grid_spacing[1], y_points)
        
        self.x = x
        self.y = y
        xx, yy = np.meshgrid(x,y)
        self.data = np.zeros((y_points,x_points,2))
        self.u = np.zeros((y_points, x_points))
        self.v = np.zeros((y_points, x_points))
        if sample is not None:
            for i in range(y_points):
                for j in range(x_points):
                    self.data[i][j] = sample(x[j],y[i])
                    self.u[i][j], self.v[i][j] = self.data[i][j]
        elif data is not None:
            self.data = data.copy()
            m, n, _ = self.data.shape
            self.u = np.zeros((m,n))
            self.v = np.zeros((m,n))
            for i in range(m):
                for j in range(n):
                    self.u[i][j], self.v[i][j] = self.data[i][j]

    def __getitem__(self, pos):
        '''Translates cartesian coord to matrice and return element

        parameters:
            - pos: cartesian coordinates of element
        '''
        x, y = pos
        return self.data[y][x]

    def plot(self, stream=False):
        if stream:
            plt.streamplot(self.x,self.y,self.u,self.v)
        else:
            plt.quiver(self.x, self.y + 0.5 * self.grid_spacing[1], self.u, np.zeros_like(self.u))
            plt.quiver(self.x + 0.5 * self.grid_spacing[0], self.y, np.zeros_like(self.v), self.v)
        plt.grid(True)
        
    def divergent(self):
        m, n, _ = self.data.shape
        
        div = np.zeros((m,n))
        for i in range(m):
            for j in range(n):
                div[i][j] = self.divergent_at_point(j,i)
        
        return ScalarGrid2(self.resolution,self.grid_spacing,data=div)

    def divergent_at_point(self, x, y):
        m, n, _ = self.data.shape
        dx, dy = self.grid_spacing
        center = self[x,y]

        uright, _ = self[x+1, y] if x != n-1 else center
        uleft, _ = self[x-1, y] if x != 0 else center
        _, vup = self[x, y+1] if y != m-1 else center
        _, vdown = self[x, y-1] if y != 0 else center
        
        return (uright - uleft) / dx + (vup - vdown) / dy
    
    def laplacian(self):
        m, n, _ = self.data.shape

        lapl = np.zeros_like(self.data)
        for i in range(m):
            for j in range(n):
                #iterate over columns which is equivalent to iterate over x axis
                lapl[i][j] = self.laplacian_at_point(j, i)
        
        return VectorGrid2(self.resolution, self.grid_spacing, data=lapl)

    def laplacian_at_point(self, x, y):
        m, n, _ = self.data.shape
        dx, dy = self.grid_spacing
        center = self[x,y]

        dright, dleft, dup, ddown = 0, 0, 0, 0
        if(x != 0):
            dleft, _ = center - self[x-1,y]
        if(x != n-1):
            dright, _ = self[x+1,y] - center        
        if(y != 0):
            _, ddown = center - self[x,y-1]
        if(y != m-1):
            _, dup = self[x,y+1] - center
        
        return [(dright - dleft) / dx**2, (dup + ddown) / dy**2]

class ScalarGrid2:
    def __init__(self, resolution, grid_spacing, sample=None, data=None):
        self.resolution = resolution
        self.grid_spacing = grid_spacing
        self.sample = sample

        x_points = int(resolution[0]/float(grid_spacing[0]))
        y_points = int(resolution[1]/float(grid_spacing[1]))
        x = np.linspace(0.5*grid_spacing[0], resolution[0] - 0.5 * grid_spacing[0], x_points)
        y = np.linspace(0.5*grid_spacing[1], resolution[1] - 0.5 * grid_spacing[1], y_points)
        
        self.x = x
        self.y = y
        self.data = np.zeros((y_points,x_points))
        self.xx, self.yy = np.meshgrid(x,y)
        if sample is not None:
            self.data = sample(self.xx,self.yy)
        elif data is not None:
            self.data = data
    
    def __repr__(self):
        return self.data

    def __str__(self):
        return str(self.data)

    def __getitem__(self, pos):
        '''Translates cartesian coord to matrice and return element
        '''
        x, y = pos
        return self.data[y][x]

    def gradient(self):
        m, n = self.data.shape

        grad = np.zeros((m,n,2))
        for i in range(1,m-1):
            for j in range(1,n-1):
                #iterate over columns which is equivalent to iterate over x axis
                grad[i][j] = self.gradient_at_point(j, i)
        
        return VectorGrid2(self.resolution, self.grid_spacing, data=grad)

    def	gradient_at_point(self, x, y):	
        m, n = self.data.shape
        dx, dy = self.grid_spacing
        center = self[x,y]

        right = self[x+1, y] if x != n-1 else center
        left = self[x-1, y] if x != 0 else center
        up = self[x, y+1] if y != m-1 else center
        down = self[x, y-1] if y != 0 else center
        
        return 0.5 * np.asarray([(center-left)/dx, (center-down)/dy])

    def laplacian(self):
        m, n = self.data.shape

        lapl = np.zeros_like(self.data)
        for i in range(m):
            for j in range(n):
                #iterate over columns which is equivalent to iterate over x axis
                lapl[i][j] = self.laplacian_at_point(j, i)
        
        return ScalarGrid2(self.resolution, self.grid_spacing, data=lapl)

    def laplacian_at_point(self, x, y):
        m, n = self.data.shape
        dx, dy = self.grid_spacing
        center = self[x,y]

        dright, dleft, dup, ddown = 0, 0, 0, 0
        if(x != 0):
            dleft = center - self[x-1,y]
        if(x != n-1):
            dright = self[x+1,y] - center
        if(y != 0):
            ddown = center - self[x,y-1]
        if(y != m-1):
            dup = self[x,y+1] - center
        
        return (dright - dleft) / dx**2 + (dup + ddown) / dy**2

def f1(x, y):
    # return x+y
    # return x**2 * np.sin(y)
    return np.sin(x) * np.sin(y)

def f2(x, y):
    # return np.asarray([1,1])
    # return np.asarray([x,y])
    return np.asarray([np.sin(x),np.sin(y)])
    # return np.asarray([-y, x])

def main():
    grids_res = (5,5)
    grids_spc = (0.1,0.1)

    grid = ScalarGrid2(grids_res, grids_spc, f1)
    grad = grid.gradient()
    lapl = grid.laplacian()
    plt.title("Scalar field gradient laplacian")
    plt.subplot(2,2,1)
    plt.title("sin(x)*sin(y)")
    plt.pcolormesh(grid.x, grid.y, grid.data, shading='auto')
    plt.colorbar()
    plt.subplot(2,2,2)
    plt.title("laplacian")
    plt.pcolormesh(lapl.x, lapl.y, lapl.data, shading='auto')
    plt.colorbar()
    plt.subplot(2,2,3)
    plt.title("gradient")
    grad.plot(True)
    plt.show()
    
    vec_grid = VectorGrid2(grids_res, grids_spc, f2)
    div = vec_grid.divergent()
    lapl = vec_grid.laplacian()
    plt.title("Vector field divergent laplacian")
    plt.subplot(2,2,1)
    plt.title("[sin(x),sin(y)]")
    vec_grid.plot(True)
    plt.subplot(2,2,2)
    plt.title("laplacian")
    lapl.plot(True)
    plt.subplot(2,2,3)
    plt.title("divergent")
    plt.pcolormesh(div.x, div.y, div.data, shading='auto')
    plt.colorbar()
    plt.show()

    

if __name__ == '__main__':
    main()