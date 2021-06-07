import numpy as np
import matplotlib.pyplot as plt
import copy

class VectorGrid2:
    '''Class representing staggered grid, u and v component are stored at face
    centers, left and bottom respectively

    atrib:
        - xx: x coordinate in relation to x velocity
        - yx: y coordinate in relation to x velocity
        - xy: x coordinate in relation to y velocity
        - yy: y coordinate in relation to y velocity
                          
                          CELL
        -------------------------------------------------                 
        |                                               |                   
        |                                               |                   
        |                                               |                   
        |                                               |                   
        |                                               |                   
        |                                               |                   
        |                                               |                   
        |      U                                        |                   
(xx,yx) |----------->                                   |                   
        |                                               |                   
        |                                               |                   
        |                                               |                   
        |                                               |                   
        |                         A                     |                   
        |                         |                     |                   
        |                         |  V                  |                   
        |                         |                     |                   
        |                         |                     |                   
        -------------------------------------------------                
                                (xy,yy)                              
                                                    
    '''
    def __init__(self, resolution, grid_spacing, sample=None, data=None):
        self.resolution = resolution
        self.grid_spacing = grid_spacing
        self.sample = sample

        x_points = int(resolution[0]/float(grid_spacing[0]))
        y_points = int(resolution[1]/float(grid_spacing[1]))
        x = np.linspace(0.5*grid_spacing[0], resolution[0] - 0.5 * grid_spacing[0], x_points)
        y = np.linspace(0.5*grid_spacing[1], resolution[1] - 0.5 * grid_spacing[1], y_points)
        
        # TODO melhorar o jeito que estão salvas as coords
        self.xx = x - 0.5*grid_spacing[0]
        self.yx = y 
        self.xy = x 
        self.yy = y - 0.5*grid_spacing[1]
        self.data = np.zeros((y_points,x_points,2))
        self.u = np.zeros((y_points, x_points))
        self.v = np.zeros((y_points, x_points))
        if sample is not None:
            for i in range(y_points):
                for j in range(x_points):
                    #iterate over columns which is equivalent to iterate over x axis
                    self.u[i][j] = sample(self.xx[j],self.yx[i])[0]
                    self.v[i][j] = sample(self.xy[j],self.yy[i])[1]
                    self.data[i][j] = self.u[i][j], self.v[i][j]
        elif data is not None:
            self.data = data.copy()
            m, n, _ = self.data.shape
            self.u = np.zeros((m,n))
            self.v = np.zeros((m,n))
            for i in range(m):
                for j in range(n):
                    self.u[i][j], self.v[i][j] = self.data[i][j]

    def __repr__(self):
        return self.data

    def __str__(self):
        return str(self.data)

    def __getitem__(self, pos):
        '''Translates cartesian coord to matrice and return element

        parameters:
            - pos: cartesian coordinates of element
        '''
        x, y = pos
        return self.data[y][x]

    def __iter__(self):
        x_points, y_points = len(self.xx), len(self.yy)
        for i in range(y_points):
            for j in range(x_points):
                yield (self.xx[j], self.yx[i], self.u[i][j]), (self.xy[j], self.yy[i], self.v[i][j])

    def interpolate(self, x, y) -> np.ndarray:
        '''Bilinear interpolation of staggered grid, not treating border case
        '''
        if(x < 0.5*self.grid_spacing[0] or y < 0.5*self.grid_spacing[1]
            or x > self.resolution[0] - 0.5*self.grid_spacing[0]
            or y > self.resolution[1] - 0.5*self.grid_spacing[1]):
            print("no border case")
            return

        # para u
        dx, dy = self.grid_spacing
        x0, y0 = 0, 0 #origin of grid
        ui = int((x-x0)/dx)
        uj = int((y-y0)/dy-0.5)
        
        xi = self.xx[ui]
        xj = self.yx[uj]
        x1, y1, q11 = xi, xj, self[(ui),(uj)][0]
        x2, _y1, q21 = xi+self.grid_spacing[0], xj, self[(ui+1),(uj)][0]
        _x1, y2, q12 = xi, xj+self.grid_spacing[1], self[(ui),(uj+1)][0]
        _x2, _y2, q22 = xi+self.grid_spacing[0], xj+self.grid_spacing[1], self[(ui+1),(uj+1)][0]

        if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
            raise ValueError('points do not form a rectangle')
        if not x1 <= x <= x2 or not y1 <= y <= y2:
            raise ValueError('(x, y) not within the rectangle')

        u = (q11 * (x2 - x) * (y2 - y) +
                q21 * (x - x1) * (y2 - y) +
                q12 * (x2 - x) * (y - y1) +
                q22 * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1) + 0.0)
        
        # para v
        dx, dy = self.grid_spacing
        x0, y0 = 0, 0 #origin of grid
        vi = int((x-x0)/dx-0.5)
        vj = int((y-y0)/dy)
        
        yi = self.xy[vi] 
        yj = self.yy[vj]
        x1, y1, q11 = yi, yj, self[(vi),(vj)][1]
        x2, _y1, q21 = yi+self.grid_spacing[0], yj, self[(vi+1),(vj)][1]
        _x1, y2, q12 = yi, yj+self.grid_spacing[1], self[(vi),(vj+1)][1]
        _x2, _y2, q22 = yi+self.grid_spacing[0], yj+self.grid_spacing[1], self[(vi+1),(vj+1)][1]

        if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
            raise ValueError('points do not form a rectangle')
        if not x1 <= x <= x2 or not y1 <= y <= y2:
            raise ValueError('(x, y) not within the rectangle')

        v = (q11 * (x2 - x) * (y2 - y) +
                q21 * (x - x1) * (y2 - y) +
                q12 * (x2 - x) * (y - y1) +
                q22 * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1) + 0.0)

        return np.asarray([u,v])

    def plot(self, ax, stream=False):
        if stream:
            '''
            Do not draw lines in upper and right boundary
            '''
            speed = np.sqrt(self.u**2 + self.v**2)
            lw = 3 * speed / speed.max()
            ax.streamplot(self.xx,self.yy,self.u,self.v, color='k', density=0.6, linewidth=lw)
        else:
            ax.quiver(self.xx, self.yx, self.u, np.zeros_like(self.u))
            ax.quiver(self.xy, self.yy, np.zeros_like(self.v), self.v)
        ax.set_xticks(self.xx)
        ax.set_yticks(self.yy)
        ax.grid(True)
        
    def divergent(self):
        m, n, _ = self.data.shape
        
        div = np.zeros((m,n))
        for i in range(m):
            for j in range(n):
                #iterate over columns which is equivalent to iterate over x axis
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

        parameters:
            - pos: cartesian coordinates of element
        '''
        x, y = pos
        return self.data[y][x]

    def interpolate(self, x, y):
        '''Bilinear interpolation of cell centered grid, not treating border case
        '''
        if(x < 0.5*self.grid_spacing[0] or y < 0.5*self.grid_spacing[1]
            or x > self.resolution[0] - 0.5*self.grid_spacing[0]
            or y > self.resolution[1] - 0.5*self.grid_spacing[1]):
            print("no border case")
            return

        dx, dy = self.grid_spacing
        x0, y0 = 0.5*self.grid_spacing[0], 0.5*self.grid_spacing[1] #origin of grid
        i = int((x-x0)//dx)
        j = int((y-y0)//dy)
        
        xi = self.x[i]
        yj = self.y[j]
        x1, y1, q11 = xi, yj, self[(i),(j)]
        x2, _y1, q21 = xi+self.grid_spacing[0], yj, self[(i+1),(j)]
        _x1, y2, q12 = xi, yj+self.grid_spacing[1], self[(i),(j+1)]
        _x2, _y2, q22 = xi+self.grid_spacing[0], yj+self.grid_spacing[1], self[(i+1),(j+1)]

        if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
            raise ValueError('points do not form a rectangle')
        if not x1 <= x <= x2 or not y1 <= y <= y2:
            raise ValueError('(x, y) not within the rectangle')

        return (q11 * (x2 - x) * (y2 - y) +
                q21 * (x - x1) * (y2 - y) +
                q12 * (x2 - x) * (y - y1) +
                q22 * (x - x1) * (y - y1)
            ) / ((x2 - x1) * (y2 - y1) + 0.0)

    def plot(self, ax, shade='auto'):
        return ax.pcolormesh(self.x, self.y, self.data, shading=shade)

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

def backtrace(vec_field, ax, x, y, timestep, sub_steps):
    w = copy.deepcopy(vec_field)
    
    ix, iy = vec_field.interpolate(x,y)
    ax.quiver(x, y, ix, iy, color='red')
    plt.pause(0.05)
    s = timestep/sub_steps
    pos_u = x, y
    for _ in range(sub_steps):
        pos_u = pos_u - s*vec_field.interpolate(pos_u[0], pos_u[1])
        ax.scatter(pos_u[0], pos_u[1], c='red')
        plt.pause(0.05)
    return w, pos_u

def f1(x, y):
    return x+y
    # return x**2 * np.sin(y)
    # return np.sin(x) * np.sin(y)

def f2(x, y):
    # return np.asarray([1,1])
    # return np.asarray([x,y])
    return np.asarray([np.sin(x),np.sin(y)])
    # return np.asarray([-y, x])

def rmse(result, ref):
    '''Compare two values with rmse
    output:
        -rse: root mean squared error of a against b
    '''
    if type(ref) == np.ndarray:
        size = ref.shape[0]
    else:
        size = 1 
    return np.sqrt(np.sum(np.square(result - ref))/size)

def test_interpolate(grid):
    '''interative plot to calculate interpolation on clicked coord
    '''
    ev = None
    def onclick(event):
        nonlocal ev
        ev = event
        if type(grid) == VectorGrid2:
            u, v = grid.interpolate(ev.xdata,ev.ydata)
            ref = f2(ev.xdata,ev.ydata)
            ax.quiver(ev.xdata,ev.ydata, u, v, color='red')
            print("RMSE: ", rmse([u,v],ref))
            print(ref, "   ", [u,v])
        elif type(grid) == ScalarGrid2:
            res = grid.interpolate(ev.xdata,ev.ydata)
            ref = f1(ev.xdata,ev.ydata)
            ax.scatter(ev.xdata,ev.ydata)
            ax.annotate("%.3f"%res, (ev.xdata,ev.ydata))
            print("RMSE: ", rmse(ref, res))
            print(ref, "   ", res)

    plt.ion()
    fig, ax = plt.subplots()
    grid.plot(ax)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show(block=True)
    return (ev.xdata, ev.ydata) if ev is not None else None
    # return (ev.x, ev.y) if ev is not None else None

def test_backtrace(sca_grid, vec_grid, fig, ax, timestep, sub_steps):
    '''interative plot to calculate interpolation on clicked coord
    '''
    ev = None
    def onclick(event):
        nonlocal ev
        ev = event
        w, pos = backtrace(vec_grid, ax, ev.xdata, ev.ydata, timestep, sub_steps) 
        value = sca_grid.interpolate(pos[0], pos[1])
        ax.scatter(ev.xdata, ev.ydata)
        ax.annotate("%.3f"%value, (ev.xdata,ev.ydata))

    plt.ion()
    
    vec_grid.plot(ax,True)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show(block=True)
    return (ev.xdata, ev.ydata) if ev is not None else None
    # return (ev.x, ev.y) if ev is not None else None

def main():
    grids_res = (10,10)
    grids_spc = (.5,.5)

    fig, ax = plt.subplots()

    ax.set_xlim(0, grids_res[0])
    ax.set_ylim(0, grids_res[1])
    sca_grid = ScalarGrid2(grids_res, grids_spc, f1)
    plot = sca_grid.plot(ax, shade='gouraud')
    fig.colorbar(plot)
    plt.show(block=False)
    plt.pause(0.1)
    # test_interpolate(sca_grid)
    # lapl = sca_grid.laplacian()
    # plt.title("Scalar field gradient laplacian")
    # plt.subplot(2,2,1)
    # plt.title("sin(x)*sin(y)")
    # plt.pcolormesh(sca_grid.x, sca_grid.y, sca_grid.data, shading='auto')
    # plt.colorbar()
    # plt.subplot(2,2,2)
    # plt.title("laplacian")
    # plt.pcolormesh(lapl.x, lapl.y, lapl.data, shading='auto')
    # plt.colorbar()
    # plt.subplot(2,2,3)
    # plt.title("gradient")
    # grad.plot(True)
    # plt.show()
    
    vec_grid = VectorGrid2(grids_res, grids_spc, f2)
    # for cell in vec_grid:
    #     print(cell)
    test_backtrace(sca_grid, vec_grid, fig, ax, 1, 10)
    # test_interpolate(vec_grid)
    # div = vec_grid.divergent()
    # lapl = vec_grid.laplacian()
    # plt.title("Vector field divergent laplacian")
    # plt.subplot(2,2,1)
    # plt.title("[sin(x),sin(y)]")
    # vec_grid.plot(True)
    # plt.subplot(2,2,2)
    # plt.title("laplacian")
    # lapl.plot(True)
    # plt.subplot(2,2,3)
    # plt.title("divergent")
    # plt.pcolormesh(div.x, div.y, div.data, shading='auto')
    # plt.colorbar()
    # plt.show()

if __name__ == '__main__':
    main()