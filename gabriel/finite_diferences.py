import numpy as np
import matplotlib.pyplot as plt
import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable

class StaggeredVectorGrid2:
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

    def __setitem__(self, pos, item):
        '''Translates cartesian coord to matrice and return element

        parameters:
            - pos: cartesian coordinates of element
            - item: item to be assigned
        '''
        x, y = pos
        self.data[y][x] = item

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
                yield (self.xx[j], self.yx[i], self.u[i][j]), (self.xy[j], self.yy[i],self.v[i][j]), i, j

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

    def plot(self, ax, stream=False, color='k'):
        if stream:
            '''
            Do not draw lines in upper and right boundary
            '''
            speed = np.sqrt(self.u**2 + self.v**2)
            lw = 3 * speed / speed.max()
            ax.streamplot(self.xx,self.yy,self.u,self.v, color=color, density=0.6, linewidth=lw)
        else:
            ax.quiver(self.xx, self.yx, self.u, np.zeros_like(self.u),color=color)
            ax.quiver(self.xy, self.yy, np.zeros_like(self.v), self.v,color=color)
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

class VectorGrid2:
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
        self.data = np.zeros((y_points,x_points,2))
        xx, yy = np.meshgrid(x,y)
        if sample is not None:
            for i in range(y_points):
                for j in range(x_points):
                    self.data[i][j] = sample(x[j],y[i])
        elif data is not None:
            self.data = data
    
    def __repr__(self):
        return self.data

    def __str__(self):
        return str(self.data)

    def __setitem__(self, pos, item):
        '''Translates cartesian coord to matrice and return element

        parameters:
            - pos: cartesian coordinates of element
            - item: item to be assigned
        '''
        x, y = pos
        self.data[y][x] = item

    def __getitem__(self, pos):
        '''Translates cartesian coord to matrice and return element

        parameters:
            - pos: cartesian coordinates of element
        '''
        x, y = pos
        return self.data[y][x]

    def __iter__(self):
        for i in range(len(self.y)):
            for j in range(len(self.x)):
                yield self.x[j], self.y[i], j, i, self[i,j]

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
            ) / ((x2 - x1) * (y2 - y1))

    def plot(self, ax, stream=False, color='black'):
        if stream:
            '''
            Do not draw lines in upper and right boundary
            '''
            speed = np.sqrt(np.sum(np.square(self.data),2))
            lw = 3 * speed / speed.max()
            ax.streamplot(self.x,self.y,self.data[:,:,0],self.data[:,:,1], color=color, density=0.6, linewidth=lw)
        else:
            ax.quiver(self.x, self.y, self.data[:,:,0],self.data[:,:,1],color=color)
        ax.set_xticks(self.x-0.5*self.grid_spacing[0])
        ax.set_yticks(self.y-0.5*self.grid_spacing[1])
        ax.grid(True)

    def gradient(self):
        #TODO make divergent
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

class ScalarGrid2:
    def __init__(self, resolution, grid_spacing, sample=None, data=None):
        self.resolution = resolution
        self.grid_spacing = grid_spacing
        self.sample = sample

        self.x_points = int(resolution[0]/float(grid_spacing[0]))
        self.y_points = int(resolution[1]/float(grid_spacing[1]))
        x = np.linspace(0.5*grid_spacing[0], resolution[0] - 0.5 * grid_spacing[0], self.x_points)
        y = np.linspace(0.5*grid_spacing[1], resolution[1] - 0.5 * grid_spacing[1], self.y_points)
        
        self.x = x
        self.y = y
        self.data = np.zeros((self.y_points,self.x_points))
        self.xx, self.yy = np.meshgrid(x,y)
        if sample is not None:
            self.data = sample(self.xx, self.yy)
        elif data is not None:
            self.data = data
    
    def __repr__(self):
        return self.data

    def __str__(self):
        return str(self.data)

    def __setitem__(self, pos, item):
        '''Translates cartesian coord to matrice and return element

        parameters:
            - pos: cartesian coordinates of element
            - item: item to be assigned
        '''
        x, y = pos
        self.data[y][x] = item

    def __getitem__(self, pos):
        '''Translates cartesian coord to matrice and return element

        parameters:
            - pos: cartesian coordinates of element
        '''
        x, y = pos
        return self.data[y][x]

    def __iter__(self):
        '''Return a quintuple with (x coord, y coord, x index, y index, value)
        '''
        for i in range(len(self.y)):
            for j in range(len(self.x)):
                yield self.x[j], self.y[i], j, i, self.data[i,j]

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

    def plot(self, ax, colorbar=False, shade='auto'):
        p = ax.pcolormesh(self.x, self.y, self.data, shading=shade)
        if colorbar: 
            divider = make_axes_locatable(ax)

            ax_cb = divider.new_horizontal(size="5%", pad=0.05)
            fig = ax.get_figure()
            fig.add_axes(ax_cb)
            
            plt.colorbar(p, cax=ax_cb)
            ax_cb.yaxis.tick_right()
            # ax_cb.yaxis.set_tick_params(labelright=False)

    def plot_surf(self, ax, colorbar=False, shade='auto'):
        p = ax.plot_surface(self.xx, self.yy, self.data,cmap='viridis')
        if colorbar: 
            divider = make_axes_locatable(ax)

            ax_cb = divider.new_horizontal(size="5%", pad=0.05)
            fig = ax.get_figure()
            fig.add_axes(ax_cb)
            
            plt.colorbar(p, cax=ax_cb)
            ax_cb.yaxis.tick_right()
            # ax_cb.yaxis.set_tick_params(labelright=False)

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

def backtrace_mass(sca_field, vec_field, timestep, sub_steps) -> ScalarGrid2:
    w = copy.deepcopy(sca_field)
    s = timestep/sub_steps
    for cell in sca_field:
        pos = cell[0], cell[1]
        i, j = cell[2], cell[3]
        if(pos[0] <= 0.5*sca_field.grid_spacing[0] or pos[1] <= 0.5*sca_field.grid_spacing[1]
            or pos[0] >= sca_field.resolution[0] - 0.5*sca_field.grid_spacing[0]
            or pos[1] >= sca_field.resolution[1] - 0.5*sca_field.grid_spacing[1]):
            continue
        for _ in range(sub_steps):
            pos = pos - s*vec_field.interpolate(pos[0], pos[1])
            
        w.data[i,j] = sca_field.interpolate(pos[0], pos[1])
    return w

def backtrace_velocity(vec_field, timestep, sub_steps) -> VectorGrid2:
    w = copy.deepcopy(vec_field)
    s = timestep/sub_steps
    for cell in vec_field:
        pos_u = cell[0][0], cell[0][1]
        pos_v = cell[1][0], cell[1][1]
        i, j = cell[2], cell[3]
        if(pos_u[0] <= 0.5*vec_field.grid_spacing[0] or pos_u[1] <= 0.5*vec_field.grid_spacing[1]
            or pos_u[0] >= vec_field.resolution[0] - vec_field.grid_spacing[0]
            or pos_u[1] >= vec_field.resolution[1] - vec_field.grid_spacing[1]):
            continue
        for _ in range(sub_steps):
            pos_u = pos_u - s*vec_field.interpolate(pos_u[0], pos_u[1])
            pos_v = pos_v - s*vec_field.interpolate(pos_v[0], pos_v[1])
        w.data[i,j] = vec_field.interpolate(pos_u[0], pos_u[1])[0], vec_field.interpolate(pos_v[0], pos_v[1])[1]
    return w

def backtrace_interactive(vec_field, ax, x, y, timestep, sub_steps):    
    ix, iy = vec_field.interpolate(x,y)
    ax.quiver(x, y, ix, iy, color='red')
    plt.pause(0.05)
    s = timestep/sub_steps
    pos_u = x, y
    for _ in range(sub_steps):
        pos_u = pos_u - s*vec_field.interpolate(pos_u[0], pos_u[1])
        ax.scatter(pos_u[0], pos_u[1], c='red')
        plt.pause(0.005)
    return pos_u

def f1(x, y):
    return x+y
    # return x**2 * np.sin(y)
    # return np.sin(x) * np.sin(y)

def f2(x, y):
    # return np.asarray([1,1])
    # return np.asarray([x,y])
    return np.asarray([np.sin(x),np.sin(y)])
    # return np.asarray([-y, x])

def poisson_exerc(x, y):
    return -5/4.0 * np.pi**2 * np.sin(np.pi*x) * np.cos(0.5*np.pi*y)

def poisson_exerc_solution(x,y):
    return np.sin(np.pi*x)*np.cos(0.5*np.pi*y)

def poisson_iter_solver(sca_grid: ScalarGrid2, x0, tol):
    max_it = 10000
    m, n = sca_grid.x_points, sca_grid.y_points
    b = np.zeros(m*n)
    dx = sca_grid.grid_spacing[0]
    for idx, cell in enumerate(sca_grid):
        x, y, idx_x, idx_y, p = cell
        if (idx_y != m-1 and idx_y != 0 and idx_x != n-1 and idx_x != 0):
            b[idx] = dx**2*p
        elif idx_y == 0:
            b[idx] = np.sin(np.pi * x)

    b = np.diag(b.reshape((m,n)))
    A = sca_grid.data
    for k in range(max_it):
        x_old  = x0.copy()
        for i in range(m):
            x0[i] = (b[i] - np.dot(A[i,:i], x0[:i]) - np.dot(A[i,(i+1):], x_old[(i+1):])) / A[i ,i]
    return A*x0

def poisson_iter_solver_2(sca_grid: ScalarGrid2, x0, tol):
    max_it = 10000
    m, n = sca_grid.x_points, sca_grid.y_points
    b = np.zeros(m*n)
    dx = sca_grid.grid_spacing[0]
    for idx, cell in enumerate(sca_grid):
        x, y, idx_x, idx_y, p = cell
        if idx_y == 0:
            x0[idx_y, idx_x] = np.sin(np.pi * x)

    A = sca_grid.data
    for _ in range(max_it):
        old = x0.copy()
        for i in range(1,n-1):
            for j in range(1,n-1):
                x0[i,j] = 0.25*(x0[i+1,j] + x0[i-1,j] + x0[i,j+1] + x0[i,j-1] - A[i,j]*dx**2)
        
        if (abs(old - x0) < tol).all():
            break
    return x0

def poisson_solver(sca_grid: ScalarGrid2):
    m, n = sca_grid.x_points, sca_grid.y_points
    mat = -4 * np.eye(m * n)+0
    b = np.zeros(m*n)
    dx = sca_grid.grid_spacing[0]
    for idx, cell in enumerate(sca_grid):
        x, y, idx_x, idx_y, p = cell
        if (idx_y != m-1 and idx_y != 0 and idx_x != n-1 and idx_x != 0):
            mat[idx,idx_y*m+idx_x+1] = 1 #right
            mat[idx,idx_y*m+idx_x-1] = 1 #left
            mat[idx,(idx_y+1)*m+idx_x] = 1 #up
            mat[idx,(idx_y-1)*m+idx_x] = 1 #down
            b[idx] = dx**2*p
        else:
            mat[idx] = np.zeros(m*n)
            mat[idx,idx] = 1

            if idx_y == 0:
                b[idx] = np.sin(np.pi * x)
            
    
    phi = np.linalg.solve(mat, b).reshape(sca_grid.data.shape)
    return 

def rmse(result, ref):
    '''Compare two values with rmse
    output:
        -rse: root mean squared error of a against b
    '''
    if type(ref) == np.ndarray:
        size = np.multiply.reduce(ref.shape)
    else:
        size = 1
    
    return np.sqrt(np.sum(np.square(result - ref))/size)

def test_interpolate(grid, fig, ax):
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
    grid.plot(ax)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show(block=True)
    return (ev.xdata, ev.ydata) if ev is not None else None
    # return (ev.x, ev.y) if ev is not None else None

def test_backtrace():
    '''interative plot to calculate interpolation on clicked coord
    '''
    grids_res = (1,1)
    grids_spc = (.031,.031)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    sca_grid = ScalarGrid2(grids_res, grids_spc, f1)
    vec_grid = VectorGrid2(grids_res, grids_spc, f2)

    timestep, sub_steps = .5, 10

    ev = None
    def onclick(event):
        nonlocal ev
        ev = event
        pos = backtrace_interactive(vec_grid, ax, ev.xdata, ev.ydata, timestep, sub_steps) 
        value = sca_grid.interpolate(pos[0], pos[1])
        ax.scatter(ev.xdata, ev.ydata)
        ax.annotate("%.3f"%value, (ev.xdata,ev.ydata))

    plt.ion()
    
    vec_grid.plot(ax,True)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    plt.show(block=True)
    return (ev.xdata, ev.ydata) if ev is not None else None
    # return (ev.x, ev.y) if ev is not None else None

def poisson_term(w3: VectorGrid2, h, rho, pa):
    w4 = copy.deepcopy(w3)
    dx = w3.resolution[0]
    poisson_coefs = np.asarray([[-4,1,1,1,0,0,1,0,0],
            [1,-4,1,0,1,0,0,1,0],
            [1,1,-4,0,0,1,0,0,1],
            [1,0,0,-4,1,1,1,0,0],
            [0,1,0,1,-4,1,0,1,0],
            [0,0,1,1,1,-4,0,0,1],
            [1,0,0,1,0,0,-4,1,1],
            [0,1,0,0,1,0,1,-4,1],
            [0,0,0,0,0,0,0,0,1]])
    u = dx * w3.u.flatten()
    v = dx * w3.v.flatten()
    b = np.zeros_like(u)
    b[0] = u[1] - u[0] + v[3] - v[0]
    b[1] = u[2] - u[1] + v[4] - v[1]
    b[2] = u[0] - u[2] + v[5] - v[2]
    b[3] = u[4] - u[3] + v[6] - v[3]
    b[4] = u[5] - u[4] + v[7] - v[4]
    b[5] = u[3] - u[5] + v[8] - v[5]
    b[6] = u[7] - u[6] + v[0] - v[6]
    b[7] = u[8] - u[7] + v[1] - v[7]
    b[8] = h/(dx*rho)*pa

    #TODO checar se u e v continuam consistentes com o campo criado
    sca_field = ScalarGrid2(w3.resolution, w3.grid_spacing, data=np.linalg.solve(poisson_coefs, b).reshape(w3.resolution))
    w4.data = w3.data - sca_field.gradient().data
    return w4, sca_field

def diffusion_term(w2: VectorGrid2, h, nu) -> VectorGrid2:
    w3 = copy.deepcopy(w2)
    dx = w2.grid_spacing[0]
    k = 4+dx**2/(h*nu)
    coefs = np.asarray([[k,-1,-1,-1,0,0,-1,0,0],
            [-1,k,-1,0,-1,0,0,-1,0],
            [-1,-1,k,0,0,-1,0,0,-1],
            [-1,0,0,k,-1,-1,-1,0,0],
            [0,-1,0,-1,k,-1,0,-1,0],
            [0,0,-1,-1,-1,k,0,0,-1],
            [-1,0,0,-1,0,0,k,-1,-1],
            [0,-1,0,0,-1,0,-1,k,-1],
            [0,0,-1,0,0,-1,-1,-1,k]])
    u = dx**2/(h*nu) * w2.u.flatten()
    v = dx**2/(h*nu) * w2.v.flatten()
    w3.u = np.linalg.solve(coefs, u).reshape(w3.resolution)
    w3.v = np.linalg.solve(coefs, v).reshape(w3.resolution)
    for i in range(w3.resolution[0]):
        for j in range(w3.resolution[1]):
            w3.data[i][j] = w3.u[i][j], w3.v[i][j]
    return w3

def print_mat(mat: np.ndarray):
    m,n = mat.shape
    for i in range(m):
        for j in range(n):
            print(mat[i,j].astype(int),end=', ')
        print()

def test_poisson_solver_iter():
    grids_res = (1,1)
    grids_spc = (.031,.031)

    # set up a figure twice as wide as it is tall
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')

    sca_grid = ScalarGrid2(grids_res, grids_spc, poisson_exerc)
    ax1.set_title("Erro")
    ax2.set_title("Nossa solução")
    ax3.set_title("Referencia")
    poisson_result = ScalarGrid2(sca_grid.resolution, sca_grid.grid_spacing,
                                data=poisson_iter_solver_2(sca_grid, np.zeros_like(sca_grid.data), 10e-5))
    poisson_result.plot_surf(ax2)
    poisson_sol = ScalarGrid2(grids_res, grids_spc, poisson_exerc_solution)
    poisson_sol.plot_surf(ax3)

    diff = poisson_result.data - poisson_sol.data
    error = ScalarGrid2(grids_res, grids_spc, data=diff)
    error.plot(ax1, True)
    print(rmse(poisson_sol.data, poisson_result.data))

def test_poisson_solver():
    grids_res = (1,1)
    grids_spc = (.031,.031)

    # set up a figure twice as wide as it is tall
    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')

    sca_grid = ScalarGrid2(grids_res, grids_spc, poisson_exerc)
    ax1.set_title("Erro")
    ax2.set_title("Nossa solução")
    ax3.set_title("Referencia")
    poisson_result = ScalarGrid2(sca_grid.resolution, sca_grid.grid_spacing, data=poisson_solver(sca_grid))
    poisson_result.plot_surf(ax2)
    poisson_sol = ScalarGrid2(grids_res, grids_spc, poisson_exerc_solution)
    poisson_sol.plot_surf(ax3)

    diff = poisson_result.data - poisson_sol.data
    error = ScalarGrid2(grids_res, grids_spc, data=diff)
    error.plot(ax1, True)
    print(rmse(poisson_sol.data, poisson_result.data))

def main():
    
    # test_poisson_solver_iter()
    test_backtrace()
    
    
    
    # ax1.set_title("Campo escalar")
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
    plt.show()
    
    # vec_grid = VectorGrid2(grids_res, grids_spc, f2)
    # vec_grid.plot(ax1)
    # plt.show()
    # w2 = backtrace_velocity(vec_grid, .7, 10)
    # w3 = diffusion_term(w2, .7, 1)
    # w4, phi = poisson_term(w3, .7, 1, 1)
    # test_interpolate(phi,fig,ax2)
    # plt.show(block=True)
    # test_interpolate(w2, fig, ax2)
    # vec_grid.plot(ax, True)    
    # for cell in vec_grid:
    #     print(cell)
    # test_backtrace(sca_grid, vec_grid, fig, ax2, .5, 10)
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