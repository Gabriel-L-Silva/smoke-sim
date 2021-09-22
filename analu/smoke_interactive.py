# -----------------------------------------------------------------------------
# From Numpy to Python
# Copyright (2017) Nicolas P. Rougier - BSD license
# More information at https://github.com/rougier/numpy-book
# -----------------------------------------------------------------------------

import numpy as np 
from smoke_solver import vel_step, dens_step


N = 64 #tamanho do grid
size = N + 2 #tamanho efetivo do grid, com ghost cells
dt = 0.1 #passo de tempo
diff = 0.0 #o que é dif??
visc = 0.0 #ignoramos esse passo
force = -9.8 #gravidade, posso adicionar, adicionar força de buoynce
source = 100.0 #porque definida assim?
dvel = False # e agora?

mouse = {"ox": 0.0, "oy": 0.0,
         "x": 0.0,  "y": 0.0,
         "button": None}  #dicionario mouse

u = np.zeros((size, size), np.float32)  # velocity
u_prev = np.zeros((size, size), np.float32)

v = np.zeros((size, size), np.float32)  # velocity
v_prev = np.zeros((size, size), np.float32)

dens = np.zeros((size, size), np.float32)  # density
dens_prev = np.zeros((size, size), np.float32)

#definição do grid
def initialization():
    global u, v, u_prev, v_prev, dens, dens_prev, size #u_prev, v_prev: campos anteriores (swap)

    u[:, :] = 0.0 #do começo ao fim do array numpy em 2 dimensoes 
    v[:, :] = 0.0
    u_prev[:, :] = 0.0
    v_prev[:, :] = 0.0
    dens[:, :] = 0.0
    dens_prev[:, :] = 0.0


def user_step(d, u, v): #chama a função user_step toda vez que o usuario executa uma ação qualquer
    global mouse

    d[:, :] = 0.0
    u[:, :] = 0.0
    v[:, :] = 0.0

    if mouse["button"] not in [1, 3]:
        #se não executou nenhuma das 3 ações, não faz nada, só retorna
        return
    if mouse["x"] is None or mouse["y"] is None:
        #se clicou fora para o programa
        return

    i = int(mouse["y"]*N) + 1 # posição no grid, é fino mas nem tanto, poderia ser mais, tentar
    j = int(mouse["x"]*N) + 1
    # i e j vão de 1 a 64, as ghost cells sao 0 4 65
    if not 0 < i < N+1 and not 0 < j < N+1:
        #se nao clicou dentro do grid também sai
        return

    #define as ações do mouse
    if mouse["button"] == 3:
        #right click: adiciona uma nova fonte aquele ponto
        d[i, j] = source
    elif mouse["button"] == 1:
        #left click: exerce uma força proporcional ao esforço???
        u[i, j] = force * (mouse["y"] - mouse["oy"])*200 # gravidade atua somente na direção y
        #pensar em como adicionar a força de empuxo
        v[i, j] =  (mouse["x"] - mouse["ox"])*200
        

#nao sei o que é ox, imaginei que fosse a origem
# a nova pos do mouse é igual a anterior? CONFUSO
    mouse["ox"] = mouse["x"]
    mouse["oy"] = mouse["y"]

# *args: passa um número aleatório de argumentos
# se fosse **kwargs: um dicionario
def update(*args):
    #o que é im?
    global im, dens, dens_prev, u, u_prev, v, v_prev, N, visc, dt, diff
    #DIF de difusao? 

    user_step(dens_prev, u_prev, v_prev) #usuario executou uma ação
    vel_step(N, u, v, u_prev, v_prev, visc, dt) #quando o tempo passa, para velocidade
    dens_step(N, dens, dens_prev, u, v, diff, dt) #quando o tempo passa atualiza densidade
    im.set_data(dens) #o que é set_data
    # im.set_clim(vmin=dens.min(), vmax=dens.max())

#funões de evento
#o id é atribuido automaticamente? ou tem algum mapa específico
def on_button_press(event):
    global mouse
    mouse["ox"] = mouse["x"] = event.xdata
    mouse["oy"] = mouse["y"] = event.ydata
    mouse["button"] = event.button


def on_button_release(event):
    global mouse
    mouse["ox"] = mouse["x"] = event.xdata
    mouse["oy"] = mouse["y"] = event.ydata
    mouse["button"] = None


def on_motion(event):
    global mouse
    mouse["x"] = event.xdata
    mouse["y"] = event.ydata


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    #criando a figura e adiocionando os eixos

    cid = fig.canvas.mpl_connect('button_press_event', on_button_press)
    cid = fig.canvas.mpl_connect('button_release_event', on_button_release)
    cid = fig.canvas.mpl_connect('motion_notify_event', on_motion)
    #mpl_connect mapeia para um inteiro ou eu que faço?
    # onde ele muda esse cid ao longo do código

    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_ylim(0, 1)
    ax.set_yticks([])

    initialization()
    #im é a imagem que vamos trabalhar em cima, im.set_data() atualiza a im a cada chamada
    im = ax.imshow(dens[1:-1, 1:-1],
                   interpolation='bicubic', extent=[0, 1, 0, 1],
                   cmap=plt.cm.gray, origin="lower", vmin=0, vmax=1)
    #[1:-1, 1:-1]: vai da primeira ate a ultima posiçao do array
    #extent: espaço que a imagem vai ocupar
    #vmin e vmax: range do color map
    animation = FuncAnimation(fig, update, interval=10, frames=800)
    #interval: delay entre os frames em milissegundos, super rápido
    #o que são os frames?
    plt.show()
